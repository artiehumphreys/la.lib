#include "lalib/scalar/type_traits.hpp"
#include <array>
#include <cassert>
#include <span>
#include <type_traits>

namespace lalib {

template <class T, std::size_t N> struct Vector {
  static_assert(N > 0, "vector size must be positive.");

  std::array<T, N> arr{};

  constexpr Vector() = default;

  template <class Iter> constexpr Vector(Iter begin, Iter end) {
    std::size_t k = 0;
    for (; begin != end && k < N; ++begin, ++k) {
      arr[k] = static_cast<T>(*begin);
    }
  }

  template <class U> constexpr Vector(std::span<const U, N> s) {
    for (std::size_t i = 0; i < N; ++i) {
      arr[i] = static_cast<T>(s[i]);
    }
  }

  constexpr T &operator[](std::size_t i) noexcept {
    assert(i >= 0 && i < N && "index out of bounds.");
    return arr[i];
  }

  constexpr const T &operator[](std::size_t i) const noexcept {
    assert(i >= 0 && i < N && "index out of bounds.");
    return arr[i];
  }

  constexpr void
  fill(const T &v) noexcept(std::is_nothrow_copy_assignable_v<T>) {
    for (T &a : arr) {
      a = v;
    }
  }

  // conditional noexcept on whether or not the resulting type can be
  // copy-assigned and default constructed w/o exceptions
  template <class U, std::size_t M, class R = std::common_type_t<T, U>>
  constexpr auto operator+(const Vector<U, M> &other) const
      noexcept(nothrow_element_v<R>) {
    static_assert(N == M, "vector dimensions must match");

    Vector<R, N> ans{};
    for (std::size_t i = 0; i < N; ++i) {
      ans[i] = static_cast<R>(arr[i]) + static_cast<R>(other.arr[i]);
    }

    return ans;
  }

  template <class U, std::size_t M>
  constexpr auto operator-(const Vector<U, M> &other) const
      noexcept(nothrow_element_v<signed_result_t<T, U>>) {
    static_assert(N == M, "vector dimensions must match");

    using R = lalib::signed_result_t<T, U>;

    Vector<R, N> neg{};
    for (std::size_t i = 0; i < N; ++i) {
      neg[i] = -static_cast<R>(other.arr[i]);
    }
    return *this + neg;
  }
};
} // namespace lalib
