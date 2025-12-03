#pragma once

#include "lalib/scalar/type_traits.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdlib>
#include <initializer_list>
#include <mdspan>
#include <span>
#include <type_traits>

namespace lalib {

template <class T, std::size_t N, std::size_t M> struct Matrix {
  static_assert(N > 0 && M > 0, "Matrix dimensions must be positive");

  std::array<T, N * M> arr{};

  constexpr Matrix() = default;

  template <class Iter> constexpr Matrix(Iter begin, Iter end) {
    // primarily runtime construction
    std::size_t i = 0;
    while (begin != end) {
      T value = static_cast<T>(*begin);
      arr[i++] = value;
      ++begin;
    }
  }

  template <class U> constexpr Matrix(std::span<const U, N * M> s) {
    std::size_t i = 0;
    for (const U &val : s) {
      arr[i++] = static_cast<T>(val);
    }
  }

  template <class U> constexpr Matrix(std::initializer_list<U> init) {
    // curly brace initialization
    std::size_t i = 0;
    for (const U &v : init) {
      arr[i++] = static_cast<T>(v);
    }
  }

  // constructor for multi-dimensional containers wrapped in an mdspan
  template <class U, class Extents, class Layout, class Accessor>
  constexpr Matrix(std::mdspan<const U, Extents, Layout, Accessor> s) {
    static_assert(Extents::rank() == 2);

    static_assert(Extents::static_extent(0) == N ||
                  Extents::static_extent(0) == std::dynamic_extent);
    static_assert(Extents::static_extent(1) == M ||
                  Extents::static_extent(1) == std::dynamic_extent);

    for (std::size_t r = 0; r < N; ++r) {
      for (std::size_t c = 0; c < M; ++c) {
        arr[r * M + c] = static_cast<T>(s(r, c));
      }
    }
  }

  static constexpr Matrix<T, N, N> identity() noexcept
    requires(N == M)
  {
    Matrix<T, N, N> ans{};
    for (std::size_t i = 0; i < N; ++i) {
      ans(i, i) = T{1};
    }

    return ans;
  }

  constexpr const T &operator[](std::size_t i) const noexcept {
    assert(i >= 0 && i < N * M);
    return arr[i];
  }

  constexpr T &operator[](std::size_t i) noexcept {
    assert(i >= 0 && i < N * M);
    return arr[i];
  }

  constexpr const T &operator()(std::size_t r, std::size_t c) const noexcept {
    assert(r < N && c < M);
    return arr[r * M + c];
  }

  constexpr T &operator()(std::size_t r, std::size_t c) noexcept {
    assert(r < N && c < M);
    return arr[r * M + c];
  }

  template <class U, std::size_t P, std::size_t Q,
            class R = std::common_type_t<T, U>>
  constexpr auto operator+(const Matrix<U, P, Q> &other) const noexcept {
    static_assert(N == P && M == Q, "matrix dimensions must match");

    Matrix<R, N, M> ans{};
    for (std::size_t i = 0; i < N * M; ++i) {
      ans[i] = static_cast<R>(arr[i]) + static_cast<R>(other.arr[i]);
    }
    return ans;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr Matrix &operator+=(const Matrix<U, P, Q> &other) noexcept {
    // in-place operations don't return a copy of the modified object (hence
    // Matrix& vs. auto)
    static_assert(N == P && M == Q, "matrix dimensions must be the same");
    for (std::size_t i = 0; i < N * M; ++i) {
      arr[i] += static_cast<T>(other.arr[i]);
    }
    return *this;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr auto operator-(const Matrix<U, P, Q> &other) const noexcept {
    static_assert(N == P && M == Q, "matrix dimensions must match");

    // force signed result in case of unsigned matrix types
    using R = lalib::signed_result_t<T, U>;

    Matrix<R, N, M> neg{};
    for (std::size_t i = 0; i < N * M; ++i) {
      neg.arr[i] = -static_cast<R>(other.arr[i]);
    }
    return (*this) + neg;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr Matrix &operator-=(const Matrix<U, P, Q> &other) noexcept {
    // do not modify matrix, even if operation leads to underflow
    static_assert(N == P && M == Q, "matrix dimensions must match");

    for (std::size_t i = 0; i < N * M; ++i) {
      arr[i] -= static_cast<T>(other.arr[i]);
    }
    return *this;
  }

  template <class U, std::size_t P, std::size_t Q,
            class R = std::common_type_t<T, U>>
  constexpr auto operator*(const Matrix<U, P, Q> &other) const noexcept {
    static_assert(M == P, "matrix dimensions must match");
    static constexpr std::size_t TILE_SIZE =
        std::max<std::size_t>(1, 64 / sizeof(R));

    Matrix<R, N, Q> ans{};

    for (std::size_t ii = 0; ii < N; ii += TILE_SIZE) {
      for (std::size_t kk = 0; kk < P; kk += TILE_SIZE) {
        for (std::size_t jj = 0; jj < Q; jj += TILE_SIZE) {

          // min is constexpr
          const std::size_t i_end = std::min(ii + TILE_SIZE, N);
          const std::size_t j_end = std::min(jj + TILE_SIZE, Q);
          const std::size_t k_end = std::min(kk + TILE_SIZE, P);

          for (std::size_t i = ii; i < i_end; ++i) {
            // i -> k -> j to ensure both matrices are accessed row-major
            for (std::size_t k = kk; k < k_end; ++k) {
              const R curr = static_cast<R>((*this)(i, k));

              for (std::size_t j = jj; j < j_end; ++j) {
                ans(i, j) += curr * static_cast<R>(other(k, j));
              }
            }
          }
        }
      }
    }
    return ans;
  }

  template <class Scalar> constexpr Matrix operator*(Scalar s) const noexcept {
    using R = decltype(std::declval<T>() * std::declval<Scalar>());
    // compile-time computation of promoted type

    Matrix<R, N, M> ans{};
    for (std::size_t i = 0; i < N * M; ++i) {
      ans[i] = static_cast<R>(arr[i] * s);
    }
    return ans;
  }

  template <class Scalar>
    requires(std::is_convertible_v<
             decltype(std::declval<T>() * std::declval<Scalar>()), T>)
  constexpr Matrix &operator*=(Scalar s) noexcept {
    // non-narrowing in-place scalar multiplication
    for (std::size_t i = 0; i < N * M; ++i) {
      arr[i] *= s;
    }
    return *this;
  }

  constexpr Matrix &transpose_inplace() noexcept
    requires(N == M)
  {
    // in-place transpose for square matrices
    // fluent interface support (returning the matrix)
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = i + 1; j < M; ++j) {
        std::swap((*this)(j, i), (*this)(i, j));
      }
    }
    return *this;
  }

  constexpr Matrix<T, M, N> transpose() const noexcept {
    Matrix<T, M, N> ans{};

    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < M; ++j) {
        ans(j, i) = (*this)(i, j);
      }
    }
    return ans;
  }

  constexpr Matrix<T, M, N> operator~() const noexcept { return transpose(); }

  constexpr auto det() const noexcept
    requires(N == M)
  {
    using R = lalib::floating_point_result_t<T, T>;

    Matrix<R, N, N> tmp = *this;
    R ans = R{1};

    for (std::size_t i = 0; i < N; ++i) {
      std::size_t pivot = i;
      for (std::size_t j = i + 1; j < N; ++j) {
        if (std::abs(tmp(j, i)) > std::abs(tmp(pivot, i))) {
          pivot = j;
        }
      }

      if (std::abs(tmp(pivot, i)) == R{0}) {
        return R{0};
      }

      if (i != pivot) {
        for (std::size_t c = 0; c < N; ++c) {
          std::swap(tmp(i, c), tmp(pivot, c));
        }
        ans = -ans;
      }

      ans *= tmp(i, i);

      for (std::size_t c = i + 1; c < N; ++c) {
        tmp(i, c) /= tmp(i, i);
      }

      for (std::size_t r = 0; r < N; ++r) {
        if (r != i && std::abs(tmp(r, i)) > R{0}) {
          const R factor = tmp(r, i);
          for (std::size_t c = i + 1; c < N; ++c) {
            tmp(r, c) -= tmp(i, c) * factor;
          }
        }
      }
    }
    return ans;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr bool operator==(const Matrix<U, P, Q> &other) const noexcept {
    if constexpr (P != N || Q != M) {
      return false;
    }

    for (std::size_t i = 0; i < N * M; ++i) {
      if (arr[i] != other.arr[i]) {
        return false;
      }
    }
    return true;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr bool operator!=(const Matrix<U, P, Q> &other) const noexcept {
    return !(*this == other);
  }

  // TODO: matrix-vector mutliplication
};
} // namespace lalib
