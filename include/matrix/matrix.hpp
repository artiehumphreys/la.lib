#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <mdspan>
#include <span>

namespace lalib {

template <class T, std::size_t N, std::size_t M> struct Matrix {
  static_assert(N > 0 && M > 0, "Matrix dimensions must be positive");

  std::array<T, N * M> arr{};

  constexpr Matrix() = default;

  template <class Iter> constexpr Matrix(Iter begin, Iter end) {
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

  template <class U, std::size_t P, std::size_t Q>
  constexpr auto operator*(Matrix<U, P, Q> other) const noexcept {
    static_assert(M == P, "matrix dimensions must match");
    static constexpr std::size_t TILE_SIZE =
        64 / std::max(sizeof(T), sizeof(U));

    using R = decltype(std::declval<T>() * std::declval<U>());

    Matrix<R, N, Q> ans{};

    for (int i = 0; i < N; i += TILE_SIZE) {
      for (int j = 0; j < Q; j += TILE_SIZE) {
        for (int k = 0; k < P; k += TILE_SIZE) {

          // min is constexpr
          const std::size_t i_end = std::min(i + TILE_SIZE, N);
          const std::size_t j_end = std::min(j + TILE_SIZE, Q);
          const std::size_t k_end = std::min(k + TILE_SIZE, P);

          for (int i0 = i; i0 < i_end; ++i0) {
            for (int j0 = j; j0 < j_end; ++j0) {
              R sum = ans(i0, j0);
              for (int k0 = k; k0 < k_end; ++k0) {
                sum += (*this)(i0, k0) * other(k0, j0);
              }
              ans(i0, j0) = sum;
            }
          }
        }
      }
    }
    return ans;
  }
};
} // namespace lalib
