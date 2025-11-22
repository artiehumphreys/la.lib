#pragma once

#include "lalib/scalar/type_traits.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <mdspan>
#include <span>
#include <type_traits>

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

  template <class U, std::size_t P, std::size_t Q,
            class R = std::common_type_t<T, U>>
  constexpr auto operator+(const Matrix<U, P, Q> &other) const noexcept {
    static_assert(N == P && M == Q, "matrix dimensions must match");

    Matrix<R, N, M> ans{};
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < M; ++j) {
        ans(i, j) = static_cast<R>((*this)(i, j) + other(i, j));
      }
    }
    return ans;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr Matrix &operator+=(const Matrix<U, P, Q> &other) noexcept {
    // in-place operations don't return a copy of the modified object (hence
    // Matrix& vs. auto)
    static_assert(N == P && M == Q, "matrix dimensions must be the same");
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < M; ++j) {
        (*this)(i, j) += static_cast<T>(other(i, j));
      }
    }
    return *this;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr auto operator-(const Matrix<U, P, Q> &other) const noexcept {
    static_assert(N == P && M == Q, "matrix dimensions must match");

    // force signed result in case of unsigned matrix types
    using R = lalib::signed_result_t<T, U>;

    Matrix<R, N, M> neg{};
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < M; ++j) {
        neg(i, j) = -static_cast<R>(other(i, j));
      }
    }
    return (*this) + neg;
  }

  template <class U, std::size_t P, std::size_t Q>
  constexpr Matrix &operator-=(const Matrix<U, P, Q> &other) noexcept {
    // do not modify matrix, even if operation leads to underflow
    static_assert(N == P && M == Q, "matrix dimensions must match");

    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < M; ++j) {
        (*this)(i, j) -= static_cast<T>(other(i, j));
      }
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

    for (std::size_t i = 0; i < N; i += TILE_SIZE) {
      for (std::size_t j = 0; j < Q; j += TILE_SIZE) {
        for (std::size_t k = 0; k < P; k += TILE_SIZE) {

          // min is constexpr
          const std::size_t i_end = std::min(i + TILE_SIZE, N);
          const std::size_t j_end = std::min(j + TILE_SIZE, Q);
          const std::size_t k_end = std::min(k + TILE_SIZE, P);

          for (std::size_t i0 = i; i0 < i_end; ++i0) {
            for (std::size_t j0 = j; j0 < j_end; ++j0) {
              R sum{};
              for (std::size_t k0 = k; k0 < k_end; ++k0) {
                // TODO: dereferencing for this many ops can accumulate overhead
                sum += static_cast<R>((*this)(i0, k0)) *
                       static_cast<R>(other(k0, j0));
              }
              ans(i0, j0) = sum;
            }
          }
        }
      }
    }
    return ans;
  }

  constexpr auto transpose_inplace() noexcept
    requires(N == M)
  {
    // in-place transpose for square matrices
    // fluent interface support (returning the matrix)
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < M; ++j) {
        std::swap((*this)(j, i), (*this)(i, j));
      }
    }
    return *this;
  }

  constexpr Matrix &transpose() const noexcept {
    Matrix<T, M, N> ans{};

    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = i + 1; j < M; ++j) {
        ans(j, i) = (*this)(i, j);
      }
    }
    return ans;
  }

  constexpr auto operator~() const noexcept { return transpose(); }

  // TODO: matrix-vector mutliplication
};
} // namespace lalib
