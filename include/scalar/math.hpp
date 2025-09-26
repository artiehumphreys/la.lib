#pragma once

#include <bit>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstdint>
#include <limits>

namespace lalib {

inline constexpr double EPS = 1e-8;

template <class T>
  requires std::floating_point<T>
constexpr T sqrt(T x) {
  const T zero{0}, one{1}, half{0.5};
  assert(x >= zero && "x must be nonnegative");

  if (x == zero || x == one) {
    return x;
  }

  auto absT = [](T t) { return t < zero ? -t : t; };

  // babylonian method: https://www.cs.utep.edu/vladik/2009/olg09-05a.pdf
  T est = (x + one) * half;
  constexpr int maxIters =
      static_cast<int>(std::numeric_limits<T>::digits) / 2 + 8;
  for (size_t i = 0; i < maxIters; ++i) {
    // TODO: scale epsilon based on input so the calculation doesn't end too
    // early / late
    const T next = (est + x / est) * half;
    T diff = absT(next - est);
    if (diff <= T{EPS} || est == next) {
      return next;
    }
    est = next;
  }

  return est;
}

template <class U>
  requires std::unsigned_integral<U>
constexpr U sqrt(U x) {
  const U one{1}, zero{0};

  if (x <= one) {
    return x;
  }

  U l = zero, r = one << ((std::bit_width(x - one) + 1) / 2);

  if (r > x) {
    r = x;
  }

  U ans = zero;
  // ensure progress when l + 1 == r
  while (l < r) {
    const U mid = l + ((r - l + one) >> one);
    // avoid division by zero
    if (mid != zero && mid <= x / mid) {
      l = mid;
    } else {
      r = mid - one;
    }
  }

  return l;
}

} // namespace lalib
