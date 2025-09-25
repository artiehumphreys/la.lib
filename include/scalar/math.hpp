#pragma once

#include <cassert>
#include <cmath>
#include <concepts>

namespace lalib {

inline constexpr double EPS = 1e-8;

template <class T>
  requires std::floating_point<T>
constexpr T sqrt(T x) {
  assert(x > 0);

  if (x == 0.0 || x == 1.0) {
    return x;
  }

  // babylonian method: https://www.cs.utep.edu/vladik/2009/olg09-05a.pdf
  double prev = 0, est = (x + 1) * 0.5;
  while (std::fabs(est - prev) > EPS) {
    // TODO: scale epsilon based on input so the calculation doesn't end too
    // early / late
    prev = est;
    est = (est + x / est) * 0.5;
  }

  return est;
}

template <class U>
  requires std::unsigned_integral<U>
constexpr U sqrt(U x) {
  // placeholder
  return U(0);
}

} // namespace lalib
