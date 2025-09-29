#include "scalar/math.hpp"
#include <cstdint>
#include <iostream>

int main() {
  // ----- constexpr checks -----
  static_assert(lalib::sqrt(4.0) == 2.0);
  static_assert(lalib::sqrt(9u) == 3u);
  static_assert(lalib::exponentiate(3, 5u) == 243);
  constexpr auto e = lalib::exponentiate(2.0, 5LL);
  static_assert(e == 32.0);

  // ----- simple runtime checks -----
  if (lalib::sqrt(0.0) != 0.0)
    return 1;
  if (lalib::sqrt(1.0) != 1.0)
    return 1;

  auto ef = lalib::exponentiate(2.0, 10LL);
  if (std::abs(ef - 1024.0) > 1e-12)
    return 1;

  return 0;
}
