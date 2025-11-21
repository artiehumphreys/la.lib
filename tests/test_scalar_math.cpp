#include "lalib/scalar/math.hpp"

int main() {
  static_assert(lalib::sqrt(4.0) == 2.0);
  static_assert(lalib::sqrt(9u) == 3u);
  static_assert(lalib::exponentiate(3, 5u) == 243);
  static_assert(lalib::exponentiate(2.0, 5u) == 32.0);

  static_assert(lalib::sqrt(0.0) == 0.0);
  static_assert(lalib::sqrt(0u) == 0u);
  static_assert(lalib::sqrt(1.0) == 1.0);

  static_assert(lalib::exponentiate(2.0, 10LL));

  return 0;
}
