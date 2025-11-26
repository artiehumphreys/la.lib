

#include <type_traits>

namespace lalib {

template <class A, class B> struct signed_result {
  using common = std::common_type_t<A, B>;

  using type = std::conditional_t<std::is_integral_v<common>,
                                  std::make_signed_t<common>, common>;
};

template <class A, class B> struct floating_point_result {
  using common = std::common_type_t<A, B>;

  using type =
      std::conditional_t<std::is_floating_point_v<common>, common, double>;
};

template <class A, class B>
using signed_result_t = typename signed_result<A, B>::type;

template <class A, class B>
using floating_point_result_t = typename floating_point_result<A, B>::type;

} // namespace lalib
