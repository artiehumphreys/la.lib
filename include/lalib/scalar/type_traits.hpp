

#include <type_traits>

namespace lalib {

template <class A, class B> struct signed_result {
  using common = std::common_type_t<A, B>;

  using type = std::conditional_t<std::is_integral_v<common>,
                                  std::make_signed_t<common>, common>;
};

template <class A, class B>
using signed_result_t = typename signed_result<A, B>::type;

} // namespace lalib
