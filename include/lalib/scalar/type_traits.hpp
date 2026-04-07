#include <type_traits>

namespace lalib {

struct scalar {};

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

// unevaluted instance of type From. Doesn't require a constructor.
// check if type To can be constructed from type From.
// braces reject narrowing.
template <class To, class From>
concept non_narrowing = requires { To{std::declval<From>()}; };

// ensure that the product type can be stored into ElemT w/o narrowing
template <class Scalar, class ElemT>
concept safe_scalar_multiply =
    std::is_arithmetic_v<Scalar> &&
    non_narrowing<decltype(std::declval<ElemT>() * std::declval<Scalar>()),
                  ElemT>;

template <class A, class B>
using signed_result_t = typename signed_result<A, B>::type;

template <class A, class B>
using floating_point_result_t = typename floating_point_result<A, B>::type;

} // namespace lalib
