#ifndef EOS_ISOTHERMAL_LINE_HPP
#define EOS_ISOTHERMAL_LINE_HPP

#include <eos/cubic_eos_traits.hpp>

namespace eos {

template <typename Eos>
class IsothermalLine {
 public:
  using Scalar = typename CubicEosTraits<Eos>::Scalar;

  /// @param[in] t Temperature
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  IsothermalLine(const Scalar& t, const Scalar& a, const Scalar& b)
      : t_{t}, a_{a}, b_{b} {}

  /// @param[in] v Volume
  Scalar pressure(const Scalar& v) const noexcept {
    return Eos::pressure_impl(t_, v, a_, b_);
  }

 private:
  Scalar t_;  /// Temperature
  Scalar a_;  /// Attraction parameter
  Scalar b_;  /// Repulsion parameter
};

}  // namespace eos

#endif  // EOS_ISOTHERMAL_LINE_HPP
