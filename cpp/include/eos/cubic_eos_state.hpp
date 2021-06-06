#ifndef EOS_CUBIC_EOS_STATE_HPP
#define EOS_CUBIC_EOS_STATE_HPP

#include <eos/cubic_eos_traits.hpp>
#include <vector>

namespace eos {

template <typename Eos>
struct CubicEosState {
 public:
  using Scalar = typename CubicEosTraits<Eos>::Scalar;

  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulstion parameter
  CubicEosState(const Scalar& a, const Scalar& b) : a_{a}, b_{b} {}

  /// @brief Computes Z-factors
  /// @tparam CubicEquationSolver Cubic equation solver
  /// @param[in] solver Cubic equation solver
  ///
  /// A cubic equation solver takes the coefficients of the cubic equation as
  /// std::array<Scalar, 3> and returns real roots as std::vector<Scalar>.
  template <typename CubicEquationSolver>
  std::vector<Scalar> zfactors(CubicEquationSolver&& solver) const {
    return solver(Eos::zfactor_cubic_eq(a_, b_));
  }

  /// @param[in] z Z-factor
  Scalar ln_fugacity_coeff(const Scalar& z) const noexcept {
    return Eos::ln_fugacity_coeff(z, a_, b_);
  }

  /// @param[in] z Z-factor
  Scalar fugacity_coeff(const Scalar& z) const noexcept {
    return Eos::fugacity_coeff(z, a_, b_);
  }

 private:
  Scalar a_;  /// Reduced attraction parameter
  Scalar b_;  /// Reduced repulsion parameter
};

}  // namespace eos

#endif  // EOS_CUBIC_EOS_STATE_HPP
