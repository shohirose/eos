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

  /// @param[in] z Z-factors
  std::vector<Scalar> ln_fugacity_coeff(
      const std::vector<Scalar>& z) const noexcept {
    std::vector<Scalar> ln_phi(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) {
      ln_phi[i] = Eos::ln_fugacity_coeff(z[i], a_, b_);
    }
    return ln_phi;
  }

  /// @param[in] z Z-factor
  Scalar fugacity_coeff(const Scalar& z) const noexcept {
    return Eos::fugacity_coeff(z, a_, b_);
  }

  /// @param[in] z Z-factors
  std::vector<Scalar> fugacity_coeff(
      const std::vector<Scalar>& z) const noexcept {
    std::vector<Scalar> phi(z.size());
    for (std::size_t i = 0; i < z.size(); ++i) {
      phi[i] = Eos::fugacity_coeff(z[i], a_, b_);
    }
    return phi;
  }

 private:
  Scalar a_;  /// Reduced attraction parameter
  Scalar b_;  /// Reduced repulsion parameter
};

}  // namespace eos

#endif  // EOS_CUBIC_EOS_STATE_HPP
