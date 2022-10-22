#ifndef EOS_SOAVE_REDLICH_KWONG_EOS_HPP
#define EOS_SOAVE_REDLICH_KWONG_EOS_HPP

#include <cmath>
#include <eos/cubic_eos_base.hpp>

namespace eos {

template <typename Scalar>
class SoaveRedlichKwongEos;

template <typename Scalar_>
struct CubicEosTraits<SoaveRedlichKwongEos<Scalar_>> {
  using Scalar = Scalar_;
  static constexpr Scalar omega_a = 0.42748;
  static constexpr Scalar omega_b = 0.08664;
};

/// @brief Soave-Redlich-Kwong Equations of State
template <typename Scalar>
class SoaveRedlichKwongEos : public CubicEosBase<SoaveRedlichKwongEos<Scalar>> {
 public:
  using Base = CubicEosBase<SoaveRedlichKwongEos<Scalar>>;

  // Static Functions

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static Scalar pressure_impl(const Scalar& t, const Scalar& v, const Scalar& a,
                              const Scalar& b) noexcept {
    return gas_constant<Scalar>() * t / (v - b) - a / (v * (v + b));
  }

  /// @brief Computes coefficients of the cubic equation of Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor
  static std::array<Scalar, 3> zfactor_cubic_eq(const Scalar& a,
                                                const Scalar& b) noexcept {
    return {-1, a - b - b * b, -a * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static Scalar ln_fugacity_coeff(const Scalar& z, const Scalar& a,
                                  const Scalar& b) noexcept {
    return z - 1 - std::log(z - b) - a / b * std::log((z + b) / z);
  }

  /// @brief Computes a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Fugacity coefficient
  static Scalar fugacity_coeff(const Scalar& z, const Scalar& a,
                               const Scalar& b) noexcept {
    return std::exp(ln_fugacity_coeff(z, a, b));
  }

  /// @brief Computes residual Helmholtz free energy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Residual Helmholtz free energy
  static Scalar residual_helmholtz_free_energy(const Scalar& z, const Scalar& a,
                                               const Scalar& b) noexcept {
    return std::log(z / (z - b)) + a / z * std::log(z / (z + b));
  }

  /// @brief Computes residual Gibbs free energy
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Residual Gibbs free energy
  static Scalar residual_gibbs_free_energy(const Scalar& z, const Scalar& a,
                                           const Scalar& b) noexcept {
    return z - 1 + residual_helmoltz_free_energy(z, a, b);
  }

  // Constructors

  SoaveRedlichKwongEos() = default;

  SoaveRedlichKwongEos(const Scalar& pc, const Scalar& tc,
                       const Scalar& omega) noexcept
      : Base{pc, tc}, omega_{omega}, m_{calc_m(omega)} {}

  SoaveRedlichKwongEos(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos(SoaveRedlichKwongEos&&) = default;

  SoaveRedlichKwongEos& operator=(const SoaveRedlichKwongEos&) = default;
  SoaveRedlichKwongEos& operator=(SoaveRedlichKwongEos&&) = default;

  // Member functions

  void set_params(const Scalar& pc, const Scalar& tc,
                  const Scalar& omega) noexcept {
    this->Base::set_params(pc, tc);
    omega_ = omega;
    m_ = calc_m(omega);
  }

  /// @param[in] tr Reduced temperature
  constexpr Scalar correction_factor(const Scalar& tr) const noexcept {
    const auto a = 1 + m_ * (1 - std::sqrt(tr));
    return a * a;
  }

 private:
  /// @brief Computes parameter \f$ m \f$ from acentric factor
  /// @param[in] omega Acentric factor
  static Scalar calc_m(const Scalar& omega) noexcept {
    return 0.48 + (1.574 - 0.176 * omega) * omega;
  }

  Scalar omega_;  /// Acentric factor
  Scalar m_;      /// Factor used to compute the temperature correction factor
};

/// @brief Makes Soave-Redlich-Kwong EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
template <typename Scalar>
inline SoaveRedlichKwongEos<Scalar> make_soave_redlich_kwong_eos(
    const Scalar& pc, const Scalar& tc, const Scalar& omega) {
  return {pc, tc, omega};
}

}  // namespace eos

#endif  // EOS_SOAVE_REDLICH_KWONG_EOS_HPP