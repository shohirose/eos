#ifndef EOS_PENG_ROBINSON_EOS_HPP
#define EOS_PENG_ROBINSON_EOS_HPP

#include <cmath>
#include <eos/cubic_eos_base.hpp>
#include <eos/mathematical_constants.hpp>

namespace eos {

template <typename Scalar>
class PengRobinsonEos;

template <typename Scalar_>
struct CubicEosTraits<PengRobinsonEos<Scalar_>> {
  using Scalar = Scalar_;
  static constexpr Scalar omega_a = 0.45724;
  static constexpr Scalar omega_b = 0.07780;
};

/// @brief Peng-Robinson Equations of State
template <typename Scalar>
class PengRobinsonEos : public CubicEosBase<PengRobinsonEos<Scalar>> {
 public:
  using Base = CubicEosBase<PengRobinsonEos<Scalar>>;

  // Static Functions

  /**
   * @brief Computes pressure at given temperature and volume.
   *
   * @param t Temperature
   * @param v Volume
   * @param a Attraction parameter
   * @param b Repulsion parameter
   * @return Scalar Pressure
   *
   * Pressure for Peng-Robinso EoS can be expressed by
   * @f[
   * P = \frac{RT}{V-b} - \frac{a}{V^2+2bV-b^2},
   * @f]
   * where @f$ P @f$ is pressure, @f$ T @f$ is temperature, @f$ V @f$ is volume,
   * @f$ R @f$ is gas constant, @f$ a @f$ is attraction parameter, and @f$ b @f$
   * is repulsion parameter.
   */
  static Scalar pressure_impl(const Scalar& t, const Scalar& v, const Scalar& a,
                              const Scalar& b) noexcept {
    constexpr auto R = gas_constant<Scalar>();
    return R * t / (v - b) - a / (v * (v + b) + b * (v - b));
  }

  /**
   * @brief Computes coefficients of the cubic equation of Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @returns Coefficients of the cubic equation of z-factor
   *
   * A cubic equation of Z-factor for Peng-Robinson EoS can be expressed by
   * @f[
   * Z^3 + (1-B)Z^2 + (3B+2)BZ + (B^2+B-A)B = 0,
   * @f]
   * where @f$ Z @f$ is Z-factor, @f$ A @f$ is reduced attraction parameter,
   * and @f$ B @f$ is reduced repulsion parameter.
   * @f[
   * A = \frac{aP}{(RT)^2}, \quad B = \frac{bP}{RT}
   * @f]
   */
  static std::array<Scalar, 3> zfactor_cubic_eq(const Scalar& a,
                                                const Scalar& b) noexcept {
    return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
  }

  /**
   * @brief Computes the natural log of a fugacity coefficient
   * @param[in] z Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @returns The natural log of a fugacity coefficient
   *
   * Fugacity coefficient @f$ \phi @f$ for Peng-Robinson EoS can be expressed by
   * @f{gather*}{
   * \ln \phi = -\ln(Z-B) + (Z-1) + I, \\
   * I = \frac{A}{2\sqrt{2}B} \ln \frac{Z-(\sqrt{2}-1)B}{Z+(\sqrt{2}+1)B}.
   * @f}
   */
  static Scalar ln_fugacity_coeff(const Scalar& z, const Scalar& a,
                                  const Scalar& b) noexcept {
    return z - 1 - std::log(z - b) - calc_q(z, a, b);
  }

  /**
   * @brief Computes a fugacity coefficient
   * @param[in] z Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @returns Fugacity coefficient
   */
  static Scalar fugacity_coeff(const Scalar& z, const Scalar& a,
                               const Scalar& b) noexcept {
    return std::exp(ln_fugacity_coeff(z, a, b));
  }

  /**
   * @brief Computes residual Helmholtz free energy
   * @param[in] z Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @returns Residual Helmholtz free energy
   *
   * Residual Helmholtz free energy @f$ F^r @f$ for Peng-Robinson EoS can be
   * expressed by
   * @f[
   * \frac{F^r}{RT} = \ln \frac{Z}{Z-B} + I.
   * @f]
   */
  static Scalar residual_helmholtz_free_energy(const Scalar& z, const Scalar& a,
                                               const Scalar& b) noexcept {
    return std::log(z / (z - b)) + calc_q(z, a, b);
  }

  /**
   * @brief Computes residual Gibbs free energy
   * @param[in] z Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   * @returns Residual Gibbs free energy
   *
   * Residual Gibbs free energy @f$ G^r @f$ for Peng-Robinson EoS can be
   * expressed by
   * @f[
   * \frac{G^r}{RT} = (Z-1) + \ln \frac{Z}{Z-B} + I.
   * @f]
   */
  static Scalar residual_gibbs_free_energy(const Scalar& z, const Scalar& a,
                                           const Scalar& b) noexcept {
    return z - 1 + residual_helmoltz_free_energy(z, a, b);
  }

  /**
   * @brief Computes residual internal energy.
   *
   * @param z Z-factor
   * @param a Attraction parameter
   * @param b Repulsion parameter
   * @param dlna_dlnt Logarithmic derivative of attraction correction factor in
   * terms of temperature
   * @return Scalar Residual internal energy
   *
   * @f[
   * \frac{U^r}{NRT} = \left( 1 - \frac{d \ln \alpha}{d \ln T}) I.
   * @f]
   */
  static Scalar residual_internal_energy(const Scalar& z, const Scalar& a,
                                         const Scalar& b,
                                         const Scalar& dlna_dlnt) noexcept {
    return (1 - dlna_dlnt) * calc_q(z, a, b);
  }

  // Constructors

  PengRobinsonEos() = default;

  /**
   * @brief Construct a new Peng Robinson Eos object
   *
   * @param pc Critical pressure
   * @param tc Critical temperature
   * @param omega Acentric factor
   */
  PengRobinsonEos(const Scalar& pc, const Scalar& tc,
                  const Scalar& omega) noexcept
      : Base{pc, tc}, omega_{omega}, m_{calc_m(omega)} {}

  PengRobinsonEos(const PengRobinsonEos&) = default;
  PengRobinsonEos(PengRobinsonEos&&) = default;

  PengRobinsonEos& operator=(const PengRobinsonEos&) = default;
  PengRobinsonEos& operator=(PengRobinsonEos&&) = default;

  // Member functions

  /**
   * @brief Set parameters
   *
   * @param pc Critical pressure
   * @param tc Critical temperature
   * @param omega Acentric factor
   */
  void set_params(const Scalar& pc, const Scalar& tc,
                  const Scalar& omega) noexcept {
    this->Base::set_params(pc, tc);
    omega_ = omega;
    m_ = calc_m(omega);
  }

  /**
   * @brief Compute correction factor
   * @param[in] tr Reduced temperature
   */
  constexpr Scalar correction_factor(const Scalar& tr) const noexcept {
    const auto a = 1 + m_ * (1 - std::sqrt(tr));
    return a * a;
  }

 private:
  /// @param[in] omega Acentric factor
  static Scalar calc_m(const Scalar& omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  /**
   * @param[in] z Z-factor
   * @param[in] a Reduced attraction parameter
   * @param[in] b Reduced repulsion parameter
   */
  static Scalar calc_q(const Scalar& z, const Scalar& a,
                       const Scalar& b) noexcept {
    constexpr auto sqrt2 = sqrt_two<Scalar>();
    constexpr auto delta1 = 1 + sqrt2;
    constexpr auto delta2 = 1 - sqrt2;
    return a / (2 * sqrt2 * b) * std::log((z + delta1 * b) / (z + delta2 * b));
  }

  Scalar omega_;  /// Acentric factor
  Scalar m_;      /// Factor used to compute the temperature correction factor
};

/// @brief Makes Peng-Robinson EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
/// @param[in] omega Acentric factor
template <typename Scalar>
inline PengRobinsonEos<Scalar> make_peng_robinson_eos(const Scalar& pc,
                                                      const Scalar& tc,
                                                      const Scalar& omega) {
  return {pc, tc, omega};
}

}  // namespace eos

#endif  // EOS_PENG_ROBINSON_EOS_HPP