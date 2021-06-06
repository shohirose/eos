#ifndef EOS_PENG_ROBINSON_EOS_HPP
#define EOS_PENG_ROBINSON_EOS_HPP

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

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static Scalar pressure_impl(const Scalar& t, const Scalar& v, const Scalar& a,
                              const Scalar& b) noexcept {
    constexpr auto R = gas_constant<Scalar>();
    return R * t / (v - b) - a / (v * (v + b) + b * (v - b));
  }

  /// @brief Computes coefficients of the cubic equation of Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor
  static std::array<Scalar, 3> zfactor_cubic_eq(const Scalar& a,
                                                const Scalar& b) noexcept {
    return {b - 1, a - (3 * b + 2) * b, (-a + b + b * b) * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static Scalar ln_fugacity_coeff(const Scalar& z, const Scalar& a,
                                  const Scalar& b) noexcept {
    return z - 1 - std::log(z - b) - calc_q(z, a, b);
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

  // Constructors

  PengRobinsonEos() = default;

  PengRobinsonEos(const Scalar& pc, const Scalar& tc,
                  const Scalar& omega) noexcept
      : Base{pc, tc}, omega_{omega}, m_{calc_m(omega)} {}

  PengRobinsonEos(const PengRobinsonEos&) = default;
  PengRobinsonEos(PengRobinsonEos&&) = default;

  PengRobinsonEos& operator=(const PengRobinsonEos&) = default;
  PengRobinsonEos& operator=(PengRobinsonEos&&) = default;

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
  /// @param[in] omega Acentric factor
  static Scalar calc_m(const Scalar& omega) noexcept {
    return 0.3796 + omega * (1.485 - omega * (0.1644 - 0.01667 * omega));
  }

  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
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