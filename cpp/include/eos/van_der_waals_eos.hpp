#ifndef EOS_VAN_DER_WAALS_EOS_HPP
#define EOS_VAN_DER_WAALS_EOS_HPP

#include <array>
#include <cmath>
#include <eos/cubic_eos_base.hpp>

namespace eos {

template <typename Scalar>
class VanDerWaalsEos;

template <typename Scalar_>
struct CubicEosTraits<VanDerWaalsEos<Scalar_>> {
  using Scalar = Scalar_;
  static constexpr Scalar omega_a = 0.421875;
  static constexpr Scalar omega_b = 0.125;
};

/// @brief Van der Waals Equations of State
template <typename Scalar>
class VanDerWaalsEos : public CubicEosBase<VanDerWaalsEos<Scalar>> {
 public:
  using Base = CubicEosBase<VanDerWaalsEos<Scalar>>;

  // Static Functions

  /// @brief Computes pressure at given temperature and volume.
  /// @param[in] t Temperature
  /// @param[in] v Volume
  /// @param[in] a Attraction parameter
  /// @param[in] b Repulsion parameter
  /// @returns Pressure
  static Scalar pressure_impl(const Scalar& t, const Scalar& v, const Scalar& a,
                              const Scalar& b) noexcept {
    return gas_constant<Scalar>() * t / (v - b) - a / (v * v);
  }

  /// @brief Computes coefficients of the cubic equation of Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns Coefficients of the cubic equation of z-factor
  static std::array<Scalar, 3> zfactor_cubic_eq(const Scalar& a,
                                                const Scalar& b) noexcept {
    return {-b - 1, a, -a * b};
  }

  /// @brief Computes the natural logarithm of a fugacity coefficient
  /// @param[in] z Z-factor
  /// @param[in] a Reduced attraction parameter
  /// @param[in] b Reduced repulsion parameter
  /// @returns The natural logarithm of a fugacity coefficient
  static Scalar ln_fugacity_coeff(const Scalar& z, const Scalar& a,
                                  const Scalar& b) noexcept {
    return -std::log(z - b) - a / z + z - 1;
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

  VanDerWaalsEos() = default;

  VanDerWaalsEos(const Scalar& pc, const Scalar& tc) noexcept : Base{pc, tc} {}

  VanDerWaalsEos(const VanDerWaalsEos&) = default;
  VanDerWaalsEos(VanDerWaalsEos&&) = default;

  VanDerWaalsEos& operator=(const VanDerWaalsEos&) = default;
  VanDerWaalsEos& operator=(VanDerWaalsEos&&) = default;

  // Member functions

  void set_params(const Scalar& pc, const Scalar& tc) noexcept {
    this->Base::set_params(pc, tc);
  }

  /// @param[in] tr Reduced temperature (unused)
  constexpr Scalar correction_factor(const Scalar&) const noexcept {
    return static_cast<Scalar>(1);
  }
};

/// @brief Makes van der Waals EoS
/// @param[in] pc Critical pressure
/// @param[in] tc Critical temperature
template <typename Scalar>
inline VanDerWaalsEos<Scalar> make_van_der_waals_eos(const Scalar& pc,
                                                     const Scalar& tc) {
  return {pc, tc};
}

}  // namespace eos

#endif  // EOS_VAN_DER_WAALS_EOS_HPP