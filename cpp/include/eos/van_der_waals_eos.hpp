#ifndef EOS_VAN_DER_WAALS_EOS_HPP
#define EOS_VAN_DER_WAALS_EOS_HPP

#include <array>
#include <cmath>
#include <eos/cubic_eos_base.hpp>

namespace eos
{

  template <typename Scalar>
  class VanDerWaalsEos;

  template <typename Scalar_>
  struct CubicEosTraits<VanDerWaalsEos<Scalar_>>
  {
    using Scalar = Scalar_;
    static constexpr Scalar omega_a = 0.421875;
    static constexpr Scalar omega_b = 0.125;
  };

  /// @brief Van der Waals Equations of State
  template <typename Scalar>
  class VanDerWaalsEos : public CubicEosBase<VanDerWaalsEos<Scalar>>
  {
  public:
    using Base = CubicEosBase<VanDerWaalsEos<Scalar>>;

    // Static Functions

    /**
     * @brief Computes pressure at given temperature and volume.
     * @param[in] t Temperature
     * @param[in] v Molar volume
     * @param[in] a Attraction parameter
     * @param[in] b Repulsion parameter
     * @returns Scalar Pressure
     *
     * @f[
     * P = \frac{RT}{v - b} - \frac{a}{v^2},
     * @f]
     * where @f$ R @f$ is gas constant, @f$ T @f$ is temperature, @f$ v @f$ is
     * molar volume, @f$ a @f$ is attraction parameter, @f$ b @f$ is repulsion
     * parameter.
     */
    static Scalar pressure_impl(const Scalar &t, const Scalar &v, const Scalar &a,
                                const Scalar &b) noexcept
    {
      return gas_constant<Scalar>() * t / (v - b) - a / (v * v);
    }

    /**
     * @brief Computes coefficients of the cubic equation of Z-factor
     * @param[in] a Reduced attraction parameter
     * @param[in] b Reduced repulsion parameter
     * @returns std::array<Scalar, 3> Coefficients of the cubic equation of z-factor
     *
     * @f[
     * Z^3 - (B+1)Z^2 + AZ -AB = 0,
     * @f]
     * where @f$ A @f$ is reduced attraction parameter, and @f$ B @f$ is reduced
     * repulsion parameter.
     * @f[
     * A = \frac{a}{(RT)^2}, \quad B = \frac{b}{RT}.
     * @f]
     */
    static std::array<Scalar, 3> zfactor_cubic_eq(const Scalar &a,
                                                  const Scalar &b) noexcept
    {
      return {-b - 1, a, -a * b};
    }

    /**
     * @brief Computes the natural logarithm of a fugacity coefficient
     * @param[in] z Z-factor
     * @param[in] a Reduced attraction parameter
     * @param[in] b Reduced repulsion parameter
     * @returns Scalar The natural log of a fugacity coefficient
     *
     * Fugacity coefficient @f$ \phi @f$ for vdW EoS can be expressed by
     * @f[
     * \ln \phi = -\frac{A}{Z} - \ln(Z-B) + Z-1,
     * @f]
     * where @f$ Z @f$ is Z-factor, @f$ A @f$ is reduced attraction parameter,
     * and @f$ B @f$ is reduced repulsion parameter.
     */
    static Scalar ln_fugacity_coeff(const Scalar &z, const Scalar &a,
                                    const Scalar &b) noexcept
    {
      return -std::log(z - b) - a / z + z - 1;
    }

    /**
     * @brief Computes a fugacity coefficient
     * @param[in] z Z-factor
     * @param[in] a Reduced attraction parameter
     * @param[in] b Reduced repulsion parameter
     * @returns Scalar Fugacity coefficient
     */
    static Scalar fugacity_coeff(const Scalar &z, const Scalar &a,
                                 const Scalar &b) noexcept
    {
      return std::exp(ln_fugacity_coeff(z, a, b));
    }

    /**
     * @brief Computes residual Helmholtz free energy
     * @param[in] z Z-factor
     * @param[in] a Reduced attraction parameter
     * @param[in] b Reduced repulsion parameter
     * @returns Scalar Residual Helmholtz free energy
     *
     * Residual Helmholtz free energy @f$ F^r @f$ for vdW EoS can be expressed by
     * @f[
     * \frac{F^r}{NRT} = \ln(\frac{Z}{Z-B}) - \frac{A}{Z}.
     * @f]
     */
    static Scalar residual_helmholtz_free_energy(const Scalar &z, const Scalar &a,
                                                 const Scalar &b) noexcept
    {
      return std::log(z / (z - b)) - a / z;
    }

    /**
     * @brief Computes residual Gibbs free energy
     * @param[in] z Z-factor
     * @param[in] a Reduced attraction parameter
     * @param[in] b Reduced repulsion parameter
     * @returns Scalar Residual Gibbs free energy
     *
     * Residual Gibbs free energy @f$ G^r @f$ for vdW EoS can be expressed by
     * @f[
     * \frac{G^r}{NRT} = (Z - 1) + \frac{F^r}{NRT}.
     * @f]
     */
    static Scalar residual_gibbs_free_energy(const Scalar &z, const Scalar &a,
                                             const Scalar &b) noexcept
    {
      return z - 1 + residual_helmoltz_free_energy(z, a, b);
    }

    /**
     * @brief Computes residual internal energy.
     *
     * @param z Z-factor
     * @param a Attraction parameter
     * @param b Repulsion parameter
     * @return Scalar Residual internal energy
     *
     * @f[
     * \frac{U^r}{NRT} = -\frac{A}{Z}.
     * @f]
     */
    static Scalar residual_internal_energy(const Scalar &z, const Scalar &a, const Scalar &b) noexcept
    {
      return -a / z;
    }

    // Constructors

    VanDerWaalsEos() = default;

    /**
     * @brief Construct a new Van Der Waals Eos object
     *
     * @param pc Critical pressure
     * @param tc Critical temperature
     */
    VanDerWaalsEos(const Scalar &pc, const Scalar &tc) noexcept : Base{pc, tc} {}

    VanDerWaalsEos(const VanDerWaalsEos &) = default;
    VanDerWaalsEos(VanDerWaalsEos &&) = default;

    VanDerWaalsEos &operator=(const VanDerWaalsEos &) = default;
    VanDerWaalsEos &operator=(VanDerWaalsEos &&) = default;

    // Member functions

    /**
     * @brief Set parameters
     *
     * @param pc Critical pressure
     * @param tc Critical temperature
     */
    void set_params(const Scalar &pc, const Scalar &tc) noexcept
    {
      this->Base::set_params(pc, tc);
    }

    /**
     * @brief Computes correction factor for attraction parameter.
     *
     * @param[in] tr Reduced temperature (unused)
     * @return Scalar Correction factor
     */
    constexpr Scalar correction_factor(const Scalar &) const noexcept
    {
      return static_cast<Scalar>(1);
    }
  };

  /**
   * @brief Makes van der Waals EoS
   * @param[in] pc Critical pressure
   * @param[in] tc Critical temperature
   * @return VanDerWaalsEos<Scalar> van der Waals EoS
   */
  template <typename Scalar>
  inline VanDerWaalsEos<Scalar> make_van_der_waals_eos(const Scalar &pc,
                                                       const Scalar &tc)
  {
    return {pc, tc};
  }

} // namespace eos

#endif // EOS_VAN_DER_WAALS_EOS_HPP