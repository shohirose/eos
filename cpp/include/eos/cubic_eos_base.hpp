#ifndef EOS_CUBIC_EOS_BASE_HPP
#define EOS_CUBIC_EOS_BASE_HPP

#include <eos/cubic_eos_state.hpp>
#include <eos/cubic_eos_traits.hpp>
#include <eos/thermodynamic_constants.hpp>

namespace eos {

template <typename Derived>
class CubicEosBase {
 public:
  using Scalar = typename CubicEosTraits<Derived>::Scalar;
  static constexpr auto omega_a = CubicEosTraits<Derived>::omega_a;
  static constexpr auto omega_b = CubicEosTraits<Derived>::omega_b;

  CubicEosBase() = default;
  CubicEosBase(const CubicEosBase&) = default;
  CubicEosBase(CubicEosBase&&) = default;

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  CubicEosBase(const Scalar& pc, const Scalar& tc) noexcept
      : pc_{pc},
        tc_{tc},
        a_{this->attraction_param(pc, tc)},
        b_{this->repulsion_param(pc, tc)} {}

  CubicEosBase& operator=(const CubicEosBase&) = default;
  CubicEosBase& operator=(CubicEosBase&&) = default;

  // Static functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static Scalar attraction_param(const Scalar& pc, const Scalar& tc) noexcept {
    constexpr auto R = gas_constant<Scalar>();
    return (omega_a * R * R) * tc * tc / pc;
  }

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  static Scalar repulsion_param(const Scalar& pc, const Scalar& tc) noexcept {
    constexpr auto R = gas_constant<Scalar>();
    return (omega_b * R) * tc / pc;
  }

  /// @brief Returns reduced attraction parameter at a given pressure and
  /// temperature without temperature correction.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static Scalar reduced_attraction_param(const Scalar& pr,
                                         const Scalar& tr) noexcept {
    return omega_a * pr / (tr * tr);
  }

  /// @brief Returns reduced repulsion parameter at a given pressure and
  /// temperature.
  /// @param[in] pr Reduced pressure
  /// @param[in] tr Reduced temperature
  static Scalar reduced_repulsion_param(const Scalar& pr,
                                        const Scalar& tr) noexcept {
    return omega_b * pr / tr;
  }

  // Member functions

  /// @param[in] pc Critical pressure
  /// @param[in] tc Critical temperature
  void set_params(const Scalar& pc, const Scalar& tc) noexcept {
    pc_ = pc;
    tc_ = tc;
    a_ = attraction_param(pc, tc);
    b_ = repulsion_param(pc, tc);
  }

  /// @brief Computes reduced pressure
  /// @param[in] p Pressure
  Scalar reduced_pressure(const Scalar& p) const noexcept { return p / pc_; }

  /// @brief Computes reduced temperature
  /// @param[in] t Temperature
  Scalar reduced_temperature(const Scalar& t) const noexcept { return t / tc_; }

  /// @brief Computes pressure at given temperature and volume
  /// @param[in] t Temperature
  /// @param[in] v Volume
  Scalar pressure(const Scalar& t, const Scalar& v) const noexcept {
    const auto tr = this->reduced_temperature(t);
    const auto a = this->derived().correction_factor(tr) * a_;
    return Derived::pressure_impl(t, v, a, b_);
  }

  /// @brief Creates isobaric-isothermal state
  /// @param[in] p Pressure
  /// @param[in] t Temperature
  CubicEosState<Derived> create_state(const Scalar& p,
                                      const Scalar& t) const noexcept {
    const auto pr = this->reduced_pressure(p);
    const auto tr = this->reduced_temperature(t);
    const auto ar = this->derived().correction_factor(tr) *
                    this->reduced_attraction_param(pr, tr);
    const auto br = this->reduced_repulsion_param(pr, tr);
    return {ar, br};
  }

 protected:
  /// @brief Get reference to derived class object
  Derived& derived() noexcept { return static_cast<Derived&>(*this); }

  /// @brief Get const reference to derived class object
  const Derived& derived() const noexcept {
    return static_cast<const Derived&>(*this);
  }

  Scalar pc_;  /// Critical pressure
  Scalar tc_;  /// Critical temperature
  Scalar a_;   /// Attraction parameter
  Scalar b_;   /// Repulsion parameter
};

}  // namespace eos

#endif  // EOS_CUBIC_EOS_BASE_HPP