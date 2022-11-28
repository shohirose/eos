#ifndef EOS_SOAVE_CORRECTION_FACTOR_HPP
#define EOS_SOAVE_CORRECTION_FACTOR_HPP

#include <cmath>
#include <eos/correction_factor_base.hpp>

namespace eos {

/**
 * @brief Temperature correction factor for attraction parameter proposed by
 * Soave (1972).
 * 
 * @tparam Scalar 
 * 
 * Correction factor @f$ \alpha @f$ is given by
 * @f[
 * \alpha(T_r) = \left[ 1 + m \left(1 - \sqrt{T_r}\right) \right]^2,
 * @f]
 * where @f$ T_r @f$ is reduced temperature, and @f$ m @f$ is a constant.
 */
template <typename Scalar>
class SoaveCorrectionFactor
    : CorrectionFactorBase<SoaveCorrectionFactor<Scalar>> {
 public:
  SoaveCorrectionFactor(const Scalar& m) : m_{m} {}
  SoaveCorrectionFactor(const SoaveCorrectionFactor&) = default;
  SoaveCorrectionFactor(SoaveCorrectionFactor&&) = default;

  SoaveCorrectionFactor& operator=(const SoaveCorrectionFactor&) = default;
  SoaveCorrectionFactor& operator=(SoaveCorrectionFactor&&) = default;

  /**
   * @brief Computes correction factor
   *
   * @param tr Reduced temperature
   * @return Scalar Correction factor
   */
  Scalar value(const Scalar& tr) const noexcept override {
    using std::sqrt;
    const auto a = 1 + m_ * (1 - sqrt(tr));
    return a * a;
  }

  /**
   * @brief Computes logatithmic derivative of correction factor
   *
   * @param tr Reduced temperature
   * @return Scalar Logarithmic derivative of correction factor
   *
   * Logarithmic derivative of correction factor is given by
   * @f[
   * \frac{d \ln \alpha}{d \ln T} = -m \sqrt{\frac{T_r}{\alpha}}
   * @f]
   */
  Scalar ln_drv(const Scalar& tr) const noexcept override {
    using std::sqrt;
    const auto sqrt_tr = sqrt(tr);
    const auto a = 1 + m_ * (1 - sqrt_tr);
    return -m * sqrt_tr / a;
  }

 private:
  Scalar m_;
};

}  // namespace eos

#endif  // EOS_SOAVE_CORRECTION_FACTOR_HPP