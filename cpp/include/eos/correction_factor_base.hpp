#ifndef EOS_CORRECTION_FACTOR_BASE_HPP
#define EOS_CORRECTION_FACTOR_BASE_HPP

#include <eos/correction_factor_traits.hpp>

namespace eos {

template <typename Derived>
class CorrectionFactorBase {
 public:
  using Scalar = CorrectionFactorTraits<Derived>::Scalar;

  Scalar value(const Scalar& tr) const noexcept {
    this->derived().value(tr);
  }

  Scalar ln_drv(const Scalar& tr) const noexcept {
    return this->derived().ln_drv(tr);
  }

 private:
  Derived& derived() noexcept { return static_cast<Derived&>(*this); }

  const Derived& derived() const noexcept {
    return static_cast<const Derived&>(*this);
  }
};

}  // namespace eos

#endif  // EOS_CORRECTION_FACTOR_BASE_HPP