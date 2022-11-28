#ifndef EOS_CORRECTION_FACTOR_TRAITS_HPP
#define EOS_CORRECTION_FACTOR_TRAITS_HPP

namespace eos {

template <typename CorrectionFactor>
class CorrectionFactorTraits {
 public:
  using Scalar = CorrectionFactor::Scalar;
};

}  // namespace eos

#endif  // EOS_CORRECTION_FACTOR_TRAITS_HPP