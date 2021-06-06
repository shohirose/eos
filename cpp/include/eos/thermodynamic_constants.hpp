#ifndef EOS_THERMODYNAMIC_CONSTANTS_HPP
#define EOS_THERMODYNAMIC_CONSTANTS_HPP

#include <type_traits>

namespace eos {

/// @{
/// @name Thermodynamic constants
///
/// The following definitions are based on the International System of Units(SI)
/// 9th edition 2019.

/// @brief Gas constant [J/mol-K]
template <typename T>
inline constexpr auto gas_constant()
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  return static_cast<T>(8.314'462'618'153'24);
}

/// @brief Avogadro constant [1/mol]
template <typename T>
inline constexpr auto avogadro_constant()
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  return static_cast<T>(6.022'140'76e23);
}

/// @brief Boltzman constant [J/K]
template <typename T>
inline constexpr auto boltzman_constant()
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  return static_cast<T>(1.380'649e-23);
}

/// @brief Planck constant [J-s]
template <typename T>
inline constexpr auto planck_consant()
    -> std::enable_if_t<std::is_floating_point_v<T>, T> {
  return static_cast<T>(6.626'070'15e-34);
}

/// @}

}  // namespace eos

#endif  // EOS_THERMODYNAMIC_CONSTANTS_HPP