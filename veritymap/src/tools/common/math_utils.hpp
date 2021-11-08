//
// Created by Andrey Bzikadze on 2/26/21.
//

#pragma once

template<typename T, typename U>
T constexpr custom_pow(T base, U exponent) {
  static_assert(std::is_integral<U>(), "exponent must be integral");
  return exponent == 0 ? 1 : base * custom_pow(base, exponent - 1);
}
