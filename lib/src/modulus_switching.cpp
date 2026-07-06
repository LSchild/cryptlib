//
// Created by leonard on 7/6/26.
//

#include <cstdint>
#include <iostream>
#include "modulus_switching.h"

long double MSErrorBT(long double input_variance, uint64_t source_modulus, uint64_t target_modulus, uint64_t expected_l0) {
    long double Qs = source_modulus;
    auto Qs2 = Qs * Qs;
    long double Qt = target_modulus;
    auto Qt2 = Qt * Qt;

    long double principal = (Qs2 * input_variance) / Qt2;

    long double hamming_mean = expected_l0;

    return principal + hamming_mean / 12.0;
}

long double EstimateModulusSwitchingVariance(long double input_variance, uint64_t source_modulus, uint64_t target_modulus, KeyDistribution kd, uint64_t expected_l0) {

    if (kd == BINARY or kd == TERNARY) {
        return MSErrorBT(input_variance, source_modulus, target_modulus, expected_l0);
    } else {
        std::cerr << "Modulus Switching Variance computation implemented for Binary/Ternary keys only for now" << std::endl;
        std::exit(1);
    }

}