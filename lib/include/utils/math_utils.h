//
// Created by leonard on 4/27/26.
//

#ifndef LARGE_FUNCTIONS_MATH_UTILS_H
#define LARGE_FUNCTIONS_MATH_UTILS_H

#include <complex>
#include <cstdint>
#include <vector>

void EvalNegacyclicAutomorphism(uint64_t* __restrict output, const uint64_t* __restrict input, uint64_t index, uint64_t n, uint64_t Q);

uint64_t IntLog2(uint64_t input);

inline int64_t unsigned_to_signed_repr(uint64_t val, uint64_t modulus) {
    const auto mod_half = modulus >> 1;
    if (val < mod_half) {
        return int64_t(val);
    } else {
        uint64_t abs_inv = modulus - val;
        return -int64_t(abs_inv);
    }
}

inline uint64_t signed_to_unsigned_repr(int64_t val, uint64_t modulus) {
    if (val < 0) {
        return modulus - uint64_t(std::abs(val));
    } else {
        return uint64_t(val);
    }
}

std::vector<long double> estimate_mean(uint64_t* seq, uint64_t n_samples, uint64_t len, uint64_t modulus);

std::vector<long double> estimate_variance(uint64_t* seq, uint64_t n_samples, uint64_t len, uint64_t modulus);

long double EstimateFailureProbability(long double variance, uint64_t modulus, uint64_t plaintext_modulus);
#endif //LARGE_FUNCTIONS_MATH_UTILS_H
