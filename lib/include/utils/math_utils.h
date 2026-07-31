//
// Created by leonard on 4/27/26.
//

#ifndef LARGE_FUNCTIONS_MATH_UTILS_H
#define LARGE_FUNCTIONS_MATH_UTILS_H

#include <cstring>
#include <complex>
#include <cstdint>
#include <vector>

void EvalNegacyclicAutomorphism(uint64_t* __restrict output, const uint64_t* __restrict input, uint64_t index, uint64_t n, uint64_t Q);

uint64_t IntLog2(uint64_t input);


inline int64_t UnsignedToSignedRepr(uint64_t val, uint64_t modulus) {
    const auto mod_half = modulus >> 1;
    if (val < mod_half) {
        return int64_t(val);
    } else {
        uint64_t abs_inv = modulus - val;
        return -int64_t(abs_inv);
    }
}

inline uint64_t SignedToUnsignedRepr(int64_t val, uint64_t modulus) {
    if (val < 0) {
        return modulus - uint64_t(std::abs(val));
    } else {
        return uint64_t(val);
    }
}



/**
 * Estimates the sample mean of *sequences* using an unbiased estimator & numerically stable Kahan summation
 * @param seq Total sequence containing the values
 * @param n_samples Number of sequences
 * @param len Sequence length
 * @param modulus Modulus Q of the ring Z_Q in which the variables live
 * @return Sequence sample mean
 */
std::vector<long double> EstimateSampleMean(uint64_t* seq, uint64_t n_samples, uint64_t len, uint64_t modulus);

/**
 * Estimates the sample variance of *sequences* using an unbiased estimator & numerically stable Kahan summation
 * @param seq Total sequence containing the values
 * @param n_samples Number of sequences
 * @param len Sequence length
 * @param modulus Modulus Q of the ring Z_Q in which the variables live
 * @return Sequence sample variance
 */
std::vector<long double> EstimateSampleVariance(uint64_t* seq, uint64_t n_samples, uint64_t len, uint64_t modulus);

/**
 * Determines the probability that a variable distributed by a (sub) Gaussian falls outside a decryptable
 * interval, i.e. that the noise term is outside the interval [-Q/2p, Q/2p] where
 * Q is the ciphertext modulus and P is the plaintext modulus
 *
 * @param variance Concrete variance
 * @param modulus Ciphertext modulus
 * @param plaintext_modulus Plaintext modulus
 * @return Probability of decryption failure
 */
long double EstimateFailureProbability(long double variance, uint64_t modulus, uint64_t plaintext_modulus);

/**
 * Given a value reverses the order of bits
 * @tparam IntType Type of value
 * @param n Concrete value
 * @return Input value with reversed bits
 */
template <class IntType>
IntType ReverseBitsGeneric(IntType n) {

#ifdef __clang__
    const bool clang_present = true;
#else
    const bool clang_present = false;
#endif

    if constexpr (clang_present) {

        if constexpr (std::is_same_v<IntType, uint8_t>) {
            return __builtin_bitreverse8(n);
        }

        if constexpr (std::is_same_v<IntType, uint16_t>) {
            return __builtin_bitreverse16(n);
        }

        if constexpr (std::is_same_v<IntType, uint32_t>) {
            return __builtin_bitreverse32(n);
        }
        if constexpr (std::is_same_v<IntType, uint8_t>) {
            return __builtin_bitreverse64(n);
        }
    }


    short bits = sizeof(IntType) * 8;
    IntType mask = ~IntType(0); // equivalent to uint32_t mask = 0b11111111111111111111111111111111;

    while (bits >>= 1) {
        mask ^= mask << (bits); // will convert mask to 0b00000000000000001111111111111111;
        n = (n & ~mask) >> bits | (n & mask) << bits; // divide and conquer
    }

    return n;

}

inline long double RelativeError(long double x, long double x_expected) {
    // avoid division by 0
    if (std::abs(x_expected) == 0.0) {
        x_expected = 1e-30;
    }

    return std::abs(x - x_expected) / std::abs(x_expected);
}

#endif //LARGE_FUNCTIONS_MATH_UTILS_H
