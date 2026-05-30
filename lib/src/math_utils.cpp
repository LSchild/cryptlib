//
// Created by leonard on 4/27/26.
//
#include <cstdint>
#include <vector>
#include <algorithm>

#include "utils/math_utils.h"
#include "hexl/hexl.hpp"

uint64_t IntLog2(uint64_t input) {

    uint64_t res = 0;
    uint64_t acc = 1;
    while (input > acc) {
        res++;
        acc <<= 1;
    }

    return res;
}

std::vector<long double> estimate_mean(uint64_t* seq, uint64_t n_samples, uint64_t sample_len, uint64_t modulus) {
    std::vector<long double> means(sample_len, 0.0);
    std::vector<long double> comp(sample_len, 0.0);

    // numerically stable sum, do not optimize
    // see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    for(uint64_t i = 0; i < n_samples; i++) {
        for(uint64_t j = 0; j < sample_len; j++) {
            long double v_ij = unsigned_to_signed_repr(seq[i * sample_len + j], modulus);
            auto y = v_ij - comp[j];
            auto t = means[j] + y;
            comp[j] = (t - means[j]) - y;
            means[j] = t;
        }
    }

    std::transform(means.begin(), means.end(), means.begin(), [n_samples](long double x) {return x / ((long double) n_samples - 1);});
    return means;
}

std::vector<long double> estimate_variance(uint64_t* seq, uint64_t n_samples, uint64_t sample_len, uint64_t modulus) {

    auto means = estimate_mean(seq, n_samples, sample_len, modulus);
    std::vector<long double> variances(sample_len, 0.0);
    std::vector<long double> comp(sample_len, 0.0);

    // numerically stable sum, do not optimize
    // see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    for(uint64_t i = 0; i < n_samples; i++) {
        for(uint64_t j = 0; j < sample_len; j++) {
            long double v_ij = unsigned_to_signed_repr(seq[i * sample_len + j], modulus);
            auto diff_square = (v_ij - means[j]) * (v_ij - means[j]);
            auto y = diff_square - comp[j];
            auto t = variances[j] + y;
            comp[j] = (t - variances[j]) - y;
            variances[j] = t;
        }
    }

    std::transform(variances.begin(), variances.end(), variances.begin(), [n_samples](long double x) {return x / ((long double) (n_samples - 1));});
    return variances;
}

void EvalNegacyclicAutomorphism(uint64_t* __restrict output, const uint64_t* __restrict input, uint64_t index, uint64_t n, uint64_t Q) {

    // todo in-place algo
    auto n_bits = IntLog2(n);
    auto n_mod_mask = n - 1;

    for(uint64_t buffer_idx = 0; buffer_idx < n; buffer_idx++) {
        auto permuted_index = buffer_idx * index;
        auto new_index = permuted_index & n_mod_mask;
        bool negate = ((permuted_index >> n_bits) & 1) == 1;
        auto read_val = input[buffer_idx];
        auto write_val = negate ? intel::hexl::SubUIntMod(0, read_val, Q) : read_val;
        output[new_index] = write_val;
    }

}