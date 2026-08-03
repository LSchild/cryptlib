//
// Created by leonard on 7/17/26.
//

#ifndef LARGE_FUNCTIONS_SAMPLE_EXTRACTION_H
#define LARGE_FUNCTIONS_SAMPLE_EXTRACTION_H

#include <cstdint>
#include <vector>

/**
 * Extracts a LWE sample from a given RLWE over the ring Z_q[X]/(X^n + 1) n= a power of 2, and for the given index
 * @tparam IntType Template param for the integer type
 * @param output_lwe Pointer to output LWE sample buffer
 * @param input_rlwe Pointer to input RLWE sample buffer
 * @param coef_index Index for which to extract the LWE sample
 * @param dim Dimension of the ring n
 * @param modulus Ring modulus
 */
template<typename IntType>
void SampleExtract(IntType *__restrict output_lwe, IntType *__restrict input_rlwe, uint64_t coef_index, uint64_t dim, IntType modulus) {
    IntType B = input_rlwe[dim + coef_index];

    // p(X) * X^{-l} -> aut_{-1}(p(X) * X^{-l}) -> p(X^{-1}) * X^{l}
    // [p_0, -p_{N - 1}, -p_{N - 2}, ... , -p_1] -> [p_{l}, p_{l-1}, p_{l - 2}, ... , p_0, -p_{N - 1}, ..., -p_{l + 1}]
    for(uint64_t i = 0; i <= coef_index; i++) {
        output_lwe[i] = input_rlwe[coef_index - i];
    }
    for(uint64_t i = coef_index + 1; i < dim; i++) {
        auto v = input_rlwe[dim - i];
        if (v != 0) [[likely]] {
            output_lwe[i] = modulus - v;
        }
    }
    output_lwe[dim] = input_rlwe[dim + coef_index];
}

/**
 * Extracts a LWE sample from a given RLWE over the ring Z_q[X]/(X^n + 1) n= a power of 2, and for the given index
 * @tparam IntType Template param for the integer type
 * @param output_lwe Output vector for LWE sample
 * @param input_rlwe Input vector for RLWE sample
 * @param coef_index Index for which to extract the LWE sample
 * @param dim Dimension of the ring n
 * @param modulus Ring modulus
 */
template<typename IntType, typename Allocator>
void SampleExtract(std::vector<IntType, Allocator>& output_lwe, std::vector<IntType, Allocator>& input_rlwe, uint64_t coef_index, uint64_t dim, IntType modulus) {
    for(uint64_t i = 0; i <= coef_index; i++) {
        output_lwe[i] = input_rlwe[coef_index - i];
    }
    for(uint64_t i = coef_index + 1; i < dim; i++) {
        auto v = input_rlwe[dim - i];
        if (v != 0) [[likely]] {
            output_lwe[i] = modulus - v;
        }
    }
}

#endif //LARGE_FUNCTIONS_SAMPLE_EXTRACTION_H
