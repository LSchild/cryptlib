//
// Created by leonard on 4/24/26.
//

#ifndef LARGE_FUNCTIONS_GADGET_DECOMP_H
#define LARGE_FUNCTIONS_GADGET_DECOMP_H


#include <cstdint>
#include <algorithm>

#include "hexl/hexl.hpp"

/**
 * Performs an UNSIGNED gadget decomposition and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Output buffer, should have size at least n * k * digits
 * @param source Input buffer, expected to have size n
 * @param n Input dimension
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring modulus
 */
template<uint32_t k>
void UnsignedGadgetDecomposeRep(uint64_t *const result, const uint64_t *const source, uint64_t n, uint64_t digits, uint64_t basis_bits, uint64_t modulus_bits) {
    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    const uint64_t mask = (1ull << basis_bits) - 1;
    uint64_t shift = modulus_bits - basis_bits;

    for(uint64_t d_i = 0; d_i < digits; d_i++) {
        uint64_t* result_row = result + d_i * (k + 1) * n;
        for (uint64_t i = 0; i < n; i++) {
            result_row[i] = (source[i] >> shift) & mask;
        }
        shift -= basis_bits;

        for(uint32_t kk = 0; kk < k; kk++) {
            std::copy(result_row, result_row + n, result_row + (kk + 1) * n);
        }
    }
}
#define UnsignedDigitDecompose UnsignedGadgetDecomposeRep<0>
#define UnsignedDigitDecomposeRep2 UnsignedGadgetDecomposeRep<1>


/**
 * Performs a SIGNED gadget decomposition and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Output buffer, should have size at least n * k * digits
 * @param source Input buffer, expected to have size n
 * @param n Input dimension
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring modulus
 */
template<uint32_t k>
void SignedGadgetDecomposeRep(uint64_t *const result, uint64_t *const source, uint64_t n, uint64_t digits, uint64_t basis_bits, uint64_t modulus, uint64_t modulus_bits) {
    // We aim to obtain a decomposition w.r.t a basis B of a vector of x such that the digits lie in [-B/2,B/2)
    // Then, note that x in [-B/2, B/2) implies x + B/2 in [0, B) i.e. a normal unsigned digit decomp
    // Next, let corr = sum_i B^i B//2 and for x' = x + corr mod Q let its digits be [x0',x1',...]
    // Then, we obtain the signed digits by setting xi = xi' - B//2 which is correct as
    // \sum B^i xi = \sum B^i (xi' - B//2) = x + corr - \sum B^i B//2 = x

    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    uint64_t corrector = digits * basis_bits == modulus_bits ? 0 : 1ull << (modulus_bits - digits * basis_bits - 1);
    for(uint64_t i = 0; i < digits; i++) {
        corrector += 1ull << (modulus_bits - basis_bits * i + basis_bits - 1);
    }
    corrector = corrector >= modulus ? corrector - modulus : corrector;
    intel::hexl::EltwiseAddMod(source, source, corrector, n, modulus);
    UnsignedGadgetDecomposeRep<k>(result, source, n, digits, basis_bits, modulus_bits);
    intel::hexl::EltwiseSubMod(result, result, 1ull << (basis_bits - 1), n, modulus);
}

#define SignedDigitDecompose SignedGadgetDecomposeRep<0>
#define SignedDigitDecomposeRep2 SignedGadgetDecomposeRep<1>

/**
 * Performs an UNSIGNED gadget decomposition with the output in NTT format, and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Result buffer of size at least n * k * digits
 * @param source Source buffer assumed of size n
 * @param ntt Pointer to ntt engine (subject to change)
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring modulus
 */
template<uint32_t k>
void UnsignedGadgetDecomposeRepNTT(uint64_t *const result, const uint64_t *const source, std::shared_ptr<intel::hexl::NTT> ntt, uint64_t digits, uint64_t basis_bits, uint64_t modulus_bits) {
    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    const uint64_t mask = (1ull << basis_bits) - 1;
    const uint64_t n = ntt->GetDegree();
    uint64_t shift = modulus_bits - basis_bits;

    for(uint64_t d_i = 0; d_i < digits; d_i++) {
        uint64_t* result_row = result + d_i * (k + 1) * n;
        for (uint64_t i = 0; i < n; i++) {
            result_row[i] = (source[i] >> shift) & mask;
        }
        shift -= basis_bits;
        ntt->ComputeForward(result_row, result_row, 1, 1);
        for(uint32_t kk = 0; kk < k; kk++) {
            std::copy(result_row, result_row + n, result_row + (kk + 1) * n);
        }
    }
}

#define UnsignedDigitDecomposeNTT UnsignedGadgetDecomposeRepNTT<0>
#define UnsignedDigitDecomposeRep2NTT UnsignedGadgetDecomposeRepNTT<1>

/**
 * Performs an SIGNED gadget decomposition with the output in NTT format, and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Result buffer of size at least n * k * digits
 * @param source Source buffer assumed of size n
 * @param ntt Pointer to ntt engine (subject to change)
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring modulus
 */
template<uint32_t k>
void SignedGadgetDecomposeRepNTT(uint64_t *const result, const uint64_t *const source, std::shared_ptr<intel::hexl::NTT> ntt,  uint64_t digits, uint64_t basis_bits,uint64_t modulus_bits) {
    // We aim to obtain a decomposition w.r.t a basis B of a vector of x such that the digits lie in [-B/2,B/2)
    // Then, note that x in [-B/2, B/2) implies x + B/2 in [0, B) i.e. a normal unsigned digit decomp
    // Next, let corr = sum_i B^i B//2 and for x' = x + corr mod Q let its digits be [x0',x1',...]
    // Then, we obtain the signed digits by setting xi = xi' - B//2 which is correct as
    // \sum B^i xi = \sum B^i (xi' - B//2) = x + corr - \sum B^i B//2 = x

    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    const uint64_t modulus = ntt->GetModulus();
    const uint64_t n = ntt->GetDegree();
    const auto last_result_row = result + (digits - 1) * (k + 1) * n;

    uint64_t corrector = digits * basis_bits == modulus_bits ? 0 : 1ull << (modulus_bits - digits * basis_bits - 1);
    for(uint64_t i = 0; i < digits; i++) {
        corrector += 1ull << (modulus_bits - basis_bits * (i + 1) + basis_bits - 1);
    }
    corrector = corrector >= modulus ? corrector - modulus : corrector;

    // TODO: revisit correctness ?
    // Should be fine though
    intel::hexl::EltwiseAddMod(last_result_row, source, corrector, n, modulus);
    const uint64_t mask = (1ull << basis_bits) - 1;
    uint64_t shift = modulus_bits - basis_bits;

    for(uint64_t d_i = 0; d_i < digits; d_i++) {
        uint64_t* result_row = result + d_i * (k + 1) * n;
        for (uint64_t i = 0; i < n; i++) {
            result_row[i] = (last_result_row[i] >> shift) & mask;
        }
        shift -= basis_bits;
        intel::hexl::EltwiseSubMod(result_row, result_row, 1ull << (basis_bits - 1), n, modulus);
        ntt->ComputeForward(result_row, result_row, 1, 1);
        for(uint32_t kk = 0; kk < k; kk++) {
            std::copy(result_row, result_row + n, result_row + (kk + 1) * n);
        }
    }
}

#define SignedDigitDecomposeNTT SignedGadgetDecomposeRepNTT<0>
#define SignedDigitDecomposeRep2NTT SignedGadgetDecomposeRepNTT<1>

#endif //LARGE_FUNCTIONS_GADGET_DECOMP_H
