//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_UTILS_H
#define LARGE_FUNCTIONS_UTILS_H

#include <chrono>
#include "setup.h"
#include "rgsw.h"

#define TIC_HD std::chrono::high_resolution_clock::now();

/**
 *  Computes acc[i] += left[i] * right[i] for i in [0..N)
 * @tparam N array dimension
 * @param acc accumulation buffer
 * @param left left buffer
 * @param right right buffer
 */
template<uint32_t N> void vector_accumulate(uint64_t* acc, uint32_t* left, uint32_t* right) {
    for (uint32_t i = 0; i < N; i++) {
        // @clang @gcc if I multiply two 32bit numbers and add it to a 64bit accumulator
        // one would think that i don't need to cast the result...
        acc[i] += left[i] * (right[i] * 1ull);
    }
}

template<uint32_t N> void vec_add_eq(uint64_t* lhs, const uint64_t* rhs) {
    for(uint32_t i = 0; i < N; i++) {
        lhs[i] += rhs[i];
    }
}

/**
 * Creates a LWE sample
 * @param sk secret key, a vector in \ZZ_q^n
 * @param msg the message in [0, q)
 * @param std the standard deviation for the error term
 * @return LWE_{sk}(m)
 */
LWECiphertext encrypt_lwe(const NativeVector& sk, const NativeInteger& msg, double std = 0);

/**
 * Creates a RLWE sample
 * @param params Ring parameters i.e. modulus & dimension
 * @param msg message polynomial in \ZZ_Q[X]/(X^N + 1) where N is a power of two
 * @param skntt secret polynomial in NTT domain
 * @return RLWE(msg)
 */
RLWECiphertext encrypt_rlwe(std::shared_ptr<RingGSWCryptoParams> &params, NativePoly& msg, NativePoly& skntt);

/**
 * Decrypts a LWE sample
 * @param sk secret key
 * @param ct LWE(m)
 * @return m
 */
NativeInteger decrypt_lwe(const NativeVector& sk, const LWECiphertext& ct);

/**
 * Creates a RGSW sample
 * @param params Ring parameters i.e. modulus & dimension
 * @param msg message polynomial in \ZZ_Q[X]/(X^N + 1) where N is a power of two
 * @param skntt secret polynomial in NTT domain
 * @return RGSW(msg)
 */
RingGSWSample encrypt_rgsw(std::shared_ptr<RingGSWCryptoParams> &params, NativePoly& msg, NativePoly& skntt);

/**
 * Given two vector in \ZZ_q^n, computes the euclidean inner product
 * @param lhs first vector
 * @param rhs second vector
 * @return \sum_i^n lhs_i * rhs_i
 */
NativeInteger DotProduct(const NativeVector& lhs, const NativeVector& rhs);

/**
 *  Modulus switch operation for a ciphertext
 * @param ct LWE_sk(delta(q, t) * m) under modulus q
 * @param target_modulus output modulus P
 * @return LWE_sk(delta(P, t) * m) under modulus P
 */
LWECiphertext ModSwitch(const LWECiphertext& ct, uint64_t target_modulus);
LWECiphertext ModSwitch(const NativeVector& A, const NativeInteger& B, uint64_t target_modulus);

/**
 * Generates a dummy database given parameters
 * @param n_records number of records it shall contain
 * @param record_size_bits the size of each record
 * @return a database D[x] where x in [0, n_records) and D[x] is a (vectorized) record of given size.
 */
std::vector<std::vector<uint32_t>> GenerateTestDatabase(uint64_t n_records, uint64_t record_size_bits);

/**
 * Method to print a record
 * @param record the vectorized record
 * @param bits_per_print how many bits should be printed before printing a whitespace
 */
void print_record_as_vec(std::vector<uint32_t>& record, uint32_t bits_per_print);

/**
 *
 * @param x integer
 * @return true if x is a power of 2
 */
bool IsPowerOf2(uint64_t x);

LWECiphertext RLWESampleExtract(lbcrypto::NativePoly &rlwe_A, lbcrypto::NativePoly &rlwe_B, uint32_t idx);

#endif //LARGE_FUNCTIONS_UTILS_H
