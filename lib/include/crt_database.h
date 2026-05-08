//
// Created by leonard on 2/17/25.
//

#ifndef CRT_DATABASE_H
#define CRT_DATABASE_H

#include <vector>
#include <cstdint>
#include <rlwe-ciphertext.h>

struct CRTDatabaseParams {

    uint64_t selector_modulus;
    uint64_t composite_ring_modulus;

    uint32_t factor_1;
    uint32_t factor_2;

    uint32_t poly_dim;

    uint32_t rlwe_prime_basis;
    uint32_t rlwe_prime_digits;

    uint32_t output_samples;
};




#endif //CRT_DATABASE_H
