//
// Created by leonard on 6/18/26.
//

#ifndef LARGE_FUNCTIONS_SOLINAS_2_32_1_M1_1_H
#define LARGE_FUNCTIONS_SOLINAS_2_32_1_M1_1_H

#include <cstdint>

struct Solinas_2_32_1_M1_1 {

    static const uint64_t modulus = 18446744069414584321ull;
    static const uint64_t K = (1ull << 32) - 1;

    static uint64_t mod_add(uint64_t x, uint64_t y) {
        auto a = x + K;
        auto t = a + y;

        if (t < a) {
            return t;
        } else {
            return t - K;
        }
    }

    static uint64_t mod_sub(uint64_t x, uint64_t y) {
        auto z = x - y;
        if (x < y) {
            return z - K;
        } else {
            return z;
        }
    }

    static uint64_t mod_mul(uint64_t x, uint64_t y) {
        __uint128_t z = x * y;
        uint64_t z_lo = z % (__uint128_t(1) << 64);
        uint64_t z_hi = z >> 64;
        uint64_t z_hi_0 = z_hi % (uint64_t(1) << 32);
        uint64_t z_hi_1 = z_hi >> 32;
        auto res = mod_sub(z_lo, z_hi_0+z_hi_1);
        return mod_add(res, z_hi_0 << 32);
    }

    static uint64_t mod_exp(uint64_t x, uint64_t e) {
        e %= modulus - 1;

        if (e == 0)
            return 1;
        if (e == 1)
            return x;

        uint64_t res = 1;
        uint64_t acc = x;
        while (e > 1) {
            if (e & 1)
                res = mod_mul(res, acc);
            acc = mod_mul(acc, acc);
            e >>= 1;
        }

        return res;
    }

};

#endif //LARGE_FUNCTIONS_SOLINAS_2_32_1_M1_1_H
