//
// Created by leonard on 8/8/26.
//

#include "backend/solinas_2_32_1_m1_1.h"
#include "utils/generic_utils.h"


Solinas_2_32_1_M1_1::Solinas_2_32_1_M1_1(uint64_t dimension) : m_ntt_dim(dimension) {
    Precompute();
}

uint64_t Solinas_2_32_1_M1_1::AddMod(uint64_t left, uint64_t right) {
    auto a = left + K;
    auto t = a + right;

    if (t < a) {
        return t;
    } else {
        return t - K;
    }
}

uint64_t Solinas_2_32_1_M1_1::SubMod(uint64_t left, uint64_t right) {
    auto z = left - right;
    if (left < right) {
        return z - K;
    } else {
        return z;
    }
}

uint64_t Solinas_2_32_1_M1_1::MulMod(uint64_t left, uint64_t right) {
    __uint128_t prod = __uint128_t(left) * right;
    return ReduceMod(prod);
}

uint64_t Solinas_2_32_1_M1_1::FMAMod(uint64_t left, uint64_t right, uint64_t summand) {
    __uint128_t prod = __uint128_t(left) * right;
    prod += summand;

    return ReduceMod(prod);
}

uint64_t Solinas_2_32_1_M1_1::PowMod(uint64_t base, uint64_t exponent) {
    exponent %= Q - 1;

    if (exponent == 0)
        return 1;
    if (exponent == 1)
        return base;

    uint64_t res = 1;
    uint64_t acc = base;
    while (exponent > 1) {
        if (exponent & 1)
            res = mod_mul(res, acc);
        acc = mod_mul(acc, acc);
        exponent >>= 1;
    }

    return res;
}

uint64_t Solinas_2_32_1_M1_1::InvMod(uint64_t left) {
    return PowMod(left, Q - 2);
}

#define UNROLL_FACTOR 4
#define UNROLL_FACTOR_SHIFT 2

void Solinas_2_32_1_M1_1::AddModElw(uint64_t *result, uint64_t *left, uint64_t *right, uint64_t len) {

    for(uint64_t i = 0; i < len; i++) {
        result[i] = AddMod(left[i], right[i]);
    }

}

void Solinas_2_32_1_M1_1::AddEqModElw(uint64_t *result, uint64_t *right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = AddMod(result[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::SubModElw(uint64_t *result, uint64_t *left, uint64_t *right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = SubMod(left[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::SubEqModElw(uint64_t *left_out, uint64_t *right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = SubMod(left_out[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::MulModElw(uint64_t *result, uint64_t *left, uint64_t *right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = MulMod(left[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::MulEqModElw(uint64_t *left_out, uint64_t *right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = MulMod(left_out[i], right[i]);
    }
}

uint64_t Solinas_2_32_1_M1_1::GetModulus() const {
    return Q;
}

uint64_t Solinas_2_32_1_M1_1::GetDimension() const {
    return m_ntt_dim;
}

void
Solinas_2_32_1_M1_1::FMAModElw(uint64_t *result, uint64_t *left, uint64_t *right, uint64_t *summand, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = FMAMod(left[i], right[i], summand[i]);
    }
}

void Solinas_2_32_1_M1_1::FMAEqModElw(uint64_t *left_out, uint64_t *right, uint64_t *summand, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = FMAMod(left_out[i], right[i], summand[i]);
    }
}

void Solinas_2_32_1_M1_1::MulSModElw(uint64_t *result, uint64_t *left, uint64_t right, uint64_t len) {

    if (right == 0) {
        ZERO_UINT64_ARR(result, len);
        return;
    }

    if (right == 1) {
        std::memcpy(result, left, sizeof(uint64_t) * len);
        return;
    }

    if (right == 2) {
        AddModElw(result, left, left, len);
        return;
    }



    for(uint64_t i = 0; i < len; i++) {
        result[i] = MulMod(left[i], right);
    }
}

void Solinas_2_32_1_M1_1::MulSEqModElw(uint64_t *left_out, uint64_t right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = MulMod(left_out[i], right);
    }
}

void Solinas_2_32_1_M1_1::Precompute() {
    // TODO
}