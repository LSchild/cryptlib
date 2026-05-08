//
// Created by leonard on 11/12/24.
//

#include <cassert>
#include "utils.h"

RLWECiphertext encrypt_rlwe(std::shared_ptr<RingGSWCryptoParams> &params, NativePoly& msg, NativePoly& skntt) {
    auto PolyP = params->GetPolyParams();
    NativeInteger Q = params->GetQ();
    auto dug = DiscreteUniformGeneratorImpl<NativeVector>();
    dug.SetModulus(Q);
    auto dgg = params->GetDgg();
    msg.SetFormat(EVALUATION);

    auto A = NativePoly(dug, PolyP, EVALUATION);
    //A *= 0;
    auto B = A * skntt;
    B += msg;
    auto E = NativePoly(dgg, PolyP, COEFFICIENT);
    E.SetFormat(EVALUATION);
    B += E;

    std::vector<NativePoly> rlwe_vec = {A, B};
    return std::make_shared<RLWECiphertextImpl>(rlwe_vec);
}

NativeInteger DotProduct(const NativeVector& lhs, const NativeVector& rhs) {
    assert(lhs.GetLength() == rhs.GetLength());
    assert(lhs.GetModulus() == rhs.GetModulus());

    const auto& q = rhs.GetModulus();
    auto n = rhs.GetLength();

    NativeInteger mu = q.ComputeMu();
    NativeInteger result = 0;
    for(int i = 0; i < n; i++) {
        auto prod_i = lhs.at(i).ModMulFast(rhs.at(i), q, mu);
        result.ModAddEq(prod_i, q);
    }

    return result;
}


LWECiphertext encrypt_lwe(const NativeVector& sk, const NativeInteger& msg, double std) {

    auto q = sk.GetModulus();
    auto n = sk.GetLength();

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    DiscreteGaussianGeneratorImpl<NativeVector> dgg(std);

    auto A = dug.GenerateVector(n, q);
    auto B = DotProduct(A, sk);
    B.ModAddEq(msg, q);
    B.ModAddEq(dgg.GenerateInteger(q), q);

    return std::make_shared<LWECiphertextImpl>(A, B);
}

RingGSWSample encrypt_rgsw(std::shared_ptr<RingGSWCryptoParams> &params, NativePoly& msg, NativePoly& skntt) {

    auto digits = params->GetDigitsG();
    auto base = params->GetBaseG();

    std::vector<RLWECiphertext> m;
    std::vector<RLWECiphertext> msm;

    NativePoly msg_poly = msg.Clone();
    msg_poly.SetFormat(EVALUATION);

    NativePoly msm_poly = -msg_poly * skntt;

    for(int i =0; i < digits; i++) {
        m.push_back(encrypt_rlwe(params, msg_poly, skntt));
        msm.push_back(encrypt_rlwe(params, msm_poly, skntt));

        msg_poly *= base;
        msm_poly *= base;
    }

    return {params, msm, m};
}

bool IsPowerOf2(uint64_t x) {
    return (x != 0) and ((x & (x - 1)) == 0);
}

LWECiphertext ModSwitch(const LWECiphertext& ct, uint64_t target_modulus) {
    return ModSwitch(ct->GetA(), ct->GetB(), target_modulus);
}

LWECiphertext ModSwitch(const NativeVector& A, const NativeInteger& B, uint64_t target_modulus) {

    auto source_modulus = A.GetModulus();
    auto A_new = A.MultiplyAndRound(target_modulus, source_modulus).Mod(target_modulus);
    auto B_new = B.MultiplyAndRound(target_modulus, source_modulus).Mod(target_modulus);
    A_new.SetModulus(target_modulus);

    return std::make_shared<LWECiphertextImpl>(A_new, B_new);
}

std::vector<std::vector<uint32_t>> GenerateTestDatabase(uint64_t n_records, uint64_t record_size_bits) {
    // change seeds as needed
    uint64_t x=123456789, y=362436069, z=521288629;
    std::vector<std::vector<uint32_t>> database(n_records);
    // I know.
    const uint32_t uint32_size_bits = 32;
    // note we assume that it's divisible by 32. Sue me.
    auto n_entries_per_record = std::max(record_size_bits / uint32_size_bits, 1ul);

    uint64_t t;
    for(uint64_t i = 0; i < n_records; i++) {
        std::vector<uint32_t> record_vector(n_entries_per_record, 0);
        for(uint64_t j = 0; j < n_entries_per_record; j++) {
            // xorshift start
            x ^= x << 16;
            x ^= x >> 5;
            x ^= x << 1;

            t = x;
            x = y;
            y = z;
            z = t ^ x ^ y;
            // xorshift end
            record_vector[j] = z;
        }
        database[i] = std::move(record_vector);
    }

    return database;
}

void print_record_as_vec(std::vector<uint32_t>& record, uint32_t bits_per_print) {
    assert((32 % bits_per_print) == 0);
    auto rounds = 32 / bits_per_print;
    auto mask = (1u << bits_per_print) - 1;
    for(auto& v : record) {
        for(uint32_t i = 0; i < rounds; i++) {
            std::cerr << ((v >> (i * bits_per_print)) & mask) << " ";
        }
    }
    std::cerr << std::endl;
}

NativeInteger decrypt_lwe(const NativeVector& sk, const LWECiphertext& ct) {
    NativeVector sk_clone = sk;
    sk_clone.SwitchModulus(ct->GetModulus());
    NativeInteger phase = ct->GetB();
    phase.ModSubEq(DotProduct(sk_clone, ct->GetA()), ct->GetModulus());
    return phase;
}

LWECiphertext RLWESampleExtract(lbcrypto::NativePoly &rlwe_A, lbcrypto::NativePoly &rlwe_B, uint32_t idx) {


    auto N = rlwe_A.GetLength();
    auto Q = rlwe_B.GetModulus();

    NativeVector A(N, Q);
    NativeInteger B = rlwe_B.at(idx);

    for (int i = 0; i <= idx; i++) {
        A[i] = rlwe_A[idx - i];
    }

    for (uint32_t i = idx + 1; i < N; i++) {
        A[i] = Q.ModSub(rlwe_A[N + idx - i], Q);
    }


    return std::make_shared<LWECiphertextImpl>(A, B);
}


