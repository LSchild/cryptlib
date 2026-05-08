//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_BLIND_ROTATOR_H
#define LARGE_FUNCTIONS_BLIND_ROTATOR_H

#include "setup.h"
#include "rgsw.h"
#include "interfaces.h"


// This is a replacement for the blind rotation implementation in openfhe that does not drop the least significant digit
class BlindRotationKey {

public:

    explicit BlindRotationKey() = default;

    BlindRotationKey(std::shared_ptr<RingGSWCryptoParams>& params, NativePoly sk_rlwe, NativeVector sk_lwe);

    RLWECiphertext BlindRotate(RLWECiphertext& in, NativeVector A);

    NativePoly m_sk;
    NativePoly m_sk_ntt;
    NativeVector m_sk_lwe;

    std::vector<uint64_t> m_acc_buffer;
    std::shared_ptr<RingGSWCryptoParams> m_params;
    std::vector<std::pair<RGSWSample, RGSWSample>> m_brk;

    std::vector<uint64_t> m_monomials_raw;
    std::vector<NativePoly> m_monomials;

};

class FastBlindRotationKey {

public:
#define ALIGN_AT 32

    using Vector = std::vector<uint64_t, intel::hexl::AlignedAllocator<uint64_t, ALIGN_AT>>;

    explicit FastBlindRotationKey() = default;

    FastBlindRotationKey(std::shared_ptr<RingGSWCryptoParams>& params, NativePoly sk_rlwe, NativeVector sk_lwe, uint32_t step_size);

    RLWECiphertext BlindRotate(RLWECiphertext& in, NativeVector A);

    RLWECiphertext BlindRotate1(RLWECiphertext& in, NativeVector A);
    RLWECiphertext BlindRotate2(RLWECiphertext& in, NativeVector A);
    RLWECiphertext BlindRotate3(RLWECiphertext& in, NativeVector A);

    void KeyGen1();

    void KeyGen2();

    void KeyGen3();

    void Step1(Vector&, uint32_t idx, uint32_t a_i);

    void Step2(Vector&, uint32_t idx, uint32_t a_i, uint32_t a_i_1);

    void Step3(Vector&, uint32_t idx, uint32_t a_i, uint32_t a_i_1, uint32_t a_i_2);

    void MakePoly2(uint64_t* buffer, uint64_t a_0, uint64_t a_1);

    void MakePoly3(uint64_t* buffer, uint64_t a_0, uint64_t a_1, uint64_t a_2);

    void ExtProd(uint64_t* rgsw, uint64_t* digits, uint64_t* mon, uint64_t* buffer);

    void RLWEVectorProd(uint64_t* __restrict rlwe_prime, uint64_t* __restrict digit_poly, uint64_t* __restrict result_buffer);

    void RLWEVectorProd2(std::array<uint64_t* __restrict, 3> ops, std::array<uint64_t* __restrict, 3> mons, const uint64_t* __restrict digit_poly, uint64_t* __restrict result_buffer);

    uint64_t m_step_size;
    uint64_t m_rgsw_count_per_step;

    uint64_t m_rlwe_size;
    uint64_t m_rlwe_prime_size;
    uint64_t m_rgsw_size;
    uint64_t m_brk_chunk_size;

    NativePoly m_sk;
    NativePoly m_sk_ntt;
    NativeVector m_sk_lwe;

    Vector m_acc_buffer;
    std::shared_ptr<RingGSWCryptoParams> m_params;

    Vector m_brk;
    Vector m_monomials_raw;
    Vector m_one_ntt;
    Vector m_poly_buffer;

    Vector m_digit_buffer;

    std::unique_ptr<intel::hexl::NTT> m_engine;

};


#endif //LARGE_FUNCTIONS_BLIND_ROTATOR_H
