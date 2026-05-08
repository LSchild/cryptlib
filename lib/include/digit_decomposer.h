//
// Created by lschild on 19/11/24.
//

#ifndef LARGE_FUNCTIONS_DIGIT_DECOMPOSER_H
#define LARGE_FUNCTIONS_DIGIT_DECOMPOSER_H

#include "setup.h"
#include "functional_bootstrap.h"
#include "lwe_key_switching.h"

struct BitDecomposerParams {
    /* LWE dimension */
    uint32_t m_n;
    /* RLWE dimension */
    uint32_t m_N;
    /* LWE modulus (blind-rotation) */
    uint64_t m_q;
    /* LWE modulus (input) */
    uint64_t m_q_large;
    uint32_t m_q_large_bits;
    /* RLWE modulus */
    uint64_t m_Q;

    /* standard deviation */
    double m_std;

    /** bases **/
    /* for blind-rotation key */
    uint32_t m_br_basis;

    /* for automorphism use during bootstrapping key */
    uint32_t m_auto_key_basis;

    /* amount of bits that (should) be removed in every iteration */
    uint32_t m_clear_bits;

    /* error width */
    uint32_t m_error_width;

};

struct DigitDecompositionParameters {
    /* LWE dimension */
    uint32_t m_n;
    /* RLWE dimension */
    uint32_t m_N;
    /* LWE modulus (blind-rotation) */
    uint64_t m_q;
    /* LWE modulus (input) */
    uint64_t m_q_large;
    uint32_t m_q_large_bits;
    /* RLWE modulus */
    uint64_t m_Q;

    /* standard deviation for RLWE/RGSW/KSK*/
    double m_std;

    /** bases **/
    /* for blind-rotation key */
    uint32_t m_br_basis;

    /* amount of bits that (should) be removed in every iteration */
    uint32_t m_clear_bits;

    /* error width */
    uint32_t m_error_width;

    /* key switch modulus */
    uint64_t m_ksk_q;

    /* key switch basis */
    uint32_t m_ksk_basis;
};

class DigitDecomposer {
public:

    DigitDecomposer() = default;

    explicit DigitDecomposer(const DigitDecompositionParameters& params);

    void SetKeys(NativePoly& ring_sk, NativeVector& lwe_sk);

    std::vector<LWECiphertext> DigitDecompose(LWECiphertext& input, int32_t max_digits = 0);

    std::vector<LWECiphertext> BitDecompose(LWECiphertext& input, int32_t max_bits = 0);

    const std::shared_ptr<RingGSWCryptoParams> GetRingParams() const;

private:

    LWECiphertext HomFloor(LWECiphertext& input);

    std::shared_ptr<LWECryptoParams> m_ksk_params;
    std::shared_ptr<RingGSWCryptoParams> m_br_params;

    FunctionalBootstrapEngine m_boot_engine;
    LWEKeySwitchingKey m_lwe_ksk;

    uint64_t m_alpha;
    uint64_t m_beta;

    NativePoly m_ring_sk;
    NativeVector m_lwe_sk;

};


class BitDecomposer {
public:

    BitDecomposer() {};

    explicit BitDecomposer(BitDecomposerParams& params);

    void SetKeys(NativePoly sk_rlwe, LWEPrivateKey& sk_lwe);

    std::vector<LWECiphertext> TightDecompose(ConstLWECiphertext& ct);

    std::vector<LWECiphertext> BitDecompose(ConstLWECiphertext& ct);

    std::vector<LWECiphertext> RelaxedDecompose(ConstLWECiphertext& ct, uint32_t padding_bits);

    std::pair<LWECiphertext, std::vector<LWECiphertext>> DoFloorAlpha(uint64_t Q, uint64_t q, NativeVector& A, NativeInteger& B);

    std::shared_ptr<RingGSWCryptoParams> m_br_params;

    FunctionalBootstrapEngine m_engine;
    BitDecomposerParams m_params;
    LWEKeySwitchingKey m_lwe_ksk;

    uint64_t m_qks;

    /* Debugging */
    NativePoly m_sk_rlwe;
    NativePoly m_sk_rlwe_ntt;
    NativeVector m_sk_lwe;
};

#endif //LARGE_FUNCTIONS_DIGIT_DECOMPOSER_H
