//
// Created by leonard on 2/17/26.
//

#ifndef DIGIT_DECOMPOSER_FAST_H
#define DIGIT_DECOMPOSER_FAST_H
#include <cstdint>

#include "operators/automorphism_evaluation.h"
#include "blind_rotator.h"
#include "setup.h"

struct FastDigitDecompositionParameters {
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
    uint32_t m_br_digits;

    /* offset */
    uint64_t m_alpha;
    uint64_t m_beta;

    /* automorphism */
    uint32_t m_auto_base;
    double m_auto_std;
};

class FastDigitDecomposer {
public:

    FastDigitDecomposer() = default;

    explicit FastDigitDecomposer(const FastDigitDecompositionParameters& params);

    void SetKeys(NativePoly& ring_sk, NativeVector& lwe_sk);

    std::vector<LWECiphertext> DigitDecompose(LWECiphertext& input);

    std::vector<LWECiphertext> PhaseDecomp(LWECiphertext& input);

    RLWECiphertext HomTrunc(RLWECiphertext& acc);

    RLWECiphertext LWE2RLWE(LWECiphertext& ct);

    const std::shared_ptr<RingGSWCryptoParams> GetRingParams() const;

private:

    std::shared_ptr<RingGSWCryptoParams> m_br_params;
    std::shared_ptr<RingGSWCryptoParams> m_auto_params;


    FastBlindRotationKey m_br_key;
    std::vector<FastAutomorphismKey> m_decimation_keys;
    //std::vector<AutomorphismKey> m_decimation_keys;

    uint64_t m_alpha;
    uint64_t m_beta;

    NativePoly m_ring_sk;
    NativeVector m_lwe_sk;

};



#endif //DIGIT_DECOMPOSER_FAST_H
