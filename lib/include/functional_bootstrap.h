//
// Created by lschild on 19/11/24.
//

#ifndef LARGE_FUNCTIONS_FUNCTIONAL_BOOTSTRAP_H
#define LARGE_FUNCTIONS_FUNCTIONAL_BOOTSTRAP_H

#include "blind_rotator.h"
#include "setup.h"
#include "operators/automorphism_evaluation.h"

class FunctionalBootstrapEngine {
public:

    explicit FunctionalBootstrapEngine() {}

    FunctionalBootstrapEngine(std::shared_ptr<RingGSWCryptoParams>& br_params,std::shared_ptr<RingGSWCryptoParams>& m_auto_params, NativePoly sk_rlwe, NativeVector sk_lwe);

    LWECiphertext Boot(ConstLWECiphertext& ct, BootFunction& F, uint64_t alpha, uint64_t output_modulus = 1);

    LWECiphertext SampleExtract(NativePoly& rlwe_A, NativePoly& rlwe_B);

    std::vector<LWECiphertext> BitDecompose(ConstLWECiphertext& ct, uint32_t n_msb);

    std::shared_ptr<RingGSWCryptoParams> m_br_params;
    std::shared_ptr<RingGSWCryptoParams> m_automorphism_params;
    FastBlindRotationKey m_br_key;

    std::vector<NativePoly> m_bit_extraction_polys;
    std::vector<AutomorphismKey> m_automorphism_keys;

    uint32_t m_pt_space = 0;
};

#endif //LARGE_FUNCTIONS_FUNCTIONAL_BOOTSTRAP_H
