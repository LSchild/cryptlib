//
// Created by leonard on 12/11/24.
//

#ifndef LARGE_FUNCTIONS_LWE_KEY_SWITCHING_H
#define LARGE_FUNCTIONS_LWE_KEY_SWITCHING_H

#include "setup.h"
#include "utils.h"

using namespace lbcrypto;

class LWEKeySwitchingKey {

public:

    explicit LWEKeySwitchingKey() = default;

    LWEKeySwitchingKey(std::shared_ptr<LWECryptoParams>& params, NativeVector sk_source, NativeVector sk_target, bool precompute_digits = false);

    LWECiphertext SwitchKey(const LWECiphertext& source);

    LWECiphertext SwitchKey(const NativeVector & A, const NativeInteger & B);

    std::shared_ptr<LWECryptoParams> m_params;
    std::vector<std::vector<std::vector<LWECiphertextImpl>>> m_key;
    /* for debugging */
    NativeVector m_sk_source;
    NativeVector m_sk_target;

    bool m_precompute_digits;
    uint32_t m_decomp_digits;
};

#endif //LARGE_FUNCTIONS_LWE_KEY_SWITCHING_H
