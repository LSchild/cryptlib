//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_RLWE_KEY_SWITCHING_H
#define LARGE_FUNCTIONS_RLWE_KEY_SWITCHING_H

#include "setup.h"
#include "utils.h"

using namespace lbcrypto;

class LWE2RLWEKey {
public:
    using el_type = uint32_t;

    LWE2RLWEKey(NativePoly& sk, uint32_t L);

    std::vector<uint64_t> SwitchKey(NativeVector& ct);

    std::vector<el_type> m_key;
    uint32_t m_basis;
    uint32_t m_digits;
    uint32_t m_dim;
    uint64_t m_modulus;
    uint32_t m_additions_before_reduction;

};


class RLWEKeySwitchingKey {

public:

    explicit RLWEKeySwitchingKey() = default;

    RLWEKeySwitchingKey(std::shared_ptr<RingGSWCryptoParams>& params, NativePoly sk_source, NativePoly sk_target);

    RLWECiphertext SwitchKey(const RLWECiphertext& source) const;

    RLWECiphertext SwitchKey(const NativePoly& A, const NativePoly& B) const;

private:

    std::shared_ptr<RingGSWCryptoParams> m_params;
    std::vector<RLWECiphertext> m_key;
    /* for debugging */
    NativePoly m_sk_source;
    NativePoly m_sk_target;
};

#endif //LARGE_FUNCTIONS_RLWE_KEY_SWITCHING_H
