//
// Created by leonard on 2/25/25.
//

#ifndef RESPONSE_COMPRESSION_H
#define RESPONSE_COMPRESSION_H

#include <cstdint>
#include "setup.h"
#include <vector>

using namespace lbcrypto;

class CompressionKey {
public:
    CompressionKey(lbcrypto::NativePoly& sk_src, lbcrypto::NativePoly& sk_dst, uint32_t basis, uint32_t digits, double std);

    lbcrypto::RLWECiphertext CompressRLWE(lbcrypto::RLWECiphertext& ct);

    lbcrypto::NativePoly Kappa(lbcrypto::NativePoly& poly);

    NativePoly KappaInv(NativePoly& poly);

    NativePoly Pack(std::vector<NativePoly> &poly);

    std::vector<lbcrypto::RLWECiphertext> m_key;

    std::shared_ptr<RingGSWCryptoParams> m_params_big;
    std::shared_ptr<RingGSWCryptoParams> m_params_small;

    uint32_t m_k;
    uint32_t m_basis;
    uint32_t m_digits;
    double m_std;

};

#endif //RESPONSE_COMPRESSION_H
