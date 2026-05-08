//
// Created by leonard on 12/11/24.
//

#include "lwe_key_switching.h"
#include "cassert"

#include <utility>

LWEKeySwitchingKey::LWEKeySwitchingKey(std::shared_ptr<LWECryptoParams> &params, NativeVector sk_source,
                                       NativeVector sk_target, bool precompute_digits) {

    m_precompute_digits = precompute_digits;
    m_params = params;
    m_sk_source = std::move(sk_source);
    m_sk_target = std::move(sk_target);

    auto modulus = params->GetqKS();

    m_sk_source.SwitchModulus(modulus);
    m_sk_target.SwitchModulus(modulus);

    auto decomp_basis = params->GetBaseKS();

    uint32_t decomp_digits = 0;
    uint64_t pows = 1;
    while (pows < modulus) {
        pows *= decomp_basis;
        decomp_digits++;
    }

    double std = m_params->GetDgg().GetStd();
    m_decomp_digits = decomp_digits;

    if (precompute_digits) {
        m_key = std::vector<std::vector<std::vector<LWECiphertextImpl>>>(m_sk_source.GetLength());
        for(int i = 0; i < m_sk_source.GetLength(); i++) {
            std::vector<std::vector<LWECiphertextImpl>> tmp(decomp_digits);
            auto sk_i = m_sk_source.at(i);
            NativeInteger val_ij = sk_i;
            for(int j = 0; j < decomp_digits; j++) {
                std::vector<LWECiphertextImpl> digits(decomp_basis - 1);
                for(int k = 1; k < decomp_basis; k++) {
                    NativeInteger val_ijk = val_ij.ModMul(k, modulus);
                    auto ct = encrypt_lwe(m_sk_target, val_ijk, std);
                    digits.at(k - 1) = std::move(*ct);
                }
                val_ij.ModMulEq(decomp_basis, modulus);
                tmp.at(j) = std::move(digits);
            }
            m_key.at(i) = std::move(tmp);
        }
    } else {
        m_key = std::vector<std::vector<std::vector<LWECiphertextImpl>>>(m_sk_source.GetLength());
        for(int i = 0; i < m_sk_source.GetLength(); i++) {
            std::vector<std::vector<LWECiphertextImpl>> tmp(decomp_digits);
            auto sk_i = m_sk_source.at(i);
            NativeInteger val_ij = sk_i;
            for(int j = 0; j < decomp_digits; j++) {
                std::vector<LWECiphertextImpl> digits(1);
                auto ct = encrypt_lwe(m_sk_target, val_ij, std);
                digits.at(0) = std::move(*ct);
                val_ij.ModMulEq(decomp_basis, modulus);
                tmp.at(j) = std::move(digits);
            }
            m_key.at(i) = std::move(tmp);
        }
    }
}

LWECiphertext LWEKeySwitchingKey::SwitchKey(const NativeVector &A, const NativeInteger &B) {

    assert(A.GetLength() == m_key.size());
    assert(A.GetModulus() == m_params->GetqKS());

    auto modulus = m_params->GetqKS();
    NativeInteger result_B = B;
    NativeVector result_A = NativeVector( m_sk_target.GetLength(), m_sk_target.GetModulus(), 0);

    auto decomp_basis = m_params->GetBaseKS();

    if (m_precompute_digits) {

        for (int i = 0; i < A.GetLength(); i++) {
            auto& m_key_i = m_key.at(i);
            uint64_t a_i = A.at(i).ConvertToInt();
            for(int j = 0; j < m_decomp_digits; j++) {
                auto dig = a_i % decomp_basis;
                if (dig) {
                    auto& ct_dig = m_key_i.at(j).at(dig - 1);
                    result_A.ModSubEq(ct_dig.GetA());
                    result_B.ModSubEq(ct_dig.GetB(), modulus);
                }
                a_i /= decomp_basis;
            }
        }

    } else {

        for (int i = 0; i < A.GetLength(); i++) {
            auto& m_key_i = m_key.at(i);
            uint64_t a_i = A.at(i).ConvertToInt();
            for(int j = 0; j < m_decomp_digits; j++) {
                auto dig = a_i % decomp_basis;
                if (dig) {
                    auto& ct_dig = m_key_i.at(j).at(0);

                    result_A.ModSubEq(ct_dig.GetA().ModMul(dig));
                    result_B.ModSubEq(ct_dig.GetB().ModMul(dig, modulus), modulus);
                }
                a_i /= decomp_basis;
            }
        }
    }

    return std::make_shared<LWECiphertextImpl>(result_A, result_B);
}

LWECiphertext LWEKeySwitchingKey::SwitchKey(const lbcrypto::LWECiphertext &source) {
    return SwitchKey(source->GetA(), source->GetB());
}