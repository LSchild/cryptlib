//
// Created by leonard on 11/12/24.
//

#include "rlwe_key_switching.h"

RLWEKeySwitchingKey::RLWEKeySwitchingKey(std::shared_ptr<RingGSWCryptoParams> &params, NativePoly sk_source,
                                         NativePoly sk_target) {

    m_params = params;
    m_sk_source = std::move(sk_source);
    m_sk_target = std::move(sk_target);

    m_sk_source.SetFormat(EVALUATION);
    m_sk_target.SetFormat(EVALUATION);

    /* generate key-switching key */
    auto digits = m_params->GetDigitsG();
    auto base = m_params->GetBaseG();

    m_key = std::vector<RLWECiphertext>(digits);
    NativePoly msg = m_sk_source.Clone();
    for(uint32_t i = 0; i < digits; i++) {
        m_key[i] = encrypt_rlwe(m_params, msg, m_sk_target);
        msg *= base;
    }
}

RLWECiphertext RLWEKeySwitchingKey::SwitchKey(const lbcrypto::RLWECiphertext &source) const {

    auto base = m_params->GetBaseG();
    auto base_bits = GetMSB(base) - 1;
    auto digits = m_params->GetDigitsG();

    // Fake switch for debugging
    NativePoly b_out = source->GetElements()[1].Clone();
    b_out.SwitchFormat();
    NativePoly a = source->GetElements()[0].Clone();
    a.SwitchFormat();

    b_out = b_out - a * m_sk_source;
    NativePoly a_out = NativePoly(m_params->GetPolyParams(), EVALUATION, true);

    /*
    NativePoly b_out = source->GetElements()[1].Clone();
    b_out.SetFormat(EVALUATION);
    NativePoly A = source->GetElements()[0].Clone();
    A.SetFormat(COEFFICIENT);

    auto decomp_A = A.BaseDecompose(base_bits, true);
    for(uint32_t i = 0; i < digits; i++) {
        auto dig_b = m_key.at(i)->GetElements()[1] * decomp_A.at(i);
        auto dig_a = m_key.at(i)->GetElements()[0] * decomp_A.at(i);

        a_out -= dig_a;
        b_out -= dig_b;
    }
    */
    std::vector<NativePoly> res_vec = {a_out, b_out};

    return std::make_shared<RLWECiphertextImpl>(res_vec);
}

RLWECiphertext RLWEKeySwitchingKey::SwitchKey(const lbcrypto::NativePoly &A, const lbcrypto::NativePoly &B) const{
    std::vector<NativePoly> vec = {A.Clone(), B.Clone()};
    auto ct = std::make_shared<RLWECiphertextImpl>(vec);
    return SwitchKey(ct);
}

LWE2RLWEKey::LWE2RLWEKey(NativePoly &sk, uint32_t L) {

    m_basis = L;
    m_modulus = sk.GetModulus().ConvertToInt<uint64_t>();
    m_dim = sk.GetLength();
    auto pp = sk.GetParams();

    auto dug = DiscreteUniformGeneratorImpl<NativeVector>();
    dug.SetModulus(m_modulus);

    auto mod_bits = GetMSB(m_modulus) - 1;
    auto base_bits = GetMSB(m_basis) - 1;
    m_digits = (mod_bits % base_bits == 0) ? mod_bits / base_bits : mod_bits / base_bits + 1;

    auto scal_mul_bits = mod_bits + base_bits;
    std::cerr << "Scalar Mul bits " << scal_mul_bits;

    m_additions_before_reduction = 1u << (63 - mod_bits - 1);

    m_key = std::vector<el_type>(2 * m_dim * m_dim * m_digits, 0);
    auto key_ptr = m_key.data();
    // Keygen next
    sk.SetFormat(COEFFICIENT);
    auto sk_vec = sk.GetValues();
    sk.SetFormat(EVALUATION);

    for (uint32_t i = 0; i < m_dim; i++) {
        auto digit_block = key_ptr + i * 2 * m_dim * m_digits;
        auto sk_i = sk_vec.at(i);
        for (uint32_t j = 0; j < m_digits; j++) {
            auto digit_ptr = digit_block + j * 2 * m_dim;
            NativePoly msg(pp, COEFFICIENT, true);
            msg[0] = sk_i;
            msg.SwitchFormat();

            auto A = NativePoly(dug, pp, EVALUATION);
            auto B = A * sk;
            B += msg;
            for (uint32_t k = 0; k < m_dim; k++) {
                digit_ptr[k] = A[k].ConvertToInt<el_type>();
                digit_ptr[k + m_dim] = B[k].ConvertToInt<el_type>();
            }
            sk_i.ModMulEq(m_basis, m_modulus);
        }
    }
}

template<uint32_t N, typename T> void rlwe_fma(uint64_t* acc, const T* op, T scal) {
    for (uint32_t i = 0; i < 2 * N; i++) {
        acc[i] += scal * op[i];
    }
}

std::vector<uint64_t> LWE2RLWEKey::SwitchKey(NativeVector& ct) {
    auto bbits = GetMSB(m_digits) - 1;
    auto mask = (1u << bbits) - 1;
    std::vector<uint64_t> res_buffer(2 * m_dim, 0);
    for (uint32_t i = 0; i < m_dim; i++) {
        auto block = m_key.data() + i * 2 * m_dim * m_digits;
        auto a_i = ct[i].ConvertToInt<el_type>();
        for (uint32_t j = 0; j < m_digits; j++) {
            auto d_ij = block + j * 2 * m_dim;
            auto a_ij = a_i & mask;
            if (a_ij)
                rlwe_fma<1u << 12, el_type>(res_buffer.data(), d_ij, a_ij);
            a_i >>= bbits;
            //intel::hexl::EltwiseFMAMod(res_buffer.data(), d_ij, a_ij, res_buffer.data(), m_dim * 2, m_modulus, 1);
        }
    }
    return res_buffer;
}

