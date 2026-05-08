//
// Created by leonard on 2/25/25.
//

#include "response_compression.h"

CompressionKey::CompressionKey(lbcrypto::NativePoly &sk_src, lbcrypto::NativePoly &sk_dst, uint32_t basis, uint32_t digits, double std) {
    m_basis = basis;
    m_digits = digits;
    m_std = std;

    m_params_small = std::make_shared<RingGSWCryptoParams>(sk_dst.GetLength(), sk_dst.GetModulus(), 2048, basis, basis, GINX, std);
    m_params_big = std::make_shared<RingGSWCryptoParams>(sk_src.GetLength(), sk_dst.GetModulus(), 2048, basis, basis, GINX, std);

    m_k = m_params_big->GetN() / m_params_small->GetN();

    /** key generation **/
    NativePoly sk_s(m_params_big->GetPolyParams(), COEFFICIENT, true);
    sk_src.SetFormat(COEFFICIENT);
    auto vec = sk_src.GetValues();
    vec.SwitchModulus(sk_dst.GetModulus());
    sk_s.SetValues(vec, COEFFICIENT);

    sk_dst.SetFormat(EVALUATION);

    auto dgg_small = m_params_small->GetDgg();
    auto dug_small = DiscreteUniformGeneratorImpl<NativeVector>();
    dug_small.SetModulus(m_params_small->GetQ());

    auto acc = sk_s;
    acc.SetFormat(EVALUATION);

    for (int i = 0; i < digits; i++) {
        std::vector<NativePoly> a_vec;
        std::vector<NativePoly> b_vec;

        for (int j = 0; j < m_k; j++) {
            a_vec.emplace_back(dug_small, m_params_small->GetPolyParams(), COEFFICIENT);
            b_vec.emplace_back(a_vec.back().Clone());
            auto& last = b_vec.back();
            last.SetFormat(EVALUATION);
            last *= sk_dst;
            auto err = NativePoly(dgg_small, m_params_small->GetPolyParams(), COEFFICIENT);
            err.SetFormat(EVALUATION);
            last += err;
            last.SwitchFormat();

        }

        auto a_packed = Pack(a_vec);
        a_packed *= 0;
        auto b_packed = Pack(b_vec);
        a_packed.SetFormat(EVALUATION);
        b_packed.SetFormat(EVALUATION);
        b_packed -= acc;

        acc *= basis;
        auto vv = {a_packed, b_packed};
        m_key.push_back(std::make_shared<RLWECiphertextImpl>(vv));

    }

}

lbcrypto::NativePoly CompressionKey::Kappa(lbcrypto::NativePoly &poly) {
    poly.SetFormat(COEFFICIENT);
    NativePoly result(m_params_big->GetPolyParams(), COEFFICIENT, true);
    for (int i = 0; i < poly.GetLength(); i++) {
        result[i * m_k] = poly[i];
    }

    return result;
}

NativePoly CompressionKey::KappaInv(NativePoly &poly) {
    poly.SetFormat(COEFFICIENT);
    NativePoly result(m_params_small->GetPolyParams(), COEFFICIENT, true);
    for (int i = 0; i < result.GetLength(); i++) {
        result[i] = poly[i * m_k];
    }

    return result;
}

NativePoly CompressionKey::Pack(std::vector<NativePoly> &poly) {
    NativePoly result(m_params_big->GetPolyParams(), COEFFICIENT, true);

    for (int i = 0; i < poly.size(); i++) {
        for (int j = 0; j < m_params_small->GetN(); j++) {
            result[i + j * m_k] = poly[i][j];
        }
    }

    return result;
}


lbcrypto::RLWECiphertext CompressionKey::CompressRLWE(lbcrypto::RLWECiphertext &ct) {
    ct->SetFormat(COEFFICIENT);
    auto A = ct->GetElements()[0];
    auto B = ct->GetElements()[1];

    auto q_in = A.GetModulus();
    auto q_out = m_params_small->GetQ();

    // mod switch A
    auto An = NativePoly(m_params_big->GetPolyParams(), COEFFICIENT, true);
    auto A_vec = A.GetValues();
    A_vec = A_vec.MultiplyAndRound(q_out, q_in).Mod(q_out);
    A_vec.SetModulus(q_out);

    An.SetValues(A_vec,COEFFICIENT);

    auto basebits = GetMSB(m_basis) - 1;
    auto A_digits = An.BaseDecompose(basebits, true);

    auto acc_A = m_key[0]->GetElements()[0] * A_digits[0];
    auto acc_B = NativePoly(m_params_big->GetPolyParams(), COEFFICIENT, true);
    auto b_vec = B.GetValues();

    b_vec = b_vec.MultiplyAndRound(q_out, q_in).Mod(q_out);
    b_vec.SetModulus(q_out);
    acc_B.SetValues(b_vec, COEFFICIENT);
    acc_B.SetFormat(EVALUATION);

    acc_B += m_key[0]->GetElements()[1] * A_digits[0];

    for (int i = 1; i < m_key.size(); i++) {
        acc_A += m_key[i]->GetElements()[0] * A_digits[i];
        acc_B += m_key[i]->GetElements()[1] * A_digits[i];
    }

    acc_A.SwitchFormat();
    acc_B.SwitchFormat();

    auto small_A = KappaInv(acc_A);
    auto small_B = KappaInv(acc_B);

    auto v = {small_A, small_B};
    return std::make_shared<RLWECiphertextImpl>(v);
}
