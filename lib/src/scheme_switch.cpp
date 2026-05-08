//
// Created by leonard on 11/12/24.
//

#include <cassert>
#include "scheme_switch.h"

SchemeSwitchKey::SchemeSwitchKey(std::shared_ptr<RingGSWCryptoParams> &params, lbcrypto::NativePoly &sk) :
        RLWEKeySwitchingKey(params, sk * sk, sk) {

}

RLWECiphertext SchemeSwitchKey::SchemeSwitch(const lbcrypto::RLWECiphertext& source) {

    auto dummy = source->GetElements()[1].Clone();
    dummy.SetValuesToZero();

    std::vector<NativePoly> source_copy_vec = {source->GetElements()[0].Clone(), dummy};
    auto source_copy = std::make_shared<RLWECiphertextImpl>(source_copy_vec);

    auto switched = SwitchKey(source_copy);
    switched->GetElements()[0] -= source->GetElements()[1];
    switched->GetElements()[0] = switched->GetElements()[0].Negate();
    switched->GetElements()[1] = switched->GetElements()[1].Negate();
    return switched;
}

RLWECiphertext SchemeSwitchKey::SchemeSwitch(const lbcrypto::NativePoly &A, const lbcrypto::NativePoly &B) {
    NativePoly zero = NativePoly(A.GetParams(), EVALUATION, true);
    auto switched = SwitchKey(A, zero);
    switched->GetElements()[0] -= B;

    return switched;
}

SchemeSwitchEngine::SchemeSwitchEngine(SchemeSwitchParameters &params) : m_params(params) {

    /* set parameters */
    m_br_params = std::make_shared<RingGSWCryptoParams>(m_params.m_N, m_params.m_Q, m_params.m_q, m_params.m_br_basis, 1, BINFHE_METHOD::GINX, m_params.m_std);
    m_squaring_params = std::make_shared<RingGSWCryptoParams>(m_params.m_N, m_params.m_Q, m_params.m_q, m_params.m_squaring_key_basis, 1, BINFHE_METHOD::GINX, m_params.m_std);
    m_automorphism_params = std::make_shared<RingGSWCryptoParams>(m_params.m_N, m_params.m_Q, m_params.m_q, m_params.m_automorphism_key_basis, 1, BINFHE_METHOD::GINX, m_params.m_std);
    m_rgsw_params = std::make_shared<RingGSWCryptoParams>(m_params.m_N, m_params.m_Q, m_params.m_q, m_params.m_rgsw_basis, 1, BINFHE_METHOD::GINX, m_params.m_std);

}

void SchemeSwitchEngine::SetKeys(lbcrypto::NativePoly &sk_rlwe, lbcrypto::LWEPrivateKey &sk_lwe) {
    m_sk_rlwe = sk_rlwe.Clone();
    m_sk_rlwe.SetFormat(COEFFICIENT);

    m_sk_rlwe_ntt = sk_rlwe.Clone();
    m_sk_rlwe_ntt.SetFormat(EVALUATION);

    m_sk_lwe = sk_lwe->GetElement();

    assert(m_sk_lwe.GetLength() == m_params.m_n);
    assert(m_sk_rlwe.GetLength() == m_params.m_N);

    /* (re)generate internal keys */
    //m_br_key = BlindRotationKey(m_br_params, m_sk_rlwe, m_sk_lwe);
    m_br_key = FastBlindRotationKey(m_br_params, m_sk_rlwe, m_sk_lwe, 2);

    m_squaring_key = SchemeSwitchKey(m_squaring_params, m_sk_rlwe_ntt);
    m_automorphism_keys = std::vector<AutomorphismKey>();
    for(int i = 2; i <= m_params.m_N; i *= 2) {
        m_automorphism_keys.emplace_back(m_automorphism_params, m_sk_rlwe, i + 1);
    }
}

RLWECiphertext SchemeSwitchEngine::CreateDigits(lbcrypto::ConstLWECiphertext &ct) {

    /* Set up polynomials */
    auto poly_params = m_br_params->GetPolyParams();
    /* lhs of the accumulator */
    NativePoly acc_A = NativePoly(poly_params, EVALUATION, true);
    /* pad poly * digit_poly will induce acc_B */
    NativePoly pad_poly = NativePoly(poly_params, COEFFICIENT, true);
    NativePoly digit_poly = NativePoly(poly_params, COEFFICIENT, true);
    /* shift poly for b value in ct */
    NativePoly shift_poly = NativePoly(poly_params, COEFFICIENT, true);

    uint64_t full_digits = std::ceil(std::log2(m_params.m_Q) / std::log2(m_params.m_rgsw_basis));
    uint32_t total_digits = full_digits - m_params.skip_first_digits;
    // TODO: change
    //std::cerr << full_digits << " " << total_digits << std::endl;
    assert(IsPowerOf2(total_digits));

    // implicit assumption: LWE error is less than chunk_size / 2 in absolute value
    auto chunk_size = m_params.m_N / total_digits;
    NativeInteger inverse_two = NativeInteger(2).ModInverse(m_params.m_Q);
    NativeInteger base = m_rgsw_params->GetBaseG();

    for(uint32_t i = m_params.skip_first_digits; i < m_rgsw_params->GetDigitsG(); i++) {
        auto value = inverse_two.ModMul(base.ModExp(i, m_params.m_Q), m_params.m_Q);
        digit_poly.at((i - m_params.skip_first_digits) * chunk_size) = m_rgsw_params->GetQ() - value;
    }

    for(uint32_t i = 0; i < chunk_size; i++)
        pad_poly.at(i) = 1;

    // implicit assumption: as we assume that |error| is <= chunk_size / 2, we force the error to be negative
    // by subtracting chunk_size / 2
    NativeInteger B_nat = ct->GetB().ModSub(chunk_size / 2, m_params.m_q);

    auto B = B_nat.ConvertToInt<uint32_t>();
    if (B >= m_params.m_N)
        shift_poly.at(B - m_params.m_N) = m_params.m_Q-1;
    else
        shift_poly.at(B) = 1;

    shift_poly.SetFormat(EVALUATION);
    pad_poly.SetFormat(EVALUATION);
    digit_poly.SetFormat(EVALUATION);
    /* acc_B's coefficients consist of L^i (rgsw basis) for a range of i, in chunks of size \chunk_size */
    auto acc_B = (shift_poly *= digit_poly) * pad_poly;

    std::vector<NativePoly> acc_vec = {acc_A, acc_B};
    auto acc = std::make_shared<RLWECiphertextImpl>(acc_vec);

    // TODO: put - in front if using the "fast" one
    auto result = m_br_key.BlindRotate(acc, -ct->GetA());

    //auto p_r = result->GetElements()[1] - result->GetElements()[0] * m_sk_rlwe_ntt;

    //p_r.SetFormat(COEFFICIENT);
    //result->SetFormat(COEFFICIENT);

    //std::cerr << p_r << std::endl;

    // transform the choice of (-1)^bit * val/2 to bit * val in the relevant slots
    for(uint32_t i = m_params.skip_first_digits; i < full_digits; i++) {
        auto value = inverse_two.ModMul(base.ModExp(i, m_params.m_Q), m_params.m_Q);
        result->GetElements()[1].at((i - m_params.skip_first_digits) * chunk_size).ModAddEq(value, m_params.m_Q);
    }

    result->SetFormat(EVALUATION);


    return result;
}

std::vector<RLWECiphertext> SchemeSwitchEngine::ExtractDigits(lbcrypto::RLWECiphertext &ct_digits) {

    uint64_t full_digits = std::ceil(std::log2(m_params.m_Q) / std::log2(m_params.m_rgsw_basis));

    uint32_t total_digits = full_digits - m_params.skip_first_digits;
    assert(IsPowerOf2(total_digits));

    // implicit assumption: LWE error is less than chunk_size / 2 in absolute value
    auto chunk_size = m_params.m_N / total_digits;
    NativePoly shift_poly(m_squaring_params->GetPolyParams(), COEFFICIENT, true);
    shift_poly.at(m_params.m_N - chunk_size) = m_rgsw_params->GetQ()-1;
    shift_poly.SetFormat(EVALUATION);

    std::vector<RLWECiphertext> output_digits;
    for(int i = 0 ; i < total_digits; i++) {
        output_digits.push_back(GetConstantTerm(ct_digits));
        ct_digits->GetElements()[0] *= shift_poly;
        ct_digits->GetElements()[1] *= shift_poly;
    }

    return output_digits;
}

RLWECiphertext SchemeSwitchEngine::GetConstantTerm(const lbcrypto::RLWECiphertext &ct) {
    auto n_automorphism = m_automorphism_keys.size();
    NativeInteger n_auto_pow = 1 << n_automorphism;
    NativeInteger premult = n_auto_pow.ModExp(m_params.m_Q - 2, m_params.m_Q);

    auto A = ct->GetElements()[0].Clone();
    auto B = ct->GetElements()[1].Clone();

    A *= premult;
    B *= premult;

    for(int i = 0; i < n_automorphism; i++) {
        auto m_key_i = m_automorphism_keys.at(n_automorphism - i - 1);
        auto rlwe_auto_i = m_key_i.AutomorphismTransform(A, B);
        A += rlwe_auto_i->GetElements()[0];
        B += rlwe_auto_i->GetElements()[1];
    }

    std::vector<NativePoly> rlwe_vec = {A, B};
    return std::make_shared<RLWECiphertextImpl>(rlwe_vec);
}

RingGSWSample SchemeSwitchEngine::SwitchToRGSW(lbcrypto::ConstLWECiphertext &ct) {

    auto digits_packed = CreateDigits(ct);
    auto digits_extracted = ExtractDigits(digits_packed);

    /* build lhs of rgsw sample */
    std::vector<RLWECiphertext> digits_mul_sk;

    for(auto& ct_dig : digits_extracted) {
        digits_mul_sk.push_back(m_squaring_key.SchemeSwitch(ct_dig));
    }

    return {m_rgsw_params, digits_mul_sk, digits_extracted};
}