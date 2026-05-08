//
// Created by leonard on 11/12/24.
//

#include <gtest/gtest.h>
#include "scheme_switch.h"
#include "rlwe-ciphertext.h"

class SchemeSwitchTests : public testing::Test {

protected:

    SchemeSwitchParameters m_params{};

    void SetUp() override {
        srand(123456);
        uint64_t Q = 72057594037641217;
        uint32_t N = 1<<11;
        uint64_t q = 2 * N;
        uint32_t n = 600;
        uint32_t L = 1 << 10;
        uint32_t L_out = 1 << 4;
        uint32_t L_auto_square = 1 << 2;
        uint32_t skip_digits = 6;
        // noise set to 0 for exact comparison
        double std = 0;
        SchemeSwitchParameters params{n, N, q, Q, std, skip_digits, L, L_auto_square, L_auto_square, L_out};
        m_params = params;

    }

};

TEST_F(SchemeSwitchTests, TestSwitchAndMux) {
    auto scheme_switch_engine = SchemeSwitchEngine(m_params);

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = genT.GenerateVector(m_params.m_n, m_params.m_q);


    auto m_sk_lwe = std::make_shared<LWEPrivateKeyImpl>(sk_source);
    auto m_sk = NativePoly(genT, scheme_switch_engine.m_br_params->GetPolyParams(), COEFFICIENT);
    auto m_sk_ntt = m_sk.Clone();
    m_sk_ntt.SetFormat(EVALUATION);

    scheme_switch_engine.SetKeys(m_sk, m_sk_lwe);

    auto gen = DiscreteUniformGeneratorImpl<NativeVector>();
    auto A = gen.GenerateVector(m_sk_lwe->GetLength(), m_params.m_q);

    // Sue me
    auto bit = 1;//random() & 1;
    NativeInteger B = bit * m_params.m_q / 2;
    for(int i = 0; i < A.GetLength(); i++) {
        B.ModAddEq(A.at(i).ModMul(m_sk_lwe->GetElement().at(i), m_params.m_q), m_params.m_q);
    }
    auto ct = std::make_shared<LWECiphertextImpl>(A, B);

    auto switched = scheme_switch_engine.SwitchToRGSW(ct);

    auto pos = NativePoly(scheme_switch_engine.m_rgsw_params->GetPolyParams(), COEFFICIENT, true);
    auto zero = pos.Clone();

    pos[0] = pos.GetModulus() / 2;
    for (int i = 1; i < 1024; i++)
        pos[i] = i;
    auto vals = {zero, pos};
    auto check = std::make_shared<RLWECiphertextImpl>(vals);

    auto prod = switched.mul(check);

    auto dec = prod->GetElements()[1] - prod->GetElements()[0] * m_sk_ntt;
    dec.SetFormat(COEFFICIENT);

    auto bit_out = dec[0].MultiplyAndRound(2, dec.GetModulus()).Mod(2);

    //std::cerr << dec << std::endl;
    //std::cerr << bit << std::endl;
    //std::cerr << bit_out << std::endl;

    EXPECT_EQ(bit, bit_out.ConvertToInt<long>());
}

TEST_F(SchemeSwitchTests, TestCreateDigits) {

    auto scheme_switch_engine = SchemeSwitchEngine(m_params);

    auto genT = TernaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = genT.GenerateVector(m_params.m_n, m_params.m_q);

    auto m_sk_lwe = std::make_shared<LWEPrivateKeyImpl>(sk_source);
    auto m_sk = NativePoly(genT, scheme_switch_engine.m_br_params->GetPolyParams(), COEFFICIENT);
    auto m_sk_ntt = m_sk.Clone();
    m_sk_ntt.SetFormat(EVALUATION);

    scheme_switch_engine.SetKeys(m_sk, m_sk_lwe);

    auto gen = DiscreteUniformGeneratorImpl<NativeVector>();
    auto A = gen.GenerateVector(m_sk_lwe->GetLength(), m_params.m_q);

    srand(time(nullptr));
    auto bit = random() & 1;
    NativeInteger B = bit * m_params.m_q / 2;
    for(int i = 0; i < A.GetLength(); i++) {
        B.ModAddEq(A.at(i).ModMul(m_sk_lwe->GetElement().at(i), m_params.m_q), m_params.m_q);
    }
    auto ct = std::make_shared<LWECiphertextImpl>(A, B);

    auto ct_digits = scheme_switch_engine.CreateDigits(ct);
    auto digit_phase = ct_digits->GetElements()[1] - ct_digits->GetElements()[0] * m_sk_ntt;
    digit_phase.SetFormat(COEFFICIENT);

    std::cerr << digit_phase << std::endl << std::flush;

    EXPECT_EQ(digit_phase.at(0).ConvertToInt(), bit == 1 ? (1 << 20) : 0);
}

TEST_F(SchemeSwitchTests, TestExtractDigits) {
    auto scheme_switch_engine = SchemeSwitchEngine(m_params);

    auto genT = TernaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = genT.GenerateVector(m_params.m_n, m_params.m_q);

    auto m_sk_lwe = std::make_shared<LWEPrivateKeyImpl>(sk_source);
    auto m_sk = NativePoly(genT, scheme_switch_engine.m_br_params->GetPolyParams(), COEFFICIENT);
    auto m_sk_ntt = m_sk.Clone();
    m_sk_ntt.SetFormat(EVALUATION);

    scheme_switch_engine.SetKeys(m_sk, m_sk_lwe);

    auto gen = DiscreteUniformGeneratorImpl<NativeVector>();
    auto A = gen.GenerateVector(m_sk_lwe->GetLength(), m_params.m_q);

    auto bit = random() & 1;
    NativeInteger B = bit * m_params.m_q / 2;
    for(int i = 0; i < A.GetLength(); i++) {
        B.ModAddEq(A.at(i).ModMul(m_sk_lwe->GetElement().at(i), m_params.m_q), m_params.m_q);
    }
    auto ct = std::make_shared<LWECiphertextImpl>(A, B);

    auto ct_digits = scheme_switch_engine.CreateDigits(ct);

    auto extracted = scheme_switch_engine.ExtractDigits(ct_digits);
    NativeInteger v0 = bit == 1 ? 1 << 20 : 0;
    (void) v0;

    for(auto& ct_i : extracted) {
        auto phase_i = ct_i->GetElements()[1] - ct_i->GetElements()[0] * m_sk_ntt;
        // If the implementation is correct, we don't actually need to switch to COEFFICIENT mode
        phase_i.SetFormat(COEFFICIENT);
        EXPECT_EQ(phase_i.at(0).ConvertToInt<uint64_t>(), v0.ConvertToInt<uint64_t>());
        v0.ModMulEq(m_params.m_rgsw_basis, m_params.m_Q);
    }
}

TEST_F(SchemeSwitchTests, TestSchemeSwitch) {
    auto scheme_switch_engine = SchemeSwitchEngine(m_params);

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = genT.GenerateVector(m_params.m_n, m_params.m_q);

    auto m_sk_lwe = std::make_shared<LWEPrivateKeyImpl>(sk_source);
    auto m_sk = NativePoly(genT, scheme_switch_engine.m_br_params->GetPolyParams(), COEFFICIENT);
    auto m_sk_ntt = m_sk.Clone();
    m_sk_ntt.SetFormat(EVALUATION);

    scheme_switch_engine.SetKeys(m_sk, m_sk_lwe);

    auto gen = DiscreteUniformGeneratorImpl<NativeVector>();
    auto A = gen.GenerateVector(m_sk_lwe->GetLength(), m_params.m_q);

    auto bit = random() & 1;
/*
    NativeInteger B = bit * m_params.m_q / 2;
    for(int i = 0; i < A.GetLength(); i++) {
        B.ModAddEq(A.at(i).ModMul(m_sk_lwe->GetElement().at(i), m_params.m_q), m_params.m_q);
    }
    auto ct = std::make_shared<LWECiphertextImpl>(A, B); */
    NativeInteger msg = bit * m_params.m_q / 2;
    auto ct = encrypt_lwe(sk_source, msg);

    auto ct_digits = scheme_switch_engine.SwitchToRGSW(ct);

    NativeInteger v0 = bit == 1 ? 1 << 24 : 0;

    for(auto& ct_i : ct_digits.m_m) {
        auto phase_i = ct_i->GetElements()[1] - ct_i->GetElements()[0] * m_sk_ntt;
        phase_i.SetFormat(COEFFICIENT);
        //std::cerr << phase_i << std::endl;
        EXPECT_EQ(phase_i.at(0).ConvertToInt<uint64_t>(), v0.ConvertToInt<uint64_t>());
        v0.ModMulEq(m_params.m_rgsw_basis, m_params.m_Q);
    }

    v0 = bit == 1 ? 1 << 24 : 0;

    for(auto& ct_i : ct_digits.m_msm) {
        auto phase_i = ct_i->GetElements()[1] - ct_i->GetElements()[0] * m_sk_ntt;
        auto current_expected_msm = v0 * m_sk_ntt;

        phase_i.SetFormat(COEFFICIENT);
        current_expected_msm.SetFormat(COEFFICIENT);

        EXPECT_TRUE(-phase_i == current_expected_msm);
        v0.ModMulEq(m_params.m_rgsw_basis, m_params.m_Q);
    }
}