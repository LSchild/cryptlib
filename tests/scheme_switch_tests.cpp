//
// Created by leonard on 11/12/24.
//

#include <gtest/gtest.h>
#include "scheme_switch.h"
#include "rlwe-ciphertext.h"
#include "operators/cggi_blind_rotation.h"
#include "operators/bmmp_blind_rotation.h"
#include "operators/lwe_to_rgsw_conversion.h"

class CGGILWE2RGSWTests : public testing::Test {

protected:

    std::vector<GenericKey> test_keys;

    std::shared_ptr<SchemeSwitchingContext> ss_context;

    void SetUp() override {

        // declare params
        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint64_t n = 100;
        double std = 0;
        uint64_t basebits = 14;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;

        auto params_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);
        auto ntt = params_bin->GetNTT();
        auto trace_params = std::make_shared<TraceEvaluationContext>(BINARY, ntt, base, digits, std);
        auto square_params = std::make_shared<RLWEConversionContext>(BINARY, ntt, base, digits, std);

        std::shared_ptr<OperatorContext<BlindRotator>> params_rot = std::reinterpret_pointer_cast<OperatorContext<BlindRotator>>(params_bin);
        ss_context = std::make_shared<SchemeSwitchingContext>(params_rot, trace_params, square_params, digits, base);

        std::vector<uint64_t> tmp_key(std::max<uint64_t>(N, n), 0);


        // gen LWE key
        for(uint32_t i = 0; i < n; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("LWE_KEY", tmp_key.data(), n);

        // gen RLWE key
        for (uint32_t i = 0; i< N; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("RLWE_KEY", tmp_key.data(), N);

        // ntt'ed RLWE key
        ntt->ForwardNTT(tmp_key.data(), tmp_key.data());
        test_keys.emplace_back("RLWE_KEY_NTT", tmp_key.data(), N);
    }

};

TEST_F(CGGILWE2RGSWTests, TestConvExactDigs) {

    std::vector<GenericKey> test(test_keys.begin(), test_keys.begin() + 2);
    auto m_converter = ss_context->ConstructOperator(test);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(ss_context->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(ss_context->GetOutputContainer(input_shape));
    auto output_params = ss_context->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;

    auto digs = ss_context->GetOutputDigits();
    auto basis = ss_context->GetOutputBasis();

    auto lwe_key = test_keys[0].GetKey();
    auto rlwe_key = test_keys[1].GetKey();
    auto rlwe_key_ntt = test_keys[2].GetKey();


    AlignedVector rgsw_out(4 * rlwe_N * digs, 0);
    AlignedVector rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedVector expected_result(2 * rlwe_N * digs, 0);
    AlignedVector lwe(lwe_n + 1,0);

    auto ntt = output_params->GetNTT();

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_key[i]) % lwe_q;
    }

    // random bit message
    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    m_converter->Convert(rgsw_out.data(), lwe.data());

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_key_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_key.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        int64_t err = int64_t(expected_result[i]) - int64_t(rgsw_phase[i]);
        EXPECT_LE(std::abs(err), 1);
    }

}

TEST_F(CGGILWE2RGSWTests, TestConvApproxDigs) {


    ss_context->SetOutputBasis(1u << 10);
    ss_context->SetOutputDigits(4);
    auto max_digits = 6;
    auto digit_difference = max_digits - 4;


    auto m_converter = ss_context->ConstructOperator(test_keys);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(ss_context->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(ss_context->GetOutputContainer(input_shape));
    auto output_params = ss_context->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;
    auto digs = ss_context->GetOutputDigits();
    auto basis = ss_context->GetOutputBasis();

    auto rlwe_key = test_keys[1].GetKey();
    auto rlwe_key_ntt = test_keys[2].GetKey();
    auto lwe_key = test_keys[0].GetKey();

    AlignedVector rgsw_out(4 * rlwe_N * digs, 0);
    AlignedVector rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedVector expected_result(2 * rlwe_N * digs, 0);
    AlignedVector lwe(lwe_n + 1);

    auto ntt = output_params->GetNTT();

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_key[i]) % lwe_q;
    }

    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    m_converter->Convert(rgsw_out.data(), lwe.data());

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_key_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digit_difference; i++) {
        scal *= basis;
    }
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_key.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        EXPECT_EQ(expected_result[i], rgsw_phase[i]);
    }
}


class BMMPLWE2RGSWTests : public testing::Test {

protected:

    AlignedVector rlwe_secret;
    AlignedVector rlwe_secret_ntt;
    AlignedVector lwe_secret_binary;

    std::shared_ptr<SchemeSwitchingContext> m_params;

    void SetUp() override {

        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint64_t n = 100;
        double std = 0;
        uint64_t basebits = 14;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;

        auto params_bin = std::make_shared<BMMPBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std, 2);

        rlwe_secret = AlignedVector(N);
        rlwe_secret_ntt = AlignedVector(N);
        lwe_secret_binary = AlignedVector(n);

        for (uint32_t i = 0; i< N; i++) {
            rlwe_secret[i] = rand() % 2;
        }

        params_bin->GetNTT()->ForwardNTT(rlwe_secret_ntt.data(), rlwe_secret.data());
        for(uint32_t i = 0; i < n; i++) {
            lwe_secret_binary[i] = rand() % 2;
        }

        auto trace_params = std::make_shared<TraceEvaluationContext>(BINARY, params_bin->GetNTT(), base, digits, std);
        auto square_params = std::make_shared<RLWEConversionContext>(BINARY, params_bin->GetNTT(), base, digits, std);

        std::shared_ptr<OperatorContext<BlindRotator>> params_rot = std::reinterpret_pointer_cast<OperatorContext<BlindRotator>>(params_bin);

        m_params = std::make_shared<SchemeSwitchingContext>(params_rot, trace_params, square_params, digits, base);



    }

};

TEST_F(BMMPLWE2RGSWTests, TestConvExactDigs) {


    std::vector<GenericKey> keys = {
            {"LWE", lwe_secret_binary.data(), lwe_secret_binary.size()},
            {"RLWE", rlwe_secret.data(), rlwe_secret.size()}
    };

    auto m_converter = m_params->ConstructOperator(keys);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(m_params->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(m_params->GetOutputContainer(input_shape));
    auto output_params = m_params->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;
    auto digs = m_params->GetOutputDigits();
    auto basis = m_params->GetOutputBasis();

    AlignedVector lwe(lwe_n + 1);

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_secret_binary[i]) % lwe_q;
    }

    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    AlignedVector rgsw_out(4 * rlwe_N * digs, 0);
    AlignedVector rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedVector expected_result(2 * rlwe_N * digs, 0);

    m_converter->Convert(rgsw_out.data(), lwe.data());

    auto ntt = output_params->GetNTT();

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_secret_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_secret.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        int64_t err = int64_t(expected_result[i]) - int64_t(rgsw_phase[i]);
        EXPECT_LE(std::abs(err), 1);
    }

}

TEST_F(BMMPLWE2RGSWTests, TestConvApproxDigs) {


    m_params->SetOutputBasis(1u << 10);
    m_params->SetOutputDigits(4);
    auto max_digits = 6;
    auto digit_difference = max_digits - 4;
    std::vector<GenericKey> keys = {
            {"LWE", lwe_secret_binary.data(), lwe_secret_binary.size()},
            {"RLWE", rlwe_secret.data(), rlwe_secret.size()}
    };


    auto m_converter = m_params->ConstructOperator(keys);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(m_params->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(m_params->GetOutputContainer(input_shape));
    auto output_params = m_params->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;
    auto digs = m_params->GetOutputDigits();
    auto basis = m_params->GetOutputBasis();

    AlignedVector lwe(lwe_n + 1);

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_secret_binary[i]) % lwe_q;
    }

    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    AlignedVector rgsw_out(4 * rlwe_N * digs, 0);
    AlignedVector rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedVector expected_result(2 * rlwe_N * digs, 0);

    m_converter->Convert(rgsw_out.data(), lwe.data());

    auto ntt = output_params->GetNTT();

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_secret_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digit_difference; i++) {
        scal *= basis;
    }
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_secret.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        EXPECT_EQ(expected_result[i], rgsw_phase[i]);
    }
}


// OLD version
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
    NativeInteger B = bit * ss_context.m_q / 2;
    for(int i = 0; i < A.GetLength(); i++) {
        B.ModAddEq(A.at(i).ModMul(m_sk_lwe->GetElement().at(i), ss_context.m_q), ss_context.m_q);
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