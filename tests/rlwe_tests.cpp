//
// Created by leonard on 11/12/24.
//

#include <gtest/gtest.h>

#include "rlwe_key_switching.h"
#include "automorphism_key.h"
#include "response_compression.h"

#include "utils.h"

class RLWETests : public testing::Test {

protected:

    std::shared_ptr<RingGSWCryptoParams> m_params;

    void SetUp() override {
        NativeInteger Q = 2251799813773313;
        uint32_t N = 1024;
        NativeInteger q = 2 * N;
        uint32_t base = 1 << 10;

        double std = 0;
        m_params = std::make_shared<RingGSWCryptoParams>(N,Q,q,base,base,BINFHE_METHOD::GINX, std);
    }

};

TEST_F(RLWETests, TestLWE2RLWE) {

    NativeInteger Q = 67239937;
    uint32_t N = 1 << 12;
    NativeInteger q = 2 * N;
    uint32_t base = 1 << 8;

    double std = 0;
    m_params = std::make_shared<RingGSWCryptoParams>(N,Q,q,base,base,BINFHE_METHOD::GINX, std);

    auto gen = TernaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);

    auto l2r = LWE2RLWEKey(sk_source, base);

    auto dug = DiscreteUniformGeneratorImpl<NativeVector>();
    dug.SetModulus(Q);

    NativeVector A = dug.GenerateVector(N);

    auto start = std::chrono::high_resolution_clock::now();
    auto res = l2r.SwitchKey(A);
    auto stop = std::chrono::high_resolution_clock::now();
    auto el = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();


    std::cerr << res[0] << " " << el << std::endl;

}

TEST_F(RLWETests, TestRLWEToRLWESwitch) {
    auto gen = TernaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);
    auto sk_target = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);

    auto switch_key = RLWEKeySwitchingKey(m_params, sk_source, sk_target);

    sk_source.SetFormat(EVALUATION);
    sk_target.SetFormat(EVALUATION);

    auto msg_source = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
    for (int i = 0; i < m_params->GetN(); i++) {
        msg_source.at(i) = i;
    }

    auto ct = encrypt_rlwe(m_params, msg_source, sk_source);

    auto ct_switched = switch_key.SwitchKey(ct);
    auto phase = ct_switched->GetElements()[1] - ct_switched->GetElements()[0] * sk_target;
    phase.SetFormat(COEFFICIENT);

    for(int i = 0; i < m_params->GetN(); i++) {
        EXPECT_EQ(phase.at(i), i);
    }

}

TEST_F(RLWETests, TestAutomorphism) {
    auto gen = TernaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);

    int auto_idx = 3;

    auto switch_key = AutomorphismKey(m_params, sk_source, auto_idx);

    sk_source.SetFormat(EVALUATION);

    auto msg_source = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
    for (int i = 0; i < m_params->GetN(); i++) {
        msg_source.at(i) = i;
    }

    auto expected_result = AutomorphismKey::ApplyAutomorphism(msg_source, auto_idx);

    auto ct = encrypt_rlwe(m_params, msg_source, sk_source);

    auto ct_switched = switch_key.AutomorphismTransform(ct);
    auto phase = ct_switched->GetElements()[1] - ct_switched->GetElements()[0] * sk_source;
    phase.SetFormat(COEFFICIENT);

    for(int i = 0; i < m_params->GetN(); i++) {
        EXPECT_EQ(phase.at(i), expected_result.at(i));
    }

}

TEST_F(RLWETests, TestRGSWMul) {
    auto gen = TernaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);
    sk_source.SetFormat(EVALUATION);

    auto msg_lhs = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
    auto msg_rhs = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);

    for (int i = 0; i < m_params->GetN(); i++) {
        msg_lhs.at(i) = i;
    }

    msg_rhs.at(0) = 3;

    auto ct_rgsw = encrypt_rgsw(m_params, msg_lhs, sk_source);
    auto ct_rlwe = encrypt_rlwe(m_params, msg_rhs, sk_source);

    auto prod = ct_rgsw.mul(ct_rlwe);

    auto phase = prod->GetElements()[1] - prod->GetElements()[0] * sk_source;
    phase.SetFormat(COEFFICIENT);

    for(int i = 0; i < m_params->GetN(); i++) {
        EXPECT_EQ(phase.at(i), 3 * i);
    }

}

TEST_F(RLWETests, BenchNewRGSW) {

    auto N = m_params->GetN();
    auto Q = m_params->GetQ().ConvertToInt();
    auto root = m_params->GetPolyParams()->GetRootOfUnity().ConvertToInt();
    auto basis = m_params->GetBaseG();
    auto digits = m_params->GetDigitsG();
    auto ntt_engine = std::make_shared<intel::hexl::NTT>(N,Q,root);

    std::vector<uint64_t> rlwe_buffer(4 * N, 0);

    auto gen = TernaryUniformGeneratorImpl<NativeVector>();
    auto sk_source = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);
    sk_source.SetFormat(EVALUATION);

    auto msg_lhs = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
    auto msg_rhs = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);

    for (int i = 0; i < m_params->GetN(); i++) {
        msg_lhs.at(i) = i+1;
    }

    msg_rhs.at(0) = 3;

    auto ct_rgsw = encrypt_rgsw(m_params, msg_lhs, sk_source);
    auto ct_rlwe = encrypt_rlwe(m_params, msg_rhs, sk_source);

    for(uint32_t i = 0; i < N; i++) {
        rlwe_buffer[i] = ct_rlwe->GetElements()[0].at(i).ConvertToInt();
        rlwe_buffer[i + N] = ct_rlwe->GetElements()[1].at(i).ConvertToInt();
    }

    auto ct_rgsw_new = RGSWSample(ntt_engine, basis,digits, ct_rgsw.m_msm, ct_rgsw.m_m);

    auto mul_def_start = TIC_HD;
    auto mul_res_default = ct_rgsw.mul(ct_rlwe);
    auto mul_def_stop = TIC_HD;


    auto mul_new_start = TIC_HD;

    auto mul_res_new = ct_rgsw_new.mul(ct_rlwe);
    auto mul_new_stop = TIC_HD;

    std::cerr << std::chrono::duration_cast<std::chrono::microseconds>(mul_def_stop-mul_def_start).count()
    << " "
    << std::chrono::duration_cast<std::chrono::microseconds>(mul_new_stop-mul_new_start).count()
    << std::endl;
    ;

    auto phase = mul_res_default->GetElements()[1] - mul_res_default->GetElements()[0] * sk_source;
    phase.SetFormat(COEFFICIENT);

    auto phase_new = mul_res_new->GetElements()[1] - mul_res_new->GetElements()[0] * sk_source;
    phase_new.SetFormat(COEFFICIENT);


    std::cerr << phase << std::endl;
    std::cerr << phase_new << std::endl;

}

TEST_F(RLWETests, TestCRT) {

    uint64_t P = 33554467;
    uint64_t Q = 33554473;
    uint64_t PQ = 1125902456980891;

    uint32_t N = 1 << 11;
    std::vector<uint64_t> random_data(2 * N * N * 2);
    std::srand(std::time(nullptr));
    std::generate(random_data.begin(), random_data.end(), std::rand);

    std::vector<uint64_t> mod_P(2 * N * N * 2);
    std::vector<uint64_t> mod_Q(2 * N * N * 2);

    intel::hexl::EltwiseReduceMod(random_data.data(), random_data.data(), 2 * N * N * 2, PQ, PQ, 1);
    auto NTTP = intel::hexl::NTT(N, P);
    auto NTTQ = intel::hexl::NTT(N, Q);

    auto start = std::chrono::high_resolution_clock::now();

    intel::hexl::EltwiseReduceMod(mod_P.data(), random_data.data(), 2 * N * N * 2, P, P, 1);
    intel::hexl::EltwiseReduceMod(mod_Q.data(), random_data.data(), 2 * N * N * 2, Q, Q, 1);

    for (int i = 0; i < N; i++) {
        NTTP.ComputeForward(mod_P.data() + i * 2 * N, mod_P.data() + i * 2 * N, 1, 1);
        NTTP.ComputeForward(mod_P.data() + i * 2 * N + N, mod_P.data() + i * 2 * N + N, 1, 1);
    }

    for (int i = 0; i < N; i++) {
        NTTQ.ComputeForward(mod_Q.data() + i * 2 * N, mod_Q.data() + i * 2 * N, 1, 1);
        NTTQ.ComputeForward(mod_Q.data() + i * 2 * N + N, mod_Q.data() + i * 2 * N + N, 1, 1);
    }

    auto stop = std::chrono::high_resolution_clock::now();

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << " TOOK " << elapsed << std::endl;

    EXPECT_EQ(mod_P[0], random_data[0] % P);
    EXPECT_EQ(mod_Q[0], random_data[0] % Q);

}

TEST_F(RLWETests, TestCompression) {

    NativeInteger Q_0 = 1073750017;
    NativeInteger Q_1 = 1073153;
    auto N_0 = 2048;
    auto N_1 = 512;
    auto basis = 4;
    auto digits = Q_0 == Q_1 ? 16 : 11;
    double std = 0;

    auto p_small = std::make_shared<ILNativeParams>(2 * N_1, Q_1);
    auto p_large = std::make_shared<ILNativeParams>(2 * N_0, Q_0);

    auto sk_source = NativePoly(p_large, COEFFICIENT, true);
    sk_source[0] = 1;
    sk_source[1234] = 1;
    sk_source[32] = 1;
    sk_source *= 1;

    auto sk_dest = NativePoly(p_small, COEFFICIENT, true);
    sk_dest[1] = 1;
    sk_dest[4] = 1;
    sk_dest *= 0;

    auto comp_key = CompressionKey(sk_source, sk_dest, basis, digits, std);

    auto dug = DiscreteUniformGeneratorImpl<NativeVector>();
    dug.SetModulus(Q_0);
    auto A = NativePoly(dug, p_large, EVALUATION);
    auto M = NativePoly(dug, p_large, COEFFICIENT);
    for (int i = 0; i < N_0; i++) {
        M[i] = (i & 3) * (Q_0 >> 2);
    }

    sk_source.SwitchFormat();
    M.SwitchFormat();
    M += A * sk_source;

    sk_source.SwitchFormat();

    auto v = {A, M};
    auto ct = std::make_shared<RLWECiphertextImpl>(v);
    auto ct_res = comp_key.CompressRLWE(ct);
    ct_res->SetFormat(EVALUATION);
    sk_dest.SetFormat(EVALUATION);
    auto phase = ct_res->GetElements()[1] - ct_res->GetElements()[0] * sk_dest;
    phase.SwitchFormat();
    std::cerr << phase.MultiplyAndRound(4, Q_1).Mod(4) << std::endl;
    EXPECT_EQ(phase[0].MultiplyAndRound(4, Q_1).Mod(4).ConvertToInt<uint64_t>(), 0);
}