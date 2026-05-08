//
// Created by leonard on 11/13/24.
//

#include "blind_rotator.h"
#include <gtest/gtest.h>
#include <chrono>

#include "utils.h"

class BlindRotationTests : public testing::Test {

protected:

    std::shared_ptr<RingGSWCryptoParams> m_params;
    std::shared_ptr<RingGSWCryptoParams> m_switch_params;
    std::shared_ptr<LWEPrivateKeyImpl> m_sk_lwe;
    NativePoly m_sk, m_sk_ntt;

    void SetUp() override {
        NativeInteger Q = 36028797018972161;
        uint32_t N = 1 << 11;
        NativeInteger q = 2 * N;
        uint32_t base = 1 << 10;
        uint32_t n = 600;

        double std = 0;
        m_params = std::make_shared<RingGSWCryptoParams>(N, Q, q,base,base,BINFHE_METHOD::GINX, std);

        auto gen = BinaryUniformGeneratorImpl<NativeVector>();
        auto sk_source = gen.GenerateVector(n, q);
        sk_source[0] = 1; sk_source[1] = 1;
        m_sk_lwe = std::make_shared<LWEPrivateKeyImpl>(sk_source);
        m_sk = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);
        m_sk_ntt = m_sk.Clone();
        m_sk_ntt.SetFormat(EVALUATION);
    }
};

TEST_F(BlindRotationTests, TestBlindRotateNoiseless) {

    auto blind_rotator = BlindRotationKey(m_params, m_sk, m_sk_lwe->GetElement());
    auto acc_A = NativePoly(m_params->GetPolyParams(), EVALUATION, true);
    auto acc_B = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
    auto shift_poly = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);

    for(int i = 0; i < acc_B.GetLength(); i++) {
        acc_B.at(i) = i;
    }

    acc_B.SetFormat(EVALUATION);

    auto gen = DiscreteUniformGeneratorImpl<NativeVector>();
    auto A = gen.GenerateVector(m_sk_lwe->GetLength(), m_params->Getq());

    srand(time(nullptr));
    auto shift = random() % m_params->Getq();

    NativeInteger B = (m_params->Getq() - shift) % m_params->Getq();

    for(int i =0; i < A.GetLength(); i++) {
        B.ModAddEq(A.at(i).ModMul(m_sk_lwe->GetElement().at(i), m_params->Getq()), m_params->Getq());
    }

    auto val = B >= m_params->GetN() ? m_params->GetQ() -1 : 1;
    auto b_idx = B.ConvertToInt() % m_params->GetN();
    shift_poly.at(b_idx) = val;
    shift_poly.SetFormat(EVALUATION);

    acc_A *= shift_poly;
    acc_B *= shift_poly;

    std::vector<NativePoly> acc_vec = {acc_A, acc_B};
    auto acc = std::make_shared<RLWECiphertextImpl>(acc_vec);
    auto ct_lwe = std::make_shared<LWECiphertextImpl>(A, B);

    auto start = std::chrono::high_resolution_clock::now();
    auto br_out = blind_rotator.BlindRotate(acc, A);

    auto stop = std::chrono::high_resolution_clock::now();
    auto ell = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    std::cerr << ell << std::endl << std::flush;
    auto phase2 = DotProduct(m_sk_lwe->GetElement(), A);

    std::cerr << "PHASE = " << phase2 << std::endl;
    auto phase = br_out->GetElements()[1] - br_out->GetElements()[0] * m_sk_ntt;
    phase.SetFormat(COEFFICIENT);

    std::cerr << phase << std::endl;
    auto p0 = phase.at(0).ConvertToInt();
    if (shift >= m_params->GetN()) {
        shift = shift - m_params->GetN();
        EXPECT_EQ(p0, m_params->GetQ().ConvertToInt() - shift);
    } else {
        EXPECT_EQ(p0, shift);
    }

}

class FastBlindRotationKeyTests : public testing::Test {
public:
    std::shared_ptr<RingGSWCryptoParams> m_params;
    std::shared_ptr<LWEPrivateKeyImpl> m_sk_lwe;
    NativePoly m_sk, m_sk_ntt;

    void SetUp() override {
        uint64_t Q = 72057594037641217;
        uint32_t N = 1<<11;
        uint64_t q = 2 * N;
        uint32_t n = 600;
        uint32_t L = 1 << 10;

        double std = 0;
        m_params = std::make_shared<RingGSWCryptoParams>(N, Q, q,L,L,BINFHE_METHOD::GINX, std);
        auto gen = BinaryUniformGeneratorImpl<NativeVector>();
        auto sk_source = gen.GenerateVector(n, q);

        m_sk_lwe = std::make_shared<LWEPrivateKeyImpl>(sk_source);
        m_sk = NativePoly(gen, m_params->GetPolyParams(), COEFFICIENT);
        m_sk_ntt = m_sk.Clone();
        m_sk_ntt.SetFormat(EVALUATION);
    }

};

TEST_F(FastBlindRotationKeyTests, TestStep1) {

    auto brk = FastBlindRotationKey(m_params, m_sk, m_sk_lwe->GetElement(), 2);
    auto acc_A = NativePoly(m_params->GetPolyParams(), EVALUATION, true);
    auto acc_B = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
    auto shift_poly = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);

    for(int i = 0; i < acc_B.GetLength(); i++) {
        acc_B.at(i) = i;
    }

    acc_B.SetFormat(EVALUATION);

    std::srand(time(nullptr));

    auto gen = DiscreteUniformGeneratorImpl<NativeVector>();
    auto A = gen.GenerateVector(m_sk_lwe->GetLength(), m_params->Getq());

    auto shift = std::rand() % m_params->Getq();

    NativeInteger B = (m_params->Getq() - shift) % m_params->Getq();

    for(int i =0; i < A.GetLength(); i++) {
        B.ModAddEq(A.at(i).ModMul(m_sk_lwe->GetElement().at(i), m_params->Getq()), m_params->Getq());
    }

    auto val = B >= m_params->GetN() ? m_params->GetQ() -1 : 1;
    auto b_idx = B.ConvertToInt() % m_params->GetN();
    shift_poly.at(b_idx) = val;
    shift_poly.SetFormat(EVALUATION);

    acc_A *= shift_poly;
    acc_B *= shift_poly;

    std::vector<NativePoly> acc_vec = {acc_A, acc_B};
    auto acc = std::make_shared<RLWECiphertextImpl>(acc_vec);
    auto ct_lwe = std::make_shared<LWECiphertextImpl>(A, B);
    volatile auto n_rounds = 10;

    std::vector<RLWECiphertext> res;
    auto start = std::chrono::high_resolution_clock::now();

    res.push_back(brk.BlindRotate(acc, -A));

    auto stop = std::chrono::high_resolution_clock::now();
    auto ell = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();

    std::cerr << "TIME = ell" << ell << std::endl << std::flush;

    auto br_out = res.back();

    auto phase = br_out->GetElements()[1] - br_out->GetElements()[0] * m_sk_ntt;
    phase.SetFormat(COEFFICIENT);

    std::cerr << phase << std::endl;
    auto p0 = phase.at(0).ConvertToInt();
    if (shift >= m_params->GetN()) {
        shift = shift - m_params->GetN();
        EXPECT_EQ(p0, m_params->GetQ().ConvertToInt() - shift);
    } else {
        EXPECT_EQ(p0, shift);
    }

}