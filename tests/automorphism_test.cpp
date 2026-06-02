//
// Created by leonard on 5/19/26.
//

#include "operators/automorphism_evaluation.h"
#include "base_crypto.h"

#include <gtest/gtest.h>

class AutoBaseTestGroup : public testing::Test {

protected:

    std::shared_ptr<AutomorphismEvaluator> eval_auto;

    AlignedVector secret_nonoise, secret_ntt;
    std::shared_ptr<AutomorphismContext> m_auto_params;


    void SetUp() override {
        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint32_t basis = 1 << 14;
        uint32_t digits = 4;
        double std = 0;

        m_auto_params = std::make_shared<AutomorphismContext>(BINARY, Q, N, basis, digits, std, 1);
        secret_nonoise = AlignedVector(N);
        secret_ntt = AlignedVector(N);

        std::srand(4325243);

        for (uint32_t i = 0; i< N; i++) {
            secret_nonoise[i] = rand() % 2;
        }

        m_auto_params->GetNTT()->ComputeForward(secret_ntt.data(),secret_nonoise.data(), 1, 1);

    }

};

TEST_F(AutoBaseTestGroup, TestAutomorphism) {

    m_auto_params->SetAutomorphismIndex(257);

    std::vector<GenericKey> key = {
            {"RLWE", secret_nonoise.data(), secret_nonoise.size()}
    };

    auto eval = m_auto_params->ConstructOperator(key);

    auto N = m_auto_params->GetDimension();
    auto e_ntt = m_auto_params->GetNTT();
    auto encryptor = RLWEEncryptor(e_ntt, m_auto_params->GetStd());

    AlignedVector msg(N, 0);
    AlignedVector ct(2 * N, 0);
    AlignedVector ct_out(2 * N, 0);

    msg[1] = m_auto_params->GetModulus() / 2;
    //msg[769] = 1;

    encryptor.MakeRLWE(ct.data(),msg.data(), secret_ntt.data(), false);
    e_ntt->ComputeInverse(ct.data(),ct.data(),1,1);
    e_ntt->ComputeInverse(ct.data() + N,ct.data() + N,1,1);

    eval->Eval(ct_out.data(), ct.data());

    std::fill(ct.begin(),ct.end(), 0);
    encryptor.PhaseRLWE(ct.data(), ct_out.data(), secret_ntt.data());

    for(uint32_t i = 0; i < N; i++) {
        EXPECT_EQ(ct[i], i == 257 ? m_auto_params->GetModulus() / 2 : 0);
    }
    std::cerr << std::endl;
}
