//
// Created by leonard on 12/11/24.
//

#include <gtest/gtest.h>
#include "openfhe.h"
#include "utils.h"
#include "lwe_key_switching.h"

class LweKeySwitchTests : public testing::Test {

protected:

    NativeVector m_sk_source, m_sk_target;
    std::shared_ptr<lbcrypto::LWECryptoParams> m_params;

    void SetUp() override {

        // Sue me
        srand(time(nullptr));
        uint64_t Q = 4294991873;
        uint32_t Nbits = 10;
        uint32_t N = 1 << Nbits;
        uint64_t q = 2 * N;
        uint32_t n = 50;
        uint32_t L = 1 << 4;
        // standard deviation to zero so we can test for equality
        double std = 0;

        m_params = std::make_shared<lbcrypto::LWECryptoParams>(n, N, q, Q, Q, std, L);

        auto genT = lbcrypto::TernaryUniformGeneratorImpl<NativeVector>();

        m_sk_source = genT.GenerateVector(n, Q);
        m_sk_target = genT.GenerateVector(N, Q);
    }

};

TEST_F(LweKeySwitchTests, TestPrecompute) {

    auto key = LWEKeySwitchingKey(m_params, m_sk_source, m_sk_target, true);

    auto val = random() % m_params->GetqKS();

    auto ct_test = encrypt_lwe(m_sk_source, val);
    auto ct_switched = key.SwitchKey(ct_test);

    auto phase = ct_switched->GetB().ModSub(DotProduct(m_sk_target, ct_switched->GetA()), m_params->GetqKS());

    EXPECT_EQ(phase, val);

}


TEST_F(LweKeySwitchTests, TestNoPrecompute) {

    auto key = LWEKeySwitchingKey(m_params, m_sk_source, m_sk_target);

    auto val = random() % m_params->GetqKS();

    auto ct_test = encrypt_lwe(m_sk_source, val);
    auto ct_switched = key.SwitchKey(ct_test);

    auto phase = ct_switched->GetB().ModSub(DotProduct(m_sk_target, ct_switched->GetA()), m_params->GetqKS());

    EXPECT_EQ(phase, val);

}
