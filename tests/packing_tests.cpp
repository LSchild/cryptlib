//
// Created by leonard on 6/2/26.
//

#include "gtest/gtest.h"

#include "operators/lwe_to_rlwe_packing.h"
#include "base_crypto.h"

class PackingTestGroup : public testing::Test {

protected:

    AlignedVector secret, secret_ntt;
    std::shared_ptr<LWEtoRLWEPackingContext> m_auto_params;

    void SetUp() override {
        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint32_t basis = 1 << 14;
        uint32_t digits = 4;
        double std = 0;

        m_auto_params = std::make_shared<LWEtoRLWEPackingContext>(BINARY, Q, N, basis, digits, std);
        secret = AlignedVector(N);
        secret_ntt = AlignedVector(N);

        std::srand(time(nullptr));

        for (uint32_t i = 0; i< N; i++) {
            secret[i] = rand() % 2;
        }
        secret[0] = 1;

        m_auto_params->GetNTT()->ComputeForward(secret_ntt.data(),secret.data(), 1, 1);
    }

};

TEST_F(PackingTestGroup, TestFullPacking) {

    auto N = m_auto_params->GetDimension();
    auto Q = m_auto_params->GetModulus();
    auto ntt = m_auto_params->GetNTT();

    std::vector<GenericKey> keys = {
            {"RLWE", secret.data(), secret.size()}
    };

    auto packer = m_auto_params->ConstructOperator(keys);
    auto enc = LWEEncryptor(Q, N, 0.0);

    AlignedVector lwe_samples((N + 1) * N, 0);
    for(uint32_t i = 0; i < N; i++) {
        enc.MakeLWE(lwe_samples.data() + (N + 1) * i, i, secret.data());
    }

    AlignedVector result(2 * N, 0);
    packer->Pack(result.data(), lwe_samples.data(), N);
    ntt->ComputeForward(result.data(), result.data(), 1, 1);
    intel::hexl::EltwiseMultMod(result.data(), result.data(), secret_ntt.data(), N, Q, 1);
    ntt->ComputeInverse(result.data(), result.data(), 1, 1);
    intel::hexl::EltwiseSubMod(result.data(), result.data() + N, result.data(), N, Q);

    for(uint32_t i = 0; i < N; i++) {
        EXPECT_EQ(result[i], i);
    }

}

TEST_F(PackingTestGroup, TestConsecutivePacking) {

    auto N = m_auto_params->GetDimension();
    auto Q = m_auto_params->GetModulus();
    auto ntt = m_auto_params->GetNTT();

    std::vector<GenericKey> keys = {
            {"RLWE", secret.data(), secret.size()}
    };

    auto packer = m_auto_params->ConstructOperator(keys);
    auto enc = LWEEncryptor(Q, N, 0.0);

    auto n_samples = 4 + (rand() % 60);
    std::cerr << n_samples << std::endl;
    AlignedVector lwe_samples((N + 1) * n_samples, 0);
    for(uint32_t i = 0; i < n_samples; i++) {
        enc.MakeLWE(lwe_samples.data() + (N + 1) * i, i, secret.data());
    }

    AlignedVector result(2 * N, 0);
    packer->PackConsecutively(result.data(), lwe_samples.data(), n_samples);
    ntt->ComputeForward(result.data(), result.data(), 1, 1);
    intel::hexl::EltwiseMultMod(result.data(), result.data(), secret_ntt.data(), N, Q, 1);
    ntt->ComputeInverse(result.data(), result.data(), 1, 1);
    intel::hexl::EltwiseSubMod(result.data(), result.data() + N, result.data(), N, Q);

    for(uint32_t i = 0; i < n_samples; i++) {
        EXPECT_EQ(result[i], i);
    }
    for(uint32_t i = n_samples; i < N; i++) {
        EXPECT_EQ(result[i], 0);
    }


}



