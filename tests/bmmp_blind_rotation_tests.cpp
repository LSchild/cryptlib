//
// Created by leonard on 4/20/26.
//
#include <gtest/gtest.h>
#include "base_crypto.h"
#include "mux_operator.h"
#include "bmmp_blind_rotator.h"

class BMMPBlindRotTestGroup : public testing::Test {

protected:

    std::shared_ptr<BMMPBlindRotationContext> ctx;
    std::shared_ptr<BMMPBlindRotator> rotatorBinary;

    AlignedVector rlwe_secret;
    AlignedVector rlwe_secret_ntt;
    AlignedVector lwe_secret_binary;
    AlignedVector lwe_secret_ternary;

    void SetUp() override {

        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint64_t n = 1024;
        double std = 0;
        uint64_t basebits = 14;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;
        uint64_t m_step_size = 2;

        ctx = std::make_shared<BMMPBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std, m_step_size);

        rlwe_secret = AlignedVector(N);
        rlwe_secret_ntt = AlignedVector(N);
        lwe_secret_binary = AlignedVector(n);

        std::srand(time(nullptr));

        for (uint32_t i = 0; i< N; i++) {
            rlwe_secret[i] = rand() % Q;
        }

        for(uint32_t i = 0; i < n; i++) {
            lwe_secret_binary[i] = rand() % 2;
        }

        auto bundle = BlindRotationKeys {lwe_secret_binary.data(), rlwe_secret.data()};
        rotatorBinary = ctx->ConstructOperator(bundle);

    }

};

TEST_F(BMMPBlindRotTestGroup, TestBinaryBlindRotate) {
    auto encryptor =  rotatorBinary->GetEncryptor();
    auto ntt= encryptor->GetNTT();
    auto Q = ctx->GetModulus();
    auto N = ctx->GetRingDimension();
    auto n = ctx->GetLWEDimension();

    ntt->ComputeForward(rlwe_secret_ntt.data(), rlwe_secret.data(), 1, 1);

    AlignedVector data(2 * N);
    AlignedVector phase(N);
    AlignedVector lwe(ctx->GetLWEDimension() + 1);

    // create bogus lwe sample
    uint64_t msg = N + (rand() % (N));
    lwe[n] = msg;
    for(uint64_t  i = 0; i < n; i++) {
        lwe[i] = rand() % (2 * N);
        lwe[n] += lwe[i] * lwe_secret_binary[i];
    }
    lwe[n] %= (2 * N);

    // Set up acc message
    for(uint64_t  i = 0; i < N; i++) {
        data[i] = i;
        data[i + N] = (Q - i) % Q;
    }

    // Set up acc
    AlignedVector rlwe_acc(2 * N);
    encryptor->MakeRLWE(rlwe_acc.data(), data.data(), rlwe_secret_ntt.data(), false);

    auto start = std::chrono::high_resolution_clock::now();
    rotatorBinary->BlindRotate(lwe.data(), rlwe_acc.data());
    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << "Took " << elapsed << std::endl;

    encryptor->PhaseRLWE(phase.data(), rlwe_acc.data(), rlwe_secret_ntt.data());

    auto start_idx = (2 * N - msg) % (2 * N);
    for(uint64_t i = 0; i < N; i++) {
        EXPECT_EQ(phase[i], data[(i + start_idx) % (2 * N)]);
    }

}
