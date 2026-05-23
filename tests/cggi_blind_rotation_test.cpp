//
// Created by leonard on 4/8/26.
//

#include <gtest/gtest.h>
#include "base_crypto.h"
#include "mux_operator.h"
#include "cggi_blind_rotator.h"

class CGGIBlindRotTestGroup : public testing::Test {

protected:

    std::shared_ptr<CGGIBlindRotator> rotatorBinary;
    std::shared_ptr<CGGIBlindRotator> rotatorTernary;

    AlignedVector rlwe_secret;
    AlignedVector rlwe_secret_ntt;
    AlignedVector lwe_secret_binary;
    AlignedVector lwe_secret_ternary;

    void SetUp() override {

        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint64_t n = 100;
        double std = 0;
        uint64_t basebits = 14;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;

        auto params_bin = CGGIBlindRotationContext(KeyDistribution::BINARY, Q, N, n, base, digits, std);
        auto params_ter = CGGIBlindRotationContext(KeyDistribution::TERNARY, Q, N, n, base, digits, std);

        rlwe_secret = AlignedVector(N);
        rlwe_secret_ntt = AlignedVector(N);
        lwe_secret_binary = AlignedVector(n);
        lwe_secret_ternary = AlignedVector(n);

        std::srand(time(nullptr));

        for (uint32_t i = 0; i< N; i++) {
            rlwe_secret[i] = rand() % Q;
        }

        for(uint32_t i = 0; i < n; i++) {
            lwe_secret_binary[i] = rand() % 2;
            lwe_secret_ternary[i] = (2 * N + (rand() % 3) - 1) % (2 * N);
        }

        rotatorBinary = std::make_shared<CGGIBlindRotator>(params_bin);
        rotatorTernary = std::make_shared<CGGIBlindRotator>(params_ter);

        rotatorBinary->KeyGen(lwe_secret_binary.data(), rlwe_secret.data());
        rotatorTernary->KeyGen(lwe_secret_ternary.data(), rlwe_secret.data());

    }

};

TEST_F(CGGIBlindRotTestGroup, TestBinaryBlindRotate) {
    auto& params = dynamic_cast<const CGGIBlindRotationContext&>(rotatorBinary->GetParams());
    auto encryptor =  rotatorBinary->GetEncryptor();
    auto ntt= encryptor->GetNTT();
    auto Q = params.GetModulus();
    auto N = params.GetRingDimension();
    auto n = params.GetLWEDimension();

    ntt->ComputeForward(rlwe_secret_ntt.data(), rlwe_secret.data(), 1, 1);

    AlignedVector data(2 * N);
    AlignedVector phase(N);
    AlignedVector lwe(params.GetLWEDimension() + 1);

    // create bogus lwe sample
    auto lwe_enc = LWEEncryptor(2 * N, n, 0.0);
    uint64_t msg = N + (rand() % (N));
    lwe_enc.MakeLWE(lwe.data(), msg, lwe_secret_binary.data());
    /*
    lwe[n] = msg;
    for(uint64_t  i = 0; i < n; i++) {
        lwe[i] = rand() % (2 * N);
        lwe[n] += lwe[i] * lwe_secret_binary[i];
    }
    lwe[n] %= (2 * N); */

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

TEST_F(CGGIBlindRotTestGroup, TestTernaryBlindRotate) {
    auto& params = dynamic_cast<const CGGIBlindRotationContext&>(rotatorTernary->GetParams());
    auto encryptor =  rotatorTernary->GetEncryptor();
    auto ntt= encryptor->GetNTT();
    auto Q = params.GetModulus();
    auto N = params.GetRingDimension();
    auto n = params.GetLWEDimension();

    ntt->ComputeForward(rlwe_secret_ntt.data(), rlwe_secret.data(), 1, 1);

    AlignedVector data(2 * N);
    AlignedVector phase(N);
    AlignedVector lwe(params.GetLWEDimension() + 1);

    // create bogus lwe sample
    uint64_t msg = N + (rand() % (N));
    lwe[n] = msg;
    for(uint64_t  i = 0; i < n; i++) {
        lwe[i] = rand() % (2 * N);
        lwe[n] += (lwe[i] * lwe_secret_ternary[i]) % (2 * N);
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
    rotatorTernary->BlindRotate(lwe.data(), rlwe_acc.data());
    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << "Took " << elapsed << std::endl;

    encryptor->PhaseRLWE(phase.data(), rlwe_acc.data(), rlwe_secret_ntt.data());

    auto start_idx = (2 * N - msg) % (2 * N);
    for(uint64_t i = 0; i < N; i++) {
        EXPECT_EQ(phase[i], data[(i + start_idx) % (2 * N)]);
    }
}