//
// Created by leonard on 4/27/26.
//

#include <gtest/gtest.h>
#include <fstream>
#include "operators/cggi_blind_rotator.h"
#include "math_utils.h"

TEST(NoiseTests, TestCGGIVariance) {

    // Setup
    uint64_t Q = 36028797018972161;
    uint32_t N = 1 << 11;
    uint64_t n = 100;
    double std = 3.19;
    uint64_t basebits = 8;
    uint64_t base = 1 << basebits;
    uint64_t digits = 7;
    uint64_t n_samples = 100;

    auto lwe_encryptor = LWEEncryptor(2 * N, n, 0.0);


    AlignedVector lwe_key(n, 0);
    AlignedVector rlwe_key(N, 0);
    AlignedVector rlwe_vector(N, 0);

    AlignedVector rlwe_message(N, 0);

    AlignedVector acc_data(N, 0);
    AlignedVector rlwe_sample(2 * N, 0);
    AlignedVector lwe_sample(n, 0);


    srand(time(nullptr));


    auto ctx_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);

    auto ref_msg = random() % (2 * N);
    AlignedVector results(N * n_samples, 0.0);

    for(uint64_t i = 0; i <n_samples; i++) {

        for(uint32_t j = 0; j < n; j++) {
            lwe_key[j] = random() % 2;
        }

        for(uint32_t j = 0; j < N; j++) {
            rlwe_key[j] = random() % 2;
            acc_data[j] = j % 2 == 0 ? 1 : Q - 1;
        }

        std::vector<GenericKey> bundle;

        bundle.emplace_back(std::string("LWE_SECRET"), lwe_key.data(), lwe_key.size());
        bundle.emplace_back(std::string("RLWE_SECRET"), rlwe_key.data(), rlwe_key.size());
        auto rotator = ctx_bin->ConstructOperator(bundle);

        auto rlwe_encryptor = rotator->GetEncryptor();
        auto ntt = rlwe_encryptor->GetNTT();

        ntt->ComputeForward(rlwe_vector.data(), rlwe_key.data(), 1,1);
        std::fill(rlwe_sample.begin(), rlwe_sample.end(), 0);
        ntt->ComputeForward(rlwe_sample.data() + N, acc_data.data(), 1, 1);

        auto sample_row = results.data() + i * N;
        lwe_encryptor.MakeLWE(lwe_sample.data(), 0, lwe_key.data());
        rotator->BlindRotate(nullptr, lwe_sample.data(), rlwe_sample.data());
        rlwe_encryptor->PhaseRLWE(rlwe_message.data(), rlwe_sample.data(), rlwe_vector.data());
        intel::hexl::EltwiseSubMod(sample_row, rlwe_message.data(), acc_data.data(), N, Q);
    }


    auto ref_variance = ctx_bin->ComputeOutputVariance();
    auto sample_variance = estimate_variance(results.data(), n_samples, N, Q);

    for(auto v_i : sample_variance) {
        EXPECT_LE(v_i, 10 * ref_variance);
    }

}