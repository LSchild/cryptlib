//
// Created by leonard on 4/27/26.
//

#include <gtest/gtest.h>
#include <fstream>
#include <random>
#include "operators/cggi_blind_rotation.h"
#include "utils/math_utils.h"

TEST(NoiseTests, TestCGGIVarianceExact) {

    uint64_t rounds = 1000;

    // Setup
    uint64_t Q = 36028797018972161;
    uint32_t N = 1 << 10;
    uint64_t n = 512;
    double std = 3.19;
    uint64_t basis_bits = 10;
    uint64_t base = 1 << basis_bits;
    uint64_t digits = 6;

    auto lwe_encryptor = LWEEncryptor(2 * N, n, 0.0);

    AlignedVector lwe_key(n, 0);
    AlignedVector rlwe_key(N, 0);
    std::vector<uint64_t> rlwe_key_ntt(N, 0);
    AlignedVector phase(N);

    AlignedVector lwe_sample(n + 1, 0);
    AlignedVector rlwe_sample(2 * N, 0);
    AlignedVector rlwe_sample_out(2 * N, 0);

    AlignedVector zero(N, 0);

    /* randomness */
    std::random_device random_device;
    std::mt19937 engine(random_device());
    std::uniform_int_distribution<uint64_t> unif_1_0(0, 1);


    auto ctx_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);

    AlignedVector output_noise(N * rounds, 0.0);

    for(uint64_t i = 0; i < rounds; i++) {

        std::fill(rlwe_sample.begin(), rlwe_sample.end(), 0);

        for(uint32_t j = 0; j < n; j++) {
            lwe_key[j] = unif_1_0(engine);
        }

        for(uint32_t j = 0; j < N; j++) {
            rlwe_key[j] = unif_1_0(engine);
        }

        std::vector<GenericKey> bundle;

        bundle.emplace_back(std::string("LWE_SECRET"), lwe_key.data(), lwe_key.size());
        bundle.emplace_back(std::string("RLWE_SECRET"), rlwe_key.data(), rlwe_key.size());

        auto rotator = ctx_bin->ConstructOperator(bundle);

        auto rlwe_encryptor = rotator->GetEncryptor();
        auto ntt = rlwe_encryptor->GetNTT();
        ntt->ForwardNTT(rlwe_key_ntt.data(), rlwe_key.data());
        rlwe_encryptor->MakeRLWE(rlwe_sample.data(), zero.data(), rlwe_key_ntt.data());
        ntt->BackwardNTT(rlwe_sample.data(), rlwe_sample.data());
        ntt->BackwardNTT(rlwe_sample.data() + N, rlwe_sample.data() + N);
        lwe_encryptor.MakeLWE(lwe_sample.data(), 0, lwe_key.data());

        rotator->BlindRotate(rlwe_sample_out.data(), lwe_sample.data(), rlwe_sample.data(), false);
        rlwe_encryptor->PhaseRLWE(phase.data(), rlwe_sample_out.data(), rlwe_key_ntt.data());
        std::copy(phase.begin(), phase.end(), output_noise.begin() + i * N);
    }

    auto expected_output_variance = ctx_bin->ComputeOutputVariance();
    auto conversion_variance = EstimateSampleVariance(output_noise.data(), rounds, N, Q);

    for(auto& var : conversion_variance) {
        auto relative_error = RelativeError(var, expected_output_variance);
        EXPECT_LE(relative_error, 1);
    }

}

TEST(NoiseTests, TestCGGIVarianceApproximate) {

    uint64_t rounds = 1000;

    // Setup
    uint64_t Q = 36028797018972161;
    uint32_t N = 1 << 10;
    uint64_t n = 512;
    double std = 3.19;
    uint64_t basis_bits = 10;
    uint64_t base = 1 << basis_bits;
    uint64_t digits = 5;

    auto lwe_encryptor = LWEEncryptor(2 * N, n, 0.0);

    AlignedVector lwe_key(n, 0);
    AlignedVector rlwe_key(N, 0);
    std::vector<uint64_t> rlwe_key_ntt(N, 0);
    AlignedVector phase(N);

    AlignedVector lwe_sample(n + 1, 0);
    AlignedVector rlwe_sample(2 * N, 0);
    AlignedVector rlwe_sample_out(2 * N, 0);

    AlignedVector zero(N, 0);

    /* randomness */
    std::random_device random_device;
    std::mt19937 engine(random_device());
    std::uniform_int_distribution<uint64_t> unif_1_0(0, 1);


    auto ctx_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);

    AlignedVector output_noise(N * rounds, 0.0);

    for(uint64_t i = 0; i < rounds; i++) {

        std::fill(rlwe_sample.begin(), rlwe_sample.end(), 0);

        for(uint32_t j = 0; j < n; j++) {
            lwe_key[j] = unif_1_0(engine);
        }

        for(uint32_t j = 0; j < N; j++) {
            rlwe_key[j] = unif_1_0(engine);
        }

        std::vector<GenericKey> bundle;

        bundle.emplace_back(std::string("LWE_SECRET"), lwe_key.data(), lwe_key.size());
        bundle.emplace_back(std::string("RLWE_SECRET"), rlwe_key.data(), rlwe_key.size());

        auto rotator = ctx_bin->ConstructOperator(bundle);

        auto rlwe_encryptor = rotator->GetEncryptor();
        auto ntt = rlwe_encryptor->GetNTT();
        ntt->ForwardNTT(rlwe_key_ntt.data(), rlwe_key.data());
        rlwe_encryptor->MakeRLWE(rlwe_sample.data(), zero.data(), rlwe_key_ntt.data());
        ntt->BackwardNTT(rlwe_sample.data(), rlwe_sample.data());
        ntt->BackwardNTT(rlwe_sample.data() + N, rlwe_sample.data() + N);
        lwe_encryptor.MakeLWE(lwe_sample.data(), 0, lwe_key.data());

        rotator->BlindRotate(rlwe_sample_out.data(), lwe_sample.data(), rlwe_sample.data(), false);
        rlwe_encryptor->PhaseRLWE(phase.data(), rlwe_sample_out.data(), rlwe_key_ntt.data());
        std::copy(phase.begin(), phase.end(), output_noise.begin() + i * N);
    }

    auto expected_output_variance = ctx_bin->ComputeOutputVariance();
    auto conversion_variance = EstimateSampleVariance(output_noise.data(), rounds, N, Q);

    for(auto& var : conversion_variance) {
        auto relative_error = RelativeError(var, expected_output_variance);
        EXPECT_LE(relative_error, 1);
    }

}