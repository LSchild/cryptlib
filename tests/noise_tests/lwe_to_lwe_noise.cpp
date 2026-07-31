//
// Created by leonard on 7/28/26.
//

#include <gtest/gtest.h>
#include <random>

#include "base_crypto.h"
#include "operators/endo_glwe_conversion.h"
#include "utils/math_utils.h"

TEST(LWEConversionNoiseTest,TestBinaryKeyExact) {

    const uint64_t rounds = 100;

    const uint64_t modulus = 1ull << 50;
    const uint64_t source_dimension = 1024;
    const uint64_t target_dimension = 512;
    const uint64_t gadget_basis = 1 << 10;
    const uint64_t gadget_digits = 5;
    const long double std = 3.19;

    /* key buffers */
    std::vector<uint64_t> source_key(source_dimension);
    std::vector<uint64_t> target_key(target_dimension);
    std::vector<uint64_t> input(source_dimension + 1);
    std::vector<uint64_t> target(target_dimension + 1);

    /* randomness */
    std::random_device random_device;
    std::mt19937 engine(random_device());
    std::uniform_int_distribution<uint64_t> unif_1_0(0, 1);

    auto conversion_context = std::make_shared<LWEConversionContext>(BINARY, modulus, source_dimension, target_dimension, gadget_basis, gadget_digits, std);
    auto expected_output_variance = conversion_context->ComputeOutputVariance();

    auto source_encryptor = std::make_shared<LWEEncryptor>(modulus, source_dimension, 0.0);
    auto target_encryptor = std::make_shared<LWEEncryptor>(modulus, target_dimension, 0.0);

    AlignedVector output_noise(rounds);
    for(uint64_t i = 0; i < rounds; i++) {

        std::fill(target.begin(), target.end(), 0);

        /* generate keys */
        for(auto& s_v : source_key) {
            s_v = unif_1_0(engine);
        }

        for(auto& s_t : target_key) {
            s_t = unif_1_0(engine);
        }

        std::vector<GenericKey> keys = {source_key, target_key};

        /* generate converter */
        auto converter = conversion_context->ConstructOperator(keys);

        /* generate source ciphertext */
        source_encryptor->MakeLWE(input.data(), 0, source_key.data());

        /* perform conversion */
        converter->Convert(target, input);

        /* decrypt and get result = pure noise */
        auto noise = target_encryptor->PhaseLWE(target.data(), target_key.data());
        output_noise[i] = noise;


    }

    auto conversion_variance = EstimateSampleVariance(output_noise.data(), rounds, 1, modulus);
    auto relative_error = RelativeError(conversion_variance[0], expected_output_variance);

    EXPECT_LE(relative_error, 1);

}

TEST(LWEConversionNoiseTest,TestBinaryKeyApproximate) {

    const uint64_t rounds = 1000;

    const uint64_t modulus = 1ull << 50;
    const uint64_t source_dimension = 1024;
    const uint64_t target_dimension = 512;
    const uint64_t gadget_basis = 1 << 10;
    const uint64_t gadget_digits = 3;
    const long double std = 3.19;

    /* key buffers */
    std::vector<uint64_t> source_key(source_dimension);
    std::vector<uint64_t> target_key(target_dimension);
    std::vector<uint64_t> input(source_dimension + 1);
    std::vector<uint64_t> target(target_dimension + 1);

    /* randomness */
    std::random_device random_device;
    std::mt19937 engine(random_device());
    std::uniform_int_distribution<uint64_t> unif_1_0(0, 1);

    auto conversion_context = std::make_shared<LWEConversionContext>(BINARY, modulus, source_dimension, target_dimension, gadget_basis, gadget_digits, std);
    auto expected_output_variance = conversion_context->ComputeOutputVariance();

    auto source_encryptor = std::make_shared<LWEEncryptor>(modulus, source_dimension, 0.0);
    auto target_encryptor = std::make_shared<LWEEncryptor>(modulus, target_dimension, 0.0);

    AlignedVector output_noise(rounds);
    for(uint64_t i = 0; i < rounds; i++) {

        std::fill(target.begin(), target.end(), 0);

        /* generate keys */
        for(auto& s_v : source_key) {
            s_v = unif_1_0(engine);
        }

        for(auto& s_t : target_key) {
            s_t = unif_1_0(engine);
        }

        std::vector<GenericKey> keys = {source_key, target_key};

        /* generate converter */
        auto converter = conversion_context->ConstructOperator(keys);

        /* generate source ciphertext */
        source_encryptor->MakeLWE(input.data(), 0, source_key.data());

        /* perform conversion */
        converter->Convert(target, input);

        /* decrypt and get result = pure noise */
        auto noise = target_encryptor->PhaseLWE(target.data(), target_key.data());
        output_noise[i] = noise;


    }

    // check relative error
    auto conversion_variance = EstimateSampleVariance(output_noise.data(), rounds, 1, modulus);
    auto relative_error = RelativeError(conversion_variance[0], expected_output_variance);

    EXPECT_LE(relative_error, 1);

    // check absolute error
    auto delta = std::powl(2.0, IntLog2(modulus)) / std::powl(gadget_basis, gadget_digits);
    auto bound = source_dimension * delta / 2.0;

    for(auto& o_i : output_noise) {
        auto signed_value = UnsignedToSignedRepr(o_i, modulus);
        EXPECT_LE(std::abs(signed_value), bound);
    }


}