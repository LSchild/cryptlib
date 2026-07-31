//
// Created by leonard on 7/29/26.
//


#include <gtest/gtest.h>
#include <random>

#include "base_crypto.h"
#include "operators/endo_glwe_conversion.h"
#include "utils/math_utils.h"

TEST(RLWEConversionNoiseTest, TestBinaryKeyExact) {

    const uint64_t rounds = 1000;

    const uint64_t modulus = 36028797018972161; // 56 bit
    const uint64_t dimension = 1 << 10;
    const uint64_t gadget_basis = 1 << 8;
    const uint64_t gadget_digits = 7;
    const long double std = 3.19;

    /* key buffers */
    std::vector<uint64_t> source_key(dimension);
    std::vector<uint64_t> target_key(dimension);
    std::vector<uint64_t> source_key_ntt(dimension);
    std::vector<uint64_t> target_key_ntt(dimension);

    std::vector<uint64_t> input(2 * dimension);
    std::vector<uint64_t> target(2 * dimension);
    std::vector<uint64_t> message(dimension);
    std::vector<uint64_t> phase(dimension);

    std::fill(message.begin(), message.end(), 0);

    /* randomness */
    std::random_device random_device;
    std::mt19937 engine(random_device());
    std::uniform_int_distribution<uint64_t> unif_1_0(0, 1);

    auto conversion_context = std::make_shared<RLWEConversionContext>(BINARY, modulus, dimension, gadget_basis, gadget_digits, std);
    auto worker = conversion_context->GetNTT();

    auto expected_output_variance = conversion_context->ComputeOutputVariance();

    auto encryptor = std::make_shared<RLWEEncryptor>(modulus, dimension, 0.0);

    AlignedVector output_noise(rounds * dimension);

    for(uint64_t i = 0; i < rounds; i++) {

        std::fill(target.begin(), target.end(), 0);

        /* generate keys */
        for(uint64_t j = 0; j < dimension; j++) {
            source_key[j] = unif_1_0(engine);
        }

        for(uint64_t j = 0; j < dimension; j++) {
            target_key[j] = unif_1_0(engine);
        }

        worker->ForwardNTT(source_key_ntt, source_key);
        worker->ForwardNTT(target_key_ntt, target_key);

        std::vector<GenericKey> keys = {source_key, target_key};

        /* generate converter */
        auto converter = conversion_context->ConstructOperator(keys);

        /* generate source ciphertext */
        encryptor->MakeRLWE(input.data(), message.data(), source_key_ntt.data());
        // TODO: Fixme
        worker->BackwardNTT(input, input);

        /* perform conversion */
        converter->Convert(target, input);

        /* decrypt and get result = pure noise */
        encryptor->PhaseRLWE(phase.data(), target.data(), target_key_ntt.data());
        std::copy(phase.begin(), phase.end(), output_noise.begin() + dimension * i);


    }

    auto conversion_variance = EstimateSampleVariance(output_noise.data(), rounds, dimension, modulus);
    for(auto& var : conversion_variance) {
        auto relative_error = RelativeError(var, expected_output_variance);
        EXPECT_LE(relative_error, 1);
    }

}

TEST(RLWEConversionNoiseTest, TestBinaryKeyApproximate) {

    const uint64_t rounds = 1000;

    const uint64_t modulus = 36028797018972161; // 56 bit
    const uint64_t dimension = 1 << 10;
    const uint64_t gadget_basis = 1 << 8;
    const uint64_t gadget_digits = 5;
    const long double std = 3.19;

    /* key buffers */
    std::vector<uint64_t> source_key(dimension);
    std::vector<uint64_t> target_key(dimension);
    std::vector<uint64_t> source_key_ntt(dimension);
    std::vector<uint64_t> target_key_ntt(dimension);

    std::vector<uint64_t> input(2 * dimension);
    std::vector<uint64_t> target(2 * dimension);
    std::vector<uint64_t> message(dimension);
    std::vector<uint64_t> phase(dimension);

    std::fill(message.begin(), message.end(), 0);

    /* randomness */
    std::random_device random_device;
    std::mt19937 engine(random_device());
    std::uniform_int_distribution<uint64_t> unif_1_0(0, 1);

    auto conversion_context = std::make_shared<RLWEConversionContext>(BINARY, modulus, dimension, gadget_basis, gadget_digits, std);
    auto worker = conversion_context->GetNTT();

    auto expected_output_variance = conversion_context->ComputeOutputVariance();

    auto encryptor = std::make_shared<RLWEEncryptor>(modulus, dimension, 0.0);

    AlignedVector output_noise(rounds * dimension);

    for(uint64_t i = 0; i < rounds; i++) {

        std::fill(target.begin(), target.end(), 0);

        uint64_t source_key_hamming_weight = 0;
        /* generate keys */
        for(uint64_t j = 0; j < dimension; j++) {
            source_key[j] = unif_1_0(engine);
        }

        for(uint64_t j = 0; j < dimension; j++) {
            target_key[j] = unif_1_0(engine);
        }

        worker->ForwardNTT(source_key_ntt, source_key);
        worker->ForwardNTT(target_key_ntt, target_key);

        std::vector<GenericKey> keys = {source_key, target_key};

        /* generate converter */
        auto converter = conversion_context->ConstructOperator(keys);

        /* generate source ciphertext */
        encryptor->MakeRLWE(input.data(), message.data(), source_key_ntt.data());
        // TODO: Fixme
        worker->BackwardNTT(input, input);

        /* perform conversion */
        converter->Convert(target, input);

        /* decrypt and get result = pure noise */
        encryptor->PhaseRLWE(phase.data(), target.data(), target_key_ntt.data());
        std::copy(phase.begin(), phase.end(), output_noise.begin() + dimension * i);
    }

    auto conversion_variance = EstimateSampleVariance(output_noise.data(), rounds, dimension, modulus);
    for(auto& var : conversion_variance) {
        auto relative_error = RelativeError(var, expected_output_variance);
        EXPECT_LE(relative_error, 1);
    }

    // check absolute error
    auto delta = std::powl(2.0, IntLog2(modulus)) / std::powl(gadget_basis, gadget_digits);
    auto bound = dimension * delta / 2.0;

    for(auto& o_i : output_noise) {
        auto signed_value = UnsignedToSignedRepr(o_i, modulus);
        EXPECT_LE(std::abs(signed_value), bound);
    }

}