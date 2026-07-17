//
// Created by leonard on 7/16/26.
//

#include "operators/sap_decomposition.h"
#include "operators/lwe_to_rgsw_conversion.h"
#include "operators/cggi_blind_rotation.h"
#include "operators/endo_glwe_conversion.h"

#include <gtest/gtest.h>

class SAPDecompositionTests : public testing::Test {

protected:

    std::vector<GenericKey> test_keys;

    std::shared_ptr<SAPDecompositionContext> sap_context;

    void SetUp() override {

        // declare params
        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint64_t n = 100;
        double std = 0;
        uint64_t basebits = 14;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;

        auto params_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);
        auto ntt = params_bin->GetNTT();
        auto trace_params = std::make_shared<TraceEvaluationContext>(BINARY, ntt, base, digits, std);
        auto lwe_conversion_params = std::make_shared<LWEConversionContext>(BINARY, Q, N, n, base, digits, std);

        std::shared_ptr<OperatorContext<BlindRotator>> params_rot = std::reinterpret_pointer_cast<OperatorContext<BlindRotator>>(params_bin);
        sap_context = std::make_shared<SAPDecompositionContext>(params_rot, trace_params, lwe_conversion_params, 32);

        std::vector<uint64_t> tmp_key(std::max<uint64_t>(N, n), 0);


        // gen sparse LWE key
        for(uint32_t i = 0; i < SAPDecompositionContext::IMPLICIT_SK_L0; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("LWE_KEY", tmp_key.data(), n);

        // gen RLWE key
        for (uint32_t i = 0; i< N; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("RLWE_KEY", tmp_key.data(), N);

        // ntt'ed RLWE key
        ntt->ComputeForward(tmp_key.data(), tmp_key.data(), 1 ,1);
        test_keys.emplace_back("RLWE_KEY_NTT", tmp_key.data(), N);
    }

};

TEST_F(SAPDecompositionTests, TestHomTrunc) {

    auto decomp = sap_context->ConstructOperator(test_keys);
    auto ntt = sap_context->GetTraceContext()->GetNTT();

    auto N = ntt->GetDegree();
    auto Q = ntt->GetModulus();
    auto radix = sap_context->GetDefaultRadix();
    auto quot_space = N / radix;

    std::random_device random_device;
    std::mt19937 gen(random_device());

    std::uniform_int_distribution<uint64_t> hi_dist(0, quot_space);
    std::uniform_int_distribution<uint64_t> lo_dist(0, radix);

    auto truncated_part = hi_dist(gen);
    auto exponent = truncated_part * radix + lo_dist(gen);

    AlignedVector rlwe_sample(2 * N);
    AlignedVector rlwe_sample_out(2 * N, 0);
    rlwe_sample[N + exponent] = 1;
    ntt->ComputeForward(rlwe_sample.data() + N, rlwe_sample.data() + N, 1, 1);

    decomp->HomTrunc(rlwe_sample_out.data(), rlwe_sample.data(), radix);

    auto ntt_key = test_keys[2];
    intel::hexl::EltwiseMultMod(rlwe_sample_out.data(), rlwe_sample_out.data(), ntt_key.GetKeyPtr(), N, Q, 1);
    intel::hexl::EltwiseSubMod(rlwe_sample_out.data() + N, rlwe_sample_out.data() + N, rlwe_sample_out.data(), N, Q);
    ntt->ComputeInverse(rlwe_sample_out.data(), rlwe_sample_out.data() + N, 1, 1);

    for(uint64_t i = 0; i < N; i++) {
        if (i != truncated_part) {
            EXPECT_EQ(rlwe_sample_out[i], 0);
        } else {
            EXPECT_EQ(rlwe_sample_out[i], 1);
        }
    }
}