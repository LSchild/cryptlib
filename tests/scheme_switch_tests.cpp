//
// Created by leonard on 11/12/24.
//

#include <gtest/gtest.h>
#include "operators/cggi_blind_rotation.h"
#include "operators/bmmp_blind_rotation.h"
#include "operators/lwe_to_rgsw_conversion.h"

class CGGILWE2RGSWTests : public testing::Test {

protected:

    std::vector<GenericKey> test_keys;

    std::shared_ptr<SchemeSwitchingContext> ss_context;

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
        auto square_params = std::make_shared<RLWEConversionContext>(BINARY, ntt, base, digits, std);

        std::shared_ptr<OperatorContext<BlindRotator>> params_rot = std::reinterpret_pointer_cast<OperatorContext<BlindRotator>>(params_bin);
        ss_context = std::make_shared<SchemeSwitchingContext>(params_rot, trace_params, square_params, digits, base);

        std::vector<uint64_t> tmp_key(std::max<uint64_t>(N, n), 0);


        // gen LWE key
        for(uint32_t i = 0; i < n; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("LWE_KEY", tmp_key.data(), n);

        // gen RLWE key
        for (uint32_t i = 0; i< N; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("RLWE_KEY", tmp_key.data(), N);

        // ntt'ed RLWE key
        ntt->ForwardNTT(tmp_key.data(), tmp_key.data());
        test_keys.emplace_back("RLWE_KEY_NTT", tmp_key.data(), N);
    }

};

TEST_F(CGGILWE2RGSWTests, TestConvExactDigs) {

    std::vector<GenericKey> test(test_keys.begin(), test_keys.begin() + 2);
    auto m_converter = ss_context->ConstructOperator(test);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(ss_context->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(ss_context->GetOutputContainer(input_shape));
    auto output_params = ss_context->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;

    auto digs = ss_context->GetOutputDigits();
    auto basis = ss_context->GetOutputBasis();

    auto lwe_key = test_keys[0].GetKey();
    auto rlwe_key = test_keys[1].GetKey();
    auto rlwe_key_ntt = test_keys[2].GetKey();


    AlignedBuffer rgsw_out(4 * rlwe_N * digs, 0);
    AlignedBuffer rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedBuffer expected_result(2 * rlwe_N * digs, 0);
    AlignedBuffer lwe(lwe_n + 1, 0);

    auto ntt = output_params->GetNTT();

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_key[i]) % lwe_q;
    }

    // random bit message
    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    m_converter->Convert(rgsw_out.data(), lwe.data());

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_key_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_key.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        int64_t err = int64_t(expected_result[i]) - int64_t(rgsw_phase[i]);
        EXPECT_LE(std::abs(err), 1);
    }

}

TEST_F(CGGILWE2RGSWTests, TestConvApproxDigs) {


    ss_context->SetOutputBasis(1u << 10);
    ss_context->SetOutputDigits(4);
    auto max_digits = 6;
    auto digit_difference = max_digits - 4;


    auto m_converter = ss_context->ConstructOperator(test_keys);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(ss_context->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(ss_context->GetOutputContainer(input_shape));
    auto output_params = ss_context->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;
    auto digs = ss_context->GetOutputDigits();
    auto basis = ss_context->GetOutputBasis();

    auto rlwe_key = test_keys[1].GetKey();
    auto rlwe_key_ntt = test_keys[2].GetKey();
    auto lwe_key = test_keys[0].GetKey();

    AlignedBuffer rgsw_out(4 * rlwe_N * digs, 0);
    AlignedBuffer rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedBuffer expected_result(2 * rlwe_N * digs, 0);
    AlignedBuffer lwe(lwe_n + 1);

    auto ntt = output_params->GetNTT();

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_key[i]) % lwe_q;
    }

    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    m_converter->Convert(rgsw_out.data(), lwe.data());

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_key_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digit_difference; i++) {
        scal *= basis;
    }
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_key.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        EXPECT_EQ(expected_result[i], rgsw_phase[i]);
    }
}


class BMMPLWE2RGSWTests : public testing::Test {

protected:

    AlignedBuffer rlwe_secret;
    AlignedBuffer rlwe_secret_ntt;
    AlignedBuffer lwe_secret_binary;

    std::shared_ptr<SchemeSwitchingContext> m_params;

    void SetUp() override {

        uint64_t Q = 36028797018972161;
        uint32_t N = 1 << 11;
        uint64_t n = 100;
        double std = 0;
        uint64_t basebits = 14;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;

        auto params_bin = std::make_shared<BMMPBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std, 2);

        rlwe_secret = AlignedBuffer(N);
        rlwe_secret_ntt = AlignedBuffer(N);
        lwe_secret_binary = AlignedBuffer(n);

        for (uint32_t i = 0; i< N; i++) {
            rlwe_secret[i] = rand() % 2;
        }

        params_bin->GetNTT()->ForwardNTT(rlwe_secret_ntt.data(), rlwe_secret.data());
        for(uint32_t i = 0; i < n; i++) {
            lwe_secret_binary[i] = rand() % 2;
        }

        auto trace_params = std::make_shared<TraceEvaluationContext>(BINARY, params_bin->GetNTT(), base, digits, std);
        auto square_params = std::make_shared<RLWEConversionContext>(BINARY, params_bin->GetNTT(), base, digits, std);

        std::shared_ptr<OperatorContext<BlindRotator>> params_rot = std::reinterpret_pointer_cast<OperatorContext<BlindRotator>>(params_bin);

        m_params = std::make_shared<SchemeSwitchingContext>(params_rot, trace_params, square_params, digits, base);



    }

};

TEST_F(BMMPLWE2RGSWTests, TestConvExactDigs) {


    std::vector<GenericKey> keys = {
            {"LWE", lwe_secret_binary.data(), lwe_secret_binary.size()},
            {"RLWE", rlwe_secret.data(), rlwe_secret.size()}
    };

    auto m_converter = m_params->ConstructOperator(keys);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(m_params->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(m_params->GetOutputContainer(input_shape));
    auto output_params = m_params->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;
    auto digs = m_params->GetOutputDigits();
    auto basis = m_params->GetOutputBasis();

    AlignedBuffer lwe(lwe_n + 1);

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_secret_binary[i]) % lwe_q;
    }

    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    AlignedBuffer rgsw_out(4 * rlwe_N * digs, 0);
    AlignedBuffer rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedBuffer expected_result(2 * rlwe_N * digs, 0);

    m_converter->Convert(rgsw_out.data(), lwe.data());

    auto ntt = output_params->GetNTT();

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_secret_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_secret.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        int64_t err = int64_t(expected_result[i]) - int64_t(rgsw_phase[i]);
        EXPECT_LE(std::abs(err), 1);
    }

}

TEST_F(BMMPLWE2RGSWTests, TestConvApproxDigs) {


    m_params->SetOutputBasis(1u << 10);
    m_params->SetOutputDigits(4);
    auto max_digits = 6;
    auto digit_difference = max_digits - 4;
    std::vector<GenericKey> keys = {
            {"LWE", lwe_secret_binary.data(), lwe_secret_binary.size()},
            {"RLWE", rlwe_secret.data(), rlwe_secret.size()}
    };


    auto m_converter = m_params->ConstructOperator(keys);

    auto input_shape = std::dynamic_pointer_cast<LWEContainerImpl>(m_params->GetInputContainer());
    auto output_shape = std::dynamic_pointer_cast<RGSWContainerImpl>(m_params->GetOutputContainer(input_shape));
    auto output_params = m_params->GetSquaringContext();

    auto lwe_n = input_shape->GetN();
    auto rlwe_N = output_shape->GetN();
    auto rlwe_Q = output_shape->GetQ();
    auto lwe_q = 2 * rlwe_N;
    auto digs = m_params->GetOutputDigits();
    auto basis = m_params->GetOutputBasis();

    AlignedBuffer lwe(lwe_n + 1);

    std::srand(time(nullptr));

    for(uint32_t i = 0; i < lwe_n; i++) {
        lwe[i] = rand() % lwe_q;
        lwe[lwe_n] = (lwe[lwe_n] + lwe[i] * lwe_secret_binary[i]) % lwe_q;
    }

    auto m = rand() % 2;
    lwe[lwe_n] += m == 1 ? lwe_q >> 1 : 0;
    lwe[lwe_n] %= lwe_q;

    AlignedBuffer rgsw_out(4 * rlwe_N * digs, 0);
    AlignedBuffer rgsw_phase(2 * rlwe_N * digs, 0);
    AlignedBuffer expected_result(2 * rlwe_N * digs, 0);

    m_converter->Convert(rgsw_out.data(), lwe.data());

    auto ntt = output_params->GetNTT();

    for(uint32_t i = 0; i < 2 * digs; i++) {
        intel::hexl::EltwiseMultMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N, rlwe_secret_ntt.data(), rlwe_N, rlwe_Q, 1);
        intel::hexl::EltwiseSubMod(rgsw_phase.data() + i * rlwe_N, rgsw_out.data() + i * 2 * rlwe_N + rlwe_N,  rgsw_phase.data() +i * rlwe_N,rlwe_N, rlwe_Q);
        ntt->BackwardNTT(rgsw_phase.data() + i * rlwe_N, rgsw_phase.data() + i * rlwe_N);
    }

    uint64_t scal = 1 * m;
    for(uint32_t  i = 0; i < digit_difference; i++) {
        scal *= basis;
    }
    for(uint32_t  i = 0; i < digs; i++) {
        expected_result[digs * rlwe_N + i * rlwe_N] = scal;
        intel::hexl::EltwiseFMAMod(expected_result.data() + i * rlwe_N, rlwe_secret.data(), rlwe_Q - scal, nullptr, rlwe_N, rlwe_Q, 1);
        scal *= basis;
    }

    for(uint32_t i = 0; i < 2 * rlwe_N * digs; i++) {
        EXPECT_EQ(expected_result[i], rgsw_phase[i]);
    }
}
