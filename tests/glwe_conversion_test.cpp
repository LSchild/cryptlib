#include <gtest/gtest.h>
#include "base_crypto.h"
#include "common_types.h"
#include "glwe_conversion.h"

//
// Created by leonard on 4/26/26.
//
class GLWEConversionTestGroup : public testing::Test {

protected:

    std::shared_ptr<RLWEEncryptor> enc_rlwe;
    std::shared_ptr<LWEEncryptor> enc_lwe_in;
    std::shared_ptr<LWEEncryptor> enc_lwe_out;
    std::shared_ptr<LWEtoLWEConverter> lwe_conv;
    std::shared_ptr<RLWEtoRLWEConverter> rlwe_conv;

    AlignedVector secret_rlwe_src, secret_rlwe_target, secret_rlwe_src_ntt, secret_rlwe_target_ntt;
    AlignedVector secret_lwe_src, secret_lwe_target;

    void SetUp() override {
        uint64_t Q = 2251799813773313;
        uint32_t N = 1024;
        uint64_t n_in = 1234;
        uint64_t n_out = 567;
        uint64_t L = 1 << 13;
        uint64_t digits = 4;
        double std = 0;

        enc_rlwe = std::make_shared<RLWEEncryptor>(Q, N, std);
        enc_lwe_in = std::make_shared<LWEEncryptor>(Q, n_in, std);
        enc_lwe_out = std::make_shared<LWEEncryptor>(Q, n_out, std);

        secret_lwe_src.resize(n_in);
        secret_lwe_target.resize(n_out);
        secret_rlwe_src.resize(N);
        secret_rlwe_target.resize(N);
        secret_rlwe_src_ntt.resize(N);
        secret_rlwe_target_ntt.resize(N);


        std::srand(432324253);

        for (uint32_t i = 0; i < N; i++) {
            secret_rlwe_src[i] = random() % 2;
            secret_rlwe_target[i] = random() % 2;
        }

        for(uint64_t i = 0; i < n_in; i++) {
            secret_lwe_src[i] = random() % 2;
        }

        for(uint64_t i = 0; i < n_out; i++) {
            secret_lwe_target[i] = random() % 2;
        }

        enc_rlwe->GetNTT()->ComputeForward(secret_rlwe_target_ntt.data(), secret_rlwe_target.data(), 1, 1);
        enc_rlwe->GetNTT()->ComputeForward(secret_rlwe_src_ntt.data(), secret_rlwe_src.data(), 1, 1);

        auto lwe_params = LWEConversionParameters(BINARY, Q, n_in, n_out, L, digits, std);
        auto rlwe_params = RLWEConversionParameters(BINARY, Q, N, L, digits, std);

        rlwe_conv = std::make_shared<RLWEtoRLWEConverter>(rlwe_params);
        rlwe_conv->KeyGen(secret_rlwe_src.data(), secret_rlwe_target.data());

        lwe_conv = std::make_shared<LWEtoLWEConverter>(lwe_params);
        lwe_conv->KeyGen(secret_lwe_src.data(), secret_lwe_target.data());

    }

};

TEST_F(GLWEConversionTestGroup, TestLWEtoLWE) {

    auto Q = enc_lwe_in->GetModulus();

    AlignedVector lwe_in(enc_lwe_in->GetDimension() + 1);
    AlignedVector lwe_out(enc_lwe_out->GetDimension() + 1);

    auto msg = random() % Q;
    enc_lwe_in->MakeLWE(lwe_in.data(), msg, secret_lwe_src.data());

    lwe_conv->Convert(lwe_out.data(), lwe_in.data());

    auto val = enc_lwe_out->PhaseLWE(lwe_out.data(), secret_lwe_target.data());

    EXPECT_EQ(val, msg);
}

TEST_F(GLWEConversionTestGroup, TestRLWEtoRLWE) {

    auto N = enc_rlwe->GetDimension();
    auto Q = enc_rlwe->GetDimension();

    AlignedVector rlwe_in(2 * N);
    AlignedVector rlwe_out(2 * N,0);
    AlignedVector msg(N, 0);
    AlignedVector res(N, 0);

    enc_rlwe->MakeRLWE(rlwe_in.data(), msg.data(), secret_rlwe_src_ntt.data());
    enc_rlwe->GetNTT()->ComputeInverse(rlwe_in.data(), rlwe_in.data(), 1, 1);
    //enc_rlwe->GetNTT()->ComputeInverse(rlwe_in.data() + N, rlwe_in.data() + N, 1, 1);


    rlwe_conv->Convert(rlwe_out.data(), rlwe_in.data());


    enc_rlwe->PhaseRLWE(res.data(), rlwe_out.data(), secret_rlwe_target_ntt.data());

    for(uint64_t i = 0; i < N; i++) {
        EXPECT_EQ(res[i], msg[i]);
    }

}