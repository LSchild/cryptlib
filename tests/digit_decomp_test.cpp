//
// Created by lschild on 20/11/24.
//

// TODO

#include <gtest/gtest.h>
#include "digit_decomposer.h"
#include "openfhe.h"

class DigitDecompTests : public testing::Test {

protected:

    BitDecomposerParams m_params;
    LWEPrivateKey m_sk_lwe;
    NativePoly m_sk, m_sk_ntt;

    void SetUp() override {

        // Sue me
        srand(time(nullptr));
        uint64_t Q = 36028797018972161;
        uint32_t Qlarge_bits = 32;
        uint64_t Qlarge = 1ull << Qlarge_bits;
        uint32_t Nbits = 11;
        uint32_t N = 1 << Nbits;
        uint64_t q = 2 * N;
        uint32_t e_bits = 7;
        uint32_t clear_bits = (Nbits + 1) - e_bits;
        uint32_t n = 1305;
        uint32_t L = 1 << 20;
        double std = 3.19;
        BitDecomposerParams params{n, N, q, Qlarge, Qlarge_bits, Q, std, L, L,clear_bits, e_bits};
        m_params = params;

    }

};

TEST_F(DigitDecompTests, TestTrivial) {

    BitDecomposer bitDecomposer(m_params);

    auto genT = TernaryUniformGeneratorImpl<NativeVector>();

    auto beta = 1u << (m_params.m_error_width - 2);

    // sparsity requirement
    auto sk_vec = genT.GenerateVector(bitDecomposer.m_params.m_n, bitDecomposer.m_params.m_q_large, beta);
    auto ring_sk_vec = genT.GenerateVector(bitDecomposer.m_params.m_N, bitDecomposer.m_params.m_Q, beta);
    //sk_vec.ModMulEq(0);
    //sk_vec.at(0) = 1;

    m_sk = NativePoly(bitDecomposer.m_br_params->GetPolyParams(), COEFFICIENT, true);
    m_sk.SetValues(ring_sk_vec, COEFFICIENT);

    auto lwe_key_q = std::make_shared<LWEPrivateKeyImpl>(sk_vec);

    bitDecomposer.SetKeys(m_sk, lwe_key_q);

    auto pt_space_bits = m_params.m_q_large_bits - m_params.m_error_width;
    auto bigT = 1u << pt_space_bits;

    auto msg_raw = random() % bigT;

    NativeInteger msg = msg_raw * m_params.m_q_large / bigT;
    auto ct = encrypt_lwe(sk_vec, msg);
    auto digits = bitDecomposer.TightDecompose(ct);

    auto alpha = 1u << m_params.m_error_width;
    auto shift = 12 - m_params.m_error_width;

    auto last_digit_pt_space = pt_space_bits % shift;
    last_digit_pt_space = last_digit_pt_space == 0 ? shift : last_digit_pt_space;

    std::cerr << "Full message is " << msg_raw << std::endl;
    for(int i = 0; i < digits.size(); i++) {
        auto ct_i = digits.at(i);
        auto c_Q = ct_i->GetModulus().ConvertToInt<int64_t>();

        sk_vec.SwitchModulus(c_Q);
        auto phase = ct_i->GetB().ModSub(DotProduct(ct_i->GetA(), sk_vec), c_Q);

        auto t_i = c_Q / alpha;
        if (i == digits.size() - 1) {
            t_i = (1 << last_digit_pt_space);
        }


        auto decoded = phase.MultiplyAndRound(t_i, c_Q).Mod(t_i);
        auto decoded_u = decoded.ConvertToInt<uint64_t>();

        auto expected_digit = msg_raw % t_i;
        int64_t error_term = (phase.ConvertToInt<int64_t>() - int64_t(expected_digit * c_Q / t_i)) % c_Q;
        error_term = error_term >= m_params.m_N ? c_Q - error_term : error_term;

        std::cerr << "Error is " << error_term << std::endl;

        std::cerr << expected_digit << " " << decoded_u << std::endl;

        EXPECT_EQ(decoded_u, expected_digit);
        msg_raw >>= shift;
    }

    /*
    for(auto& ct_i : digits) {

        auto c_Q = ct_i->GetModulus();
        auto sk_Q = sk_source_Q;
        sk_Q.SwitchModulus(c_Q);

        auto phase = ct_i->GetB().ModSub(DotProduct(ct_i->GetA(), sk_Q), c_Q);

        auto t_i = c_Q / alpha;
        auto decoded = phase.MultiplyAndRound(t_i, c_Q).Mod(t_i);

        std::cerr << "Phase result " << phase << " " << decoded << " modulus is " << c_Q << std::endl;
    }
    */
}


TEST_F(DigitDecompTests, TestBitDecomp) {

    BitDecomposer bitDecomposer(m_params);

    auto genT = TernaryUniformGeneratorImpl<NativeVector>();

    auto beta = 1u << (m_params.m_error_width - 2);

    // sparsity requirement
    auto sk_vec = genT.GenerateVector(bitDecomposer.m_params.m_N, bitDecomposer.m_params.m_Q, beta);
    //sk_vec.ModMulEq(0);
    //sk_vec.at(0) = 1;

    m_sk = NativePoly(bitDecomposer.m_br_params->GetPolyParams(), COEFFICIENT, true);
    m_sk.SetValues(sk_vec, COEFFICIENT);

    //m_sk.SetValuesToZero();
    //m_sk.at(0) = 1;

    auto sk_source_big_q = m_sk.GetValues();
    sk_source_big_q.SwitchModulus(m_params.m_q_large);

    auto sk_source_q = m_sk.GetValues();
    sk_source_q.SwitchModulus(m_params.m_q);

    auto sk_source_Q = m_sk.GetValues();
    sk_source_Q.SwitchModulus(m_params.m_Q);

    auto lwe_key_q = std::make_shared<LWEPrivateKeyImpl>(sk_source_q);

    bitDecomposer.SetKeys(m_sk, lwe_key_q);

    auto pt_space_bits = m_params.m_q_large_bits - m_params.m_error_width;
    auto bigT = 1u << pt_space_bits;

    auto msg_raw = random() % bigT;

    NativeInteger msg = msg_raw * m_params.m_q_large / bigT;
    auto ct = encrypt_lwe(sk_source_big_q, msg);
    auto digits = bitDecomposer.BitDecompose(ct);

    for(int i = 0; i < pt_space_bits; i++) {
        auto ct_i = digits.at(i);
        auto c_Q = ct_i->GetModulus();
        auto sk_Q = sk_source_Q;
        sk_Q.SwitchModulus(c_Q);

        auto phase = ct_i->GetB().ModSub(DotProduct(ct_i->GetA(), sk_Q), c_Q);

        auto decoded = phase.MultiplyAndRound(2, c_Q).Mod(2);
        auto decoded_u = decoded.ConvertToInt<uint64_t>();

        auto expected_digit = msg_raw & 1;
        EXPECT_EQ(decoded_u, expected_digit);

        msg_raw >>= 1;
    }

}

TEST(DigitDecompositionTests, TestDigitDecomp) {

    uint32_t q_bits = 48;
    uint32_t error_bits = 8;
    uint64_t alpha = 1ull << error_bits;
    uint64_t msg_bits = (q_bits - error_bits);
    uint64_t digit_bits = 12 - error_bits;

    uint32_t T = 1u << digit_bits;
    double std = 0;
    uint32_t n = 3 * 1024;
    // parameter offer roughly 120bit security
    DigitDecompositionParameters params {
        n, 2048, 4096, 1ull << q_bits, q_bits, 72057594037641217, std, 1 << 12,
        (uint32_t)digit_bits, error_bits, 1ull << 42, 32
    };

    auto sk_hamming_weight = 1u << (params.m_error_width - 1);
    auto decomposer = DigitDecomposer(params);

    auto genT = TernaryUniformGeneratorImpl<NativeVector>();
    NativeVector sk_lwe = genT.GenerateVector(params.m_n, params.m_q_large, sk_hamming_weight);
    NativePoly sk_rlwe = NativePoly(genT, decomposer.GetRingParams()->GetPolyParams(), COEFFICIENT);

    for (uint32_t i = 0; i < sk_lwe.GetLength(); i++) {
        if (sk_lwe[i] >= 1) {
            sk_lwe[i] = 1;
        }
    }

    for (uint32_t i = 0; i < sk_rlwe.GetLength(); i++) {
        if (sk_rlwe[i] >= 1) {
            sk_rlwe[i] = 1;
        }
    }

    decomposer.SetKeys(sk_rlwe, sk_lwe);

    srand(time(nullptr));

    auto msg = random() % (1ull << msg_bits);
    auto msg_encoded = msg * alpha;
    std::vector<uint64_t> msg_digits;
    for (int shift = 0; shift < msg_bits; shift+=digit_bits) {
        auto msg_digit = (msg >> shift) & (T - 1);
        msg_digits.push_back(msg_digit);
    }

    auto ct_large = encrypt_lwe(sk_lwe, msg_encoded, std);

    auto tic = std::chrono::high_resolution_clock::now();

    std::vector<LWECiphertext> digits_ct = decomposer.DigitDecompose(ct_large);

    auto toc = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic).count();
    std::cerr << "Took " << elapsed << std::endl;

    EXPECT_EQ(digits_ct.size(), msg_digits.size());

    /*
    std::vector<uint64_t> recovered_digits;
    for (auto& ct : digits_ct) {
        auto phase = decrypt_lwe(sk_lwe, ct);
        auto phase_decoded = phase.MultiplyAndRound(T, ct->GetModulus()).Mod(T);
        std::cerr << phase << " " << phase.ModSub(phase_decoded * alpha, ct->GetModulus()) << std::endl;
        recovered_digits.push_back(phase_decoded.ConvertToInt<uint64_t>());
    }

    std::cerr << "Message was " << msg << std::endl;

    for (int i = 0; i < msg_digits.size(); i++) {
        EXPECT_EQ(msg_digits[i], recovered_digits[i]);
    }
    */

}

TEST(DigitDecompositionTests, TestBitDecomp) {

    uint32_t msg_bits = 25;
    uint32_t digit_bits = 1;
    uint32_t T = 1u << digit_bits;
    uint64_t alpha = 1u << (32 - 25);
    double std = 3.19;

    DigitDecompositionParameters params {
        1800, 2048, 4096, 1ull << 32, 32, 72057594037641217, std, 1 << 12,
        5, 7, 1ull << 42, 8
    };

    auto sk_hamming_weight = 16;
    auto decomposer = DigitDecomposer(params);

    auto genT = TernaryUniformGeneratorImpl<NativeVector>();
    NativeVector sk_lwe = genT.GenerateVector(params.m_n, params.m_q_large, sk_hamming_weight);
    NativePoly sk_rlwe = NativePoly(decomposer.GetRingParams()->GetPolyParams(), COEFFICIENT, true);

    for (uint32_t i = 0; i < 1800; i++) {
        if (sk_lwe[i].ConvertToInt<uint64_t>() > 1)
            sk_lwe[i] = 1;
    }

    decomposer.SetKeys(sk_rlwe, sk_lwe);

    srand(time(nullptr));

    auto msg = random() % (1ull << msg_bits);
    auto msg_encoded = msg * alpha;
    std::vector<uint64_t> msg_digits;
    for (int shift = 0; shift < msg_bits; shift+=digit_bits) {
        auto msg_digit = (msg >> shift) & (T - 1);
        msg_digits.push_back(msg_digit);
    }

    std::cerr << "Message is " << msg << std::endl;

    auto ct_large = encrypt_lwe(sk_lwe, msg_encoded, std);

    std::vector<LWECiphertext> digits_ct = decomposer.BitDecompose(ct_large);

    EXPECT_EQ(digits_ct.size(), msg_digits.size());

    std::vector<uint64_t> recovered_digits;
    for (auto& ct : digits_ct) {
        auto phase = decrypt_lwe(sk_rlwe.GetValues(), ct);
        auto phase_decoded = phase.MultiplyAndRound(T, ct->GetModulus()).Mod(T);
        recovered_digits.push_back(phase_decoded.ConvertToInt<uint64_t>());
    }

    //std::cerr << "Message was " << msg << std::endl;

    for (int i = 0; i < msg_digits.size(); i++) {
        // recovered digits is MSB to LSB
        // msg_digits is LSB to MSB to we compare the ends
        std::cerr << i << " " << msg_digits[i] << " " << recovered_digits[msg_digits.size() - i - 1] << std::endl;
        EXPECT_EQ(msg_digits[i], recovered_digits[msg_digits.size() - i - 1]);
    }

}