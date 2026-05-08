//
// Created by leonard on 2/17/26.
//
#include <gtest/gtest.h>
#include "digit_decomposer_fast.h"
#include "openfhe.h"


TEST(FastDigitDecompositionTests, TestPhaseDecomp) {

    uint32_t msg_bits = 25;
    uint32_t digit_bits = 5;
    uint32_t T = 1u << digit_bits;
    uint64_t N = 1 << 11;
    uint64_t alpha = 32;
    uint64_t alpha_bits = 5;
    uint64_t beta = alpha;
    uint32_t lwe_n = 1600;
    double std = 3.19;
    uint64_t q_large = 1ull << 32;

    // parameter offer roughly 120bit security
    FastDigitDecompositionParameters params {
        lwe_n, 2048, 4096, q_large, 32, 36028797018972161, std, 1 << 10,
        5, alpha, beta, 32
    };

    auto sk_hamming_weight = beta;
    auto decomposer = FastDigitDecomposer(params);

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();

    NativeVector huge_key_vec = genT.GenerateVector(lwe_n, q_large);

    // there's no nicer way to set the hamming weight...
    auto one_ctr = 0;
    for (uint32_t i = 0; i < lwe_n; i++) {
        if (huge_key_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                huge_key_vec[i] = 0;
            }

        }
    }

    NativeVector sk_lwe = huge_key_vec;
    NativeVector sk2N = sk_lwe;
    sk2N.SetModulus(N);
    NativePoly sk_rlwe = NativePoly(genT, decomposer.GetRingParams()->GetPolyParams(), COEFFICIENT);

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
    auto ct_digits = decomposer.PhaseDecomp(ct_large);

    std::vector<uint64_t> phases;
    auto phase_large = decrypt_lwe(sk_lwe, ct_large);

    for (auto& ct : ct_digits) {
        auto p_i = decrypt_lwe(sk2N, ct);
        phases.push_back(p_i.ConvertToInt<uint64_t>());
    }

    uint64_t acc = 0;
    for (uint32_t i = 0; i < ct_digits.size(); i++) {
        uint64_t tmp = uint64_t(phases[i]) << (alpha_bits * i);
        acc += tmp;
    }
    acc %= q_large;

    EXPECT_EQ(acc, phase_large.ConvertToInt<uint64_t>());
}

TEST(FastDigitDecompositionTests, TestHomTrunc) {

    uint32_t msg_bits = 25;
    uint32_t digit_bits = 5;
    uint32_t T = 1u << digit_bits;
    uint64_t N = 1 << 11;
    uint64_t alpha = 32;
    uint64_t alpha_bits = 5;
    uint64_t beta = alpha;
    uint32_t lwe_n = 1600;
    double std = 0;
    uint64_t q_large = 1ull << 32;

    // parameter offer roughly 120bit security
    FastDigitDecompositionParameters params {
        lwe_n, 2048, 4096, q_large, 32, 36028797018972161, std, 1 << 10,
        5, alpha, beta, 1 << 10
    };

    auto sk_hamming_weight = beta;
    auto decomposer = FastDigitDecomposer(params);

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();

    NativeVector huge_key_vec = genT.GenerateVector(lwe_n, q_large);

    // there's no nicer way to set the hamming weight...
    auto one_ctr = 0;
    for (uint32_t i = 0; i < lwe_n; i++) {
        if (huge_key_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                huge_key_vec[i] = 0;
            }

        }
    }

    NativeVector sk_lwe = huge_key_vec;
    NativeVector sk2N = sk_lwe;
    sk2N.SetModulus(N);
    NativePoly sk_rlwe = NativePoly(genT, decomposer.GetRingParams()->GetPolyParams(), COEFFICIENT);

    decomposer.SetKeys(sk_rlwe, sk_lwe);

    srand(time(nullptr));


    auto ring_params = decomposer.GetRingParams();
    NativePoly msg(ring_params->GetPolyParams(), COEFFICIENT, true);


    for (uint32_t i = 0; i < N; i+=alpha) {
        msg[i] = i;
    }

    msg.SetFormat(EVALUATION);

    NativePoly skNtt = sk_rlwe.Clone();
    skNtt.SetFormat(EVALUATION);
    auto acc = encrypt_rlwe(ring_params, msg, skNtt);
    auto trunc = decomposer.HomTrunc(acc);
    trunc->SetFormat(EVALUATION);
    std::cerr << "WTF" << std::endl;
    auto trunc_phase = trunc->GetElements()[1] - trunc->GetElements()[0] * skNtt;
    trunc_phase.SetFormat(COEFFICIENT);

    std::cerr << trunc_phase << std::endl;
}

TEST(FastDigitDecompositionTests, TestSampleExtract) {

    uint32_t msg_bits = 25;
    uint32_t digit_bits = 5;
    uint32_t T = 1u << digit_bits;
    uint64_t N = 1 << 11;
    uint64_t alpha = 32;
    uint64_t alpha_bits = 5;
    uint64_t beta = alpha;
    uint32_t lwe_n = 1600;
    double std = 0;
    uint64_t q_large = 1ull << 32;

    // parameter offer roughly 120bit security
    FastDigitDecompositionParameters params {
        lwe_n, 2048, 4096, q_large, 32, 36028797018972161, std, 1 << 10,
        5, alpha, beta, 1 << 10
    };

    auto sk_hamming_weight = beta;
    auto decomposer = FastDigitDecomposer(params);

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();

    NativeVector huge_key_vec = genT.GenerateVector(lwe_n, q_large);

    // there's no nicer way to set the hamming weight...
    auto one_ctr = 0;
    for (uint32_t i = 0; i < lwe_n; i++) {
        if (huge_key_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                huge_key_vec[i] = 0;
            }

        }
    }

    NativeVector sk_lwe = huge_key_vec;
    NativeVector sk2N = sk_lwe;
    sk2N.SetModulus(N);
    NativePoly sk_rlwe = NativePoly(genT, decomposer.GetRingParams()->GetPolyParams(), COEFFICIENT);

    decomposer.SetKeys(sk_rlwe, sk_lwe);

    srand(time(nullptr));


    auto ring_params = decomposer.GetRingParams();
    NativePoly msg(ring_params->GetPolyParams(), COEFFICIENT, true);


    for (uint32_t i = 0; i < N; i++) {
        msg[i] = i;
    }

    msg.SetFormat(EVALUATION);

    NativePoly skNtt = sk_rlwe.Clone();
    skNtt.SetFormat(COEFFICIENT);
    NativeVector sk_vec = skNtt.GetValues();
    skNtt.SetFormat(EVALUATION);
    auto acc = encrypt_rlwe(ring_params, msg, skNtt);
    acc->SetFormat(COEFFICIENT);
    std::vector<uint64_t> phase_slots;
    for (uint32_t i =0; i < N; i++ ) {
        auto ct_i = RLWESampleExtract(acc->GetElements()[0], acc->GetElements()[1], i);
        auto p_i = decrypt_lwe(sk_vec, ct_i).ConvertToInt<uint64_t>();
        EXPECT_EQ(p_i, i);
    }

}

TEST(FastDigitDecompositionTests, TestSwitchToRLWE) {

    uint32_t msg_bits = 25;
    uint32_t digit_bits = 5;
    uint32_t T = 1u << digit_bits;
    uint64_t N = 1 << 11;
    uint64_t alpha = 32;
    uint64_t alpha_bits = 5;
    uint64_t beta = alpha;
    uint32_t lwe_n = 1600;
    double std = 0;
    uint64_t q_large = 1ull << 32;

    // parameter offer roughly 120bit security
    FastDigitDecompositionParameters params {
        lwe_n, 2048, 4096, q_large, 32, 36028797018972161, std, 1 << 10,
        5, alpha, beta, 1 << 10
    };

    auto sk_hamming_weight = beta;
    auto decomposer = FastDigitDecomposer(params);

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();

    NativeVector huge_key_vec = genT.GenerateVector(lwe_n, q_large);

    // there's no nicer way to set the hamming weight...
    auto one_ctr = 0;
    for (uint32_t i = 0; i < lwe_n; i++) {
        if (huge_key_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                huge_key_vec[i] = 0;
            }

        }
    }

    NativeVector sk_lwe = huge_key_vec;
    NativeVector sk2N = sk_lwe;
    sk2N.SetModulus(N);
    NativePoly sk_rlwe = NativePoly(genT, decomposer.GetRingParams()->GetPolyParams(), COEFFICIENT);
    NativeVector sk_vec = sk_rlwe.GetValues();
    decomposer.SetKeys(sk_rlwe, sk_lwe);

    srand(time(nullptr));


    auto ring_params = decomposer.GetRingParams();
    NativePoly msg(ring_params->GetPolyParams(), COEFFICIENT, true);

    auto ct = encrypt_lwe(sk_vec, 123456);
    auto sw = decomposer.LWE2RLWE(ct);
    sw->SetFormat(EVALUATION);
    sk_rlwe.SetFormat(EVALUATION);
    auto phase = sw->GetElements()[1] - sw->GetElements()[0] * sk_rlwe;
    phase.SetFormat(COEFFICIENT);
    std::cerr << phase << std::endl;
}


TEST(FastDigitDecompositionTests, TestDigitDecomp) {

    uint32_t msg_bits = 22;
    uint32_t digit_bits = 4;
    uint64_t N = 1 << 11;
    uint64_t alpha = 1 << digit_bits;
    uint32_t T = alpha;
    uint64_t beta = 110;
    uint32_t lwe_n = 3 * 1024;
    double std = 0;
    uint32_t q_bits = 48;
    uint64_t q_large = 1ull << q_bits;
    // TODO: Packing lmao
    // parameter offer roughly 120bit security
    FastDigitDecompositionParameters params {
        lwe_n, 2048, 4096, q_large, q_bits, 72057594037641217, std, 1 << 12,
        5, alpha, beta, 1 << 10
    };

    auto sk_hamming_weight = beta;
    auto decomposer = FastDigitDecomposer(params);

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();

    NativeVector huge_key_vec = genT.GenerateVector(lwe_n, q_large);

    // there's no nicer way to set the hamming weight...
    auto one_ctr = 0;
    for (uint32_t i = 0; i < lwe_n; i++) {
        if (huge_key_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                huge_key_vec[i] = 0;
            }

        }
    }
    NativeVector ring_sk_vec = genT.GenerateVector(N, params.m_Q);
    one_ctr = 0;
    for (uint32_t i = 0; i < params.m_N; i++) {
        if (ring_sk_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                ring_sk_vec[i] = 0;
            }

        }
    }

    NativeVector sk_lwe = huge_key_vec;

    NativePoly sk_rlwe = NativePoly( decomposer.GetRingParams()->GetPolyParams(), COEFFICIENT, true);
    sk_rlwe.SetValues(ring_sk_vec, COEFFICIENT);

    NativeVector sk2N = sk_lwe;
    sk2N.SetModulus(2*N);
    NativeVector skQ = sk_rlwe.GetValues();

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

    auto ct_digits = decomposer.DigitDecompose(ct_large);
    auto toc = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(toc-tic).count();
    std::cerr << "Took " << elapsed << std::endl;

    std::vector<uint64_t> phases;
    auto phase_large = decrypt_lwe(sk_lwe, ct_large);

    std::cerr << "TOTAL PHASE = " << phase_large << std::endl;
    auto p0 = decrypt_lwe(sk2N, ct_digits[0]);
    std::cerr << "Phase 0 " << p0 << std::endl;

    for (uint32_t i = 1; i < ct_digits.size(); i++) {
        auto p_i = decrypt_lwe(skQ, ct_digits[i]);
        std::cerr << "i = " << i << " " << p_i << " " << p_i.MultiplyAndRound(T, params.m_Q) << std::endl;
    }

}