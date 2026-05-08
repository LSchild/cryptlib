//
// Created by leonard on 12/11/24.
//
#include <gtest/gtest.h>
#include "giant_lut_evaluator.h"
#include "scheme_switch.h"
#include "digit_decomposer.h"
#include "utils.h"
#include <chrono>

TEST(GiantLUTTest, TestDBWithPolyPrecompRaw) {

    std::srand(0);
    // set up parameters

    // variance
    double std = 0; // for testing

    // moduli
    uint64_t ring_modulus = 36028797018972161;
    uint64_t small_lwe_modulus = 1 << 12;
    uint32_t large_lwe_modulus_bits = 32;
    uint64_t large_lwe_modulus = 1ull << large_lwe_modulus_bits;
    uint64_t key_switch_modulus = 1 << 30;
    uint32_t large_lwe_n = 1305;

    // (r)lwe dimensions
    uint32_t N_bits = std::round(std::log2((long double) small_lwe_modulus)) - 1;
    uint32_t N = 1u << N_bits;
    uint32_t lwe_n = 600;

    // bound values
    uint32_t error_width = 7;
    uint32_t alpha = 1 << error_width;
    uint32_t beta = 1 << (error_width - 2);
    uint32_t clear_bits = (N_bits + 1) - error_width;

    // bases
    uint32_t decomp_br_basis = 1 << 12;
    uint32_t decomp_auto_basis = 1 << 4;

    uint32_t sswitch_br_basis = 1 << 10;
    uint32_t sswitch_auto_basis = 1 << 4;
    uint32_t sswitch_squaring_basis = 1 << 8;
    // TODO: revisit scheme switching
    uint32_t sswitch_rgsw_basis = 1 << 2;
    uint32_t ksk_basis = 1 << 4;

    uint32_t skip_n_first_digits = 12;


    BitDecomposerParams decomposerParams{lwe_n, N, small_lwe_modulus,
                                         large_lwe_modulus, large_lwe_modulus_bits, ring_modulus,std,decomp_br_basis,
                                         decomp_auto_basis, clear_bits, error_width};

    SchemeSwitchParameters switchParameters{lwe_n, N, small_lwe_modulus, ring_modulus, std, skip_n_first_digits, sswitch_br_basis, sswitch_squaring_basis, sswitch_auto_basis, sswitch_rgsw_basis};

    auto ksk_params = std::make_shared<lbcrypto::LWECryptoParams>(lwe_n, N, small_lwe_modulus, ring_modulus, key_switch_modulus, std, ksk_basis);

    // generate keys


    /*
    auto genT = TernaryUniformGeneratorImpl<NativeVector>();
    NativeVector huge_key_vec = genT.GenerateVector(large_lwe_n, large_lwe_modulus, beta);
    NativeVector ring_sk_vec = genT.GenerateVector(N, ring_modulus, beta);
    */

    auto genT = BinaryUniformGeneratorImpl<NativeVector>();
    NativeVector huge_key_vec = genT.GenerateVector(large_lwe_n, large_lwe_modulus);
    NativeVector ring_sk_vec = genT.GenerateVector(N, ring_modulus);



    NativeInteger rmod = ring_modulus;
    auto pp = std::make_shared<ILNativeParams>(2 * N, rmod);
    NativePoly sk_out(pp, COEFFICIENT, true);

    sk_out.SetValues(ring_sk_vec, COEFFICIENT);

    auto huge_lwe_key = std::make_shared<LWEPrivateKeyImpl>(huge_key_vec);

    auto small_lwe_key_vec = genT.GenerateVector(lwe_n, small_lwe_modulus);
    auto small_lwe_key = std::make_shared<LWEPrivateKeyImpl>(small_lwe_key_vec);

    // construct lut evaluator
    HugeLUTConfig cfg{decomposerParams, huge_lwe_key, switchParameters, small_lwe_key, sk_out, ksk_params};

    auto evaluator = GiantLutEvaluator(cfg);

    // create ciphertext
    auto pt_space_bits = large_lwe_modulus_bits - error_width;
    auto pt_space = 1ull << pt_space_bits;

    auto message = random() % pt_space;
    NativeInteger value = message * alpha;

    auto ct_in = encrypt_lwe(huge_key_vec, value);

    auto database = GenerateTestDatabase(pt_space, 64);
    auto database_precomp = evaluator.SetupRecordsRaw(database, 4, 2);

    auto start = std::chrono::high_resolution_clock::now();
    auto rlwe_res = evaluator.EvalLUTPolyMulRaw(ct_in, database_precomp, 4, 2);
    auto stop = std::chrono::high_resolution_clock ::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();

    std::cerr << "Took " << elapsed << std::endl;

    auto target_entries = database.at(message);

    auto entry_pt_space = 1 << (4 * 2);

    // check correctness
    sk_out.SetFormat(EVALUATION);
    auto phase = rlwe_res->GetElements()[1] - rlwe_res->GetElements()[0] * sk_out;
    phase.SetFormat(COEFFICIENT);
    std::cerr << " Got " << phase.MultiplyAndRound(entry_pt_space, ring_modulus).Mod(entry_pt_space) << std::endl;
    std::cerr << "Expected : ";
    print_record_as_vec(target_entries, 2 * 4);

}

TEST(GiantLUTTest, TestDBCRT) {

    std::srand(0);
    // set up parameters

    // variance
    double std = 3.19; // for testing

    // moduli
    uint64_t ring_modulus = 36028797018972161;
    uint64_t small_lwe_modulus = 1 << 12;
    uint32_t large_lwe_modulus_bits = 32;
    uint64_t large_lwe_modulus = 1ull << large_lwe_modulus_bits;
    uint64_t key_switch_modulus = 1 << 30;
    uint32_t large_lwe_n = 1304;

    // (r)lwe dimensions
    uint32_t N_bits = std::round(std::log2((long double) small_lwe_modulus)) - 1;
    uint32_t N = 1u << N_bits;
    uint32_t lwe_n = 600;

    // bound values
    uint32_t error_width = 7;
    uint32_t alpha = 1 << error_width;
    uint32_t beta = 1 << (error_width - 2);
    uint32_t clear_bits = (N_bits + 1) - error_width;

    // bases
    uint32_t decomp_br_basis = 1 << 10;
    uint32_t decomp_auto_basis = 1 << 4;

    uint32_t sswitch_br_basis = 1 << 10;
    uint32_t sswitch_auto_basis = 1 << 4;
    uint32_t sswitch_squaring_basis = 1 << 8;
    // TODO: revisit scheme switching
    uint32_t sswitch_rgsw_basis = 1 << 2;
    uint32_t ksk_basis = 1 << 3;

    // crt
    uint64_t P = 268496897;
    uint64_t Q = 268460033;

    // set up parameter sets
    uint32_t skip_n_first_digits = 12;

    BitDecomposerParams decomposerParams{lwe_n, N, small_lwe_modulus,
                                         large_lwe_modulus, large_lwe_modulus_bits, ring_modulus,std,decomp_br_basis,
                                         decomp_auto_basis, clear_bits, error_width};

    SchemeSwitchParameters switchParameters{lwe_n, N, small_lwe_modulus, ring_modulus, std, skip_n_first_digits, sswitch_br_basis, sswitch_squaring_basis, sswitch_auto_basis, sswitch_rgsw_basis};

    auto ksk_params = std::make_shared<lbcrypto::LWECryptoParams>(lwe_n, N, small_lwe_modulus, ring_modulus, key_switch_modulus, std, ksk_basis);
    auto genT = BinaryUniformGeneratorImpl<NativeVector>();
    NativeVector huge_key_vec = genT.GenerateVector(large_lwe_n, large_lwe_modulus);

    // there's no nicer way to set the hamming weight...
    auto one_ctr = 0;
    for (uint32_t i = 0; i < large_lwe_n; i++) {
        if (huge_key_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                huge_key_vec[i] = 0;
            }

        }
    }
    NativeVector ring_sk_vec = genT.GenerateVector(N, ring_modulus);
    one_ctr = 0;
    for (uint32_t i = 0; i < large_lwe_n; i++) {
        if (ring_sk_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                ring_sk_vec[i] = 0;
            }

        }
    }

    NativeInteger rmod = ring_modulus;
    auto pp = std::make_shared<ILNativeParams>(2 * N, rmod);
    NativePoly sk_out(pp, COEFFICIENT, true);
    sk_out.SetValues(ring_sk_vec, COEFFICIENT);

    auto huge_lwe_key = std::make_shared<LWEPrivateKeyImpl>(huge_key_vec);
    auto small_lwe_key_vec = genT.GenerateVector(lwe_n, small_lwe_modulus);
    auto small_lwe_key = std::make_shared<LWEPrivateKeyImpl>(small_lwe_key_vec);

    // construct lut evaluator
    HugeLUTConfig cfg{decomposerParams, huge_lwe_key, switchParameters, small_lwe_key, sk_out, ksk_params};

    auto evaluator = GiantLutEvaluator(cfg);

    // create ciphertext
    auto pt_space_bits = large_lwe_modulus_bits - error_width;
    auto pt_space = 1ull << pt_space_bits;

    auto message = random() % pt_space;
    NativeInteger value = message * alpha;

    auto ct_in = encrypt_lwe(huge_key_vec, value);

    auto database = GenerateTestDatabase(pt_space, 64);
    auto database_precomp = evaluator.SetupRecordsCRT(database, 4, 2, P, Q);

    auto start = std::chrono::high_resolution_clock::now();
    auto rlwe_res = evaluator.QueryCRT(ct_in, database_precomp, 4, 2, P, Q);
    auto stop = std::chrono::high_resolution_clock ::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();

    std::cerr << "Took " << elapsed << std::endl;

    auto target_entries = database.at(message);
    auto entry_pt_space = 1 << (4 * 2);
    // check correctness
    sk_out.SetFormat(EVALUATION);
    auto phase = rlwe_res->GetElements()[1] - rlwe_res->GetElements()[0] * sk_out;
    phase.SetFormat(COEFFICIENT);
    std::cerr << phase << std::endl;
    std::cerr << " Got " << phase.MultiplyAndRound(entry_pt_space, ring_modulus).Mod(entry_pt_space) << std::endl;
    std::cerr << "Expected : ";
    print_record_as_vec(target_entries, 2 * 4);

}