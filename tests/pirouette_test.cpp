//
// Created by leonard on 3/10/25.
//

#include "pirouette.h"
#include <gtest/gtest.h>

TEST(PIROUETTE_TEST, TESTPRIME) {

    auto n_entries_bits = 25;
    auto record_size_bits = 64;

    auto default_config = PirouettePIR::GenerateDefaultConfig(PRIME);
    default_config->SetMetaParams(11, n_entries_bits - 11, 11, n_entries_bits);

    auto db = GenerateTestDatabase(1ull << n_entries_bits, record_size_bits);
    default_config->PrepareDB(n_entries_bits, record_size_bits, db);

    auto query_idx = 1234;
    auto query = default_config->CreateQuery(query_idx, 0);

    std::vector<uint32_t> record = db.at(query_idx);
    db.clear();

    auto rlwe_res = default_config->Query(query);

    auto entry_pt_space = 1 << (4 * 2);
    auto sk_out = default_config->GetOutputKey();
    // check correctness
    sk_out.SetFormat(EVALUATION);
    auto phase = rlwe_res->GetElements()[1] - rlwe_res->GetElements()[0] * sk_out;
    phase.SetFormat(COEFFICIENT);
    std::cerr << phase << std::endl;
    std::cerr << " Got " << phase.MultiplyAndRound(entry_pt_space, phase.GetModulus()).Mod(entry_pt_space) << std::endl;
    std::cerr << "Expected : ";
    print_record_as_vec(record, 2 * 4);

}

TEST(PIROUETTE_TEST, TESTCRT) {

    auto n_entries_bits = 25;
    auto record_size_bits = 64;

    auto default_config = PirouettePIR::GenerateDefaultConfig(CRT);
    default_config->SetMetaParams(11, n_entries_bits - 11, 11, n_entries_bits);

    auto db = GenerateTestDatabase(1ull << n_entries_bits, record_size_bits);
    default_config->PrepareDB(n_entries_bits, record_size_bits, db);

    auto query_idx = 1234;
    auto query = default_config->CreateQuery(query_idx, 0);

    std::vector<uint32_t> record = db.at(query_idx);
    db.clear();

    auto rlwe_res = default_config->Query(query);
    auto entry_pt_space = 1 << (4 * 2);
    auto sk_out = default_config->GetOutputKey();
    // check correctness


    sk_out.SetFormat(EVALUATION);
    auto phase = rlwe_res->GetElements()[1] - rlwe_res->GetElements()[0] * sk_out;
    phase.SetFormat(COEFFICIENT);
    std::cerr << phase << std::endl;
    std::cerr << " Got " << phase.MultiplyAndRound(entry_pt_space, phase.GetModulus()).Mod(entry_pt_space) << std::endl;
    std::cerr << "Expected : ";
    print_record_as_vec(record, 2 * 4);

}

TEST(PIROUETTE_TEST,TEST_DECOMP_INTO_SWITCH) {

    auto n_entries_bits = 25;
    auto record_size_bits = 64;

    srand(time(nullptr));

    auto default_config = PirouettePIR::GenerateDefaultConfig(CRT);
    default_config->SetMetaParams(11, n_entries_bits - 11, 11, n_entries_bits);

    auto db = GenerateTestDatabase(1ull << n_entries_bits, record_size_bits);
    default_config->PrepareDB(n_entries_bits, record_size_bits, db);

    auto query_idx = rand() % (1 << 25);
    auto query = default_config->CreateQuery(query_idx, 0);

    NativeVector A(1800, 1ull << 32);
    NativeInteger B = query.second;
    std::mt19937_64 stream(query.first);
    std::uniform_int_distribution<> distr(0, 1ull << 32);

    for (int i = 0; i < 1800; i++) {
        A[i] = distr(stream);
    }

    auto query2 = std::make_shared<LWECiphertextImpl>(A, B);

    auto bits = default_config->m_digit_decomposer.BitDecompose(query2);
    auto switched = default_config->BuildSelectors(bits);

    std::vector<uint8_t> qbits;
    for (uint32_t i = 0; i < 25; i++) {
        qbits.push_back((query_idx >> i) & 1);
    }

    NativePoly m_sk_ntt = default_config->m_params.m_scheme_switch_rlwe_key;
    m_sk_ntt.SetFormat(EVALUATION);
    auto basis = 1 << 4;
    uint64_t Q = 72057594037641217;

    for (uint32_t i = 0; i < 25; i++) {

        auto bit = qbits[25 - 1 - i];
        auto ct_digits = switched[i];

        NativeInteger v0 = bit == 1 ? 1 << 24 : 0;

        for(auto& ct_i : ct_digits.m_m) {
            auto phase_i = ct_i->GetElements()[1] - ct_i->GetElements()[0] * m_sk_ntt;
            phase_i.SetFormat(COEFFICIENT);
            //std::cerr << phase_i << std::endl;
            EXPECT_EQ(phase_i.at(0).ConvertToInt<uint64_t>(), v0.ConvertToInt<uint64_t>());
            v0.ModMulEq(basis, Q);
        }

        v0 = bit == 1 ? 1 << 24 : 0;

        for(auto& ct_i : ct_digits.m_msm) {
            auto phase_i = ct_i->GetElements()[1] - ct_i->GetElements()[0] * m_sk_ntt;
            auto current_expected_msm = v0 * m_sk_ntt;

            phase_i.SetFormat(COEFFICIENT);
            current_expected_msm.SetFormat(COEFFICIENT);

            EXPECT_TRUE(-phase_i == current_expected_msm);
            v0.ModMulEq(basis, Q);
        }



    }

}

TEST(PIROUETTE_TEST,TEST_DECOMP_INTO_SELECTOR) {

    auto n_entries_bits = 25;
    auto record_size_bits = 64;

    srand(time(nullptr));

    auto default_config = PirouettePIR::GenerateDefaultConfig(CRT);
    default_config->SetMetaParams(11, n_entries_bits - 11, 11, n_entries_bits);

    auto db = GenerateTestDatabase(1ull << n_entries_bits, record_size_bits);
    default_config->PrepareDB(n_entries_bits, record_size_bits, db);

    auto query_idx = rand() % (1 << 25);
    auto query = default_config->CreateQuery(query_idx, 0);

    NativeVector A(1800, 1ull << 32);
    NativeInteger B = query.second;
    std::mt19937_64 stream(query.first);
    std::uniform_int_distribution<> distr(0, 1ull << 32);

    for (int i = 0; i < 1800; i++) {
        A[i] = distr(stream);
    }

    auto query2 = std::make_shared<LWECiphertextImpl>(A, B);

    auto bits = default_config->m_digit_decomposer.BitDecompose(query2);
    auto switched = default_config->BuildSelectors(bits);

    std::vector<uint8_t> qbits;
    for (uint32_t i = 0; i < 25; i++) {
        qbits.push_back((query_idx >> i) & 1);
    }

    NativePoly m_sk_ntt = default_config->m_params.m_scheme_switch_rlwe_key;
    m_sk_ntt.SetFormat(EVALUATION);
    auto basis = 1 << 4;
    uint64_t Q = 72057594037641217;

    auto rlwep_digits = 2;
    auto rlwep_bits = 4;


    std::vector<NativePoly> rlwe_prime_scales;
    for(uint32_t i = 0; i < rlwep_bits; i++) {
        auto pol_i = m_sk_ntt.Clone();
        pol_i.SetValuesToZero();
        pol_i.SetFormat(COEFFICIENT);

        pol_i[0] = Q >> ((rlwep_digits - i) * rlwep_bits);
        pol_i.SetFormat(EVALUATION);
        rlwe_prime_scales.push_back(pol_i);
    }

    auto n0 = m_sk_ntt.Clone();
    n0.SetValuesToZero();
    n0.SetFormat(COEFFICIENT);
    n0[0] = 1;
    n0.SwitchFormat();

    auto kbit_sel = default_config->ConstructKBitSelector(switched, rlwe_prime_scales);
    auto first_block = query_idx >> (25 - 11);

    for (uint32_t i = 0; i < kbit_sel.size(); i++) {
        auto& rlwep = kbit_sel[i];
        auto phase = rlwep[0]->GetElements()[1] - rlwep[0]->GetElements()[0] * m_sk_ntt;
        phase.SetFormat(COEFFICIENT);
        auto flag = phase[0].MultiplyAndRound(1 << (rlwep_bits * rlwep_digits), Q).ConvertToInt<uint64_t>();
        EXPECT_EQ(flag, i == first_block ? 1 : 0);
        std::cerr << phase << std::endl;
    }

}