//
// Created by leonard on 12/10/24.
//

#include <cassert>
#include <algorithm>
#include <chrono>
#include "giant_lut_evaluator.h"
#include <hexl/hexl.hpp>


GiantLutEvaluator::GiantLutEvaluator(HugeLUTConfig& cfg) {

    m_bitdecomp_engine = BitDecomposer(cfg.m_decomp_params);
    m_switch_engine = SchemeSwitchEngine(cfg.m_switch_params);

    NativeVector sk_large_Q_vec = cfg.m_scheme_switch_rlwe_key.GetValues();

    sk_large_Q_vec.SwitchModulus(m_bitdecomp_engine.m_params.m_Q);
    NativePoly sk_large_Q_poly(m_bitdecomp_engine.m_br_params->GetPolyParams(), COEFFICIENT, true);
    sk_large_Q_poly.SetValues(sk_large_Q_vec, COEFFICIENT);

    m_bitdecomp_engine.SetKeys(sk_large_Q_poly, cfg.m_huge_lwe_key);
    m_switch_engine.SetKeys(cfg.m_scheme_switch_rlwe_key,cfg.m_scheme_switch_lwe_key);

    m_lwe_ksk = LWEKeySwitchingKey(cfg.m_ksk_params, sk_large_Q_vec, cfg.m_scheme_switch_lwe_key->GetElement());

    auto N = cfg.m_scheme_switch_rlwe_key.GetLength();
    uint32_t offset = N / 2;
    while (offset >= 1) {
        NativePoly rot_i(m_switch_engine.m_br_params->GetPolyParams(), COEFFICIENT, true);
        rot_i.at(N - offset) = cfg.m_scheme_switch_rlwe_key.GetModulus() - 1;
        rot_i.at(0) += cfg.m_scheme_switch_rlwe_key.GetModulus() - 1;
        rot_i.SetFormat(EVALUATION);
        m_rotation_polys.push_back(rot_i);
        offset >>= 1;
    }

}

std::vector<RingGSWSample> GiantLutEvaluator::BuildSelectors(std::vector<LWECiphertext> &bits) {

    auto Q_in = m_bitdecomp_engine.m_params.m_Q;
    uint64_t Q_ksk = m_lwe_ksk.m_params->GetqKS().ConvertToInt();
    uint64_t q = m_switch_engine.m_params.m_q;

    auto keyks = m_lwe_ksk.m_sk_target;
    keyks.SwitchModulus(q);

    std::reverse(bits.begin(), bits.end());
    std::vector<RingGSWSample> selectors;

    for(auto& bit_i : bits) {
        assert(bit_i->GetModulus() == Q_in);
        // switch to qks
        auto bit_i_ks = ModSwitch(bit_i, Q_ksk);
        auto switched_bit_i = m_lwe_ksk.SwitchKey(bit_i_ks);
        auto bit_i_q = ModSwitch(switched_bit_i, q);

        selectors.push_back(m_switch_engine.SwitchToRGSW(bit_i_q));

    }

    return selectors;
}


RLWECiphertext GiantLutEvaluator::PerformPhase2(std::vector<RLWECiphertext>& current_choices, std::vector<RingGSWSample>& selectors, uint32_t first_idx, uint32_t last_idx) {
    while (current_choices.size() > 1) {
        std::vector<RLWECiphertext> selected_cts;
        auto &selector_i = selectors.at(first_idx);
        auto new_sel = RGSWSample::convert(selector_i);
        auto offset = current_choices.size() >> 1;
        for (long i = 0; i < offset; i++) {
            auto& case_0 = current_choices.at(i);
            auto& case_1 = current_choices.at(i + offset);
            auto op = {case_1->GetElements()[0] - case_0->GetElements()[0],
                       case_1->GetElements()[1] - case_0->GetElements()[1]};

            auto op_ct = std::make_shared<RLWECiphertextImpl>(op);

            auto selected = new_sel.mul(op_ct);

            selected->GetElements()[0] += case_0->GetElements()[0];
            selected->GetElements()[1] += case_0->GetElements()[1];

            selected_cts.push_back(selected);
        }
        if (first_idx >= last_idx) {
            std::cerr << "Exceeding range that should be used during reduction (Not always an error)." << std::endl;
        }
        first_idx++;
        current_choices = std::move(selected_cts);
    }

    return current_choices.at(0);

}

std::vector<std::vector<RLWECiphertext>>
GiantLutEvaluator::BuildRLWEPrimeSelectors(std::vector<RingGSWSample> &selectors, uint32_t log2N,
                                           std::vector<NativePoly> &target_polys) {

    // Note: we assume that the first bit in selectors contains the "most significant bit"

    std::vector<std::vector<RLWECiphertext>> results;
    NativePoly zero(target_polys.at(0).GetParams(), EVALUATION, true);

    std::vector<std::vector<RLWECiphertext>> current_choices;
    std::vector<std::vector<RLWECiphertext>> new_choices;

    std::vector<RLWECiphertext> entry_left;
    std::vector<RLWECiphertext> entry_right;

    // TODO [low priority]: Get rid of me

    std::vector<RGSWSample> new_selectors;
    for(auto& sel : selectors) {
        new_selectors.push_back(RGSWSample::convert(sel));
    }

    for(auto& poly : target_polys) {

        auto rlwe_ct_acc_args = {zero.Clone(), poly.Clone()};
        auto rlwe_ct_acc = std::make_shared<RLWECiphertextImpl>(rlwe_ct_acc_args);
        auto acc0_right = new_selectors.at(0).mul(rlwe_ct_acc);
        auto& acc0_left = rlwe_ct_acc;
        acc0_left->GetElements()[0] -= acc0_right->GetElements()[0];
        acc0_left->GetElements()[1] -= acc0_right->GetElements()[1];

        entry_left.push_back(acc0_left);
        entry_right.push_back(acc0_right);

    }

    current_choices = {entry_left, entry_right};

    for(int i = 1; i < log2N; i++) {
        auto& current_selector = new_selectors.at(i);

        for(auto& entry : current_choices) {

            std::vector<RLWECiphertext> new_left;
            std::vector<RLWECiphertext> new_right;

            for(const auto& rlwect : entry) {
                auto ct_right = current_selector.mul(rlwect);
                rlwect->GetElements()[0] -= ct_right->GetElements()[0];
                rlwect->GetElements()[1] -= ct_right->GetElements()[1];
                new_left.push_back(rlwect);
                new_right.push_back(ct_right);
            }

            new_choices.push_back(new_left);
            new_choices.push_back(new_right);

        }


        current_choices = std::move(new_choices);
    }

    return current_choices;
}

std::vector<uint64_t> GiantLutEvaluator::SetupRLWEPrimeRaw(std::vector<std::vector<RLWECiphertext>> &cts) {

    auto digits = cts.at(0).size();
    auto N = cts.size();

    std::vector<uint64_t> samples_continuous(2 * N * N * digits);

    auto samples_raw = samples_continuous.data();
    auto offset = 0;
    for(auto& rlwe_prime : cts) {
        for(auto& rlwe : rlwe_prime) {
            auto& A = rlwe->GetElements()[0];
            auto& B = rlwe->GetElements()[1];

            for(long k = 0; k < N; k++) {
                samples_raw[offset + k] = A.at(k).ConvertToInt<uint64_t>();
                samples_raw[offset + k + N] = B.at(k).ConvertToInt<uint64_t>();
            }

            offset += 2u * N;
        }
    }

    return samples_continuous;
}


template<uint32_t N> void vec_mul(uint64_t* out, uint64_t* left, uint64_t* right) {
    for(uint32_t idx = 0; idx < N; idx++) {
        out[idx] = left[idx] * right[idx];
    }
}

template<uint32_t N> void vec_mul2(uint64_t* out, uint32_t* left, uint32_t* right) {
    for(uint32_t idx = 0; idx < N; idx++) {
        out[idx] = left[idx] * right[idx];
    }
}


std::vector<RLWECiphertext>
GiantLutEvaluator::PerformPhase1Raw(std::vector<uint64_t> &rlwe_prime, std::vector<uint64_t> &records, uint32_t digits,
                                    uint32_t columns) {

    assert(digits == 1 || digits == 2 || digits == 4 || digits == 8);

    // Load required parameters
    auto params = m_switch_engine.m_br_params;
    auto modulus = params->GetQ().ConvertToInt<uint64_t>();
    auto poly_dim = params->GetN();

    // Size of database column
    uint64_t column_length = poly_dim * poly_dim * digits;

    // Setup buffers
    auto rlwe_prime_data = rlwe_prime.data();

    // buffer in which we store the result(s) of poly-poly products
    std::vector<uint64_t> column_buffer(2 * poly_dim, 0);
    auto col_A = column_buffer.data();
    auto col_B = col_A + poly_dim;

    // temporary storage location
    auto final_column_loc = records.data();
    std::vector<uint64_t> prod_buffer(2 * poly_dim, 0);
    auto prod_buffer_A = prod_buffer.data();
    auto prod_buffer_B = prod_buffer_A + poly_dim;

    std::vector<RLWECiphertext> result;

    // During Phase 1 we need to perform a large amount of modular additions
    // In practice log(modulus) < 64, so we can delay the modular reduction
    // Here we determine the number of additions we can do before it becomes necessary
    auto Qbits = GetMSB(modulus) - 1;
    auto free_bits = 62 - Qbits;
    auto n_additions_before_reduction =  1u << free_bits;

    for(uint64_t kk = 0; kk < columns; kk++) {
        auto current_column = records.data() + kk * column_length;

        intel::hexl::EltwiseMultMod(col_A, current_column, rlwe_prime_data, poly_dim, modulus, 1);
        intel::hexl::EltwiseMultMod(col_B, current_column, rlwe_prime_data + poly_dim, poly_dim, modulus, 1);

        uint32_t ll;
        for(ll = 1; ll < digits * poly_dim; ll++) {
            // A component of current RLWE' sample
            auto rlwe_A = rlwe_prime_data + ll * 2 * poly_dim;
            // B component
            auto rlwe_B = rlwe_A + poly_dim;

            // database entry
            auto poly = current_column + ll * poly_dim;

            intel::hexl::EltwiseMultMod(prod_buffer_A, poly, rlwe_A, poly_dim, modulus, 1);
            intel::hexl::EltwiseMultMod(prod_buffer_B, poly, rlwe_B, poly_dim, modulus, 1);

            vec_add_eq<1u << 11>(col_A, prod_buffer_A);
            vec_add_eq<1u << 11>(col_B, prod_buffer_B);

            if (ll % n_additions_before_reduction == 0) {
                intel::hexl::EltwiseReduceMod(col_A, col_A, poly_dim, modulus, modulus, 1);
                intel::hexl::EltwiseReduceMod(col_B, col_B, poly_dim, modulus, modulus, 1);
            }

        }

        intel::hexl::EltwiseReduceMod(col_A, col_A, poly_dim, modulus, modulus, 1);
        intel::hexl::EltwiseReduceMod(col_B, col_B, poly_dim, modulus, modulus, 1);

        std::copy(col_A, col_A + poly_dim, final_column_loc);
        std::copy(col_B, col_B + poly_dim, final_column_loc + poly_dim);
        final_column_loc += 2 * poly_dim;

    }

    // desired result is now in the columns

    NativeVector A_vec(params->GetN(), params->GetQ());
    NativeVector B_vec(params->GetN(), params->GetQ());

    for(int i = 0; i < columns; i++) {

        auto ct_block = records.data() + i * 2 * poly_dim;
        for(long k = 0; k < poly_dim; k++) {
            A_vec[k] = ct_block[k];
            B_vec[k] = ct_block[k + poly_dim];
        }

        NativePoly A(params->GetPolyParams(), EVALUATION, false);
        NativePoly B(params->GetPolyParams(), EVALUATION, false);

        A.SetValues(A_vec, EVALUATION);
        B.SetValues(B_vec, EVALUATION);

        auto ct_vec = {A, B};
        result.push_back(std::make_shared<RLWECiphertextImpl>(ct_vec));
    }

    return result;
}

std::vector<uint64_t>
GiantLutEvaluator::SetupRecordsRaw(std::vector<std::vector<uint32_t>> &database, uint32_t basebits, uint32_t digits) {

    std::cerr << "[!] Preparing database..." << std::endl;

    // assume that we have *at least* 1 record per poly in phase 2
    // then the number of records per poly is N / (record_slots_needed) = N / (total_record_size_bits / entry_size_bits)

    auto poly_params = m_switch_engine.m_br_params->GetPolyParams();
    uint32_t rgsw_N = m_switch_engine.m_params.m_N;

    auto total_record_size_bits = database.at(0).size() * sizeof(uint32_t) * 8;
    auto entry_size_bits = basebits * digits;
    auto records_per_polynomial = (rgsw_N * entry_size_bits) / total_record_size_bits;
    auto rounds_per_entry = (sizeof(uint32_t) * 8) / entry_size_bits;

    auto dim_1_step_size = database.size() / rgsw_N;
    auto dim_2_step_size = dim_1_step_size / records_per_polynomial;

    std::vector<uint64_t> buffer(rgsw_N, 0);
    auto mask = (1u << entry_size_bits) - 1;
    auto mask_small = (1u << basebits)  -1;

    std::vector<uint64_t> records_continuous(digits * rgsw_N * rgsw_N * dim_2_step_size);
    auto record_ptr = records_continuous.data();
    auto memory_large_step_size = rgsw_N * rgsw_N * digits;

    for(uint64_t block_idx = 0 ; block_idx < dim_2_step_size; block_idx++) {
        auto current_memory_block = record_ptr + block_idx * memory_large_step_size;

        for(uint64_t i = 0; i < rgsw_N; i++) {
            auto current_rlwe_prime_in_memory = current_memory_block + i * rgsw_N * digits;
            auto database_index = i * dim_1_step_size + block_idx * records_per_polynomial;

            auto buffer_index = 0;
            for(uint64_t k = 0; k < records_per_polynomial; k++) {
                auto full_index = database_index + k;
                auto& entry = database[full_index];
                for(auto& v : entry) {
                    for(uint32_t kk = 0; kk < rounds_per_entry; kk++) {
                        buffer[buffer_index++] = (v >> (kk * entry_size_bits)) & mask;
                    }

                }
            }

            // next, we decompose w.r.t a basis
            std::vector<NativePoly> dim2_col(digits);
            for(uint32_t k = 0; k < digits; k++) {
                auto digit_location = current_rlwe_prime_in_memory + k * rgsw_N;

                auto poly_entry = NativePoly(poly_params, COEFFICIENT, true);
                for(long kk = 0 ; kk < rgsw_N; kk++) {
                    poly_entry[kk] = buffer[kk] & mask_small;
                    buffer[kk] >>= basebits;
                }
                poly_entry.SetFormat(EVALUATION);
                // copy polynomial

                for(long kk = 0; kk < rgsw_N; kk++) {
                    digit_location[kk] = poly_entry.at(kk).ConvertToInt<uint64_t>();
                }
                //std::memcpy(digit_location + rgsw_N, digit_location, rgsw_N * sizeof(uint64_t));
            }
        }

    }

    return records_continuous;
}


RLWECiphertext
GiantLutEvaluator::EvalLUTPolyMulRaw(lbcrypto::ConstLWECiphertext &ct_idx, std::vector<uint64_t> &records,
                                     uint32_t basebits, uint32_t digits) {

    auto sk = m_switch_engine.m_sk_rlwe_ntt;

    auto start = std::chrono::high_resolution_clock::now();

    // LWE sample to bits
    auto bits = m_bitdecomp_engine.BitDecompose(ct_idx);
    std::cerr << "[!] Finished decomposition " << std::endl;

    // convert LWE bits into RGSW bits
    auto selectors = BuildSelectors(bits);
    auto stop = std::chrono::high_resolution_clock::now();

    std::cerr << "[!] Finished building selectors" << std::endl;

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << "[!] BitDecomposition & LWE-RGSW conversion took " << elapsed << std::endl;

    auto test = NativePoly(m_switch_engine.m_br_params->GetPolyParams(), COEFFICIENT, true);
    test.at(0) = m_switch_engine.m_params.m_Q >> 8;
    test.SetFormat(EVALUATION);
    auto zero = NativePoly(test.GetParams(), EVALUATION, true);
    auto tt_vec = {zero, test};
    auto tt = std::make_shared<RLWECiphertextImpl>(tt_vec);

    uint64_t rgsw_N = m_switch_engine.m_params.m_N;
    uint32_t log2N = GetMSB(rgsw_N) - 1;

    auto rgsw_Q = m_switch_engine.m_params.m_Q;
    auto poly_params = m_switch_engine.m_br_params->GetPolyParams();

    std::vector<NativePoly> rlwe_prime_scales;
    for(uint32_t i = 0; i < digits; i++) {
        auto pol_i = NativePoly(poly_params, COEFFICIENT, true);
        pol_i[0] = rgsw_Q >> ((digits - i) * basebits);
        pol_i.SetFormat(EVALUATION);
        rlwe_prime_scales.push_back(pol_i);
    }

    auto start_sel = std::chrono::high_resolution_clock::now();

    // Here we build C_i = RLWE'(v_i) v_i = 1 for exactly one i, v_i = 0 else
    auto rlwe_prime_selectors = BuildRLWEPrimeSelectors(selectors, log2N, rlwe_prime_scales);

    auto stop_sel = std::chrono::high_resolution_clock::now();
    auto elapsed_sel = std::chrono::duration_cast<std::chrono::milliseconds>(stop_sel-start_sel).count();
    std::cerr << "[!] RLWE' selector construction took " << elapsed_sel << std::endl;

    start_sel = std::chrono::high_resolution_clock::now();

    uint64_t rs = records.size();
    uint64_t n_cols = rs / (digits * rgsw_N * rgsw_N);

    auto prepared_samples = SetupRLWEPrimeRaw(rlwe_prime_selectors);

    // perform phase 1 i.e. \sum_i C_i p_i, where p_i encodes the database
    auto current_choices = PerformPhase1Raw(prepared_samples, records, digits, n_cols);

    stop_sel = std::chrono::high_resolution_clock::now();
    elapsed_sel = std::chrono::duration_cast<std::chrono::milliseconds>(stop_sel-start_sel).count();
    std::cerr << "[!] Phase 1 took " << elapsed_sel << std::endl;

    // we can now proceed with the remaining stuff
    auto choice_size_log2 = GetMSB(current_choices.size()) - 1;
    auto selector_idx = log2N;

    auto start_phase_2 = std::chrono::high_resolution_clock::now();

    // Perform phase 2 i.e. SELECT
    auto last_ct = PerformPhase2(current_choices, selectors, selector_idx, selector_idx + choice_size_log2);
    auto stop_phase_2 = std::chrono::high_resolution_clock::now();
    auto phase_2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop_phase_2-start_phase_2).count();
    std::cerr << "[!] Phase 2 took " << phase_2 << std::endl;

    selector_idx += choice_size_log2;

    auto rot_offset = rgsw_N >> 1;

    auto start_phase_3 = std::chrono::high_resolution_clock::now();

    // below, we perform phase 3
    for(long i = selector_idx; i < selectors.size(); i++) {

        auto selector_i = RGSWSample::convert(selectors.at(i));
        selector_i.cmux(last_ct, 2 * rgsw_N - rot_offset);
        rot_offset >>= 1;

    }

    auto stop_phase_3 = std::chrono::high_resolution_clock::now();
    auto phase_3 = std::chrono::duration_cast<std::chrono::milliseconds>(stop_phase_3-start_phase_3).count();
    std::cerr << "[!] Phase 3 took " << phase_3 << std::endl;

    // Here comes output/response compression

    return last_ct;
}

void print_poly_signed(NativePoly& poly) {
    auto mod = poly.GetModulus().ConvertToInt<int64_t>();
    auto len = poly.GetLength();
    std::cerr << (poly.GetFormat() == COEFFICIENT ? "COEF" : "EVAL") << " [ ";
    for (long i = 0; i < len; i++) {
        auto v_i = poly[i].ConvertToInt<int64_t>();
        v_i = v_i > (mod / 2) ? mod - v_i : v_i;
        std::cerr << v_i;
        if (i != (len - 1))
            std::cerr << " ";
    }
    std::cerr << " ]" << std::endl;
}


std::vector<uint32_t> GiantLutEvaluator::SetupRLWEPrimeCRT(std::vector<std::vector<RLWECiphertext> > &cts, uint64_t P, uint64_t Q) {
    auto digits = cts.at(0).size();
    auto N = cts.at(0).at(0)->GetElements()[0].GetLength();

    std::vector<uint32_t> samples_continuous(2 * 2 * N * N * digits);
    std::vector<uint64_t> hexl_in_buffer(6 * N);

    auto ntt_engine_P = intel::hexl::NTT(N, P);
    auto ntt_engine_Q = intel::hexl::NTT(N, Q);

    uint64_t offset = 0;

    uint64_t PQ = P * Q;
    auto old_mod = cts[0][0]->GetElements()[0].GetModulus();

    for(auto& rlwe_prime : cts) {
        for(auto& rlwe : rlwe_prime) {
            auto& A = rlwe->GetElements()[0];
            auto& B = rlwe->GetElements()[1];

            A.SetFormat(COEFFICIENT);
            B.SetFormat(COEFFICIENT);

            auto A_vec = A.GetValues();
            auto B_vec = B.GetValues();

            A_vec.SetModulus(PQ);
            B_vec.SetModulus(PQ);

            A_vec.MultiplyAndRoundEq(PQ, old_mod).ModEq(PQ);
            B_vec.MultiplyAndRoundEq(PQ, old_mod).ModEq(PQ);

            // hexl buffer [A || B || 0 || 0 || 0 || 0]
            for(long k = 0; k < N; k++) {
                hexl_in_buffer[k] = A_vec.at(k).ConvertToInt<uint64_t>();
                hexl_in_buffer[k + N] = B_vec.at(k).ConvertToInt<uint64_t>();
            }

            // revisit in case of bugs
            // hexl buffer [A || B || Ap || Bp || Aq || Bq]
            intel::hexl::EltwiseReduceMod(hexl_in_buffer.data() + 2 * N, hexl_in_buffer.data(), 2 * N, P, P, 1);
            intel::hexl::EltwiseReduceMod(hexl_in_buffer.data() + 4 * N, hexl_in_buffer.data(), 2 * N, Q, Q, 1);

            ntt_engine_P.ComputeForward(hexl_in_buffer.data() + 2 * N, hexl_in_buffer.data() + 2 * N, 1, 1);
            ntt_engine_P.ComputeForward(hexl_in_buffer.data() + 3 * N, hexl_in_buffer.data() + 3 * N, 1, 1);
            ntt_engine_Q.ComputeForward(hexl_in_buffer.data() + 4 * N, hexl_in_buffer.data() + 4 * N, 1, 1);
            ntt_engine_Q.ComputeForward(hexl_in_buffer.data() + 5 * N, hexl_in_buffer.data() + 5 * N, 1, 1);

            // should be correct but maybe revisit
            //std::copy(hexl_in_buffer.data() + 2 * N, hexl_in_buffer.data() + 6 * N, samples_continuous.data() + offset);
            for (uint32_t k = 0; k < 4 * N; k++) {
                samples_continuous[offset + k] = hexl_in_buffer[k + 2 * N];
            }

            offset += 4u * N;
        }
    }

    return samples_continuous;
}

std::vector<uint32_t> GiantLutEvaluator::SetupRecordsCRT(std::vector<std::vector<uint32_t> > &database, uint32_t basebits, uint32_t digits, uint32_t P, uint32_t Q) {
    std::cerr << "[!] Preparing database..." << std::endl;

    // assume that we have *at least* 1 record per poly in phase 2
    // then the number of records per poly is N / (record_slots_needed) = N / (total_record_size_bits / entry_size_bits)

    auto poly_params = m_switch_engine.m_br_params->GetPolyParams();
    uint32_t rgsw_N = m_switch_engine.m_params.m_N;

    auto total_record_size_bits = database.at(0).size() * sizeof(uint32_t) * 8;
    auto entry_size_bits = basebits * digits;
    auto records_per_polynomial = (rgsw_N * entry_size_bits) / total_record_size_bits;
    auto rounds_per_entry = (sizeof(uint32_t) * 8) / entry_size_bits;

    auto dim_1_step_size = database.size() / rgsw_N;
    auto dim_2_step_size = dim_1_step_size / records_per_polynomial;

    std::vector<uint64_t> buffer(rgsw_N, 0);
    auto mask = (1u << entry_size_bits) - 1;
    auto mask_small = (1u << basebits)  -1;

    std::vector<uint32_t> records_continuous(2 * digits * rgsw_N * rgsw_N * dim_2_step_size);
    auto record_ptr = records_continuous.data();
    auto memory_large_step_size = 2 * rgsw_N * rgsw_N * digits;


    std::vector<uint64_t> hexl_buffer(4 * rgsw_N);
    auto ntt_engine_P = intel::hexl::NTT(rgsw_N, P);
    auto ntt_engine_Q = intel::hexl::NTT(rgsw_N, Q);

    for(uint64_t block_idx = 0 ; block_idx < dim_2_step_size; block_idx++) {
        auto current_memory_block = record_ptr + block_idx * memory_large_step_size;

        for(uint64_t i = 0; i < rgsw_N; i++) {
            auto current_rlwe_prime_in_memory = current_memory_block + i * 2 * rgsw_N * digits;
            auto database_index = i * dim_1_step_size + block_idx * records_per_polynomial;

            auto buffer_index = 0;
            for(uint64_t k = 0; k < records_per_polynomial; k++) {
                auto full_index = database_index + k;
                auto& entry = database[full_index];
                for(auto& v : entry) {
                    for(uint32_t kk = 0; kk < rounds_per_entry; kk++) {
                        buffer[buffer_index++] = (v >> (kk * entry_size_bits)) & mask;
                    }

                }
            }

            // next, we decompose w.r.t a basis
            //std::vector<NativePoly> dim2_col(digits);
            for(uint32_t k = 0; k < digits; k++) {
                auto digit_location = current_rlwe_prime_in_memory + k * 2 * rgsw_N;
                for(long kk = 0 ; kk < rgsw_N; kk++) {
                    hexl_buffer[kk] = buffer[kk] & mask_small;
                    buffer[kk] >>= basebits;
                }
                // TODO [in case of bugs] revisit input/output mod factor
                // In theory, an additional modular reduction would be necessary, but in practice
                // the record entries will be so small be don't care
                ntt_engine_P.ComputeForward(hexl_buffer.data() + rgsw_N, hexl_buffer.data(), 1, 1);
                ntt_engine_Q.ComputeForward(hexl_buffer.data() + rgsw_N * 2, hexl_buffer.data(), 1, 1);

                std::copy(hexl_buffer.data() + rgsw_N, hexl_buffer.data() + rgsw_N * 3, digit_location);
            }
        }

    }

    return records_continuous;
}

std::vector<RLWECiphertext> GiantLutEvaluator::PerformPhase1CRT(std::vector<uint32_t> &rlwe_prime, std::vector<uint32_t> &records, uint32_t digits, uint32_t columns, uint64_t P, uint64_t Q) {

    auto poly_params = m_switch_engine.m_br_params->GetPolyParams();
    auto output_modulus = poly_params->GetModulus();
    auto P_bits = GetMSB(P) - 1;
    auto Q_bits = GetMSB(Q) - 1;
    auto poly_dim = m_switch_engine.m_params.m_N;

    /* crt params */
    auto alpha = intel::hexl::InverseMod(P % Q, Q);
    alpha = intel::hexl::MultiplyMod(alpha, P, P * Q);

    auto beta = intel::hexl::InverseMod(Q % P, P);
    beta = intel::hexl::MultiplyMod(beta, Q, P * Q);

    assert(poly_dim == (1u << 11));

    auto ntt_engine_P = intel::hexl::NTT(poly_dim, P);
    auto ntt_engine_Q = intel::hexl::NTT(poly_dim, Q);

    auto mod_p_product_bits = 2 * P_bits - 1;
    auto mod_q_product_bits = 2 * Q_bits - 1;
    auto PQ = P * Q;

    auto number_of_additions_before_reduction = 1ull << (63 - std::max(mod_p_product_bits, mod_q_product_bits));
    assert(digits == 1 || digits == 2 || digits == 4 || digits == 8);

    std::array<uint64_t, 6 * (1u << 11)> column_buffer;

    auto column_buffer_AP = column_buffer.data();
    auto column_buffer_BP = column_buffer_AP + poly_dim;
    auto column_buffer_AQ = column_buffer_BP + poly_dim;
    auto column_buffer_BQ = column_buffer_AQ + poly_dim;

    auto rlwe_prime_data = rlwe_prime.data();

    std::vector<RLWECiphertext> results;
    for (uint32_t col_idx = 0; col_idx < columns; col_idx++) {
        auto record_column_ptr = records.data() + col_idx * poly_dim * digits * 2 * poly_dim;

        std::fill(column_buffer.begin(), column_buffer.end(), 0);

        auto acc_ctr = 0;
        for (uint32_t selector_idx = 0; selector_idx < poly_dim; selector_idx++) {
            auto rlwe_block_ptr = rlwe_prime_data + selector_idx * 2 * digits * 2 * poly_dim;

            auto record_block_ptr = record_column_ptr + selector_idx * 2 * digits * poly_dim;

            for (uint32_t digit_idx = 0; digit_idx < digits; digit_idx++) {

                auto record_digit_P = record_block_ptr + digit_idx * 2 * poly_dim;
                auto record_digit_Q = record_digit_P + poly_dim;

                auto rlwe_AP = rlwe_block_ptr + digit_idx * 4 * poly_dim;
                auto rlwe_BP = rlwe_AP + poly_dim;
                auto rlwe_AQ = rlwe_BP + poly_dim;
                auto rlwe_BQ = rlwe_AQ + poly_dim;

                vector_accumulate<1u << 11>(column_buffer_AP, rlwe_AP, record_digit_P);
                vector_accumulate<1u << 11>(column_buffer_BP, rlwe_BP, record_digit_P);
                vector_accumulate<1u << 11>(column_buffer_AQ, rlwe_AQ, record_digit_Q);
                vector_accumulate<1u << 11>(column_buffer_BQ, rlwe_BQ, record_digit_Q);

                acc_ctr++;

            }

            if ((acc_ctr >= number_of_additions_before_reduction) or (selector_idx == poly_dim - 1)) {
                // Note, usually we'd reduce mod P or mod Q once necessary
                // BUT, the results of the coefficient wise product are approx. as big as P * Q
                // So reducing mod P (resp. Q) is pointless
                // and instead of reducing for P and Q, we just do it once for P * Q

                intel::hexl::EltwiseReduceMod(column_buffer.data(), column_buffer.data(), 4 * poly_dim, PQ, PQ, 1);
                acc_ctr = 0;
            }
        }

        // we are done with the column, so we reconstruct
        // note: we need to do the correct modular reduction before the intt since the roots of unity will be different in general
        intel::hexl::EltwiseReduceMod(column_buffer.data(), column_buffer.data(), 2 * poly_dim, P, P, 1);
        intel::hexl::EltwiseReduceMod(column_buffer.data() + 2 * poly_dim, column_buffer.data() + 2 * poly_dim, 2 * poly_dim, Q, Q, 1);

        ntt_engine_P.ComputeInverse(column_buffer_AP, column_buffer_AP, 1, 1);
        ntt_engine_P.ComputeInverse(column_buffer_BP, column_buffer_BP, 1, 1);
        ntt_engine_Q.ComputeInverse(column_buffer_AQ, column_buffer_AQ, 1, 1);
        ntt_engine_Q.ComputeInverse(column_buffer_BQ, column_buffer_BQ, 1, 1);

        NativeVector A_PQ(poly_dim, P * Q);
        NativeVector B_PQ(poly_dim, P * Q);

        std::vector<uint64_t> acc(2 * poly_dim, 0);

        intel::hexl::EltwiseFMAMod(acc.data(), column_buffer_AP, beta, nullptr, poly_dim, P * Q, 1);
        intel::hexl::EltwiseFMAMod(acc.data() + poly_dim, column_buffer_BP, beta, nullptr, poly_dim, P * Q, 1);
        intel::hexl::EltwiseFMAMod(acc.data(), column_buffer_AQ, alpha, acc.data(), poly_dim, P * Q, 1);
        intel::hexl::EltwiseFMAMod(acc.data() + poly_dim, column_buffer_BQ, alpha, acc.data() + poly_dim,poly_dim, P * Q, 1);

        for(uint32_t i = 0; i < poly_dim; i++) {
            A_PQ[i] = acc[i];
            B_PQ[i] = acc[i + poly_dim];
        }

        A_PQ = A_PQ.MultiplyAndRound(output_modulus, P * Q).Mod(output_modulus);
        B_PQ = B_PQ.MultiplyAndRound(output_modulus, P * Q).Mod(output_modulus);
        A_PQ.SetModulus(output_modulus);
        B_PQ.SetModulus(output_modulus);

        auto a_poly = NativePoly(poly_params, COEFFICIENT, true);
        a_poly.SetValues(A_PQ, COEFFICIENT);

        auto b_poly = NativePoly(poly_params, COEFFICIENT, true);
        b_poly.SetValues(B_PQ, COEFFICIENT);

        a_poly.SwitchFormat();
        b_poly.SwitchFormat();

        auto p_vec = {a_poly, b_poly};

        results.push_back(std::make_shared<RLWECiphertextImpl>(p_vec));
    }

    return results;
}

RLWECiphertext GiantLutEvaluator::QueryCRT(lbcrypto::ConstLWECiphertext &ct_idx, std::vector<uint32_t> &records, uint32_t basebits, uint32_t digits, uint32_t P, uint32_t Q) {

    // TODO: There's a bug in the CRT stuff. Find it.
    auto sk = m_switch_engine.m_sk_rlwe_ntt;
    auto start = std::chrono::high_resolution_clock::now();

    // LWE sample to bits
    auto bits = m_bitdecomp_engine.BitDecompose(ct_idx);
    std::cerr << "[!] Finished decomposition " << std::endl;

    // convert LWE bits into RGSW bits
    auto selectors = BuildSelectors(bits);
    auto stop = std::chrono::high_resolution_clock::now();

    std::cerr << "[!] Finished building selectors" << std::endl;

    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << "[!] BitDecomposition & LWE-RGSW conversion took " << elapsed << std::endl;

    auto test = NativePoly(m_switch_engine.m_br_params->GetPolyParams(), COEFFICIENT, true);
    test.at(0) = m_switch_engine.m_params.m_Q >> 8;
    test.SetFormat(EVALUATION);
    auto zero = NativePoly(test.GetParams(), EVALUATION, true);
    auto tt_vec = {zero, test};
    auto tt = std::make_shared<RLWECiphertextImpl>(tt_vec);

    uint64_t rgsw_N = m_switch_engine.m_params.m_N;
    uint32_t log2N = GetMSB(rgsw_N) - 1;

    auto rgsw_Q = m_switch_engine.m_params.m_Q;
    auto poly_params = m_switch_engine.m_br_params->GetPolyParams();

    std::vector<NativePoly> rlwe_prime_scales;
    for(uint32_t i = 0; i < digits; i++) {
        auto pol_i = NativePoly(poly_params, COEFFICIENT, true);
        pol_i[0] = rgsw_Q >> ((digits - i) * basebits);
        pol_i.SetFormat(EVALUATION);
        rlwe_prime_scales.push_back(pol_i);
    }

    auto start_sel = std::chrono::high_resolution_clock::now();

    // Here we build C_i = RLWE'(v_i) v_i = 1 for exactly one i, v_i = 0 else
    auto rlwe_prime_selectors = BuildRLWEPrimeSelectors(selectors, log2N, rlwe_prime_scales);
    auto stop_sel = std::chrono::high_resolution_clock::now();
    auto elapsed_sel = std::chrono::duration_cast<std::chrono::milliseconds>(stop_sel-start_sel).count();
    std::cerr << "[!] RLWE' selector construction took " << elapsed_sel << std::endl;

    start_sel = std::chrono::high_resolution_clock::now();

    uint64_t rs = records.size();
    uint64_t n_cols = rs / (2 * digits * rgsw_N * rgsw_N);

    auto prepared_samples = SetupRLWEPrimeCRT(rlwe_prime_selectors, P, Q);

    // perform phase 1 i.e. \sum_i C_i p_i, where p_i encodes the database
    auto current_choices = PerformPhase1CRT(prepared_samples, records, digits, n_cols, P, Q);


    stop_sel = std::chrono::high_resolution_clock::now();
    elapsed_sel = std::chrono::duration_cast<std::chrono::milliseconds>(stop_sel-start_sel).count();
    std::cerr << "[!] Phase 1 took " << elapsed_sel << std::endl;

    // we can now proceed with the remaining stuff
    auto choice_size_log2 = GetMSB(current_choices.size()) - 1;
    auto selector_idx = log2N;

    auto start_phase_2 = std::chrono::high_resolution_clock::now();

    // Perform phase 2 i.e. SELECT
    auto last_ct = PerformPhase2(current_choices, selectors, selector_idx, selector_idx + choice_size_log2);
    auto stop_phase_2 = std::chrono::high_resolution_clock::now();
    auto phase_2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop_phase_2-start_phase_2).count();
    std::cerr << "[!] Phase 2 took " << phase_2 << std::endl;

    selector_idx += choice_size_log2;

    auto rot_offset = rgsw_N >> 1;

    auto start_phase_3 = std::chrono::high_resolution_clock::now();

    // below, we perform phase 3
    for(long i = selector_idx; i < selectors.size(); i++) {

        auto selector_i = RGSWSample::convert(selectors.at(i));
        selector_i.cmux(last_ct, 2 * rgsw_N - rot_offset);
        rot_offset >>= 1;
    }

    auto stop_phase_3 = std::chrono::high_resolution_clock::now();
    auto phase_3 = std::chrono::duration_cast<std::chrono::milliseconds>(stop_phase_3-start_phase_3).count();
    std::cerr << "[!] Phase 3 took " << phase_3 << std::endl;

    // Here comes output/response compression

    return last_ct;
}
