//
// Created by lschild on 19/11/24.
//

#include "digit_decomposer.h"

#include <cassert>

BitDecomposer::BitDecomposer(BitDecomposerParams &params) : m_params(params) {
    std::cerr << __FILE__ << "BitDecomposer has been depreciated !!!" << std::endl;
    m_br_params = std::make_shared<RingGSWCryptoParams>(m_params.m_N, m_params.m_Q, m_params.m_q, m_params.m_br_basis, 1, BINFHE_METHOD::GINX, 3.19);
}

void BitDecomposer::SetKeys(lbcrypto::NativePoly sk_rlwe, lbcrypto::LWEPrivateKey &sk_lwe) {
    auto m_auto_params = std::make_shared<RingGSWCryptoParams>(m_params.m_N, m_params.m_Q, m_params.m_q, m_params.m_auto_key_basis, 1, BINFHE_METHOD::GINX, m_params.m_std);

    sk_rlwe.SetFormat(COEFFICIENT);
    m_sk_rlwe = sk_rlwe;
    m_sk_rlwe_ntt = m_sk_rlwe.Clone();
    m_sk_rlwe_ntt.SetFormat(EVALUATION);
    m_sk_lwe = sk_lwe->GetElement();

    m_engine = FunctionalBootstrapEngine(m_br_params,m_auto_params,  sk_rlwe, sk_lwe->GetElement());

    m_qks = 1ull << 42;
    auto ksk_params = std::make_shared<lbcrypto::LWECryptoParams>(sk_lwe->GetLength(), 512, 1024, m_qks, m_qks, 3.19, 16);
    m_lwe_ksk = LWEKeySwitchingKey(ksk_params, m_sk_rlwe.GetValues(), sk_lwe->GetElement());

}

std::vector<LWECiphertext> BitDecomposer::TightDecompose(lbcrypto::ConstLWECiphertext &ct) {

    auto current_ct_modulus = m_params.m_q_large;
    auto digit_modulus = m_params.m_q;

    auto ct_clone = std::make_shared<LWECiphertextImpl>(ct->GetA(), ct->GetB());

    std::vector<LWECiphertext> results;

    auto alpha = 1u << m_params.m_error_width;

    while (current_ct_modulus > digit_modulus) {

        // Get the digit
        auto digit_A = ct_clone->GetA().Mod(digit_modulus);
        digit_A.SetModulus(digit_modulus);
        auto digit_B = ct_clone->GetB().Mod(digit_modulus);
        auto digit_ct = std::make_shared<LWECiphertextImpl>(digit_A, digit_B);

        results.push_back(digit_ct);

        // clear LSBs of input CT
        auto AA = ct_clone->GetA();
        auto BB = ct_clone->GetB();
        auto floor_res = DoFloorAlpha(current_ct_modulus, digit_modulus, AA, BB);

        // Now we can modulus switch
        auto floored = floor_res.first;

        auto new_A = floored->GetA();
        auto new_B = floored->GetB();

        auto new_ct_modulus = (alpha * current_ct_modulus) / digit_modulus;

        if (new_ct_modulus < digit_modulus) {
            new_ct_modulus = digit_modulus;
        }

        new_A.MultiplyAndRoundEq(new_ct_modulus, current_ct_modulus).Mod(new_ct_modulus);
        new_B.MultiplyAndRoundEq(new_ct_modulus, current_ct_modulus).Mod(new_ct_modulus);

        new_A.ModEq(new_ct_modulus);
        new_B.ModEq(new_ct_modulus);

        new_A.SetModulus(new_ct_modulus);

        current_ct_modulus = new_ct_modulus;

        /*
        current_sk = m_sk_lwe;
        current_sk.SwitchModulus(current_ct_modulus);

        current_phase = new_B.ModSub(DotProduct(current_sk, new_A), current_ct_modulus);
        std::cerr << "Floored phase is " << current_phase << " After MS" << std::endl; */

        ct_clone = std::make_shared<LWECiphertextImpl>(new_A, new_B);
    }

    results.push_back(ct_clone);

    return results;
}

void log(NativeVector& A, NativeInteger& B, NativeVector& sk) {
    NativeVector skQQ = sk;
    skQQ.SwitchModulus(A.GetModulus());
    std::cerr << B.ModSub(DotProduct(skQQ, A), A.GetModulus()) << " ";
}


std::pair<LWECiphertext, std::vector<LWECiphertext>> BitDecomposer::DoFloorAlpha(uint64_t Q, uint64_t q, NativeVector &A, NativeInteger &B) {

    auto blind_rotate_Q = m_params.m_Q;
    // - 2 as we assume that alpha >= 4 * beta
    auto beta = 1u << (m_params.m_error_width - 2);
    //auto sk_Q = m_sk_lwe;

    B.ModAddEq(beta, Q);

    // TLDR, the authors of LMP21 casually forgot that when they write -q/4 for f0 that the minus applies w.r.t Q, not q, not ring Q...
    BootFunction F0 = [Q](long m, long q, long ring_Q) {return (m < q/2 ? Q - q/4 : q/4);};
    BootFunction F1 = [](long m, long q, long ring_Q) {return (m < q/2 ? m : (q + q/2 - m) % q);};

    auto Aq = A.Mod(q);
    auto Bq = B.Mod(q);

    Aq.SetModulus(q);

    auto ctq1 = std::make_shared<LWECiphertextImpl>(Aq, Bq);

    //// clear msb
    auto f0_res = m_engine.Boot(ctq1, F0, 1, Q);

    auto a_before_ksk = f0_res->GetA().MultiplyAndRound(m_qks, blind_rotate_Q).Mod(m_qks);
    auto b_before_ksk = f0_res->GetB().MultiplyAndRound(m_qks, blind_rotate_Q).Mod(m_qks);
    a_before_ksk.SetModulus(m_qks);

    auto ct_after_ks = m_lwe_ksk.SwitchKey(a_before_ksk, b_before_ksk);

    auto A0 = ct_after_ks->GetA().MultiplyAndRound(Q, m_qks).Mod(Q);
    auto B0 = ct_after_ks->GetB().MultiplyAndRound(Q, m_qks).Mod(Q);
    A0.SetModulus(Q);

    // Phase(f0_res) is \pm q / 4, adding q / 4 makes it 0 or q / 2

    B0.ModAddEq(q / 4, Q);

    // Subtract to clear msb
    A.ModSubEq(A0);
    B.ModSubEq(B0, Q);

    /*
    std::cerr << "MSB clear | Phase is = ";
    log(A, B, sk_Q);
    std::cerr << std::endl;
    */

    // add beta again to make error positive again
    B.ModAddEq(beta, Q);

    // next create new corrector term with using ct with known MSB
    Aq = A.Mod(q);
    Bq = B.Mod(q);
    Aq.SetModulus(q);
    ctq1 = std::make_shared<LWECiphertextImpl>(Aq, Bq);

    auto f1_res = m_engine.Boot(ctq1, F1, 1, Q);

    a_before_ksk = f1_res->GetA().MultiplyAndRound(m_qks, blind_rotate_Q).Mod(m_qks);
    b_before_ksk = f1_res->GetB().MultiplyAndRound(m_qks, blind_rotate_Q).Mod(m_qks);
    a_before_ksk.SetModulus(m_qks);

    ct_after_ks = m_lwe_ksk.SwitchKey(a_before_ksk, b_before_ksk);

    auto A1 = ct_after_ks->GetA().MultiplyAndRound(Q, m_qks).Mod(Q);
    auto B1 = ct_after_ks->GetB().MultiplyAndRound(Q, m_qks).Mod(Q);
    A1.SetModulus(Q);

    // clear lsbs of input LWE
    A.ModSubEq(A1);
    B.ModSubEq(B1, Q);

    /*
    std::cerr << "Corrected | ";
    log(A, B, m_sk_lwe);
    std::cerr << std::endl;
    */

    auto res_vec = std::vector<LWECiphertext>();
    auto res_ct = std::make_shared<LWECiphertextImpl>(A, B);

    return std::make_pair(res_ct, res_vec);
}

std::vector<LWECiphertext> BitDecomposer::BitDecompose(lbcrypto::ConstLWECiphertext &ct) {

    auto large_digits = TightDecompose(ct);
    auto bits_per_digit = m_params.m_clear_bits;

    std::vector<LWECiphertext> bits;
    auto pt_space_bits = m_params.m_q_large_bits - m_params.m_error_width;
    auto clear_bits = m_params.m_clear_bits;
    auto bits_rem = pt_space_bits % clear_bits;
    auto last_digit_offset = (clear_bits - bits_rem) % clear_bits;

    int idx = 0;
    //auto sk = m_sk_lwe;
    //sk.SwitchModulus(large_digits.at(0)->GetModulus());
    for(auto & ct_i : large_digits) {

        //auto phase_i = ct_i->GetB().ModSub(DotProduct(sk, ct_i->GetA()), ct_i->GetA().GetModulus());
        //std::cerr << idx << "| " << phase_i << std::endl;

        auto ct_i_bits = m_engine.BitDecompose(ct_i, bits_per_digit);

        if  ((idx != large_digits.size() - 1) or last_digit_offset == 0) {
            bits.insert(bits.end(), ct_i_bits.begin(), ct_i_bits.end());
        } else {
            bits.insert(bits.end(), ct_i_bits.begin() + last_digit_offset, ct_i_bits.end());
        }
        idx++;
    }

    return bits;
}


DigitDecomposer::DigitDecomposer(const DigitDecompositionParameters &params) {

    m_br_params = std::make_shared<RingGSWCryptoParams>(params.m_N, params.m_Q, params.m_q, params.m_br_basis, 1, BINFHE_METHOD::GINX, params.m_std);
    m_ksk_params = std::make_shared<LWECryptoParams>(params.m_n, params.m_N, params.m_q, params.m_Q, params.m_ksk_q, params.m_std, params.m_ksk_basis);

    m_alpha = 1u << params.m_error_width;
    m_beta = 1u << (params.m_error_width - 2);
}

void DigitDecomposer::SetKeys(NativePoly &ring_sk, NativeVector &lwe_sk) {

    ring_sk.SetFormat(COEFFICIENT);
    m_ring_sk = ring_sk;
    m_lwe_sk = lwe_sk;

    assert(m_ring_sk.GetLength() == m_ksk_params->GetN());
    assert(m_lwe_sk.GetLength() == m_ksk_params->Getn());

    m_boot_engine = FunctionalBootstrapEngine(m_br_params,m_br_params, m_ring_sk, m_lwe_sk);
    m_lwe_ksk = LWEKeySwitchingKey(m_ksk_params, ring_sk.GetValues(), m_lwe_sk);

}

LWECiphertext DigitDecomposer::HomFloor(LWECiphertext &input) {

    auto small_modulus = m_br_params->Getq().ConvertToInt<uint64_t>();
    auto input_modulus = input->GetModulus().ConvertToInt<uint64_t>();
    auto keyswitch_modulus = m_ksk_params->GetqKS().ConvertToInt<uint64_t>();

    auto skk = m_lwe_sk;
    skk.SwitchModulus(input_modulus);
    ///auto ppp = input->GetB().ModSub(DotProduct(skk, input->GetA()), input_modulus);
    //std::cerr << "DIGIT-FLOOR-IN " << ppp << std::endl;;

    // TLDR, the authors of LMP21 casually forgot that when they write -q/4 for f0 that the minus applies w.r.t Q, not q, not ring Q...
    BootFunction F0 = [input_modulus](long m, long q, long ring_Q) {return (m < q/2 ? input_modulus - q/4 : q/4);};
    BootFunction F1 = [](long m, long q, long ring_Q) {return (m < q/2 ? m : (q + q/2 - m) % q);};

    // new
    input->GetB().ModAddEq(m_beta, input_modulus);

    auto input_A_q = input->GetA().Mod(small_modulus);
    auto input_B_q = input->GetB().Mod(small_modulus);


    input_A_q.SetModulus(small_modulus);
    // old
    //input_B_q.ModAddEq(m_beta, input_modulus);

    auto ctq1 = std::make_shared<LWECiphertextImpl>(input_A_q, input_B_q);

    auto ct_msb = m_boot_engine.Boot(ctq1, F0, 1, input_modulus);

    auto ct_msb_qks = ModSwitch(ct_msb, keyswitch_modulus);
    auto ct_ksk_ks = m_lwe_ksk.SwitchKey(ct_msb_qks);
    auto ct_ksk_Q = ModSwitch(ct_ksk_ks, input_modulus);

    // old
    //ct_ksk_Q->GetB().ModAddEq(small_modulus / 4, input_modulus);

    input->GetA().ModSubEq(ct_ksk_Q->GetA());
    input->GetB().ModSubEq(ct_ksk_Q->GetB(), input_modulus);

    ///ppp = ct_ksk_Q->GetB().ModSub(DotProduct(skk, ct_ksk_Q->GetA()), input_modulus);
    //std::cerr << "DIGIT-FLOOR-FIRST-BOOT " << ppp << std::endl;;

    //ppp = input->GetB().ModSub(DotProduct(skk, input->GetA()), input_modulus);
    //std::cerr << "DIGIT-FLOOR-FIRST " << ppp << std::endl;;

    input->GetB().ModAddEq(m_beta, input_modulus);
    input->GetB().ModSubEq(small_modulus / 4, input_modulus);

    input_A_q = input->GetA().Mod(small_modulus);
    input_B_q = input->GetB().Mod(small_modulus);

    input_A_q.SetModulus(small_modulus);
    ctq1 = std::make_shared<LWECiphertextImpl>(input_A_q, input_B_q);

    auto ct_remainder = m_boot_engine.Boot(ctq1, F1, 1, input_modulus);
    auto ct_remainder_qks = ModSwitch(ct_remainder, keyswitch_modulus);
    auto ct_remainder_ks = m_lwe_ksk.SwitchKey(ct_remainder_qks);
    auto ct_remainder_Q = ModSwitch(ct_remainder_ks, input_modulus);

    //ppp = ct_remainder_Q->GetB().ModSub(DotProduct(skk, ct_remainder_Q->GetA()), input_modulus);
    //std::cerr << "DIGIT-FLOOR-LAST-BOOT " << ppp << std::endl;;

    input->GetA().ModSubEq(ct_remainder_Q->GetA());
    input->GetB().ModSubEq(ct_remainder_Q->GetB(), input_modulus);

    //ppp = input->GetB().ModSub(DotProduct(skk, input->GetA()), input_modulus);
    //std::cerr << "DIGIT-FLOOR-LAST " << ppp << std::endl;;

    return input;
}

std::vector<LWECiphertext> DigitDecomposer::DigitDecompose(LWECiphertext &input, int32_t max_digits) {
    auto current_modulus = input->GetModulus().ConvertToInt<uint64_t>();
    auto small_modulus = m_br_params->Getq().ConvertToInt<uint64_t>();

    std::vector<LWECiphertext> results;
    results.reserve(64);

    while (true) {
        std::cerr << results.size() << std::endl;
        if (current_modulus <= small_modulus) {
            results.push_back(input);
            break;
        }

        auto digit_A = input->GetA().Mod(small_modulus);
        digit_A.SetModulus(small_modulus);
        auto digit_B = input->GetB().Mod(small_modulus);
        auto digit_ct = std::make_shared<LWECiphertextImpl>(digit_A, digit_B);
        results.push_back(digit_ct);

        if ((max_digits != 0) and (results.size() == max_digits))
            break;

        HomFloor(input);



        auto new_modulus = (m_alpha * current_modulus) / small_modulus;
        auto input_ms = ModSwitch(input, new_modulus);

        current_modulus = new_modulus;
        input = input_ms;
    }


    return results;
}

std::vector<LWECiphertext> DigitDecomposer::BitDecompose(LWECiphertext &input, int32_t max_bits) {

    auto small_modulus = m_br_params->Getq().ConvertToInt<uint64_t>();
    auto full_bits = GetMSB(small_modulus / m_alpha) - 1;
    auto first_bits = 0;

    auto needed_digits = 0;
    if (max_bits != 0) {
        needed_digits = std::ceil(double(max_bits) / double(full_bits));
    }

    auto large_digits = DigitDecompose(input, needed_digits);

    /*
    for (auto& ct : large_digits) {
        auto skk = m_lwe_sk;
        skk.SwitchModulus(ct->GetModulus());
        auto ppp = ct->GetB().ModSub(DotProduct(ct->GetA(), skk), ct->GetModulus());
        std::cerr << "DIGIT "  << ppp << " " << ct->GetModulus() << std::endl;
    }

    auto rlwesk = m_ring_sk;
    rlwesk.SetFormat(COEFFICIENT);
    for (auto& ct : large_digits) {
        auto deced = m_boot_engine.BitDecompose(ct, full_bits);
        for (auto& bitd : deced) {
            auto skk = rlwesk.GetValues();
            skk.SwitchModulus(bitd->GetModulus());
            auto ppp = bitd->GetB().ModSub(DotProduct(bitd->GetA(), skk), bitd->GetModulus());
            std::cerr << ppp << " " << bitd->GetModulus() << std::endl;
        }
    } */

    auto ms_digit = large_digits[large_digits.size() - 1];

    if (ms_digit->GetModulus() != small_modulus) {
        first_bits = full_bits - GetMSB(small_modulus / (ms_digit->GetModulus().ConvertToInt<uint64_t>())) - 1;
        ms_digit->SetModulus(small_modulus);
    }

    if (max_bits != 0 and (max_bits % full_bits) != 0) {
        first_bits = full_bits - (max_bits % full_bits);
    }

    // output is MSB to LSB so we to add them to the vector top to bottom since we want a MSB to LSB order
    auto top_bits = m_boot_engine.BitDecompose(ms_digit, full_bits);

    auto all_bits = std::vector<LWECiphertext>(top_bits.rbegin(), top_bits.rend() - first_bits);

    for (int32_t i = large_digits.size() - 2; i >= 0; i--) {
        auto ct_i_bits = m_boot_engine.BitDecompose(large_digits[i], full_bits);
        all_bits.insert(all_bits.end(), ct_i_bits.rbegin(), ct_i_bits.rend());
    }

    return all_bits;
}

const std::shared_ptr<RingGSWCryptoParams> DigitDecomposer::GetRingParams() const {
    return m_br_params;
}

