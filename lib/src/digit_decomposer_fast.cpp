//
// Created by leonard on 2/17/26.
//

#include "digit_decomposer_fast.h"

#include <cassert>

FastDigitDecomposer::FastDigitDecomposer(const FastDigitDecompositionParameters &params) {

    m_br_params = std::make_shared<RingGSWCryptoParams>(params.m_N, params.m_Q, params.m_q, params.m_br_basis, 1, BINFHE_METHOD::GINX, params.m_std);
    m_auto_params = std::make_shared<RingGSWCryptoParams>(params.m_N, params.m_Q, params.m_q, params.m_auto_base, 1, BINFHE_METHOD::GINX, 0);
    m_alpha = params.m_alpha;
    m_beta = params.m_beta;
    assert(m_alpha <= 2 * m_beta);
}

void FastDigitDecomposer::SetKeys(NativePoly &ring_sk, NativeVector &lwe_sk) {
    ring_sk.SetFormat(COEFFICIENT);
    m_ring_sk = ring_sk;
    m_lwe_sk = lwe_sk;
    auto ntt_engine = std::make_shared<intel::hexl::NTT>(m_auto_params->GetN(), m_auto_params->GetQ().ConvertToInt<uint64_t>());
    auto L_bits = intel::hexl::Log2(m_auto_params->GetBaseG());

    assert(m_ring_sk.GetLength() == m_br_params->GetN());

    m_br_key = FastBlindRotationKey(m_br_params, m_ring_sk, m_lwe_sk, 2);
    m_decimation_keys = std::vector<FastAutomorphismKey>();
    //m_decimation_keys = std::vector<AutomorphismKey>();

    uint32_t idx = 2;
    while (idx <= m_br_params->GetN()) {
        uint32_t auto_idx = idx + 1;
        m_decimation_keys.emplace_back(ntt_engine, L_bits, m_ring_sk, auto_idx);
        //m_decimation_keys.emplace_back(m_auto_params, m_ring_sk, auto_idx);

        idx *= 2;
    }
}

int find_and_printnz(NativePoly& p_i) {
    for (int i = 0; i < p_i.GetLength(); i++) {
        if (p_i[i].ConvertToInt<uint64_t>() != 0) {
            return i;
        }
    }
}

std::vector<LWECiphertext> FastDigitDecomposer::DigitDecompose(LWECiphertext &input) {

    auto lwe_q = input->GetModulus();
    auto N = m_br_params->GetN();
    auto t = m_alpha;
    auto alpha_bits = intel::hexl::Log2(m_alpha);
    NativeInteger scale = m_auto_params->GetQ().DividedBy(t);

    NativePoly acc_A(m_br_params->GetPolyParams(), EVALUATION, true);
    NativePoly acc_B(m_br_params->GetPolyParams(), COEFFICIENT, true);
    acc_B[0]= scale;
    auto digs = {acc_A, acc_B};
    RLWECiphertext acc = std::make_shared<RLWECiphertextImpl>(digs);
    acc->SetFormat(EVALUATION);

    auto phases = PhaseDecomp(input);
    NativePoly extractor(m_br_params->GetPolyParams(), COEFFICIENT, true);
    for (uint32_t i = 0; i < N; i++) {
        extractor[N - i - 1] = (t - i - 1) % t;
    }
    extractor.SetFormat(EVALUATION);

    LWECiphertext ct0 = std::make_shared<LWECiphertextImpl>(phases[0]->GetA(), phases[0]->GetB());
    auto zero_scale = N / m_alpha;
    ct0->GetA().ModMulEq(zero_scale);
    ct0->GetB().ModMulEq(zero_scale, N);

    std::vector<LWECiphertext> digits = {ct0};
    NativePoly b_poly(m_br_params->GetPolyParams(), COEFFICIENT, true);
    b_poly[phases[0]->GetB().ConvertToInt<uint64_t>()] = 1;
    b_poly.SwitchFormat();
    acc->GetElements()[0] *= b_poly;
    acc->GetElements()[1] *= b_poly;
    auto tk = m_lwe_sk;
    tk.SetModulus(2 * N);
    phases[0]->SetModulus(2 * N);
    //auto check = decrypt_lwe(tk, phases[0]);
    //std::cerr << "Phase 0 = " << check << std::endl;
    acc = m_br_key.BlindRotate(acc, -phases[0]->GetA());

    NativePoly skNTT = m_ring_sk;
    skNTT.SetFormat(EVALUATION);
    for (uint32_t i = 1; i < phases.size() - 1; i++) {
        //auto p_i = acc->GetElements()[1] - acc->GetElements()[0] * skNTT;
        //p_i.SwitchFormat();
        //uint32_t idx = find_and_printnz(p_i);
        //std::cerr << i << " " << idx << " " << idx % m_alpha << " " << (idx >> alpha_bits) << p_i << std::endl;
        acc = HomTrunc(acc);
        NativePoly b_poly(m_br_params->GetPolyParams(), COEFFICIENT, true);
        b_poly[phases[i]->GetB().ConvertToInt<uint64_t>()] = 1;
        b_poly.SwitchFormat();
        acc->GetElements()[0] *= b_poly;
        acc->GetElements()[1] *= b_poly;
        acc = m_br_key.BlindRotate(acc, -phases[i]->GetA());
        auto acc_i_A = acc->GetElements()[0] * extractor;
        acc_i_A.SwitchFormat();
        auto acc_i_B = acc->GetElements()[1] * extractor;
        acc_i_B.SwitchFormat();
        digits.push_back(RLWESampleExtract(acc_i_A, acc_i_B, 0));
    }

    return digits;
}


std::vector<LWECiphertext> FastDigitDecomposer::PhaseDecomp(LWECiphertext &input) {

    NativeInteger offset_lo = m_alpha * m_beta;
    NativeInteger offset_hi = m_beta;

    uint64_t current_modulus = input->GetModulus().ConvertToInt<uint64_t>();
    auto N = m_br_params->GetN();

    auto current_A = input->GetA();
    auto current_B = input->GetB();
    auto shift_bits = intel::hexl::Log2(m_alpha);

    NativeVector tmp(current_A.GetLength(), current_A.GetModulus());

    std::vector<LWECiphertext> results;

    while (current_modulus > m_alpha) {
        auto digit_A = current_A.Mod(m_alpha);
        auto digit_B = current_B.Mod(m_alpha);

        digit_A.SetModulus(2*N);
        digit_B.ModAddEq(offset_lo, 2*N);

        results.push_back(std::make_shared<LWECiphertextImpl>(digit_A, digit_B));

        current_modulus >>= shift_bits;

        for (int i = 0; i < current_A.GetLength(); i++) {
            tmp[i] = current_A[i].ConvertToInt<uint64_t>() >> shift_bits;
        }
        current_A = tmp;
        current_A.SetModulus(current_modulus);

        current_B = (current_B.ConvertToInt<uint64_t>() >> shift_bits);
        current_B.ModSubEq(offset_hi, current_modulus);
    }

    current_A.SetModulus(2 * N);
    results.push_back(std::make_shared<LWECiphertextImpl>(current_A, current_B));

    return results;
}

RLWECiphertext FastDigitDecomposer::HomTrunc(RLWECiphertext &acc) {
    auto N = m_auto_params->GetN();
    auto alpha = m_alpha;
    // TODO
    auto rat = m_beta >> 1;

    // TODO can precompute all this
    NativePoly mon(m_auto_params->GetPolyParams(), COEFFICIENT, true);
    mon[1] = 1;
    mon.SwitchFormat();
    NativePoly shift = NativePoly(m_auto_params->GetPolyParams(), COEFFICIENT, true)
    ;
    shift[0]=  1;
    shift.SetFormat(EVALUATION);

    std::vector<LWECiphertext> coefs(N);
    acc->SetFormat(COEFFICIENT);
    auto elems = acc->GetElements();
    for (uint32_t i = 0; i < N; i++) {
        coefs[i] = RLWESampleExtract(elems[0], elems[1], i);
    }

    NativePoly acc_A(m_auto_params->GetPolyParams(), EVALUATION, true);
    NativePoly acc_B(m_auto_params->GetPolyParams(), EVALUATION, true);

    for (uint32_t i = 0; i < rat; i++) {
        auto base_i = coefs[i * alpha];
        auto base_A = base_i->GetA();
        auto base_B = base_i->GetB();
        for (uint32_t j = 1; j < alpha; j++) {
            auto& acc_j = coefs[i * alpha + j];
            base_A.ModAddEq(acc_j->GetA());
            base_B.ModAddEq(acc_j->GetB(), base_A.GetModulus());
        }
        auto lwe = std::make_shared<LWECiphertextImpl>(base_A, base_B);
        auto switched = LWE2RLWE(lwe);
        // TODO get rid of NTT
        switched->SetFormat(EVALUATION);
        acc_A += shift * switched->GetElements()[0];
        acc_B += shift * switched->GetElements()[1];

        shift *= mon;
    }

    auto digs = {acc_A, acc_B};

    return std::make_shared<RLWECiphertextImpl>(digs);

}

RLWECiphertext FastDigitDecomposer::LWE2RLWE(LWECiphertext &ct) {
    NativePoly A(m_auto_params->GetPolyParams(), COEFFICIENT, true);
    NativePoly B(m_auto_params->GetPolyParams(), COEFFICIENT, true);
    A.SetValues(ct->GetA(), COEFFICIENT);
    // TODO: remove extra NTT
    A.SwitchFormat();
    A = A.Transpose();
    A.SwitchFormat();
    B[0] = ct->GetB();


    auto n_automorphism = m_decimation_keys.size();
    NativeInteger n_auto_pow = 1 << n_automorphism;
    NativeInteger premult = n_auto_pow.ModExp(m_auto_params->GetQ() - 2, m_auto_params->GetQ());

    A *= premult;
    B *= premult;

    m_ring_sk.SetFormat(EVALUATION);

    for(int i = 0; i < n_automorphism; i++) {
        auto m_key_i = m_decimation_keys.at(n_automorphism - i - 1);
        auto rlwe_auto_i = m_key_i.AutomorphismTransform(A, B);
        //rlwe_auto_i->SetFormat(COEFFICIENT);
        A+=rlwe_auto_i->GetElements()[0];
        B+=rlwe_auto_i->GetElements()[1];
    }


    std::vector<NativePoly> rlwe_vec = {A, B};
    return std::make_shared<RLWECiphertextImpl>(rlwe_vec);
}


const std::shared_ptr<RingGSWCryptoParams> FastDigitDecomposer::GetRingParams() const {
    return m_br_params;
}
