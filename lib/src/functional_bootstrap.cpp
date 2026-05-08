//
// Created by lschild on 19/11/24.
//

#include "functional_bootstrap.h"
#include "utils.h"

FunctionalBootstrapEngine::FunctionalBootstrapEngine(std::shared_ptr<RingGSWCryptoParams> &params,std::shared_ptr<RingGSWCryptoParams> &auto_params,
                                                     lbcrypto::NativePoly sk_rlwe, NativeVector sk_lwe) : m_br_params(params), m_automorphism_params(auto_params), m_br_key(params, sk_rlwe, sk_lwe, 2) {
    m_automorphism_keys = std::vector<AutomorphismKey>();
    for(int i = 2; i <= params->GetN(); i *= 2) {// TODO: check whether sk is ever modified
        m_automorphism_keys.emplace_back(m_automorphism_params, m_br_key.m_sk_ntt, i + 1);
    }
}


LWECiphertext FunctionalBootstrapEngine::Boot(ConstLWECiphertext& ct, BootFunction& F, uint64_t alpha, uint64_t output_modulus) {
    auto q = m_br_params->Getq();
    auto Q = m_br_params->GetQ();

    uint32_t N = m_br_params->GetN();
    uint64_t hi = output_modulus != 1 ? Q.ConvertToInt<uint64_t>() : 1;
    uint64_t lo = output_modulus != 1 ? output_modulus : 1;

    NativePoly rot_poly(m_br_params->GetPolyParams(), COEFFICIENT, true);
    NativePoly zero(m_br_params->GetPolyParams(), EVALUATION, true);

    // optionally shift message to make error positive
    // Sometimes, the desired function takes the error into account (alpha = 1), in which case we don't do anything
    uint32_t B_neg_shifted = q.ModSub(ct->GetB() + (alpha > 1 ? alpha / 4 : 0), q).ConvertToInt();


    for(long i = 0; i < N; i++) {
        auto index = (i + B_neg_shifted) % (2 * N);
        auto decoded_value = NativeInteger(i).DivideAndRound(alpha).Mod(q);
        NativeInteger f_val = F(decoded_value.ConvertToInt(), q.ConvertToInt(), Q.ConvertToInt());
        f_val.MultiplyAndRoundEq(hi, lo);
        f_val.ModEq(Q);

        rot_poly[index % N] = index >= N ? NativeInteger(Q).ModSub(f_val, Q) : f_val;

    }

    rot_poly.SetFormat(EVALUATION);
    auto vec = {zero, rot_poly};
    auto rlwect = std::make_shared<RLWECiphertextImpl>(vec);

    // Todo : put - if it is "not a fast br key
    m_br_key.BlindRotate(rlwect, ct->GetA());
    //m_br_key.BlindRotate(rlwect, -ct->GetA());

    auto A_T = rlwect->GetElements()[0].Transpose();
    A_T.SetFormat(COEFFICIENT);
    rlwect->GetElements()[1].SetFormat(COEFFICIENT);
    auto B = rlwect->GetElements()[1].at(0);

    return std::make_shared<LWECiphertextImpl>(A_T.GetValues(), B);
}


LWECiphertext
FunctionalBootstrapEngine::SampleExtract(lbcrypto::NativePoly &rlwe_A, lbcrypto::NativePoly &rlwe_B) {
    auto N = rlwe_A.GetLength();
    auto Q = rlwe_B.GetModulus();

    NativeVector A(N, Q);
    NativeInteger B = rlwe_B.at(0);

    A[0] = rlwe_A[0];
    for(uint32_t idx = 1; idx < N; idx++) {
        A[idx] = Q.ModSub(rlwe_A[N - idx], Q);
    }

    return std::make_shared<LWECiphertextImpl>(A, B);
}


std::vector<LWECiphertext> FunctionalBootstrapEngine::BitDecompose(lbcrypto::ConstLWECiphertext &ct, uint32_t n_msb) {
    // TODO [Low priority]: Investigate: do we actually need to do it here or can we do it during the Digit-Decomposition
    // For q = N
    // simply set the modulus to 2N, without modulus switching and blindrotate.
    // With alpha = Q/2, this implies there are no signflips

    // For q = 2N (THIS IS THE ONE WE USE)
    // set acc = Q/4 * \sum_i^{2N / t} X^i
    // To get the msb, multiply by p = \sum_i X^{i * 2N / t} and extract constant term as usual
    // For the other bits, multiply 2 * acc by whatever you need

    // lwe modulus
    auto q = ct->GetModulus();

    // rlwe modulus
    auto Q = m_br_params->GetQ().ConvertToInt();

    // polynomial dimension
    auto N = m_br_params->GetN();

    // plaintext space
    auto T = 1 << n_msb;

    if (m_pt_space == 0 or m_pt_space != T) {
        m_bit_extraction_polys.clear();

        // m_bit_extraction_polys start at msb
        NativePoly msb_poly(m_br_params->GetPolyParams(), COEFFICIENT, true);
        for(uint32_t i = 0; i < T / 2; i++) {
            msb_poly[i * 2 * N / T] = 1;
        }
        msb_poly.SwitchFormat();
        m_bit_extraction_polys.push_back(msb_poly);

        for(uint32_t msb_index = 0; msb_index < n_msb - 1; msb_index++) {
            NativePoly bit_poly(m_br_params->GetPolyParams(), COEFFICIENT, true);
            for(uint32_t i = 0; i < T / 2; i++) {
                bit_poly[i * 2 * N / T] = (i >> (n_msb - 2 - msb_index)) & 1;
            }
            bit_poly.SwitchFormat();
            m_bit_extraction_polys.push_back(bit_poly);
        }

        m_pt_space = T;
    }

    NativePoly acc(m_br_params->GetPolyParams(), COEFFICIENT, true);
    NativePoly zero(m_br_params->GetPolyParams(), EVALUATION, true);

    // we subtract 2 * N / (2 * T) = N / (T) to force the error to be negative
    auto b_uint = ct->GetB().ModAdd(N / (T), q).ConvertToInt();

    for(uint32_t i = 0; i < (2 * N) / T; i++) {

        // shifting the coefficient index corresponds to multiplying by X^b later, but we save 1 NTT
        auto idx = (2 * N + i - b_uint) % (2 * N);
        acc[idx % N] = idx >= N ? Q / 4 : NativeInteger(Q).ModSub(Q / 4, Q);
    }

    acc.SwitchFormat();
    auto ct_vec = {zero, acc};
    auto ct_acc = std::make_shared<RLWECiphertextImpl>(ct_vec);

    // Note: we do -A, as we aim to compute X^{<a,s>} but the br key would compute X^{- <a,s>}
    // Todo : put - if it is "not a fast br key
    auto br_result = m_br_key.BlindRotate(ct_acc, ct->GetA());
    auto br_A = br_result->GetElements()[0];
    auto br_B = br_result->GetElements()[1];

    std::vector<LWECiphertext> result;


    // recover MSB first
    auto& msb_poly = m_bit_extraction_polys.at(0);
    auto msb_A = br_A * msb_poly;
    auto msb_B = br_B * msb_poly;
    msb_A.SwitchFormat();
    msb_B.SwitchFormat();

    auto ct_msb_bit = SampleExtract(msb_A, msb_B);
    // Add Q / 4 since the msb is encoded w.r.t the sign flip i.e. MSB == 0 => -Q/4, else Q/4
    // adding Q/4 maps it to 0 or Q/2
    ct_msb_bit->GetB().ModAddEq(Q / 4, Q);

    // we no longer care about the sign flip, only about nonzero entries
    br_A *= 2;
    br_B *= 2;

    for(uint32_t msb_index = m_bit_extraction_polys.size() - 1; msb_index >= 1; msb_index--) {
        auto bit_A = br_A * m_bit_extraction_polys.at(msb_index);
        auto bit_B = br_B * m_bit_extraction_polys.at(msb_index);
        bit_A.SwitchFormat();
        bit_B.SwitchFormat();

        result.push_back(SampleExtract(bit_A, bit_B));
    }

    result.push_back(ct_msb_bit);

    /*
    NativeInteger phase = ct->GetB();
    NativeInteger dot = DotProduct(ct->GetA(), m_br_key.m_sk_lwe);
    phase.ModSubEq(dot, q);

    NativeVector ring_sk_vec = m_br_key.m_sk.GetValues();
    ring_sk_vec.SwitchModulus(Q);

    auto decoded = phase.MultiplyAndRound(T, q).Mod(T).ConvertToInt<uint32_t>();
    for(int i = 0; i < n_msb; i++) {

        auto ct_i = result[i];
        NativeInteger pi = ct_i->GetB().ModSub(DotProduct(ring_sk_vec, ct_i->GetA()), Q);
        pi.MultiplyAndRoundEq(2, Q);

        auto bit_i = (decoded >> i) & 1;
        std::cerr << "i = " << i << " bit = " << bit_i << "| Got " << pi.ConvertToInt() % 2 << std::endl;
    } */

    return result;
}