//
// Created by leonard on 11/12/24.
//

#include "utils.h"
#include "blind_rotator.h"

#include <cassert>

BlindRotationKey::BlindRotationKey(std::shared_ptr<RingGSWCryptoParams> &params, lbcrypto::NativePoly sk_rlwe,
                                   NativeVector sk_lwe) : m_params(params) {
    sk_rlwe.SetFormat(COEFFICIENT);
    m_sk = sk_rlwe.Clone();
    sk_rlwe.SetFormat(EVALUATION);
    auto digits = params->GetDigitsG();
    auto base = params->GetBaseG();


    m_sk_ntt = sk_rlwe.Clone();
    m_sk_lwe = sk_lwe;
    m_sk_lwe.SwitchModulus(params->Getq());

    auto N = params->GetN();
    m_acc_buffer.resize(2 * N);
    auto modulus = params->GetQ().ConvertToInt();
    auto root = params->GetPolyParams()->GetRootOfUnity().ConvertToInt();

    auto engine = std::make_shared<intel::hexl::NTT>(N, modulus, root);

    for(int i = 0; i < sk_lwe.GetLength(); i++) {
        int bit_0 = 0;
        int bit_1 = 0;

        if (sk_lwe.at(i) == 1) {
            bit_0 = 1;
        }
        if (sk_lwe.at(i) > 1) {
            bit_1 = 1;
        }

        NativePoly b0_poly = NativePoly(params->GetPolyParams(), COEFFICIENT, true);
        b0_poly.at(0) = bit_0;
        b0_poly.SetFormat(EVALUATION);

        NativePoly b1_poly = NativePoly(params->GetPolyParams(), COEFFICIENT, true);
        b1_poly.at(0) = bit_1;
        b1_poly.SetFormat(EVALUATION);

        std::vector<RLWECiphertext> m0;
        std::vector<RLWECiphertext> msm0;

        std::vector<RLWECiphertext> m1;
        std::vector<RLWECiphertext> msm1;

        NativePoly b0_poly_sk = -sk_rlwe * b0_poly.Clone();
        NativePoly b1_poly_sk = -sk_rlwe * b1_poly.Clone();

        for(int j = 0; j < digits; j++) {
            auto ct_0 = encrypt_rlwe(params, b0_poly, sk_rlwe);
            auto ct_1 = encrypt_rlwe(params, b1_poly, sk_rlwe);
            m0.push_back(ct_0);
            m1.push_back(ct_1);

            auto ct_ms0 = encrypt_rlwe(params, b0_poly_sk, sk_rlwe);
            auto ct_ms1 = encrypt_rlwe(params, b1_poly_sk, sk_rlwe);

            msm0.push_back(ct_ms0);
            msm1.push_back(ct_ms1);

            b0_poly *= base;
            b1_poly *= base;
            b0_poly_sk *= base;
            b1_poly_sk *= base;

        }

        auto b0_rgsw = RGSWSample(engine, base, digits, msm0, m0);
        auto b1_rgsw = RGSWSample(engine, base, digits, msm1, m1);

        m_brk.emplace_back(b0_rgsw, b1_rgsw);
    }

    auto Q = m_params->GetQ();
    auto q = m_params->Getq().ConvertToInt();

    m_monomials_raw.resize(q * N);
    for(long idx = 0; idx < q; idx++) {
        auto block_start = m_monomials_raw.data() + idx * N;
        NativePoly monA(m_params->GetPolyParams(), COEFFICIENT, true);
        monA.at(idx % N) = idx >= N ? (Q  - 1) : 1;
        monA.at(0).ModSubEq(1, Q);
        monA.SetFormat(EVALUATION);

        for(uint32_t k = 0; k < N; k++) {
            block_start[k] = monA[k].ConvertToInt();
        }
    }

}

RLWECiphertext BlindRotationKey::BlindRotate(lbcrypto::RLWECiphertext &in, NativeVector A) {
    auto N = m_params->GetN();

    // copy to buffer
    for(uint32_t i = 0; i < m_params->GetN(); i++) {
        m_acc_buffer[i] = in->GetElements()[0][i].ConvertToInt();
        m_acc_buffer[i + N] = in->GetElements()[1][i].ConvertToInt();
    }

    for(int i = 0; i < A.GetLength(); i++) {
        auto a_idx = A.at(i).ConvertToInt<uint32_t>();

        // TODO: as weird as it is, it might make sense to sort the LWE sample first (the  vector at least)
        // TODO: for n \approx q duplicates become quite likely, so we could re-use the monomial and it's nice for cache
        if (a_idx == 0) [[unlikely]] {
        } else {
            auto a_neg = (2 * N - a_idx) % (2 * N);
            //m_brk.at(i).first.cmux(m_acc_buffer.data(), a_neg);
            m_brk.at(i).first.ternary_mux(m_acc_buffer.data(), m_monomials_raw.data() + a_neg * N, m_monomials_raw.data() + a_idx * N, m_brk.at(i).second);
        }
    }


    NativeVector W(N, in->GetElements()[0].GetModulus());
    for(uint32_t i = 0; i < N; i++)
        W[i] = m_acc_buffer[i];

    in->GetElements()[0].SetValues(W, EVALUATION);

    for(uint32_t i = 0; i < N; i++)
        W[i] = m_acc_buffer[i + N];

    in->GetElements()[1].SetValues(W, EVALUATION);

    return in;
}

FastBlindRotationKey::FastBlindRotationKey(std::shared_ptr<RingGSWCryptoParams> &params, NativePoly sk_rlwe, NativeVector sk_lwe, uint32_t step_size)
: m_params(params) {
    m_step_size = step_size;

    sk_rlwe.SetFormat(COEFFICIENT);
    m_sk = sk_rlwe.Clone();
    sk_rlwe.SetFormat(EVALUATION);

    assert(sk_lwe.GetLength() % 3 != 2 or step_size != 3);

    m_sk_ntt = sk_rlwe.Clone();
    m_sk_lwe = sk_lwe;
    m_sk_lwe.SwitchModulus(params->Getq());

    auto Q = m_params->GetQ().ConvertToInt<uint64_t>();
    auto N = m_params->GetN();
    auto digits = m_params->GetDigitsG();

    m_rlwe_size = 2 * N;
    m_rlwe_prime_size = m_rlwe_size * digits;
    m_rgsw_size = m_rlwe_prime_size * 2;

    auto root = m_params->GetPolyParams()->GetRootOfUnity().ConvertToInt<uint64_t>();

    m_engine = std::make_unique<intel::hexl::NTT>(N, Q, root);
    m_rgsw_count_per_step = (1 << (step_size)) - 1;
    m_digit_buffer.resize(m_rgsw_count_per_step * m_rgsw_size);

    m_one_ntt.resize(N);
    m_one_ntt[0] = 1;
    m_engine->ComputeForward(m_one_ntt.data(), m_one_ntt.data(), 1, 1);

    m_poly_buffer.resize(m_rgsw_count_per_step * N);

    switch (m_step_size) {
        case 2: {
            KeyGen2();
            break;
        }
        case 3 : {
            KeyGen3();
            break;
        }
        default: {
            m_step_size = 1;
            KeyGen1();
            break;
        }
    }

}

void FastBlindRotationKey::KeyGen1() {
    auto lwe_n = m_sk_lwe.GetLength();
    auto q = m_params->Getq();
    auto Q = m_engine->GetModulus();
    auto N = m_engine->GetDegree();
    m_monomials_raw.resize(2 * N * N);

    m_brk_chunk_size = m_rgsw_size;

    m_brk.resize(m_brk_chunk_size * lwe_n);
    m_acc_buffer.resize(2 * m_rlwe_size);

    for(long idx = 0; idx < q; idx++) {
        auto block_start = m_monomials_raw.data() + idx * N;
        NativePoly monA(m_params->GetPolyParams(), COEFFICIENT, true);
        monA.at(idx % N) = idx >= N ? (Q  - 1) : 1;
        monA.at(0).ModSubEq(1, Q);
        monA.SetFormat(EVALUATION);

        for(uint32_t k = 0; k < N; k++) {
            block_start[k] = monA[k].ConvertToInt();
        }
    }

    auto brk_start = m_brk.data();

    for (uint32_t i = 0 ; i < lwe_n; i++) {
        auto v_0 = m_sk_lwe.at(i).ConvertToInt<uint64_t>();

        NativePoly p0(m_params->GetPolyParams(), COEFFICIENT, true);
        p0[0] = v_0;

        p0.SwitchFormat();

        auto c0 = encrypt_rgsw(m_params, p0, m_sk_ntt);
        auto C0 = RGSWSample::convert(c0);
        std::copy(C0.m_key.begin(), C0.m_key.end(), brk_start);

        brk_start += m_brk_chunk_size;
    }
}


void FastBlindRotationKey::KeyGen2() {

    auto lwe_n = m_sk_lwe.GetLength();
    auto tail = lwe_n % 2;

    m_brk_chunk_size = m_rgsw_size * 3;
    auto brk_size = m_brk_chunk_size * (lwe_n >> 1) + tail * m_rgsw_size;

    m_brk.resize(brk_size);
    m_acc_buffer.resize(4 * m_rlwe_size);

    auto q = m_params->Getq();
    auto Q = m_engine->GetModulus();
    auto N = m_engine->GetDegree();
    m_monomials_raw.resize(4 * N * N);

    for(long idx = 0; idx < q; idx++) {
        auto block_start = m_monomials_raw.data() + idx * N;
        NativePoly monA(m_params->GetPolyParams(), COEFFICIENT, true);
        monA.at(idx % N) = idx >= N ? (Q  - 1) : 1;
        monA.SetFormat(EVALUATION);

        for(uint32_t k = 0; k < N; k++) {
            block_start[k] = monA[k].ConvertToInt();
        }
    }

    for(long idx = 0; idx < q; idx++) {
        auto block_start = m_monomials_raw.data() + (idx+2*N) * N;
        NativePoly monA(m_params->GetPolyParams(), COEFFICIENT, true);
        monA.at(idx % N) = idx >= N ? (Q  - 1) : 1;
        monA.at(0).ModSubEq(1, Q);
        monA.SetFormat(EVALUATION);

        for(uint32_t k = 0; k < N; k++) {
            block_start[k] = monA[k].ConvertToInt();
        }
    }

    auto brk_start = m_brk.data();

    for (int i = 0; i < lwe_n; i+=2) {

        auto v_0 = m_sk_lwe.at(i).ConvertToInt<uint64_t>();
        auto v_1 = m_sk_lwe.at(i + 1).ConvertToInt<uint64_t>();
        auto v_01 = v_0 * v_1;

        NativePoly p0(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p1(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p01(m_params->GetPolyParams(), COEFFICIENT, true);

        p0[0] = v_0;
        p1[0] = v_1;
        p01[0] = v_01;

        p0.SwitchFormat();
        p1.SwitchFormat();
        p01.SwitchFormat();

        auto c0 = encrypt_rgsw(m_params, p0, m_sk_ntt);
        auto c1 = encrypt_rgsw(m_params, p1, m_sk_ntt);
        auto c01 = encrypt_rgsw(m_params, p01, m_sk_ntt);

        auto C0 = RGSWSample::convert(c0);
        auto C1 = RGSWSample::convert(c1);
        auto C01 = RGSWSample::convert(c01);

        std::copy(C0.m_key.begin(), C0.m_key.end(), brk_start);
        std::copy(C1.m_key.begin(), C1.m_key.end(), brk_start + m_rgsw_size);
        std::copy(C01.m_key.begin(), C01.m_key.end(), brk_start + m_rgsw_size * 2);

        brk_start += m_brk_chunk_size;

    }

    if (tail == 1) {

        auto v_0 = m_sk_lwe.at(lwe_n - 1).ConvertToInt<uint64_t>();
        NativePoly p0(m_params->GetPolyParams(), COEFFICIENT, true);
        p0[0] = v_0;
        p0.SwitchFormat();

        auto c0 = encrypt_rgsw(m_params, p0, m_sk_ntt);
        auto C0 = RGSWSample::convert(c0);
        std::copy(C0.m_key.begin(), C0.m_key.end(), brk_start);

    }

}

void FastBlindRotationKey::KeyGen3() {
    auto lwe_n = m_sk_lwe.GetLength();
    auto tail = lwe_n % 3;

    m_brk_chunk_size = m_rgsw_size * 7;
    auto brk_size = m_brk_chunk_size * (lwe_n - tail) + tail * m_rgsw_size;

    m_brk.resize(brk_size);
    m_acc_buffer.resize(7 * m_rlwe_size);

    auto brk_start = m_brk.data();

    for (int i = 0; i < lwe_n; i+=3) {

        auto v_0 = m_sk_lwe.at(i).ConvertToInt<uint64_t>();
        auto v_1 = m_sk_lwe.at(i + 1).ConvertToInt<uint64_t>();
        auto v_2 = m_sk_lwe.at(i + 2).ConvertToInt<uint64_t>();

        auto v_01 = v_0 * v_1;
        auto v_02 = v_0 * v_2;
        auto v_12 = v_1 * v_2;
        auto v_012 = v_0 * v_1 * v_2;

        NativePoly p0(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p1(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p2(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p01(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p02(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p12(m_params->GetPolyParams(), COEFFICIENT, true);
        NativePoly p012(m_params->GetPolyParams(), COEFFICIENT, true);

        p0[0] = v_0;
        p1[0] = v_1;
        p2[0] = v_2;

        p01[0] = v_01;
        p02[0] = v_02;
        p12[0] = v_12;

        p012[0] = v_012;

        p0.SwitchFormat();
        p1.SwitchFormat();
        p2.SwitchFormat();

        p01.SwitchFormat();
        p12.SwitchFormat();
        p02.SwitchFormat();

        p012.SwitchFormat();

        auto c0 = encrypt_rgsw(m_params, p0, m_sk_ntt);
        auto c1 = encrypt_rgsw(m_params, p1, m_sk_ntt);
        auto c2 = encrypt_rgsw(m_params, p2, m_sk_ntt);

        auto c01 = encrypt_rgsw(m_params, p01, m_sk_ntt);
        auto c12 = encrypt_rgsw(m_params, p12, m_sk_ntt);
        auto c02 = encrypt_rgsw(m_params, p02, m_sk_ntt);
        auto c012 = encrypt_rgsw(m_params, p012, m_sk_ntt);

        auto C0 = RGSWSample::convert(c0);
        auto C1 = RGSWSample::convert(c1);
        auto C2 = RGSWSample::convert(c2);

        auto C01 = RGSWSample::convert(c01);
        auto C12 = RGSWSample::convert(c12);
        auto C02 = RGSWSample::convert(c02);

        auto C012 = RGSWSample::convert(c012);

        std::copy(C0.m_key.begin(), C0.m_key.end(), brk_start);
        std::copy(C1.m_key.begin(), C1.m_key.end(), brk_start + m_rgsw_size);
        std::copy(C2.m_key.begin(), C2.m_key.end(), brk_start + m_rgsw_size * 2);

        std::copy(C01.m_key.begin(), C0.m_key.end(), brk_start + m_rgsw_size * 3);
        std::copy(C02.m_key.begin(), C02.m_key.end(), brk_start + m_rgsw_size * 4);
        std::copy(C12.m_key.begin(), C12.m_key.end(), brk_start + m_rgsw_size * 5);

        std::copy(C012.m_key.begin(), C012.m_key.end(), brk_start + 7 * m_rgsw_size);

        brk_start += m_brk_chunk_size;

    }

    for (uint32_t tail_i = 0; tail_i < tail; tail_i++) {
        auto v_0 = m_sk_lwe.at(lwe_n - 1 - tail_i).ConvertToInt<uint64_t>();
        NativePoly p0(m_params->GetPolyParams(), COEFFICIENT, true);
        p0[0] = v_0;
        p0.SwitchFormat();

        auto c0 = encrypt_rgsw(m_params, p0, m_sk_ntt);
        auto C0 = RGSWSample::convert(c0);
        std::copy(C0.m_key.begin(), C0.m_key.end(), brk_start);
        brk_start += m_brk_chunk_size;
    }

}

void FastBlindRotationKey::MakePoly2(uint64_t* buffer, uint64_t a_0, uint64_t a_1) {
    // TODO: check whether that's faster than just multiplying.
    // It should be because we we do 3 * N additions versus N multiplications but we need to check

    // returns (X^a0 - 1) * (X^a1 - 1)
    auto N = m_engine->GetDegree();
    auto modulus = m_engine->GetModulus();
    auto prod = (a_0 + a_1) % (2 * N);

    auto offset_1 = prod * N;
    auto offset_2 = a_0 * N;
    auto offset_3 = a_1 * N;

    auto data_ptr = buffer;
    auto mon_ptr = m_monomials_raw.data();

    intel::hexl::EltwiseAddMod(data_ptr, mon_ptr + offset_2, mon_ptr + offset_3, N, modulus);
    intel::hexl::EltwiseSubMod(data_ptr, mon_ptr + offset_1, data_ptr, N, modulus);
    intel::hexl::EltwiseAddMod(data_ptr, m_one_ntt.data(), data_ptr, N, modulus);

}

void FastBlindRotationKey::MakePoly3(uint64_t* buffer, uint64_t a_0, uint64_t a_1, uint64_t a_2) {
    // TODO: check whether that's faster than just multiplying.
    // It should be because we we do 3 * N additions versus N multiplications but we need to check

    // returns (X^a0 - 1) * (X^a1 - 1)
    auto N = m_engine->GetDegree();
    auto modulus = m_engine->GetModulus();

    auto prod_01 = (a_0 * a_1) % (2 * N);
    auto prod_02 = (a_0 * a_2) % (2 * N);
    auto prod_12 = (a_1 * a_2) % (2 * N);
    auto prod_012 = (prod_01 * a_2) % (2 * N);

    auto offset_0 = a_0 * N;
    auto offset_1 = a_1 * N;
    auto offset_2 = a_2 * N;

    auto offset_01 = prod_01 * N;
    auto offset_12 = prod_12 * N;
    auto offset_02 = prod_02 * N;

    auto offset_012 = prod_012 * N;

    auto data_ptr = buffer;
    auto mon_ptr = m_monomials_raw.data();

    intel::hexl::EltwiseAddMod(data_ptr, mon_ptr + offset_01, mon_ptr + offset_12, N, modulus);
    intel::hexl::EltwiseAddMod(data_ptr, mon_ptr + offset_02, data_ptr, N, modulus);
    intel::hexl::EltwiseSubMod(data_ptr, mon_ptr + offset_012, data_ptr, N, modulus);

    intel::hexl::EltwiseAddMod(data_ptr, mon_ptr + offset_0, data_ptr, N, modulus);
    intel::hexl::EltwiseAddMod(data_ptr, mon_ptr + offset_1, data_ptr, N, modulus);
    intel::hexl::EltwiseAddMod(data_ptr, mon_ptr + offset_2, data_ptr, N, modulus);

    intel::hexl::EltwiseSubMod(data_ptr, data_ptr, m_one_ntt.data(), N, modulus);
}

void FastBlindRotationKey::Step1(Vector& buffer,uint32_t idx, uint32_t a_i) {

    const auto N = m_engine->GetDegree();
    auto modulus = m_engine->GetModulus();

    // we have enough space for 4 polys
    m_engine->ComputeInverse(m_digit_buffer.data(), buffer.data(), 1, 1);
    RLWEVectorProd(m_brk.data() + idx * m_brk_chunk_size, m_digit_buffer.data(), m_acc_buffer.data());

    m_engine->ComputeInverse(m_digit_buffer.data(), buffer.data() + N, 1, 1);
    RLWEVectorProd(m_brk.data() + idx * m_brk_chunk_size + m_rlwe_prime_size, m_digit_buffer.data(), m_acc_buffer.data() + 2 * N);

    for (uint32_t i = 0; i < 2 * N; i++) {
        m_acc_buffer[i] += m_acc_buffer[i + 2 * N];
    }
    auto mon_ptr = m_monomials_raw.data() + a_i * N;

    intel::hexl::EltwiseReduceMod(m_acc_buffer.data(),m_acc_buffer.data(),N,modulus,modulus,1);
    intel::hexl::EltwiseMultMod(m_acc_buffer.data(), mon_ptr, m_acc_buffer.data(),N, modulus, 1);
    intel::hexl::EltwiseMultMod(m_acc_buffer.data() + N, mon_ptr, m_acc_buffer.data() + N,N, modulus, 1);

    intel::hexl::EltwiseAddMod(buffer.data(), buffer.data(), m_acc_buffer.data(),2 * N, modulus);
}

void FastBlindRotationKey::Step2(Vector& buffer,uint32_t idx, uint32_t a_i, uint32_t a_i1) {

    auto N = m_engine->GetDegree();
    auto rgsw_ptr = m_brk.data() + (idx >> 1) * m_brk_chunk_size;

    auto mon_ptr_0 = m_monomials_raw.data() + 2 * N * N + N * a_i;
    auto mon_ptr_1 = m_monomials_raw.data() + 2 * N * N + N * a_i1;
    std::copy(mon_ptr_0, mon_ptr_0 + N, m_poly_buffer.data());
    std::copy(mon_ptr_1, mon_ptr_1 + N, m_poly_buffer.data() + N);
    MakePoly2(m_poly_buffer.data() + 2 * N, a_i, a_i1);

    ExtProd(rgsw_ptr, nullptr, m_poly_buffer.data(), buffer.data());

}

void FastBlindRotationKey::Step3(Vector& buffer,uint32_t idx, uint32_t a_i, uint32_t a_i1, uint32_t a_i2) {

    auto basis = m_params->GetBaseG();
    auto basis_bits = GetMSB(basis) - 1;
    auto mask = basis - 1;
    auto digits = m_params->GetDigitsG();

    auto N = m_engine->GetDegree();
    auto modulus = m_engine->GetModulus();
    auto acc_ptr = m_acc_buffer.data();
    auto digit_ptr = m_digit_buffer.data();

    m_engine->ComputeInverse(acc_ptr, buffer.data(), 1, 1);
    m_engine->ComputeInverse(acc_ptr + N, buffer.data() + N, 1, 1);

    // digit decomp
    for (uint32_t i = 0; i < digits; i++) {
        auto shift = (digits - 1 - i) * basis_bits;
        for (uint32_t j = 0; j < 2 * N; j++) {
            m_digit_buffer[i * 2 * N + j] = (acc_ptr[j] >> (shift)) & mask;
        }
    }

    for (uint32_t i = 0; i < 2 * digits; i++) {
        m_engine->ComputeForward(m_digit_buffer.data() + i * N, m_digit_buffer.data() + i * N, 1, 1);
    }

    auto rgsw_ptr = m_brk.data() + (idx / 3) * m_brk_chunk_size;
    auto rgsw0 = rgsw_ptr;
    auto rgsw1 = rgsw_ptr + m_rgsw_size;
    auto rgsw2 = rgsw_ptr + m_rgsw_size * 2;
    auto rgsw01 = rgsw_ptr + m_rgsw_size * 3;
    auto rgsw02 = rgsw_ptr + m_rgsw_size * 4;
    auto rgsw12 = rgsw_ptr + m_rgsw_size * 5;
    auto rgsw012 = rgsw_ptr + m_rgsw_size * 6;

    // First the "normal" ones
    auto mon_ptr = m_monomials_raw.data() + N * a_i;
    ExtProd(rgsw0, digit_ptr, mon_ptr, acc_ptr);
    intel::hexl::EltwiseAddMod(buffer.data(), acc_ptr, buffer.data(), 2 * N, modulus);

    mon_ptr = m_monomials_raw.data() + N * a_i1;
    ExtProd(rgsw1, digit_ptr, mon_ptr, acc_ptr);
    intel::hexl::EltwiseAddMod(buffer.data(), acc_ptr, buffer.data(), 2 * N, modulus);

    mon_ptr = m_monomials_raw.data() + N * a_i2;
    ExtProd(rgsw2, digit_ptr, mon_ptr, acc_ptr);
    intel::hexl::EltwiseAddMod(buffer.data(), acc_ptr, buffer.data(), 2 * N, modulus);

    // Next, the 2 variable ones
    mon_ptr = m_poly_buffer.data();

    MakePoly2(mon_ptr, a_i, a_i1);
    ExtProd(rgsw01, digit_ptr, mon_ptr, acc_ptr);
    intel::hexl::EltwiseAddMod(buffer.data(), acc_ptr, buffer.data(), 2 * N, modulus);

    MakePoly2(mon_ptr, a_i, a_i2);
    ExtProd(rgsw02, digit_ptr, mon_ptr, acc_ptr);
    intel::hexl::EltwiseAddMod(buffer.data(), acc_ptr, buffer.data(), 2 * N, modulus);

    MakePoly2(mon_ptr, a_i1, a_i2);
    ExtProd(rgsw12, digit_ptr, mon_ptr, acc_ptr);
    intel::hexl::EltwiseAddMod(buffer.data(), acc_ptr, buffer.data(), 2 * N, modulus);

    // Finally, the 3 variable one
    MakePoly3(mon_ptr, a_i, a_i1, a_i2);
    ExtProd(rgsw012, digit_ptr, mon_ptr, acc_ptr);
    intel::hexl::EltwiseAddMod(buffer.data(), acc_ptr, buffer.data(), 2 * N, modulus);
}

void FastBlindRotationKey::ExtProd(uint64_t* __restrict rgsw_block_start, uint64_t * __restrict digit_ptr, uint64_t * __restrict mon_ptr, uint64_t *acc_ptr) {
    // Note to future self and others:
    // This function copies a lot, which seems bad.
    // Previous implementations didn't e.g. now \dgg_ptr is copied twice instead of re-using the data
    // but this lead to MUCH worse performance
    // Reasons unclear, conjecture: better branch prediction, all data is single use afterwards so can be evicted from cache, ...
    // VTune also was confused...

    auto N = m_engine->GetDegree();
    auto modulus = m_engine->GetModulus();
    auto digits = m_params->GetDigitsG();
    auto basis = m_params->GetBaseG();
    auto basis_bits = GetMSB(basis) - 1;
    auto mask = basis - 1;

    const auto dgg_ptr = m_digit_buffer.data();
    auto acc_intt_ptr = m_acc_buffer.data();

    const auto rgsw_size = m_rgsw_size;
    const auto rlwe_prime_size = m_rlwe_prime_size;

    const auto& ntt_engine = m_engine;

    ntt_engine->ComputeInverse(acc_intt_ptr, acc_ptr, 1, 1);
    ntt_engine->ComputeInverse(acc_intt_ptr + N, acc_ptr + N, 1, 1);

    // digit decomp
    for (uint32_t i = 0; i < digits; i++) {
        auto shift = (digits - 1 - i) * basis_bits;
        auto digit_i_a_start = dgg_ptr + i * 2 * N;
        auto digit_i_b_start = dgg_ptr + i * 2 * N + rlwe_prime_size;

        for (uint32_t j = 0; j < N; j++) {
            digit_i_a_start[j] = (acc_intt_ptr[j] >> (shift)) & mask;
            digit_i_b_start[j] = (acc_intt_ptr[j + N] >> shift) & mask;
        }
        ntt_engine->ComputeForward(digit_i_a_start, digit_i_a_start, 1, 1);
        std::copy(digit_i_a_start, digit_i_a_start + N, digit_i_a_start + N);
        ntt_engine->ComputeForward(digit_i_b_start, digit_i_b_start, 1, 1);
        std::copy(digit_i_b_start, digit_i_b_start + N, digit_i_b_start + N);
    }

    // now repeat
    // currently this assumes step_size = 2
    std::copy(dgg_ptr, dgg_ptr + rgsw_size, dgg_ptr + rgsw_size);
    std::copy(dgg_ptr, dgg_ptr + rgsw_size, dgg_ptr + rgsw_size * 2);

    intel::hexl::EltwiseMultMod(dgg_ptr, rgsw_block_start, dgg_ptr, 3 * rgsw_size, modulus, 1);

    // reduce
    for (uint32_t rgsw_block_idx = 0; rgsw_block_idx < 3; rgsw_block_idx++) {
        auto rgsw_start = dgg_ptr + rgsw_block_idx * rgsw_size;
        for (uint32_t rlwe_prime_block_idx = 0; rlwe_prime_block_idx < 2; rlwe_prime_block_idx++) {
            auto accumulate_start = rgsw_start + rlwe_prime_block_idx * rlwe_prime_size;
            auto digit_start = accumulate_start;
            for (uint32_t digit_idx = 1; digit_idx < digits; digit_idx++) {
                auto digit_i_start = digit_start + digit_idx * 2 * N;
                for (uint32_t idx = 0; idx < 2 * N; idx ++) {
                    accumulate_start[idx] += digit_i_start[idx];
                }
            }
        }
    }

    // multiply by the monomials
    for (uint32_t rlwe_prime_idx  = 0; rlwe_prime_idx < 6; rlwe_prime_idx++) {
        auto rlwe_ptr = dgg_ptr + rlwe_prime_idx * rlwe_prime_size;
        auto mon_idx = rlwe_prime_idx >> 1;
        // do a final modular reduction
        intel::hexl::EltwiseReduceMod(rlwe_ptr, rlwe_ptr, 2 * N, modulus, modulus, 1);
        intel::hexl::EltwiseMultMod(rlwe_ptr, rlwe_ptr, mon_ptr + N * mon_idx, N, modulus, 1);
        intel::hexl::EltwiseMultMod(rlwe_ptr + N, rlwe_ptr + N, mon_ptr + N * mon_idx, N, modulus, 1);
        intel::hexl::EltwiseAddMod(acc_ptr, rlwe_ptr, acc_ptr, 2 * N, modulus);
    }

    /*

    // uses first 4 * N slots of buffer
    // writes result into first 2 * N
    auto rlwe_m_ptr = rgsw + m_rlwe_prime_size;
    auto rlwe_msm_ptr = rgsw;
    auto digits = ss_context->GetDigitsG();

    auto actual_acc = acc_ptr + 2 * N;
    std::fill(acc_ptr, acc_ptr + 4 * N, 0);

    for (uint32_t i = 0; i < digits; i++) {
        auto A_i = digit_ptr + i * 2 * N;
        auto B_i = A_i + N;

        intel::hexl::EltwiseMultMod(acc_ptr, A_i, rlwe_msm_ptr, N, modulus, 1);
        intel::hexl::EltwiseMultMod(acc_ptr + N, A_i, rlwe_msm_ptr + N, N, modulus, 1);

        for (uint32_t j = 0; j < 2 * N; j++) {
            actual_acc[j] += acc_ptr[j];
        }

        intel::hexl::EltwiseMultMod(acc_ptr, B_i, rlwe_m_ptr, N, modulus, 1);
        intel::hexl::EltwiseMultMod(acc_ptr + N, B_i, rlwe_m_ptr + N, N, modulus, 1);


        for (uint32_t j = 0; j < 2 * N; j++) {
            actual_acc[j] += acc_ptr[j];
        }

        rlwe_msm_ptr += 2 * N;
        rlwe_m_ptr += 2 * N;
    }

    intel::hexl::EltwiseReduceMod(actual_acc, actual_acc, 2 * N, modulus, modulus, 1);

    // monomial prod
    intel::hexl::EltwiseMultMod(actual_acc, actual_acc, mon_ptr, N, modulus,1 );
    intel::hexl::EltwiseMultMod(actual_acc + N, actual_acc + N, mon_ptr, N, modulus,1 );
    */

}

void FastBlindRotationKey::RLWEVectorProd(uint64_t *rlwe_prime, uint64_t *digit_polys, uint64_t *result_buffer) {

    const auto N = m_engine->GetDegree();
    const auto modulus = m_engine->GetModulus();
    const auto digits = m_params->GetDigitsG();
    auto shift = GetMSB(m_params->GetBaseG()) - 1;
    auto current_shift = (digits - 1) * shift;
    auto mask = (1u << shift) - 1;
    std::array<uint64_t, 1 << 11> tmp;
    std::copy(digit_polys, digit_polys + N, tmp.data());

    std::fill(result_buffer, result_buffer + 2 *N ,0);
    for (uint32_t i = 0; i < digits; i++) {
        for (uint32_t j = 0; j < N; j++) {
            auto v = (tmp[j] >> current_shift) & mask;
            digit_polys[i * 2 * N + j] = v;
        }
        m_engine->ComputeForward(digit_polys + i * 2 * N, digit_polys + i * 2 * N, 1, 1);
        std::copy(digit_polys + i * 2 * N, digit_polys + i * 2 * N + N, digit_polys + i * 2 * N + N);
        current_shift -= shift;
    }

    intel::hexl::EltwiseMultMod(digit_polys, digit_polys, rlwe_prime, 2 * N * digits, modulus, 1);
    for (int32_t i = 0; i < digits - 1; i++) {
        intel::hexl::EltwiseAddMod(digit_polys + (i + 1) * 2 * N, digit_polys + i * 2 * N, digit_polys + (i + 1) * 2 * N, 2 * N, modulus);
    }
    std::copy(digit_polys + (digits - 1) * 2 * N, digit_polys + digits * 2 * N, result_buffer);
}


void FastBlindRotationKey::RLWEVectorProd2(std::array<uint64_t* __restrict, 3> ops, std::array<uint64_t* __restrict, 3> mons, const uint64_t* __restrict digit_polys, uint64_t* __restrict result_buffer) {
    //assume result buffer has size

    const auto N = m_engine->GetDegree();
    const auto modulus = m_engine->GetModulus();
    const auto digits = m_params->GetDigitsG();
    auto shift = GetMSB(m_params->GetBaseG()) - 1;
    auto current_shift = (digits - 1) * shift;
    auto mask = (1u << shift) - 1;

    auto dgg_ptr = m_digit_buffer.data();

    const auto BLOCK_0 = dgg_ptr;
    const auto BLOCK_1 = dgg_ptr + 2 * N * digits;
    const auto BLOCK_2 = dgg_ptr + 4 * N * digits;

    for (uint32_t i = 0; i < digits; i++) {
        auto offset_0 = i * 2 * N;


        for (uint32_t j = 0; j < N; j++) {
            auto v = (digit_polys[j] >> current_shift) & mask;
            BLOCK_0[j + offset_0] = v;
        }
        m_engine->ComputeForward(BLOCK_0 + offset_0, BLOCK_0 + offset_0, 1, 1);

        std::copy(BLOCK_0 + offset_0, BLOCK_0 + offset_0 + N, BLOCK_0 + offset_0 + N);

        current_shift -= shift;
    }

    std::copy(BLOCK_0, BLOCK_1, BLOCK_1);
    std::copy(BLOCK_0, BLOCK_1, BLOCK_2);

    //m_engine->ComputeInverse(BLOCK_0 + N, BLOCK_0 + N, 1, 1);
    // TODO: Fuse memory regions
    //auto block_ptr = ops[0];
    intel::hexl::EltwiseMultMod(BLOCK_0, ops[0], BLOCK_0, 2 * N * digits, modulus, 1);
    intel::hexl::EltwiseMultMod(BLOCK_1, ops[1], BLOCK_1, 2 * N * digits, modulus, 1);
    intel::hexl::EltwiseMultMod(BLOCK_2, ops[2], BLOCK_2, 2 * N * digits, modulus, 1);



    for (uint32_t i = 1; i < digits; i++) {
        intel::hexl::EltwiseAddMod(BLOCK_0, BLOCK_0, BLOCK_0 + i * 2 * N, 2 * N, modulus);
        intel::hexl::EltwiseAddMod(BLOCK_1, BLOCK_1, BLOCK_1 + i * 2 * N, 2 * N, modulus);
        intel::hexl::EltwiseAddMod(BLOCK_2, BLOCK_2, BLOCK_2 + i * 2 * N, 2 * N, modulus);
    }

    intel::hexl::EltwiseMultMod(BLOCK_0, BLOCK_0, mons[0], N, modulus,1);
    intel::hexl::EltwiseMultMod(BLOCK_0 + N, BLOCK_0 + N, mons[0], N, modulus,1);
    intel::hexl::EltwiseAddMod(result_buffer, BLOCK_0, result_buffer, 2 * N, modulus);

    intel::hexl::EltwiseMultMod(BLOCK_1, BLOCK_1, mons[1], N, modulus,1);
    intel::hexl::EltwiseMultMod(BLOCK_1 + N, BLOCK_1 + N, mons[1], N, modulus,1);
    intel::hexl::EltwiseAddMod(result_buffer, BLOCK_1, result_buffer, 2 * N, modulus);

    intel::hexl::EltwiseMultMod(BLOCK_2, BLOCK_2, mons[2], N, modulus,1);
    intel::hexl::EltwiseMultMod(BLOCK_2 + N, BLOCK_2 + N, mons[2], N, modulus,1);
    intel::hexl::EltwiseAddMod(result_buffer, BLOCK_2, result_buffer, 2 * N, modulus);



    /*
    for (uint32_t i = 0; i < 3; i++) {
        for (uint32_t d_i = 0; d_i < digits; d_i ++) {
            intel::hexl::EltwiseMultMod(m_digit_buffer_2.data() + d_i * 2 * N, mons[i], m_digit_buffer.data() + d_i * N, N, modulus, 1);
            std::copy(m_digit_buffer_2.data() + d_i * 2 * N, m_digit_buffer_2.data() + d_i * 2 * N + N, m_digit_buffer_2.data() + d_i * 2 * N + N);
        }
        intel::hexl::EltwiseMultMod(m_digit_buffer_2.data(), m_digit_buffer_2.data(), ops[i], 2 * N * digits, modulus, 1);
        for (uint32_t d_i = 0; d_i < digits; d_i++) {
            intel::hexl::EltwiseAddMod(result_buffer, result_buffer, m_digit_buffer_2.data() + d_i * 2 * N, 2 * N, modulus);
        }
    }*/

    //AccVec<1<<11,3,4>(ops,mons,{m_digit_buffer_2.data(), m_digit_buffer.data()}, result_buffer, modulus, digits);

    /*
    m_engine->ComputeInverse(result_buffer + N, result_buffer + N, 1, 1);
    for (uint32_t i = 0; i < 2048; i += 64) {
        std::cerr << i << " | ";
        for (uint32_t j = 0; j < 64; j++) {
            std::cerr << result_buffer[N + i + j] << " ";
        }
        std::cerr << std::endl;
    } */
}



RLWECiphertext FastBlindRotationKey::BlindRotate1(RLWECiphertext &in, NativeVector A) {

    auto N = m_engine->GetDegree();

    Vector buffer(2 * N);

    for (uint32_t i = 0; i < N; i++) {
        buffer[i] = in->GetElements()[0][i].ConvertToInt<uint64_t>();
        buffer[i + N] = in->GetElements()[1][i].ConvertToInt<uint64_t>();
    }

    for (uint32_t i = 0; i < A.GetLength(); i++) {
        Step1(buffer, i, A[i].ConvertToInt<uint64_t>());
    }

    for (uint32_t i = 0; i < N; i++) {
        in->GetElements()[0][i] = buffer[i];
        in->GetElements()[1][i] = buffer[i + N];
    }

    return in;
}

RLWECiphertext FastBlindRotationKey::BlindRotate2(RLWECiphertext &in, NativeVector A) {
    auto N = m_engine->GetDegree();

    Vector buffer(2 * N);

    auto lwe_n = A.GetLength();
    for (uint32_t i = 0; i < N; i++) {
        buffer[i] = in->GetElements()[0][i].ConvertToInt<uint64_t>();
        buffer[i + N] = in->GetElements()[1][i].ConvertToInt<uint64_t>();
    }

    for (uint32_t i = 0; i < lwe_n; i+=2) {
        Step2(buffer, i, A[i].ConvertToInt<uint64_t>(), A[i + 1].ConvertToInt<uint64_t>());
    }

    if (lwe_n & 1) {
        Step1(buffer, lwe_n - 1, A[lwe_n - 1].ConvertToInt<uint64_t>());
    }

    for (uint32_t i = 0; i < N; i++) {
        in->GetElements()[0][i] = buffer[i];
        in->GetElements()[1][i] = buffer[i + N];
    }

    return in;
}

RLWECiphertext FastBlindRotationKey::BlindRotate3(RLWECiphertext &in, NativeVector A) {
    auto N = m_engine->GetDegree();

    Vector buffer(2 * N);

    auto lwe_n = A.GetLength();
    for (uint32_t i = 0; i < N; i++) {
        buffer[i] = in->GetElements()[0][i].ConvertToInt<uint64_t>();
        buffer[i + N] = in->GetElements()[1][i].ConvertToInt<uint64_t>();
    }

    for (uint32_t i = 0; i < lwe_n; i+=3) {
        Step3(buffer, i, A[i].ConvertToInt<uint64_t>(), A[i + 1].ConvertToInt<uint64_t>(), A[i + 2].ConvertToInt<uint64_t>());
    }

    if (lwe_n % 3 != 0) {
        Step1(buffer, lwe_n - 2, A[lwe_n - 2].ConvertToInt<uint64_t>());
    }

    if (lwe_n % 3 == 2) {
        std::cerr << "NOPE" << std::endl;
        std::exit(-1);
    }

    for (uint32_t i = 0; i < N; i++) {
        in->GetElements()[0][i] = buffer[i];
        in->GetElements()[1][i] = buffer[i + N];
    }

    return in;
}

RLWECiphertext FastBlindRotationKey::BlindRotate(RLWECiphertext &in, NativeVector A) {
    if (m_step_size == 1) {
        return BlindRotate1(in, A);
    }

    if (m_step_size == 2) {
        return BlindRotate2(in, A);
    }

    return BlindRotate3(in, A);
}
