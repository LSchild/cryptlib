//
// Created by leonard on 11/12/24.
//

#include <cassert>
#include <memory>
#include "rgsw.h"
#include "hexl/hexl.hpp"

// Long story
#define MAX_DIGITS 16
#define DIGIT_T uint64_t
alignas(4096) DIGIT_T DECOMP_BUFFER [(1u << 11) * 2 * MAX_DIGITS];
DIGIT_T dig_buffer[MAX_DIGITS * 2];

RingGSWSample::RingGSWSample(const std::shared_ptr<RingGSWCryptoParams>& params,
                             std::vector<RLWECiphertext> &msm, std::vector<RLWECiphertext> &m) : m_params(params), m_msm(msm), m_m(m) {
    assert(m_m.size() == m_msm.size());
}

RingGSWSample::RingGSWSample(const RingGSWSample& other) {
    m_m.insert(m_m.begin(), other.m_m.begin(), other.m_m.end());
    m_msm.insert(m_msm.begin(), other.m_msm.begin(), other.m_msm.end());
    m_params = other.m_params;
}

RLWECiphertext RingGSWSample::mul(const lbcrypto::RLWECiphertext &rhs) const{

    // TODO: do it with hexl
    // TODO fixme, we're using the incorrect bassi
    auto dig = m_params->GetDigitsG();
    auto base = m_params->GetBaseG();
    auto basebits = GetMSB(base) - 1;
    uint32_t mask = (1 << (basebits)) - 1;

    auto actual_digit_count = m_m.size();
    auto offset = dig - actual_digit_count;

    NativePoly A = rhs->GetElements()[0];
    NativePoly B = rhs->GetElements()[1];

    A.SetFormat(COEFFICIENT);
    B.SetFormat(COEFFICIENT);

    assert(actual_digit_count <= MAX_DIGITS);
    //assert(A.GetModulus() == ss_context->GetQ());
    //assert(B.GetModulus() == ss_context->GetQ());

    auto N = A.GetLength();
    std::vector<NativePoly> digs(2 * dig, {m_params->GetPolyParams(), COEFFICIENT, true});

    auto acc_B = NativePoly(m_params->GetPolyParams(), EVALUATION, true);
    auto acc_A = NativePoly(m_params->GetPolyParams(), EVALUATION, true);

    // Current fastest version
    for(long i = 0; i < N; i++) {
        auto val_A_i = A.at(i).ConvertToInt<uint64_t>();
        auto val_B_i = B.at(i).ConvertToInt<uint64_t>();;
        for(short j = 0; j < MAX_DIGITS; j++) {
            dig_buffer[2 * j] = (val_A_i >> (j * basebits)) & mask;
            dig_buffer[2 * j + 1] = (val_B_i >> (j * basebits)) & mask;
        }
        std::memcpy(DECOMP_BUFFER + i * 2 * MAX_DIGITS, dig_buffer, 2 * MAX_DIGITS * sizeof(DIGIT_T) );
    }

    for(long j = 0; j < dig; j++) {
        auto& dig_A_j = digs.at(2 * j);
        auto& dig_B_j = digs.at(2 * j + 1);
        //auto scal = NativeInteger(base).Exp(j);
        for(long i = 0; i < N; i++) {
            // NativePoly::at is an actual bottleneck wtf
            dig_A_j.at(i) = DECOMP_BUFFER[i * 2 * MAX_DIGITS + 2 * j];
            dig_B_j.at(i) = DECOMP_BUFFER[i * 2 * MAX_DIGITS + 2 * j + 1];
        }
        // TODO: remove me (sanity check)
        //check_A += dig_A_j * scal;
        //check_B += dig_B_j * scal;
    }

    // TODO: remove me
    //  ssert(check_A == A);
    // assert(check_B == B);

    // EXP END

    // TODO: enable parallelization later
    //#pragma omp parallel for
    for(long i = 2u * offset; i < digs.size(); i++) {
        digs.at(i).SetFormat(EVALUATION);
    }

    //auto decomp_A = A.BaseDecompose(basebits, true);
    //auto decomp_B = B.BaseDecompose(basebits, true);

    //assert(actual_digit_count % 2 == 0);
    for(int i = 0; i < actual_digit_count; i++) {

        auto idx_base = 2 * (i + offset);
        acc_A += m_msm.at(i)->GetElements()[0] * digs.at(idx_base);
        acc_B += digs.at(idx_base) *= m_msm.at(i)->GetElements()[1];

        acc_A += m_m.at(i)->GetElements()[0] * digs.at(idx_base + 1);
        acc_B += digs.at(idx_base + 1) *= m_m.at(i)->GetElements()[1];

    }

    std::vector<NativePoly> res_vec = {acc_A, acc_B};
    return std::make_shared<RLWECiphertextImpl>(res_vec);
}

void RingGSWSample::mul_inplace(lbcrypto::RLWECiphertext &rhs) const {

    auto dig = m_params->GetDigitsG();
    auto base = m_params->GetBaseG();
    auto basebits = GetMSB(base) - 1;
    uint32_t mask = (1 << (basebits)) - 1;

    auto actual_digit_count = m_m.size();
    auto offset = dig - actual_digit_count;

    NativePoly A = rhs->GetElements()[0];
    NativePoly B = rhs->GetElements()[1];

    A.SetFormat(COEFFICIENT);
    B.SetFormat(COEFFICIENT);

    auto N = A.GetLength();
    std::vector<uint64_t> full_ct(2 * N, 0);
    for(long i =0; i < N; i++) {
        full_ct[i] = A.at(i).ConvertToInt<uint64_t>();
        full_ct[i + N] = B.at(i).ConvertToInt<uint64_t>();
    }
    auto Q = m_params->GetQ().ConvertToInt<uint64_t>();
    auto root = m_params->GetPolyParams()->GetRootOfUnity().ConvertToInt<uint64_t>();
    auto nttengine = intel::hexl::NTT(N, Q, root);
    for(long j = 0; j < dig; j++) {
        auto write_offset = DECOMP_BUFFER + j * 2 * N;
        for(long k = 0; k < 2 * N; k++) {
            write_offset[k] = (full_ct[k] >> (j * basebits)) & mask;
        }
    }

    std::vector<NativePoly> digs;
    for(long i = offset; i < dig; i++) {
        auto ct_a_dig = DECOMP_BUFFER + i * 2 * N;
        auto ct_b_dig = DECOMP_BUFFER + i * 2 * N + N;
        nttengine.ComputeForward(ct_a_dig, ct_a_dig, 1, 1);
        nttengine.ComputeForward(ct_b_dig, ct_b_dig, 1, 1);

        NativeVector Adig(N, Q);
        NativeVector Bdig(N, Q);
        for(long l = 0; l < N; l++) {
            Adig.at(l) = ct_a_dig[l];
            Bdig.at(l) = ct_b_dig[l];
        }
        NativePoly x(m_params->GetPolyParams());
        NativePoly y(m_params->GetPolyParams());
        x.SetValues(Adig, EVALUATION);
        y.SetValues(Bdig, EVALUATION);
        digs.push_back(x);
        digs.push_back(y);
    }

    auto& ret_A = rhs->GetElements()[0];
    auto& ret_B = rhs->GetElements()[1];

    // do first digit manually
    ret_A = m_m.at(0)->GetElements()[0] * digs.at(2 * (0) + 1);
    ret_B = digs.at(2 * (0 + offset) + 1) *= m_m.at(0)->GetElements()[1];

    ret_A += m_msm.at(0)->GetElements()[0] * digs.at(2 * (0));
    ret_B += digs.at(2 * (0)) *= m_msm.at(0)->GetElements()[1];

    for(int i = 1; i < actual_digit_count; i++) {

        ret_A += m_m.at(i)->GetElements()[0] * digs.at(2 * (i) + 1);
        ret_B += digs.at(2 * (i) + 1) *= m_m.at(i)->GetElements()[1];

        ret_A += m_msm.at(i)->GetElements()[0] * digs.at(2 * (i));
        ret_B += digs.at(2 * (i)) *= m_msm.at(i)->GetElements()[1];

    }
}

RLWECiphertext RingGSWSample::poly_mul_with_precomp(std::vector<NativePoly> &precomp_digs) {

    // if we multiply by polynomials rgsw <=> rlwe'

    auto dig = m_params->GetDigitsG();

    auto acc_B = NativePoly(m_params->GetPolyParams(), EVALUATION, true);
    auto acc_A = NativePoly(m_params->GetPolyParams(), EVALUATION, true);

    auto actual_digit_count = m_m.size();
    auto offset = dig - actual_digit_count;

    for(int i = 0; i < actual_digit_count; i++) {

        auto& decomp = precomp_digs.at(i + offset);

        acc_A += m_m.at(i)->GetElements()[0] * decomp;
        acc_B += decomp *= m_m.at(i)->GetElements()[1];
    }

    std::vector<NativePoly> res_vec = {acc_A, acc_B};
    return std::make_shared<RLWECiphertextImpl>(res_vec);
}

void RingGSWSample::cmux(RLWECiphertext &in, uint32_t idx) const{

    NativePoly acc_A = in->GetElements()[0].Clone();
    NativePoly acc_B = in->GetElements()[1].Clone();

    auto N = acc_B.GetLength();

    NativePoly monA(m_params->GetPolyParams(), COEFFICIENT, true);
    monA.at(idx % N) = idx >= N ? (acc_B.GetModulus() - 1 ) : 1;
    monA.at(0) += acc_B.GetModulus() -1;
    monA.SetFormat(EVALUATION);

    in->GetElements()[0] *= monA;
    in->GetElements()[1] *= monA;

    //mul_inplace(in);
    //in->GetElements()[0] += acc_A;
    //in->GetElements()[1] += acc_B;

    auto res = mul(in);

    in->GetElements()[0] = (acc_A += res->GetElements()[0]);
    in->GetElements()[1] = (acc_B += res->GetElements()[1]);

}


void RingGSWSample::cmux(RLWECiphertext &in, NativePoly& mon) const{

    NativePoly acc_A = in->GetElements()[0].Clone();
    NativePoly acc_B = in->GetElements()[1].Clone();

    in->GetElements()[0] *= mon;
    in->GetElements()[1] *= mon;

    auto res = mul(in);

    in->GetElements()[0] = (acc_A += res->GetElements()[0]);
    in->GetElements()[1] = (acc_B += res->GetElements()[1]);

}

std::vector<RLWECiphertext>
RingGSWSample::MultiCMUX(std::vector<RLWECiphertext> &choices) {
    // TODO: Does not work atm but also isn't used.
    // setup
    auto poly_params = m_params->GetPolyParams();
    auto Q = m_params->GetQ().ConvertToInt<uint64_t>();
    auto N = m_params->GetN();
    auto root = poly_params->GetRootOfUnity().ConvertToInt<uint64_t>();
    auto digs = m_params->GetDigitsG();
    auto ntt_engine = intel::hexl::NTT(Q, N, root);

    // bases
    auto offset = choices.size() / 2;
    auto dig = m_params->GetDigitsG();
    (void) dig;
    auto base = m_params->GetBaseG();
    auto basebits = GetMSB(base) - 1;
    uint32_t mask = (1 << (basebits)) - 1;

    // allocate buffers
    auto ct_buffer = new uint64_t[2 * N * choices.size()];
    auto lhs_buffer = ct_buffer;
    auto rhs_buffer = ct_buffer + choices.size() * N;
    auto rgsw_ct_buffer = new uint64_t[4 * 2 * N];
    auto digit_buffer = new uint64_t[digs * 2 * N * offset];

    // rearrange memory first
    for(long i =0; i < offset; i++) {
        auto ct_hi = choices.at(i + offset);
        auto ct_lo = choices.at(i);

        for(long l = 0; l < N; l++) {
            lhs_buffer[i * 2 * N + l] = ct_hi->GetElements()[0].at(l).ConvertToInt();
            lhs_buffer[i * 2 * N + N + l] = ct_hi->GetElements()[1].at(l).ConvertToInt();
            rhs_buffer[i * 2 * N + l] = ct_lo->GetElements()[0].at(l).ConvertToInt();
            rhs_buffer[i * 2 * N + N + l] = ct_lo->GetElements()[1].at(l).ConvertToInt();
        }
    }

    // compute difference
    intel::hexl::EltwiseSubMod(lhs_buffer, lhs_buffer, rhs_buffer, N * choices.size(), Q);

    // base decompose
    for(long l = 0; l < offset; l++) {
        auto write_ct_offset = digit_buffer + l * 2 * N * digs;
        for(long k = 0; k < digs; k++) {
            auto dig_ct_offset = write_ct_offset + k * 2 * N;
            for(long j = 0; j < 2 * N; j++) {
                dig_ct_offset[j] = (lhs_buffer[l * 2 * N + j] >> (k * basebits)) & mask;
            }
        }
    }

    // do the NTTs
    for(long l = 0; l < 2 * offset * digs; l++) {
        ntt_engine.ComputeForward(digit_buffer + l * N, digit_buffer + l * N, 1, 1);
    }

    delete[] rgsw_ct_buffer;
    delete[] ct_buffer;
    delete[] digit_buffer;

    return {};
}

RingGSWSample RingGSWSample::create_one() {

    auto dig = m_params->GetDigitsG();
    auto base = m_params->GetBaseG();
    auto Q = m_params->GetQ();

    auto actual_digit_count = m_m.size();
    auto offset = dig - actual_digit_count;

    std::vector<RLWECiphertext> lhs(actual_digit_count);
    std::vector<RLWECiphertext> rhs(actual_digit_count);

    NativeInteger first_power = NativeInteger(base).ModExp(offset, Q);
    auto zero = NativePoly(m_params->GetPolyParams(), EVALUATION, true);

    auto left = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);
    auto right = NativePoly(m_params->GetPolyParams(), COEFFICIENT, true);

    left[0] = first_power;
    right[0] = first_power;

    left.SwitchFormat();
    right.SwitchFormat();

    for(int i = 0; i < actual_digit_count; i++) {

        auto ct_left_args = {left.Clone(), zero.Clone()};
        auto ct_right_args = {zero.Clone(), right.Clone()};

        auto ct_left = std::make_shared<RLWECiphertextImpl>(ct_left_args);
        auto ct_right = std::make_shared<RLWECiphertextImpl>(ct_right_args);

        lhs[i] = ct_left;
        rhs[i] = ct_right;

        left *= base;
        right *= base;
    }

    return {m_params, lhs, rhs};
}

RingGSWSample RingGSWSample::flip_bit(RingGSWSample &one) {

    RingGSWSample cloned(one);

    for(int i = 0; i < m_m.size(); i++) {
        cloned.m_m[i]->GetElements()[0] -= m_m[i]->GetElements()[0];
        cloned.m_m[i]->GetElements()[1] -= m_m[i]->GetElements()[1];
        cloned.m_msm[i]->GetElements()[0] -= m_msm[i]->GetElements()[0];
        cloned.m_msm[i]->GetElements()[1] -= m_msm[i]->GetElements()[1];
    }

    return cloned;
}

RGSWSample::RGSWSample(std::shared_ptr<intel::hexl::NTT> &ntt_engine, uint32_t L,
                       uint32_t digits, const uint64_t *msm, const uint64_t *m) : m_engine(ntt_engine), m_basis(L), m_digits(digits) {

    // NOTE: EXPECTED ORDER OF DIGITS IS HI TO LO

    m_basis_bits = GetMSB(L) - 1;
    m_key.resize(ntt_engine->GetDegree() * 4 * digits);
    m_scratch_space.resize(ntt_engine->GetDegree() * 4);
    m_mask = (1 << m_basis_bits) - 1;

    double modulus_bits = GetMSB(ntt_engine->GetModulus()) - 1;
    double ceiled_digits = std::ceil(modulus_bits / m_basis_bits);
    m_first_shift = m_basis_bits * uint64_t(ceiled_digits - 1);

    std::copy(m, m + 2 * digits * ntt_engine->GetDegree(), m_key.data());
    std::copy(msm, msm + 2 * digits * ntt_engine->GetDegree(), m_key.data());
}

RGSWSample::RGSWSample(std::shared_ptr<intel::hexl::NTT> &ntt_engine, uint32_t L, uint32_t digits,
                       const std::vector<RLWECiphertext> &msm, const std::vector<RLWECiphertext> &m) : m_engine(ntt_engine), m_basis(L), m_digits(digits) {

    // NOTE: EXPECTED ORDER OF DIGITS IS HI TO LO
    // TODO: Maybe use general scratch space instead of for each RGSW Sample?

    m_basis_bits = GetMSB(L) - 1;
    m_key.resize(ntt_engine->GetDegree() * 4 * digits);
    m_scratch_space.resize(ntt_engine->GetDegree() * 4);
    m_mask = (1 << m_basis_bits) - 1;

    double modulus_bits = GetMSB(ntt_engine->GetModulus()) - 1;
    double ceiled_digits = std::ceil(modulus_bits / m_basis_bits);
    m_first_shift = m_basis_bits * uint64_t(ceiled_digits - 1);

    auto data = m_key.data();
    auto N = ntt_engine->GetDegree();
    for(uint32_t i = 0; i < digits; i++) {
        auto digit_m = msm.at(digits - i - 1);
        auto A = digit_m->GetElements()[0];
        auto B = digit_m->GetElements()[1];
        for(uint64_t j = 0; j < N; j++) {
            data[i * 2 * N + j] = A[j].ConvertToInt();
            data[i * 2 * N + j + N] = B[j].ConvertToInt();
        }
    }

    data = m_key.data() + 2 * N * digits;
    for(uint32_t i = 0; i < digits; i++) {
        auto digit_m = m.at(digits - i - 1);
        auto A = digit_m->GetElements()[0];
        auto B = digit_m->GetElements()[1];
        for(uint64_t j = 0; j < N; j++) {
            data[i * 2 * N + j] = A[j].ConvertToInt();
            data[i * 2 * N + j + N] = B[j].ConvertToInt();
        }
    }
}

void RGSWSample::mul_dir(uint64_t *res, const uint64_t* poly, bool mul_left, bool add_to_result) {

    // Uses scratch_space [0, 2N)

    auto N = m_engine->GetDegree();
    auto Q = m_engine->GetModulus();

    uint64_t* data;
    uint64_t* scratch_space = m_scratch_space.data();

    if (mul_left) {
        data = m_key.data() ;
    } else {
        data = m_key.data() + 2 * m_digits * N;
    }

    auto current_shift = m_first_shift;
    auto loop_start = 0;

    if (add_to_result) {
        loop_start = 0;
    } else {

        loop_start = 1;

        for (uint64_t x_i = 0; x_i < N; x_i++) {
            scratch_space[x_i] = (poly[x_i] >> current_shift) & m_mask;
        }

        m_engine->ComputeForward(scratch_space, scratch_space, 1, 1);
        intel::hexl::EltwiseMultMod(res, data, scratch_space, N, Q, 1);
        intel::hexl::EltwiseMultMod(res + N, data + N, scratch_space, N, Q, 1);

        current_shift -= m_basis_bits;
        data += 2 * N;

    }

    for(uint32_t d_i = loop_start; d_i < m_digits; d_i++) {

        for (uint64_t x_i = 0; x_i < N; x_i++) {
            scratch_space[x_i] = (poly[x_i] >> current_shift) & m_mask;
        }

        // Apply NTT to current digit
        m_engine->ComputeForward(scratch_space, scratch_space, 1, 1);
        // multiply by A component of RLWE' sample
        intel::hexl::EltwiseMultMod(scratch_space + N, data, scratch_space, N, Q, 1);
        // Add A * d_i to A component of accumulator
        intel::hexl::EltwiseAddMod(res, res, scratch_space + N, N, Q);
        // multiply digit with B component
        intel::hexl::EltwiseMultMod(scratch_space, data + N, scratch_space, N, Q, 1);
        // add to B component
        intel::hexl::EltwiseAddMod(res + N, scratch_space, res + N, N, Q);

        current_shift -= m_basis_bits;
        data += 2 * N;
    }

}

void RGSWSample::mul_poly(uint64_t *__restrict result, const uint64_t *__restrict poly, bool add_to_result) {

    mul_dir(result, poly, false, add_to_result);

}

void RGSWSample::mul(uint64_t *__restrict result, const uint64_t *__restrict rhs, bool add_to_result) {

    // uses scratch_space [0, 2N)

    mul_dir(result, rhs, true, add_to_result);
    mul_dir(result, rhs + m_engine->GetDegree(), false, true);

}

RLWECiphertext RGSWSample::mul_poly(const lbcrypto::NativePoly &poly) {

    auto N = m_engine->GetDegree();
    std::vector<uint64_t> buffer(N * 3);
    for(uint64_t i = 0; i < N; i++) {
        buffer[i] = poly[i].ConvertToInt();
    }

    mul_poly(buffer.data() + N, buffer.data());

    auto vec = NativeVector(N, poly.GetModulus());

    auto A = NativePoly(poly.GetParams(), EVALUATION, false);
    auto B = NativePoly(poly.GetParams(), EVALUATION, false);

    for(uint64_t i = 0; i < N; i++) {
        vec[i] = buffer[i + N];
    }

    A.SetValues(vec, EVALUATION);

    for(uint64_t i = 0; i < N; i++) {
        vec[i] = buffer[i + 2 * N];
    }

    B.SetValues(vec, EVALUATION);

    auto ct_vec = {A, B};
    return std::make_shared<RLWECiphertextImpl>(ct_vec);

}

RLWECiphertext RGSWSample::mul(const lbcrypto::RLWECiphertext &rhs) {
    auto N = m_engine->GetDegree();

    std::vector<uint64_t> buffer(N * 4);
    for(uint64_t i = 0; i < N; i++) {
        buffer[i] = rhs->GetElements()[0][i].ConvertToInt();
        buffer[i + N] = rhs->GetElements()[1][i].ConvertToInt();
    }

    if (rhs->GetElements()[0].GetFormat() != COEFFICIENT) {
        m_engine->ComputeInverse(buffer.data(), buffer.data(), 1, 1);
        m_engine->ComputeInverse(buffer.data() + N, buffer.data() + N, 1, 1);
    }

    mul(buffer.data() + 2 * N, buffer.data());

    auto vec = NativeVector(N, m_engine->GetModulus());

    auto A = NativePoly(rhs->GetElements()[0].GetParams(), EVALUATION, false);
    auto B = NativePoly(rhs->GetElements()[0].GetParams(), EVALUATION, false);

    for(uint64_t i = 0; i < N; i++) {
        vec[i] = buffer[i + 2 * N];
    }

    A.SetValues(vec, EVALUATION);

    for(uint64_t i = 0; i < N; i++) {
        vec[i] = buffer[i + 3 * N];
    }

    B.SetValues(vec, EVALUATION);

    auto ct_vec = {A, B};
    return std::make_shared<RLWECiphertextImpl>(ct_vec);
}

void RGSWSample::cmux(uint64_t *in_out, uint64_t *mon) {

    auto N = m_engine->GetDegree();
    auto Q = m_engine->GetModulus();

    intel::hexl::EltwiseMultMod(m_scratch_space.data() + 2 * N, in_out, mon, N, Q, 1);
    intel::hexl::EltwiseMultMod(m_scratch_space.data() + 3 * N, in_out + N, mon, N, Q, 1);

    m_engine->ComputeInverse(m_scratch_space.data() + 2 * N, m_scratch_space.data() + 2 * N, 1, 1);
    m_engine->ComputeInverse(m_scratch_space.data() + 3 * N, m_scratch_space.data() + 3 * N, 1, 1);

    mul(in_out, m_scratch_space.data() + 2 * N, true);

}

void RGSWSample::cmux(uint64_t *in_out, uint32_t idx) {

    auto N = m_engine->GetDegree();
    auto Q = m_engine->GetModulus();
    for(uint32_t i = 0; i < N; i++) {
        m_scratch_space[i + 2 * N] = 0;
    }

    m_scratch_space[2 * N] = Q - 1;
    m_scratch_space[2 * N + (idx % N)] += idx >= N ? Q - 1 : 1;

    m_engine->ComputeForward(m_scratch_space.data() + 3 * N, m_scratch_space.data() + 2 * N, 1, 1);

    intel::hexl::EltwiseMultMod(m_scratch_space.data() + 2 * N, in_out, m_scratch_space.data() + 3 * N, N, Q, 1);
    intel::hexl::EltwiseMultMod(m_scratch_space.data() + 3 * N, in_out + N, m_scratch_space.data() + 3 * N, N, Q, 1);

    mul(in_out, m_scratch_space.data() + 2 * N, true);

}

void RGSWSample::cmux(lbcrypto::RLWECiphertext &in, lbcrypto::NativePoly &mon) {
}

void RGSWSample::cmux(RLWECiphertext &in, uint32_t idx) {

    auto tmpA = in->GetElements()[0].Clone();
    auto tmpB = in->GetElements()[1].Clone();

    auto N = m_engine->GetDegree();
    NativePoly mon(in->GetElements()[0].GetParams(), COEFFICIENT, true);
    if (idx >= N) {
        idx -= N;
        mon[idx] = m_engine->GetModulus() - 1;
    } else {
        mon[idx] = 1;
    }
    mon[0] += m_engine->GetModulus() - 1;
    mon.SwitchFormat();

    in->GetElements()[0] *= mon;
    in->GetElements()[1] *= mon;

    auto prod = mul(in);

    in->GetElements()[0] = prod->GetElements()[0] + tmpA;
    in->GetElements()[1] = prod->GetElements()[1] + tmpB;

}


RGSWSample RGSWSample::convert(const RingGSWSample &old_version) {
    auto old_params = old_version.m_params;
    uint32_t digit_count = old_version.m_m.size();
    auto basis = old_params->GetBaseG();
    auto N = old_params->GetN();
    auto Q = old_params->GetQ().ConvertToInt();
    auto root = old_params->GetPolyParams()->GetRootOfUnity().ConvertToInt();

    auto ntt_engine = std::make_shared<intel::hexl::NTT>(N,Q,root);

    return {ntt_engine, basis, digit_count, old_version.m_msm, old_version.m_m};
}

void RGSWSample::ternary_mux(uint64_t *in_out, uint64_t *mon_pos, uint64_t *mon_neg, RGSWSample &high_bit) {

    auto N = m_engine->GetDegree();
    auto modulus = m_engine->GetModulus();
    auto scratch_space = m_scratch_space.data();

    // During Phase 1 we need to perform a large amount of modular additions
    // In practice log(modulus) < 64, so we can delay the modular reduction
    // Here we determine the number of additions we can do before it becomes necessary
    auto Qbits = GetMSB(modulus) - 1;
    auto free_bits = 62 - Qbits;
    auto n_additions_before_reduction =  1u << free_bits;
    const uint32_t tile_size = 64;

    assert(m_digits < n_additions_before_reduction);

    auto intt_buffer = scratch_space;
    auto digit_buffer = scratch_space + N;
    auto rlwe_prod_buffer_A = scratch_space + 2 * N;
    auto rlwe_prod_buffer_B = scratch_space + 3 * N;

    // Important: keep buffers in continuous memory
    auto accum_space = high_bit.m_scratch_space.data();
    auto lo_accumulator_A = accum_space;
    auto lo_accumulator_B = accum_space + N;
    auto hi_accumulator_A = accum_space + 2 * N;
    auto hi_accumulator_B = accum_space + 3 * N;

    std::fill(high_bit.m_scratch_space.begin(), high_bit.m_scratch_space.end(), 0);

    auto lo_bit_key_left = m_key.data();
    auto lo_bit_key_right = m_key.data() + 2 * N * m_digits;

    auto hi_bit_key_left = high_bit.m_key.data();
    auto hi_bit_key_right = high_bit.m_key.data() + 2 * N * m_digits;

    m_engine->ComputeInverse(intt_buffer, in_out + N, 1, 1);

    auto current_shift = m_first_shift;
    for(uint32_t d_i = 0; d_i < m_digits; d_i++) {

        for (uint64_t x_i = 0; x_i < N; x_i++) {
            digit_buffer[x_i] = (intt_buffer[x_i] >> current_shift) & m_mask;
        }

        m_engine->ComputeForward(digit_buffer, digit_buffer, 1, 1);

        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_A, lo_bit_key_right, digit_buffer, N, modulus, 1);
        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_B, lo_bit_key_right + N, digit_buffer, N, modulus, 1);

        for(uint32_t i = 0; i < N; i+= tile_size) {
            for(uint32_t k = 0; k < tile_size; k++) {
                lo_accumulator_A[i + k] += rlwe_prod_buffer_A[i + k];
            }

            for(uint32_t k = 0; k < tile_size; k++) {
                lo_accumulator_B[i + k] += rlwe_prod_buffer_B[i + k];
            }
        }

        //intel::hexl::EltwiseAddMod(lo_accumulator_A, rlwe_prod_buffer_A, lo_accumulator_A, N, modulus);
        //intel::hexl::EltwiseAddMod(lo_accumulator_B, rlwe_prod_buffer_B, lo_accumulator_B, N, modulus);

        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_A, hi_bit_key_right, digit_buffer, N, modulus, 1);
        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_B, hi_bit_key_right + N, digit_buffer, N, modulus, 1);

        //intel::hexl::EltwiseAddMod(hi_accumulator_A, rlwe_prod_buffer_A, hi_accumulator_A, N, modulus);
        //intel::hexl::EltwiseAddMod(hi_accumulator_B, rlwe_prod_buffer_B, hi_accumulator_B, N, modulus);

        for(uint32_t i = 0; i < N; i+= tile_size) {
            for(uint32_t k = 0; k < tile_size; k++) {
                hi_accumulator_A[i + k] += rlwe_prod_buffer_A[i + k];
            }

            for(uint32_t k = 0; k < tile_size; k++) {
                hi_accumulator_B[i + k] += rlwe_prod_buffer_B[i + k];
            }
        }

        lo_bit_key_right += 2 * N;
        hi_bit_key_right += 2 * N;
        current_shift -= m_basis_bits;


    }

    intel::hexl::EltwiseReduceMod(accum_space, accum_space, N, modulus, modulus, 1);

    intel::hexl::EltwiseMultMod(lo_accumulator_A, lo_accumulator_A, mon_pos, N, modulus, 1);
    intel::hexl::EltwiseMultMod(lo_accumulator_B, lo_accumulator_B, mon_pos, N, modulus, 1);

    intel::hexl::EltwiseMultMod(hi_accumulator_A, hi_accumulator_A, mon_neg, N, modulus, 1);
    intel::hexl::EltwiseMultMod(hi_accumulator_B, hi_accumulator_B, mon_neg, N, modulus, 1);

    m_engine->ComputeInverse(intt_buffer, in_out, 1, 1);

    intel::hexl::EltwiseAddMod(in_out, accum_space, in_out, 2 * N, modulus);
    intel::hexl::EltwiseAddMod(in_out, accum_space + 2 * N, in_out, 2 * N, modulus);

    std::fill(high_bit.m_scratch_space.begin(), high_bit.m_scratch_space.end(), 0);

    current_shift = m_first_shift;
    for(uint32_t d_i = 0; d_i < m_digits; d_i++) {

        for (uint64_t x_i = 0; x_i < N; x_i++) {
            digit_buffer[x_i] = (intt_buffer[x_i] >> current_shift) & m_mask;
        }

        m_engine->ComputeForward(digit_buffer, digit_buffer, 1, 1);

        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_A, lo_bit_key_left, digit_buffer, N, modulus, 1);
        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_B, lo_bit_key_left + N, digit_buffer, N, modulus, 1);

        for(uint32_t i = 0; i < N; i+= tile_size) {
            for(uint32_t k = 0; k < tile_size; k++) {
                lo_accumulator_A[i + k] += rlwe_prod_buffer_A[i + k];
            }

            for(uint32_t k = 0; k < tile_size; k++) {
                lo_accumulator_B[i + k] += rlwe_prod_buffer_B[i + k];
            }
        }

        //intel::hexl::EltwiseAddMod(lo_accumulator_A, rlwe_prod_buffer_A, lo_accumulator_A, N, modulus);
        //intel::hexl::EltwiseAddMod(lo_accumulator_B, rlwe_prod_buffer_B, lo_accumulator_B, N, modulus);

        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_A, hi_bit_key_left, digit_buffer, N, modulus, 1);
        intel::hexl::EltwiseMultMod(rlwe_prod_buffer_B, hi_bit_key_left + N, digit_buffer, N, modulus, 1);

        for(uint32_t i = 0; i < N; i+= tile_size) {
            for(uint32_t k = 0; k < tile_size; k++) {
                hi_accumulator_A[i + k] += rlwe_prod_buffer_A[i + k];
            }

            for(uint32_t k = 0; k < tile_size; k++) {
                hi_accumulator_B[i + k] += rlwe_prod_buffer_B[i + k];
            }
        }

        //intel::hexl::EltwiseAddMod(hi_accumulator_A, rlwe_prod_buffer_A, hi_accumulator_A, N, modulus);
        //intel::hexl::EltwiseAddMod(hi_accumulator_B, rlwe_prod_buffer_B, hi_accumulator_B, N, modulus);

        lo_bit_key_left += 2 * N;
        hi_bit_key_left += 2 * N;
        current_shift -= m_basis_bits;

    }

    intel::hexl::EltwiseReduceMod(accum_space, accum_space, N, modulus, modulus, 1);

    intel::hexl::EltwiseMultMod(lo_accumulator_A, lo_accumulator_A, mon_pos, N, modulus, 1);
    intel::hexl::EltwiseMultMod(lo_accumulator_B, lo_accumulator_B, mon_pos, N, modulus, 1);

    intel::hexl::EltwiseMultMod(hi_accumulator_A, hi_accumulator_A, mon_neg, N, modulus, 1);
    intel::hexl::EltwiseMultMod(hi_accumulator_B, hi_accumulator_B, mon_neg, N, modulus, 1);

    intel::hexl::EltwiseAddMod(in_out, accum_space, in_out, 2 * N, modulus);
    intel::hexl::EltwiseAddMod(in_out, accum_space + 2 * N, in_out, 2 * N, modulus);

}
