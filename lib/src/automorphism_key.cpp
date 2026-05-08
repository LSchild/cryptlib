//
// Created by leonard on 11/12/24.
//

#include <cassert>
#include "automorphism_key.h"

AutomorphismParameters::AutomorphismParameters(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N,
                                               uint64_t basis, uint64_t digits, double std,
                                               uint32_t automorphism_index) : RLWEConversionParameters(source_key_distribution,modulus,N,basis,digits,std), m_automorphism_index(automorphism_index) {}

AutomorphismParameters::AutomorphismParameters(KeyDistribution source_key_distribution,
                                               std::shared_ptr<intel::hexl::NTT> ntt, uint64_t basis, uint64_t digits,
                                               double std, uint32_t automorphism_index) : RLWEConversionParameters(source_key_distribution, ntt, basis, digits, std), m_automorphism_index(automorphism_index) {}

AutomorphismParameters::AutomorphismParameters(AutomorphismParameters &other) : RLWEConversionParameters(other) {
    m_automorphism_index = other.GetAutomorphismIndex();
}


void AutomorphismParameters::SetAutomorphismIndex(uint32_t idx) {
    m_automorphism_index = idx;
}

uint32_t AutomorphismParameters::GetAutomorphismIndex() const {
    return m_automorphism_index;
}


AutomorphismEvaluator::AutomorphismEvaluator(AutomorphismParameters &params) : RLWEtoRLWEConverter(params), m_automorphism_params(params) {}

void AutomorphismEvaluator::SetParams(AutomorphismParameters params) {
    RLWEtoRLWEConverter::SetParams(params);
    m_automorphism_params = params;
}

void AutomorphismEvaluator::KeyGen(const uint64_t *const key) {
    m_auto_buffer.resize(2 * m_params.GetDimension());
    std::vector<uint64_t> auto_key(m_params.GetDimension(), 0);
    ApplyAutomorphism(auto_key.data(), key);
    RLWEtoRLWEConverter::KeyGen(key, auto_key.data());
}

void AutomorphismEvaluator::KeyGen(const std::vector<uint64_t> &key) {
    m_auto_buffer.resize(2 * m_params.GetDimension());
    std::vector<uint64_t> auto_key(m_params.GetDimension(), 0);
    ApplyAutomorphism(auto_key, key);
    RLWEtoRLWEConverter::KeyGen(auto_key, key);
}

void AutomorphismEvaluator::Eval(uint64_t *output, const uint64_t *const input) {
    auto N = m_params.GetDimension();
    ApplyAutomorphism(m_auto_buffer.data(), input);
    ApplyAutomorphism(m_auto_buffer.data() + N , input + N);
    RLWEtoRLWEConverter::Convert(output, m_auto_buffer.data());
}

void AutomorphismEvaluator::Eval(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    auto N = m_params.GetDimension();
    ApplyAutomorphism(m_auto_buffer.data(), input.data());
    ApplyAutomorphism(m_auto_buffer.data() + N , input.data() + N);
    RLWEtoRLWEConverter::Convert(output.data(), m_auto_buffer.data());
}

const AutomorphismParameters &AutomorphismEvaluator::GetParams() const {
    return m_automorphism_params;
}

void AutomorphismEvaluator::ApplyAutomorphism(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    ApplyAutomorphism(output.data(), input.data());
}

void AutomorphismEvaluator::ApplyAutomorphism(uint64_t *output, const uint64_t *const input) {
    // TODO benchmark
    auto N = m_params.GetDimension();
    auto dim2_mask = 2 * N - 1;
    auto dim_mask = N - 1;
    auto Q = m_params.GetModulus();
    const uint32_t stride = 4;

    for(uint64_t i = 0; i < N; i+=stride) {
        for(uint32_t j = 0; j < stride; j++) {
            auto idx = i + j;
            auto flip_sign = (idx & dim2_mask) >= N;
            auto new_idx = idx & dim_mask;
            output[new_idx] = input[idx];
            if (flip_sign) {
                output[new_idx] = intel::hexl::SubUIntMod(0, output[new_idx], Q);
            }
        }
    }

}

// OLD VERSION
NativePoly AutomorphismKey::ApplyAutomorphism(const lbcrypto::NativePoly &A, uint32_t automorphism_idx) {
    NativePoly A2 = A.Clone();
    A2.SetFormat(COEFFICIENT);

    const auto& params = A2.GetParams();
    auto N = A.GetLength();
    auto Q = params->GetModulus();
    NativePoly res(params, COEFFICIENT, true);

    res.at(0) = A2.at(0);
    for(int i = 1; i < N; i++) {
        auto new_idx = (i * automorphism_idx) % (2 * N);
        if (new_idx >= N){
            new_idx -= N;
            res.at(new_idx) = (Q - A2.at(i)).Mod(Q);
        } else {
            res.at(new_idx) = A2.at(i) ;
        }
    }

    return res;
}

AutomorphismKey::AutomorphismKey(std::shared_ptr<RingGSWCryptoParams> &params, lbcrypto::NativePoly sk,
                                 uint32_t automorphism_idx) :
        RLWEKeySwitchingKey(params,  ApplyAutomorphism(sk, automorphism_idx), sk),
        m_automorphism_idx(automorphism_idx) {
    assert(std::gcd(automorphism_idx, 2 * params->GetN()) == 1 && "Automorphism index given does not induce actual automorphism over ring.");
}

RLWECiphertext AutomorphismKey::AutomorphismTransform(const lbcrypto::RLWECiphertext &source) const {
    auto A_T = ApplyAutomorphism(source->GetElements()[0], m_automorphism_idx);
    auto B_T = ApplyAutomorphism(source->GetElements()[1], m_automorphism_idx);
    return SwitchKey(A_T, B_T);
}

RLWECiphertext AutomorphismKey::AutomorphismTransform(const lbcrypto::NativePoly &A, const lbcrypto::NativePoly &B) {
    auto A_T = ApplyAutomorphism(A, m_automorphism_idx);
    auto B_T = ApplyAutomorphism(B, m_automorphism_idx);
    return SwitchKey(A_T, B_T);
}

FastAutomorphismKey::FastAutomorphismKey(std::shared_ptr<intel::hexl::NTT> &ntt_engine, uint32_t L_bits, NativePoly sk, uint32_t automorphism_idx) : ntt_engine(ntt_engine), m_automorphism_idx(automorphism_idx) {
    sk.SetFormat(COEFFICIENT);
    auto sk_auto = AutomorphismKey::ApplyAutomorphism(sk, automorphism_idx);
    sk_auto.SwitchFormat();
    std::vector<RLWECiphertext> m0;
    std::vector<RLWECiphertext> m1;

    auto pp = sk.GetParams();
    auto zero = NativePoly(pp, EVALUATION, true);
    NativeInteger basis = 1ull << L_bits;

    auto log2Q = intel::hexl::Log2(ntt_engine->GetModulus());
    uint32_t offset = 0;

    uint32_t digits = 0;
    while (offset < log2Q) {
        std::vector<NativePoly> rlwe0 = {zero.Clone(), zero.Clone()};
        m0.emplace_back(std::make_shared<RLWECiphertextImpl>(rlwe0));
        // TODO use "safe" encryption
        std::vector<NativePoly> rlwe1 = {zero.Clone(), sk_auto.Clone()};
        m1.emplace_back(std::make_shared<RLWECiphertextImpl>(rlwe1));
        sk_auto *= basis;
        offset += L_bits;
        digits++;
    }
    m_key = RGSWSample(ntt_engine, 1ull << L_bits,digits , m0,m1);
}

RLWECiphertext FastAutomorphismKey::AutomorphismTransform(const NativePoly &A_in, const NativePoly &B_in) {

    auto AA = AutomorphismKey::ApplyAutomorphism(A_in, m_automorphism_idx);
    auto BB = AutomorphismKey::ApplyAutomorphism(B_in, m_automorphism_idx);

    auto N = A_in.GetLength();
    std::vector<uint64_t> buffer(3 *  N);
    for (uint32_t i = 0; i < N; i++) {
        buffer[i + 2 * N] = AA[i].ConvertToInt<uint64_t>();
    }
    m_key.mul_dir(buffer.data(), buffer.data() + 2 * N);
    auto vec = NativeVector(N, A_in.GetModulus());

    auto A = NativePoly(A_in.GetParams(), EVALUATION, false);
    auto B = NativePoly(A_in.GetParams(), EVALUATION, false);

    for(uint64_t i = 0; i < N; i++) {
        vec[i] = buffer[i];
    }

    A.SetValues(vec, EVALUATION);
    A.SwitchFormat();

    for(uint64_t i = 0; i < N; i++) {
        vec[i] = buffer[i + N];
    }

    B.SetValues(vec, EVALUATION);
    B.SwitchFormat();

    auto ct_vec = {A.Negate(), BB - B};
    return std::make_shared<RLWECiphertextImpl>(ct_vec);
}


