//
// Created by leonard on 3/24/26.
//
#include "base_crypto.h"

#include "common_types.h"
#include "openfhe.h"

RLWEEncryptor::RLWEEncryptor(uint64_t modulus, uint32_t ring_dimension, double std) {
    m_std = std;
    m_ntt = std::make_shared<intel::hexl::NTT>(ring_dimension, modulus);
}

RLWEEncryptor::RLWEEncryptor(uint64_t modulus, uint32_t ring_dimension, uint64_t root_of_unity, double std) {
    m_std = std;
    m_ntt = std::make_shared<intel::hexl::NTT>(ring_dimension, modulus, root_of_unity);
}

RLWEEncryptor::RLWEEncryptor(std::shared_ptr<intel::hexl::NTT> ntt, double std) {
    m_ntt = ntt;
    m_std = std;
}

void RLWEEncryptor::MakeRLWE(uint64_t *result, uint64_t *msg, uint64_t *secret_ntt, bool msg_is_ntt) {
    // TODO: for now keep using the OPENFHE samplers so we don't mess up
    auto N = m_ntt->GetDegree();
    auto Q = m_ntt->GetModulus();

    // TODO: double check
    AlignedVector tmp = AlignedVector(N);
    lbcrypto::DiscreteUniformGeneratorImpl<NativeVector> sampler_unif(Q);
    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> sampler_gauss(m_std);

    // uniformly random = is already in ntt format ;)
    auto A = sampler_unif.GenerateVector(N);
    auto E = sampler_gauss.GenerateVector(N, Q);

    for (uint32_t i = 0; i < m_ntt->GetDegree(); i++) {
        result[i] = A[i].ConvertToInt<uint64_t>();
        tmp[i] = E[i].ConvertToInt<uint64_t>();
    }


    intel::hexl::EltwiseMultMod(result + N, result, secret_ntt, N, Q, 1);

    if (msg_is_ntt) {
        m_ntt->ComputeForward(tmp.data(),tmp.data(),1,1);
        intel::hexl::EltwiseAddMod(tmp.data(), tmp.data(), msg, N, Q);
    } else {
        intel::hexl::EltwiseAddMod(tmp.data(), tmp.data(), msg, N, Q);
        m_ntt->ComputeForward(tmp.data(),tmp.data(),1,1);
    }

    intel::hexl::EltwiseAddMod(result + N, result + N, tmp.data(), N, Q);
}

void RLWEEncryptor::MakeRGSW(uint64_t *result, uint64_t *msg, uint64_t *secret_ntt, uint64_t rgsw_basis, uint64_t rgsw_digits, bool msg_is_ntt) {

    auto N = m_ntt->GetDegree();
    auto Q = m_ntt->GetModulus();
    AlignedVector tmp = AlignedVector(3 * N);
    AlignedVector basis = AlignedVector(2 * N);
    std::fill(tmp.begin(), tmp.end(), 0);

    uint64_t max_digits = std::ceil(std::log2((long double)Q)/std::log2(rgsw_basis));

    // make tmp = [msg, -s * msg]
    if (msg_is_ntt) {
        intel::hexl::EltwiseMultMod(tmp.data(), msg, secret_ntt, N, Q, 1);
        intel::hexl::EltwiseSubMod(tmp.data() + N, tmp.data() + N, tmp.data(), N, Q);
        std::copy(msg, msg + N, tmp.data());
    } else {
        // tmp = [ NTT(msg) | 0 | 0]
        m_ntt->ComputeForward(tmp.data(), msg, 1, 1);
        // tmp = [ NTT(msg) | 0 | NTT(secret * msg)]
        intel::hexl::EltwiseMultMod(tmp.data() + 2 * N, tmp.data(), secret_ntt, N, Q, 1);
        // tmp = [ NTT(msg) | -NTT(secret * msg) | NTT(secret * msg)]
        intel::hexl::EltwiseSubMod(tmp.data() + N, tmp.data() + N, tmp.data() + 2 * N, N, Q);
    }


    for (uint64_t d_i = 0; d_i < max_digits - rgsw_digits; d_i++) {
        intel::hexl::EltwiseFMAMod(tmp.data(), tmp.data(), rgsw_basis, nullptr, 2 * N, Q, 1);
    }

    for (uint64_t d_i = 0; d_i < rgsw_digits; d_i++) {
        uint64_t true_offset = (rgsw_digits - d_i - 1);
        MakeRLWE(result + true_offset * 2 * N, tmp.data() + N, secret_ntt, true);
        MakeRLWE(result + true_offset * 2 * N + 2 * N * rgsw_digits, tmp.data(), secret_ntt, true);
        intel::hexl::EltwiseFMAMod(tmp.data(), tmp.data(), rgsw_basis, nullptr, 2 * N, Q, 1);
    }

}

void RLWEEncryptor::PhaseRLWE(uint64_t *result, uint64_t *rlwe, uint64_t *secret_ntt) {
    auto N = m_ntt->GetDegree();
    auto Q = m_ntt->GetModulus();
    AlignedVector tmp = AlignedVector(N);

    intel::hexl::EltwiseMultMod(tmp.data(), rlwe, secret_ntt, N, Q, 1);
    std::copy(rlwe + N, rlwe + 2 * N, result);
    intel::hexl::EltwiseSubMod(result, result, tmp.data(), N, Q);
    m_ntt->ComputeInverse(result, result, 1, 1);
}

void RLWEEncryptor::PhaseRGSW(uint64_t *result, uint64_t *rgsw, uint64_t *secret_ntt, uint64_t rgsw_digits) {
    auto N = m_ntt->GetDegree();

    for (uint64_t i = 0; i < 2 * rgsw_digits; i++) {
        PhaseRLWE(result + i * 2 * N, rgsw + i * 2 * N, secret_ntt);
    }
}

uint64_t RLWEEncryptor::GetDimension() {
    return m_ntt->GetDegree();
}

uint64_t RLWEEncryptor::GetModulus() {
    return m_ntt->GetModulus();
}

double RLWEEncryptor::GetStd() {
    return m_std;
}

std::shared_ptr<intel::hexl::NTT> RLWEEncryptor::GetNTT() {
    return m_ntt;
}


LWEEncryptor::LWEEncryptor(uint64_t modulus, uint64_t lwe_dimension, double std) : m_modulus(modulus), m_dimension(lwe_dimension), m_std(std) {}

void LWEEncryptor::MakeLWE(uint64_t *result, uint64_t msg, const uint64_t *const secret) const {

    lbcrypto::DiscreteUniformGeneratorImpl<NativeVector> sampler_unif(m_modulus);
    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> sampler_gauss(m_std);
    auto vec = sampler_unif.GenerateVector(m_dimension);

    auto b = sampler_gauss.GenerateInteger(m_modulus).ConvertToInt<uint64_t>();

    auto precon = intel::hexl::MultiplyFactor(1, 64, m_modulus).BarrettFactor();

    for (uint64_t i = 0; i < m_dimension; i++) {
        auto a_i = vec[i].ConvertToInt<uint64_t>();
        auto prod = intel::hexl::MultiplyMod(a_i, secret[i], precon, m_modulus);
        result[i] = a_i;
        b = intel::hexl::AddUIntMod(b, prod, m_modulus);
    }
    result[m_dimension] = intel::hexl::AddUIntMod(b, msg, m_modulus);
}

void LWEEncryptor::PhaseLWE(uint64_t *result, uint64_t *lwe_vec, const uint64_t *const secret) const {
    uint64_t  acc = 0;
    auto precon = intel::hexl::MultiplyFactor(1, 64, m_modulus).BarrettFactor();

    for (uint64_t i = 0; i < m_dimension; i++) {
        auto prod = intel::hexl::MultiplyMod(lwe_vec[i], secret[i], precon, m_modulus);
        acc = intel::hexl::AddUIntMod(acc, prod, m_modulus);
    }

    *result = intel::hexl::SubUIntMod(lwe_vec[m_dimension], acc, m_modulus);
}

uint64_t LWEEncryptor::PhaseLWE(uint64_t *lwe_vec,const uint64_t *const secret) const {
    uint64_t res;
    PhaseLWE(&res, lwe_vec, secret);
    return res;
}

double LWEEncryptor::GetStd() {
    return m_std;
}

uint64_t LWEEncryptor::GetModulus() {
    return m_modulus;
}

uint64_t LWEEncryptor::GetDimension() {
    return m_dimension;
}
