//
// Created by leonard on 4/26/26.
//

#include <cassert>
#include "base_crypto.h"
#include "glwe_conversion.h"
#include "speed_utils.h"
#include "math_utils.h"
#include "gadget_decomp.h"

RLWEConversionParameters::RLWEConversionParameters(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N,
                                                   uint64_t basis, uint64_t digits, double std) {
    SetSourceKeyDistribution(source_key_distribution);
    SetModulus(modulus);
    SetDimension(N);
    SetGadgetBasis(basis);
    SetGadgetDigits(digits);
    SetStd(std);

    auto ntt = std::make_shared<intel::hexl::NTT>(N, modulus);

    SetNTT(ntt);
}

RLWEConversionParameters::RLWEConversionParameters(KeyDistribution source_key_distribution,
                                                   std::shared_ptr<intel::hexl::NTT> ntt, uint64_t basis, uint64_t digits,
                                                   double std) {
    SetSourceKeyDistribution(source_key_distribution);
    SetNTT(ntt);
    SetModulus(ntt->GetModulus());
    SetDimension(ntt->GetDegree());
    SetGadgetBasis(basis);
    SetGadgetDigits(digits);
    SetStd(std);
}

RLWEConversionParameters::RLWEConversionParameters(RLWEConversionParameters &other) {
    SetSourceKeyDistribution(other.GetSourceKeyDistribution());
    SetNTT(other.GetNTT());
    SetModulus(other.GetModulus());
    SetDimension(other.GetDimension());
    SetGadgetBasis(other.GetGadgetBasis());
    SetGadgetDigits(other.GetGadgetDigits());
    SetStd(other.GetStd());
}


long double RLWEConversionParameters::ComputeOutputVariance(long double input_variance = 0) const {

    long double delta = (m_modulus >> (m_basis_log2 * m_digits));
    delta /= 12.0;

    long double var_nrm2_key = m_source_distribution == BINARY ? double(m_N >> 1) / 2 : double((3 * m_N) >> 1) * 2.0 / 3.0;
    long double fresh_var = m_std * m_std;
    long double rlwe_p_prod_var = m_digits * (m_basis * m_basis + 2) * fresh_var * m_N;
    rlwe_p_prod_var /= 12.0;

    long double var = rlwe_p_prod_var + delta * (1 + var_nrm2_key);
    return input_variance + var;
}

KeyDistribution RLWEConversionParameters::GetSourceKeyDistribution() const {
    return m_source_distribution;
}

uint64_t RLWEConversionParameters::GetModulus() const {
    return m_modulus;
}
uint64_t RLWEConversionParameters::GetDimension() const {
    return m_N;
}

uint64_t RLWEConversionParameters::GetGadgetBasis() const {
    return m_basis;
}

uint64_t RLWEConversionParameters::GetGadgetBasisLog2() const {
    return m_basis_log2;
}

uint64_t RLWEConversionParameters::GetGadgetDigits() const {
    return m_digits;
}

double RLWEConversionParameters::GetStd() const {
    return m_std;
}

std::shared_ptr<intel::hexl::NTT> RLWEConversionParameters::GetNTT() const {
    return m_ntt;
}

void RLWEConversionParameters::SetSourceKeyDistribution(KeyDistribution distribution) {
    m_source_distribution = distribution;
}

void RLWEConversionParameters::SetStd(double std) {
    m_std = std;
}

void RLWEConversionParameters::SetDimension(uint64_t input_dimension) {
    m_N = input_dimension;
}

void RLWEConversionParameters::SetNTT(std::shared_ptr<intel::hexl::NTT> ntt) {
    m_ntt = std::move(ntt);
}

void RLWEConversionParameters::SetModulus(uint64_t modulus) {
    m_modulus = modulus;
}

void RLWEConversionParameters::SetGadgetDigits(uint64_t digits) {
    m_digits = digits;
}

void RLWEConversionParameters::SetGadgetBasis(uint64_t basis) {
    m_basis = basis;
    m_basis_log2 = IntLog2(m_basis);
}


RLWEtoRLWEConverter::RLWEtoRLWEConverter(RLWEConversionParameters &params) : m_params(params) {
    auto N = params.GetDimension();
    auto digits = params.GetGadgetDigits();

    m_acc.resize(2 * N * digits);
    m_ksk.resize(2 * N * digits);
    m_params_set = true;
}

void RLWEtoRLWEConverter::SetParams(RLWEConversionParameters &params) {
    m_params = params;
    auto N = m_params.GetDimension();
    auto digits = m_params.GetGadgetDigits();

    m_acc.resize(2 * N * digits);
    m_ksk.resize(2 * N * digits);
    m_params_set = true;

}

void RLWEtoRLWEConverter::KeyGen(const std::vector<uint64_t> &source_key, const std::vector<uint64_t> &target_key) {
    assert(m_params_set);
    m_keys_generated = true;
    KeyGen(source_key.data(), target_key.data());
}

void RLWEtoRLWEConverter::KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) {
    assert(m_params_set);
    m_keys_generated = true;
    // assumption: both keys in non-ntt form

    auto N = m_params.GetDimension();
    auto modulus = m_params.GetModulus();
    auto modulus_bits = IntLog2(m_params.GetModulus());
    auto basis_bits = m_params.GetGadgetBasisLog2();
    auto digits = m_params.GetGadgetDigits();
    auto ntt = m_params.GetNTT();
    RLWEEncryptor encryptor(ntt, m_params.GetStd());

    AlignedVector target_key_ntt(N, 0);
    AlignedVector gadget_entry(N, 0);
    ntt->ComputeForward(target_key_ntt.data(), target_key, 1, 1);

    auto ksk = m_ksk.data();
    auto g_ij = 1ull << (modulus_bits - basis_bits);

    for(uint64_t i = 0; i < digits; i++) {
        auto chunk = ksk + i * 2 * N;
        intel::hexl::EltwiseFMAMod(gadget_entry.data(), source_key, g_ij, nullptr, N, modulus, 1);
        encryptor.MakeRLWE(chunk, gadget_entry.data(), target_key_ntt.data());
        g_ij >>= basis_bits;
    }

}

void RLWEtoRLWEConverter::Convert(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    Convert(output.data(), input.data());
}

void RLWEtoRLWEConverter::Convert(uint64_t *output, const uint64_t *const input) {
    auto modulus = m_params.GetModulus();
    auto N = m_params.GetDimension();
    auto ntt = m_params.GetNTT();
    auto digits = m_params.GetGadgetDigits();
    auto basis_bits = m_params.GetGadgetBasisLog2();
    auto modulus_bits = IntLog2(modulus);
    auto acc = m_acc.data();
    ZERO_UINT64_ARR(acc, m_acc.size());

    SignedDigitDecomposeRep2NTT(acc, input, ntt, digits, basis_bits, modulus_bits);

    intel::hexl::EltwiseMultMod(acc, acc, m_ksk.data(), 2 * N * digits, modulus, 1);

    for(uint64_t i = 1; i < digits; i++) {
        intel::hexl::EltwiseAddMod(acc, acc, acc + i * 2 * N, 2 * N, modulus);
    }

    intel::hexl::EltwiseSubMod(output, output, acc, 2 * N, modulus);

    // DANGER: Correct provided b @input + N is in NTT format
    intel::hexl::EltwiseAddMod(output + N, output + N, input + N, N, modulus);
}

Container RLWEtoRLWEConverter::GetInputContainer() const {
    return std::make_shared<RLWEContainerImpl>(m_params.GetModulus(), m_params.GetDimension(), 0.0);
}

Container RLWEtoRLWEConverter::GetOutputContainer() const {
    return std::make_shared<RLWEContainerImpl>(m_params.GetModulus(), m_params.GetDimension(), m_params.ComputeOutputVariance());
}

const RLWEConversionParameters &RLWEtoRLWEConverter::GetParams() const {
    return m_params;
}