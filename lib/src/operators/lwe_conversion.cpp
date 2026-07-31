//
// Created by leonard on 4/22/26.
//

#include <cassert>
#include <iostream>
#include "static/gadget_decomp.h"
#include "operators/endo_glwe_conversion.h"
#include "utils/generic_utils.h"
#include "utils/math_utils.h"
#include "base_crypto.h"


LWEConversionContext::LWEConversionContext(KeyDistribution source_distribution, uint64_t modulus, uint64_t source_dimension, uint64_t target_dimension, uint64_t basis,
                                           uint64_t digits, double std) : m_source_distribution(source_distribution), m_modulus(modulus),
                                             m_source_dimension(source_dimension), m_target_dimension(target_dimension), m_basis(basis),
                                             m_basis_log2(IntLog2(basis)), m_digits(digits), m_std(std)
                                             {

                                             }

LWEConversionContext::LWEConversionContext(LWEConversionContext &other) {
    SetModulus(other.GetModulus());
    SetSourceDimension(other.GetSourceDimension());
    SetTargetDimension(other.GetTargetDimension());
    SetStd(other.GetStd());
    SetGadgetBasis(other.GetGadgetBasis());
    SetGadgetDigits(other.GetGadgetDigits());
}

Container LWEConversionContext::GetInputContainer() const {
    return std::make_shared<LWEContainerImpl>(m_source_dimension, m_modulus, 0.0);
}

Container LWEConversionContext::GetOutputContainer(Container input) const {
    auto input_lwe = std::dynamic_pointer_cast<RLWEContainerImpl>(input);
    auto o_var = ComputeOutputVariance(input_lwe->GetVariance());
    return std::make_shared<LWEContainerImpl>(m_target_dimension, m_modulus, o_var);
}

OperatorID LWEConversionContext::GetOperatorID() const {
    return CONV_LWE_LWE;
}

std::unique_ptr<LWEtoLWEConverter> LWEConversionContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    auto op = std::unique_ptr<LWEtoLWEConverter>(new LWEtoLWEConverter(this->shared_from_this()));
    op->KeyGen(keys[0].GetKeyPtr(),keys[1].GetKeyPtr());
    
    return std::move(op);
}

uint64_t LWEConversionContext::GetModulus() const {
    return m_modulus;
}

uint64_t LWEConversionContext::GetSourceDimension() const {
    return m_source_dimension;
}

uint64_t LWEConversionContext::GetTargetDimension() const {
    return m_target_dimension;
}

uint64_t LWEConversionContext::GetGadgetBasis() const {
    return m_basis;
}

uint64_t LWEConversionContext::GetGadgetBasisLog2() const {
    return m_basis_log2;
}

uint64_t LWEConversionContext::GetGadgetDigits() const {
    return m_digits;
}

double LWEConversionContext::GetStd() const {
    return m_std;
}

KeyDistribution LWEConversionContext::GetSourceKeyDistribution() const {
    return m_source_distribution;
}

void LWEConversionContext::SetModulus(uint64_t modulus) {
    m_modulus = modulus;
}

void LWEConversionContext::SetSourceDimension(uint64_t input_dimension) {
    m_source_dimension = input_dimension;
}

void LWEConversionContext::SetTargetDimension(uint64_t output_dimension) {
    m_target_dimension = output_dimension;
}

void LWEConversionContext::SetGadgetBasis(uint64_t basis) {
    m_basis = basis;
    m_basis_log2 = IntLog2(basis);

}

void LWEConversionContext::SetStd(double std) {
    m_std = std;
}

void LWEConversionContext::SetGadgetDigits(uint64_t digits) {
    m_digits = digits;
}

void LWEConversionContext::SetSourceKeyDistribution(KeyDistribution distribution) {
    m_source_distribution = distribution;
}

long double LWEConversionContext::ComputeOutputVariance(long double input_variance) const {

    auto modulus_bits = IntLog2(m_modulus);
    long double delta = std::powl(2.0, modulus_bits) / std::powl(m_basis, m_digits);
    if (delta <= 1) {
        delta = 0.0;
    }

    delta *= delta / 12.0;

    long double message_norm = 1;
    if (m_source_distribution == GAUSSIAN) {
        std::cerr << "Gaussian keys are not supported yet. Setting key element norm to 5" << std::endl;
        message_norm = 5 * 5; // norm square
    }

    auto expected_hamming_weight = m_source_dimension;
    if (m_source_distribution == BINARY) {
        expected_hamming_weight = m_source_dimension / 2;
    }
    if (m_source_distribution == TERNARY) {
        expected_hamming_weight = 2 * m_source_dimension / 3;
    }

    long double encryption_var = m_std * m_std;
    long double lwe_prime_mul = m_digits * m_basis * m_basis * encryption_var;
    lwe_prime_mul /= 12.0;


    long double var = m_source_dimension * lwe_prime_mul + expected_hamming_weight * delta * message_norm;

    return input_variance + var;
}

LWEtoLWEConverter::LWEtoLWEConverter(std::shared_ptr<const LWEConversionContext> params) : m_params(params), m_params_set(true) {
    auto dim_in = params->GetSourceDimension();
    auto dim_out = params->GetTargetDimension();
    auto digits = params->GetGadgetDigits();
    m_ksk.resize(dim_in * (dim_out + 1) * digits);
    m_acc.resize(dim_out + 1);
    m_params_set = true;
}


void LWEtoLWEConverter::KeyGen(const std::vector<uint64_t> &source_key, const std::vector<uint64_t> &target_key) {
    KeyGen(source_key.data(), target_key.data());
}

void LWEtoLWEConverter::Convert(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    Convert(output.data(), input.data());
}



const std::shared_ptr<const LWEConversionContext> LWEtoLWEConverter::GetContext() const {
    return m_params;
}



void LWEtoLWEConverter::KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) {

    assert(m_params_set);

    const auto encryptor = LWEEncryptor(m_params->GetModulus(), m_params->GetTargetDimension(), m_params->GetStd());
    const auto modulus = m_params->GetModulus();
    const auto modulus_bits = IntLog2(m_params->GetModulus());
    const auto digits = m_params->GetGadgetDigits();
    const auto basis_bits = m_params->GetGadgetBasisLog2();
    const auto ksk = m_ksk.data();

    auto precon = intel::hexl::MultiplyFactor(1, 64, modulus).BarrettFactor();

    auto n_out = m_params->GetTargetDimension();
    auto n_in = m_params->GetSourceDimension();

    for(uint32_t i = 0; i < n_in; i++) {
        auto chunk = ksk + i * digits * (n_out + 1);
        // Note digits from top to bottom
        auto g_ij = 1ull << (modulus_bits - basis_bits);
        for(uint64_t j = 0; j < digits; j++, chunk += (n_out + 1)) {
            auto entry = intel::hexl::MultiplyMod(g_ij, source_key[i], precon, modulus);
            encryptor.MakeLWE(chunk, entry, target_key);
            g_ij >>= basis_bits;
        }
    }

    m_keys_generated = true;

}

void LWEtoLWEConverter::Convert(uint64_t *const output, const uint64_t *const input) {
    assert(m_keys_generated);
    const auto n_in = m_params->GetSourceDimension();
    const auto n_out = m_params->GetTargetDimension();
    const auto modulus = m_params->GetModulus();
    const auto modulus_bits = IntLog2(modulus);
    const auto basis_bits = m_params->GetGadgetBasisLog2();
    const auto mask = (1ull << basis_bits) - 1;
    const auto digits = m_params->GetGadgetDigits();
    const auto ksk = m_ksk.data();


    ZERO_UINT64_ARR(m_acc.data(), n_out + 1);
    const auto acc = m_acc.data();
    output[n_out] = input[n_in];

    auto [corrector, least_digit_correction] = GetCorrectorForSignedToUnsigned(modulus, modulus_bits, basis_bits, digits);

    uint64_t max_digits = (modulus_bits / basis_bits) + (modulus_bits % basis_bits == 0 ? 0 : 1);

    for(uint64_t i = 0; i < n_in; i++) {
        auto chunk = ksk + i * digits * (n_out + 1);
        auto shift = modulus_bits - basis_bits;
        auto a_i =  intel::hexl::AddUIntMod(input[i], corrector, modulus);

        for(uint64_t j = 0; j < digits; j++, chunk += n_out + 1) {
            auto d_ij = (a_i >> shift) & mask;
            auto delta = j != (max_digits - 1) ? 1ull << (basis_bits - 1) : least_digit_correction;
            d_ij = intel::hexl::SubUIntMod(d_ij, delta, modulus);
            intel::hexl::EltwiseFMAMod(acc, chunk, d_ij, acc, n_out + 1, modulus, 1);
            shift -= basis_bits;
        }
    }
    intel::hexl::EltwiseSubMod(output, output, acc, n_out + 1, modulus);
}
