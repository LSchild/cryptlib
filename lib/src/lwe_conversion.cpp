//
// Created by leonard on 4/22/26.
//

#include <cassert>
#include "operators/endo_glwe_conversion.h"
#include "utils/math_utils.h"
#include "base_crypto.h"
#include "utils/speed_utils.h"

LWEConversionParameters::LWEConversionParameters(KeyDistribution source_distribution, uint64_t modulus, uint64_t source_dimension, uint64_t target_dimension, uint64_t basis,
                                                 uint64_t digits, double std) : m_source_distribution(source_distribution), m_modulus(modulus),
                                             m_source_dimension(source_dimension), m_target_dimension(target_dimension), m_basis(basis),
                                             m_basis_log2(IntLog2(basis)), m_digits(digits), m_std(std)
                                             {

                                             }

LWEConversionParameters::LWEConversionParameters(LWEConversionParameters &other) {
    SetModulus(other.GetModulus());
    SetSourceDimension(other.GetSourceDimension());
    SetTargetDimension(other.GetTargetDimension());
    SetStd(other.GetStd());
    SetGadgetBasis(other.GetGadgetBasis());
    SetGadgetDigits(other.GetGadgetDigits());
}

Container LWEConversionParameters::GetInputContainer() const {
    return std::make_shared<LWEContainerImpl>(m_source_dimension, m_modulus, 0.0);
}

Container LWEConversionParameters::GetOutputContainer(Container input) const {
    auto input_lwe = std::dynamic_pointer_cast<RLWEContainerImpl>(input);
    auto o_var = ComputeOutputVariance(input_lwe->GetVariance());
    return std::make_shared<LWEContainerImpl>(m_target_dimension, m_modulus, o_var);
}

OperatorID LWEConversionParameters::GetOperatorID() const {
    return CONV_LWE_LWE;
}

std::unique_ptr<LWEtoLWEConverter> LWEConversionParameters::ConstructOperator(const std::vector<GenericKey> &keys) const {
    auto op = std::unique_ptr<LWEtoLWEConverter>(new LWEtoLWEConverter(this->shared_from_this()));
    op->KeyGen(keys[0].GetKeyPtr(),keys[1].GetKeyPtr());
    
    return std::move(op);
}

uint64_t LWEConversionParameters::GetModulus() const {
    return m_modulus;
}

uint64_t LWEConversionParameters::GetSourceDimension() const {
    return m_source_dimension;
}

uint64_t LWEConversionParameters::GetTargetDimension() const {
    return m_target_dimension;
}

uint64_t LWEConversionParameters::GetGadgetBasis() const {
    return m_basis;
}

uint64_t LWEConversionParameters::GetGadgetBasisLog2() const {
    return m_basis_log2;
}

uint64_t LWEConversionParameters::GetGadgetDigits() const {
    return m_digits;
}

double LWEConversionParameters::GetStd() const {
    return m_std;
}

KeyDistribution LWEConversionParameters::GetSourceKeyDistribution() const {
    return m_source_distribution;
}

void LWEConversionParameters::SetModulus(uint64_t modulus) {
    m_modulus = modulus;
}

void LWEConversionParameters::SetSourceDimension(uint64_t input_dimension) {
    m_source_dimension = input_dimension;
}

void LWEConversionParameters::SetTargetDimension(uint64_t output_dimension) {
    m_target_dimension = output_dimension;
}

void LWEConversionParameters::SetGadgetBasis(uint64_t basis) {
    m_basis = basis;
    m_basis_log2 = IntLog2(basis);

}

void LWEConversionParameters::SetStd(double std) {
    m_std = std;
}

void LWEConversionParameters::SetGadgetDigits(uint64_t digits) {
    m_digits = digits;
}

void LWEConversionParameters::SetSourceKeyDistribution(KeyDistribution distribution) {
    m_source_distribution = distribution;
}

long double LWEConversionParameters::ComputeOutputVariance(long double input_variance) const {

    long double delta = (m_modulus >> (m_basis_log2 * m_digits));
    delta /= 12.0;
    long double var_nrm2_key = m_source_distribution == BINARY ? double(m_source_dimension >> 1) / 2 : double((3 * m_source_dimension) >> 1) * 2.0 / 3.0;

    long double enc_var = m_std * m_std;
    long double lwe_prime_mul = m_digits * m_basis * m_basis * enc_var;
    lwe_prime_mul /= 12.0;
    long double var = m_source_dimension * (lwe_prime_mul + delta * var_nrm2_key);

    return var;
}

LWEtoLWEConverter::LWEtoLWEConverter(std::shared_ptr<const LWEConversionParameters> params) : m_params(params), m_params_set(true) {
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



const std::shared_ptr<const LWEConversionParameters> LWEtoLWEConverter::GetContext() const {
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

    for(uint64_t i = 0; i < n_in; i++) {
        auto chunk = ksk + i * digits * (n_out + 1);
        auto shift = modulus_bits - basis_bits;
        auto a_i = input[i];

        for(uint64_t j = 0; j < digits; j++, chunk += n_out + 1) {
            auto d_ij = (a_i >> shift) & mask;
            intel::hexl::EltwiseFMAMod(acc, chunk, d_ij, acc, n_out + 1, modulus, 1);
            shift -= basis_bits;
        }
    }
    intel::hexl::EltwiseSubMod(output, output, acc, n_out + 1, modulus);
}
