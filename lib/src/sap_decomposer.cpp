//
// Created by leonard on 6/3/26.
//

#include "operators/sap_decomposition.h"

std::shared_ptr<LWEConversionParameters> SAPDecompositionContext::GetLWEConversionContext() const {
    return m_conversion_context;
}

std::shared_ptr<OperatorContext<BlindRotator>> SAPDecompositionContext::GetBlindRotationContext() const {
    return m_rotation_context;
}

std::shared_ptr<LWEtoRLWEPackingContext> SAPDecompositionContext::GetPackingContext() const {
    return m_packing_context;
}

Container SAPDecompositionContext::GetInputContainer() const {
    auto dummy = m_conversion_context->GetInputContainer();
    auto out = m_conversion_context->GetOutputContainer(dummy);
    auto out_lwe = std::dynamic_pointer_cast<LWEContainerImpl>(out);

    auto in_n = out_lwe->GetN();
    auto in_q = UINT64_MAX;

    return std::make_shared<LWEContainerImpl>(in_n, in_q, 0.0);
}

long double SAPDecompositionContext::ComputeOutputVariance(long double input_variance) const {

    auto N = std::dynamic_pointer_cast<LWEContainerImpl>(m_conversion_context->GetInputContainer())->GetN();
    auto iter_ratio = N / m_default_radix;

    long double out_var = 0.0;

    for(uint32_t i = 0; i < m_restart_iter; i++) {
        auto tmp_var = m_rotation_context->ComputeOutputVariance(out_var);
        auto cluster_var = tmp_var * iter_ratio;
        out_var = m_conversion_context->ComputeOutputVariance(cluster_var);
    }

    return out_var;
}

Container SAPDecompositionContext::GetOutputContainer(Container container) const {
    auto dummy_rlwe_in = m_rotation_context->GetInputContainer();
    auto dummy_rlwe_out = std::dynamic_pointer_cast<RLWEContainerImpl>(m_rotation_context->GetOutputContainer(dummy_rlwe_in));

    auto o_var = ComputeOutputVariance(0.0);

    return std::make_shared<LWEContainerImpl>(dummy_rlwe_out->GetN(), dummy_rlwe_out->GetQ(), o_var);
}

const uint64_t SAPDecompositionContext::GetDefaultRadix() const {
    return m_default_radix;
}

const uint64_t SAPDecompositionContext::GetRestartIteration() const {
    return m_restart_iter;
}

OperatorID SAPDecompositionContext::GetOperatorID() const {
    return DECOMP_SAP;
}

void SAPDecompositionContext::SetDefaultRadix(uint64_t new_radix) {
    m_default_radix = new_radix;
    // TODO update restart iter
}

void SAPDecompositionContext::SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context) {
    m_rotation_context = new_rot_context;
}

void SAPDecompositionContext::SetPackingContext(std::shared_ptr<LWEtoRLWEPackingContext> new_trace_context) {
    m_packing_context = new_trace_context;
}

void SAPDecompositionContext::SetLWEConversionContext(std::shared_ptr<LWEConversionParameters> new_lwe_context) {
    m_conversion_context = new_lwe_context;
}

std::unique_ptr<SAPDecomposer> SAPDecompositionContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    // TODO
}

void SAPDecomposer::HomTrunc(uint64_t *__restrict output, uint64_t *__restrict input, uint64_t radix) {
    // TODO
}

void SAPDecomposer::ResetAccumulator(uint64_t *acc) {
    // TODO
}

void SAPDecomposer::Decompose(uint64_t *output, const uint64_t *const input, uint64_t radix) {
    // TODO
}

void SAPDecomposer::Decompose(std::vector<uint64_t> &output, const std::vector<uint64_t> &input, uint64_t radix) {
    Decompose(output.data(), input.data(), radix);
}

const std::shared_ptr<SAPDecompositionContext> SAPDecomposer::GetContext() const {
    return m_context;
}

SAPDecomposer::SAPDecomposer(std::shared_ptr<SAPDecompositionContext> ctx, std::unique_ptr<BlindRotator> rotator,
                             std::unique_ptr<LWEtoRLWEPacker> packer, std::unique_ptr<LWEtoLWEConverter> conv, uint64_t max_radix,
                             uint64_t lwe_sk_hamming_weight, uint64_t reset_period) : m_context(std::move(ctx)),
                             m_rotator(std::move(rotator)), m_packer(std::move(packer)), m_lwe_converter(std::move(conv)),
                             m_max_radix(max_radix), m_beta(lwe_sk_hamming_weight), m_restart_iteration(reset_period) {

}