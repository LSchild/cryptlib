//
// Created by leonard on 7/27/26.
//

#include "operators/ks_functional_bootstrap.h"

KSFunctionEvaluationContext::KSFunctionEvaluationContext(
        std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
        std::shared_ptr<LWEConversionContext> converter) : m_rotation_context(std::move(blind_rotation_context)), m_lwe_to_lwe_converter(std::move(converter)) {

}

std::shared_ptr<OperatorContext<BlindRotator>> KSFunctionEvaluationContext::GetBlindRotationContext() const {
    return m_rotation_context;
}

std::shared_ptr<LWEConversionContext> KSFunctionEvaluationContext::GetLWEConversionContext() const {
    return m_lwe_to_lwe_converter;
}

OperatorID KSFunctionEvaluationContext::GetOperatorID() const {
    return FUNC_BOOT_KS;
}

Container KSFunctionEvaluationContext::GetInputContainer() const {

    auto br_input_gen = m_rotation_context->GetInputContainer();
    auto br_input = std::dynamic_pointer_cast<TupleContainerImpl>(br_input_gen);

    return br_input->GetElem(0);
}

Container KSFunctionEvaluationContext::GetOutputContainer(Container input) const {

    auto br_input_gen = m_rotation_context->GetInputContainer();

    auto br_input = std::dynamic_pointer_cast<TupleContainerImpl>(br_input_gen);
    auto lwe_input_gen = br_input->GetElem(0);
    auto lwe_input = std::dynamic_pointer_cast<LWEContainerImpl>(lwe_input_gen);
    auto rlwe_input_gen = br_input->GetElem(1);
    auto rlwe_input = std::dynamic_pointer_cast<RLWEContainerImpl>(rlwe_input_gen);

    // todo
    uint64_t output_modulus = m_lwe_to_lwe_converter->GetModulus();

    auto lwe_conversion_in = m_lwe_to_lwe_converter->GetInputContainer();
    auto lwe_conversion_out = m_lwe_to_lwe_converter->GetOutputContainer()



}