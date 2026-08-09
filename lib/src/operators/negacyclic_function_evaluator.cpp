//
// Created by leonard on 7/27/26.
//

#include <iostream>
#include "operators/negacyclic_function_evaluator.h"
#include "utils/generic_utils.h"
#include "static/sample_extraction.h"

NegacyclicFunctionEvaluationContext::NegacyclicFunctionEvaluationContext(
        std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
        std::shared_ptr<LWEConversionContext> converter) : m_rotation_context(std::move(blind_rotation_context)), m_lwe_to_lwe_converter(std::move(converter)) {

}

std::shared_ptr<OperatorContext<BlindRotator>> NegacyclicFunctionEvaluationContext::GetBlindRotationContext() const {
    return m_rotation_context;
}

std::shared_ptr<LWEConversionContext> NegacyclicFunctionEvaluationContext::GetLWEConversionContext() const {
    return m_lwe_to_lwe_converter;
}

OperatorID NegacyclicFunctionEvaluationContext::GetOperatorID() const {
    return NEGACYCLIC_FUNC_BOOT;
}

Container NegacyclicFunctionEvaluationContext::GetInputContainer() const {

    // todo: revisit
    auto br_input_gen = m_rotation_context->GetInputContainer();
    auto br_input = std::dynamic_pointer_cast<TupleContainerImpl>(br_input_gen);

    return br_input->GetElem(0);
}

Container NegacyclicFunctionEvaluationContext::GetOutputContainer(Container input) const {

    auto br_input_gen = m_rotation_context->GetInputContainer();

    auto br_input = std::dynamic_pointer_cast<TupleContainerImpl>(br_input_gen);
    auto rlwe_input_gen = br_input->GetElem(1);
    auto rlwe_input = std::dynamic_pointer_cast<RLWEContainerImpl>(rlwe_input_gen);

    std::vector<Container> all_containers;

    if (input->GetLabel() == Tuple) {
        auto tuple = std::dynamic_pointer_cast<TupleContainerImpl>(input);
        for(uint64_t i = 0; i < tuple->GetNumElems(); i++) {
            all_containers.push_back(tuple->GetElem(i));
        }
    } else {
        all_containers.push_back(input);
    }


    bool output_as_ring_key = false;
    long double acc_var = 0.0;

    for(auto& container : all_containers) {

        if (container->GetLabel() == Flag) {
            auto flag_container = std::dynamic_pointer_cast<FlagContainerImpl<bool>>(container);
            if (flag_container->GetName() == OutputKeyTypeLabels[OutputKeyType::BlindRotationKey]) {
                output_as_ring_key = true;
                continue;
            }
        }

        if (container->GetLabel() == RLWE) {
            auto rlwe_acc_container = std::dynamic_pointer_cast<RLWEContainerImpl>(container);
            acc_var = rlwe_acc_container->GetVariance();
            continue;
        }

        if (container->GetLabel() == LWE) {
            continue;
        }

        std::cerr << "Invalid input container type. Exiting..." << std::endl;
        std::exit(1);
    }

    auto post_blind_rotation_var = m_rotation_context->ComputeOutputVariance(acc_var);
    if (output_as_ring_key) {
        return std::make_shared<LWEContainerImpl>(rlwe_input->GetN(), rlwe_input->GetQ(), post_blind_rotation_var);
    } else {

        auto current_modulus = rlwe_input->GetQ();
        auto conv_modulus = m_lwe_to_lwe_converter->GetModulus();

        if (current_modulus != conv_modulus) {
            std::cerr << "Warning: Conversion Q differs from ring Q. The output variance computation does not take this into account yet" << std::endl;
        }

        auto lwe_conv_input_container = std::make_shared<LWEContainerImpl>(rlwe_input->GetN(), conv_modulus, post_blind_rotation_var);
        return m_lwe_to_lwe_converter->GetOutputContainer(lwe_conv_input_container);
    }

}

long double NegacyclicFunctionEvaluationContext::ComputeOutputVariance(long double input_variance) const {
    auto post_br_variance = m_rotation_context->ComputeOutputVariance(input_variance);
    auto post_switch_variance = m_lwe_to_lwe_converter->ComputeOutputVariance(post_br_variance);

    return post_switch_variance;
}

std::unique_ptr<FunctionEvaluator>
NegacyclicFunctionEvaluationContext::ConstructOperator(const std::vector<GenericKey> &keys) const {

    auto blind_rotator = m_rotation_context->ConstructOperator(keys);

    std::vector<GenericKey> switch_keys;
    switch_keys.push_back(keys[1]);
    switch_keys.push_back(keys[0]);

    auto converter = m_lwe_to_lwe_converter->ConstructOperator(switch_keys);

    return std::unique_ptr<NegacyclicFunctionEvaluator>(new NegacyclicFunctionEvaluator(shared_from_this(), std::move(blind_rotator), std::move(converter)));
}

NegacyclicFunctionEvaluator::NegacyclicFunctionEvaluator(std::shared_ptr<const NegacyclicFunctionEvaluationContext> ctx,
                                                         std::unique_ptr<BlindRotator> rotator,
                                                         std::unique_ptr<LWEtoLWEConverter> conv) :
                                         m_context(std::move(ctx)), m_rotator(std::move(rotator)), m_converter(std::move(conv)) {
    m_acc_buffer.resize(4 * conv->GetContext()->GetSourceDimension());
    std::fill(m_acc_buffer.begin(), m_acc_buffer.end(), 0);
}

void NegacyclicFunctionEvaluator::EvalFunc(uint64_t *__restrict result, const uint64_t *__restrict input_lwe,
                                           std::function<uint64_t(uint64_t)> &function, OutputKeyType output) {
    auto rot_ctx = m_rotator->GetContext();
    auto br_output = std::dynamic_pointer_cast<RLWEContainerImpl>(rot_ctx->GetOutputContainer(rot_ctx->GetInputContainer()));
    auto N = br_output->GetN();

    AlignedBuffer tmp_rlwe(2 * N);
    std::fill(tmp_rlwe.begin(), tmp_rlwe.end(), 0);

    for(uint64_t i = 0; i < N; i++) {
        tmp_rlwe[i + N] = function(i);
    }

    EvalFunc(result, input_lwe, tmp_rlwe.data(), output);
}


void NegacyclicFunctionEvaluator::EvalFunc(uint64_t * RESTRICTED result, const uint64_t* RESTRICTED input, uint64_t *RESTRICTED lut_rlwe, OutputKeyType output) {
    auto rot_ctx = m_rotator->GetContext();
    auto br_output = std::dynamic_pointer_cast<RLWEContainerImpl>(rot_ctx->GetOutputContainer(rot_ctx->GetInputContainer()));
    auto N = br_output->GetN();
    auto Q = br_output->GetQ();

    m_rotator->BlindRotate(m_acc_buffer.data(), input, lut_rlwe, true);
    SampleExtract(m_acc_buffer.data() + 2 * N, m_acc_buffer.data(), 0, N, Q);

    if (output == OutputKeyType::BlindRotationKey) {
        std::copy(m_acc_buffer.data() + 2 * N, m_acc_buffer.data() + 3 * N + 1, result);
    } else {
        m_converter->Convert(result, m_acc_buffer.data() + 2 * N);
    }
}

std::shared_ptr<const NegacyclicFunctionEvaluationContext> NegacyclicFunctionEvaluator::GetContext() const {
    return m_context;
}

void NegacyclicFunctionEvaluator::EvalFunc(std::vector<uint64_t> &result, const std::vector<uint64_t>& input, std::vector<uint64_t> &lut_rlwe, OutputKeyType output) {
    EvalFunc(result.data(), input.data(), lut_rlwe.data(), output);
}

void NegacyclicFunctionEvaluator::EvalFunc(std::vector<uint64_t> &result, const std::vector<uint64_t>& input, std::function<uint64_t(uint64_t)> &function,
                                           OutputKeyType output) {
    EvalFunc(result.data(), input.data(), function, output);
}

void NegacyclicFunctionEvaluator::Finalize(uint64_t *RESTRICTED result, const uint64_t *RESTRICTED input) {
    // todo Q switch
    m_converter->Convert(result, input);
}

void NegacyclicFunctionEvaluator::Finalize(std::vector<uint64_t> &result, const std::vector<uint64_t> &input) {
    Finalize(result.data(), input.data());
}



