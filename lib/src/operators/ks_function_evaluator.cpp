//
// Created by leonard on 7/31/26.
//

#include <numeric>
#include "operators/ks_function_evaluator.h"

std::shared_ptr<OperatorContext<BlindRotator>> KSFunctionEvaluationContext::GetBlindRotationContext() const {
    return m_blind_rotation_context;
}

std::shared_ptr<LWEConversionContext> KSFunctionEvaluationContext::GetLWEConversionContext() const {
    return m_conversion_context;
}


KSFunctionEvaluator::KSFunctionEvaluator(std::shared_ptr<const KSFunctionEvaluationContext>& context,
                                         std::unique_ptr<BlindRotator> rotator,
                                         std::unique_ptr<LWEtoLWEConverter> converter, bool estimate_padding,
                                         uint64_t padding_width) : m_context(context), m_rotator(std::move(rotator)), m_converter(std::move(converter)), m_estimate_padding(estimate_padding), m_padding_width(padding_width) {}

uint64_t KSFunctionEvaluator::EstimatePadding(const uint64_t *poly) const {

   auto in = m_context->GetBlindRotationContext()->GetInputContainer();
   auto out = std::dynamic_pointer_cast<RLWEContainerImpl>(m_context->GetOutputContainer(in));

   auto N = out->GetN();

   std::vector<uint64_t> block_sizes;
   block_sizes.reserve(N);
   block_sizes.emplace_back(1);

   auto current_block_index = 0;

   for(uint64_t i = 1; i < N; i++) {
       // tee-hee I like to be """clever"""
       if (poly[i] == poly[i - 1])
           block_sizes[current_block_index]++;
       else {
           ++current_block_index;
           block_sizes.emplace_back(1);
       }
   }

   auto padding_width = block_sizes[0];
   for(uint64_t i = 1; i < block_sizes.size(); i++) {
       padding_width = std::gcd(padding_width, block_sizes[i]);
   }

    return padding_width;
}

uint64_t KSFunctionEvaluator::EstimatePadding(std::vector<uint64_t> &poly) const {
    return EstimatePadding(poly.data());
}

std::shared_ptr<const KSFunctionEvaluationContext> KSFunctionEvaluator::GetContext() const {
    return m_context;
}

uint64_t KSFunctionEvaluator::GetEstimatePadding() const {
    return m_estimate_padding;
}

uint64_t KSFunctionEvaluator::GetPaddingWidth() const {
    return m_padding_width;
}