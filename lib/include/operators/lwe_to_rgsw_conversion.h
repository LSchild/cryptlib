//
// Created by leonard on 4/28/26.
//

#ifndef LARGE_FUNCTIONS_LWE_TO_RGSW_CONVERSION_H
#define LARGE_FUNCTIONS_LWE_TO_RGSW_CONVERSION_H

#include <utility>

#include "static/container_types.h"
#include "interfaces/blindrotation_operator.h"

#include "operators/trace_evaluation.h"
#include "operators/endo_glwe_conversion.h"

#include "static/common_types.h"

#include "automorphism_evaluation.h"
#include "utils/math_utils.h"
#include "utils/generic_utils.h"

struct LWEtoRGSWConverter;

struct SchemeSwitchingContext :
        OperatorContext<LWEtoRGSWConverter>,
        public std::enable_shared_from_this<SchemeSwitchingContext> {

    SchemeSwitchingContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
                           std::shared_ptr<TraceEvaluationContext> trace_params,
                           std::shared_ptr<RLWEConversionContext> squaring_params, uint64_t output_digits, uint64_t basis) : m_rotation_context(std::move(blind_rotation_context)),
                                                                                                                             m_trace_params(std::move(trace_params)),
                                                                                                                             m_squaring_params(std::move(squaring_params)),
                                                                                                                             m_output_digits(output_digits), m_output_basis(basis),
                                                                                                                             m_output_basis_log2(IntLog2(basis)){


    }

    long double ComputeOutputVariance(long double input_variance = 0) const override;

    [[nodiscard]] const std::shared_ptr<OperatorContext<BlindRotator>> GetBlindRotationContext() const;

    [[nodiscard]] const std::shared_ptr<TraceEvaluationContext> GetTraceContext() const;

    [[nodiscard]] const std::shared_ptr<RLWEConversionContext> GetSquaringContext() const;

    [[nodiscard]] uint64_t GetOutputDigits() const;

    [[nodiscard]] uint64_t GetOutputBasis() const;

    [[nodiscard]] uint64_t GetOutputBasisLog2() const;

    void SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context);

    void SetTraceContext(std::shared_ptr<TraceEvaluationContext> new_trace_context);

    void SetSquaringParameters(std::shared_ptr<RLWEConversionContext> new_squaring_params);

    void SetOutputDigits(uint64_t new_digits);

    void SetOutputBasis(uint64_t new_basis);

    [[nodiscard]] Container GetInputContainer() const override;

    [[nodiscard]] Container GetOutputContainer(Container container) const override;

    std::unique_ptr<LWEtoRGSWConverter> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    [[nodiscard]] OperatorID GetOperatorID() const override;

    std::shared_ptr<OperatorContext<BlindRotator>> m_rotation_context;
    std::shared_ptr<TraceEvaluationContext> m_trace_params;
    std::shared_ptr<RLWEConversionContext> m_squaring_params;

    uint64_t m_output_digits;
    uint64_t m_output_basis;
    uint64_t m_output_basis_log2;

};

struct LWEtoRGSWConverter : SchemeConverter<SchemeSwitchingContext> {

    friend struct SchemeSwitchingContext;

    LWEtoRGSWConverter(std::shared_ptr<const SchemeSwitchingContext> context,
                       std::unique_ptr<BlindRotator> rotator,
                       std::unique_ptr<TraceEvaluator> trace_eval,
                       std::unique_ptr<RLWEtoRLWEConverter> squaring_eval);

    [[nodiscard]] const std::shared_ptr<const SchemeSwitchingContext> GetContext() const override;

    void Convert(uint64_t* output, const uint64_t* const input) override;

    void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) override;

    bool m_params_set = false;
    bool m_keys_generated = false;

    std::shared_ptr<const SchemeSwitchingContext> m_params;
    std::unique_ptr<BlindRotator> m_rotator;
    std::unique_ptr<TraceEvaluator> m_trace_eval;
    std::unique_ptr<RLWEtoRLWEConverter> m_square_eval;

    AlignedVector m_acc;
    AlignedVector m_extraction_buffer;
};

#endif //LARGE_FUNCTIONS_LWE_TO_RGSW_CONVERSION_H
