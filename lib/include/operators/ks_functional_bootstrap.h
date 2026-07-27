//
// Created by leonard on 7/27/26.
//

#ifndef TOOTHPASTE_KS_FUNCTIONAL_BOOTSTRAP_H
#define TOOTHPASTE_KS_FUNCTIONAL_BOOTSTRAP_H

#include "interfaces/functional_bootstrap.h"
#include "interfaces/operator_context.h"
#include "interfaces/blindrotation_operator.h"

#include "operators/endo_glwe_conversion.h"

struct KSFunctionEvaluator;

struct KSFunctionEvaluationContext : OperatorContext<FunctionEvaluator>, public std::enable_shared_from_this<KSFunctionEvaluationContext> {

    KSFunctionEvaluationContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
                                std::shared_ptr<LWEConversionContext> converter);

    [[nodiscard]] std::shared_ptr<OperatorContext<BlindRotator>> GetBlindRotationContext() const;

    [[nodiscard]] std::shared_ptr<LWEConversionContext> GetLWEConversionContext() const;

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    [[nodiscard]] std::unique_ptr<FunctionEvaluator> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    /** Describes the input LWE sample
     *
     * @return A Container describing the expected input
     */
    [[nodiscard]] Container GetInputContainer() const override;

    /** Describes the output LWE sample. Note that several input containers may be present
     * - LWEContainer (ignored)
     * - RLWEContainer, if the bootstrapping accumulator is provided
     * - FlagContainer, to specify which output type is expected
     *
     * @param input Input container description.
     * @return LWEContainer
     */
    [[nodiscard]] Container GetOutputContainer(Container input) const override;

    /** Returns an ID describing the Operator
     *
     * @return The ID
     */
    [[nodiscard]] OperatorID GetOperatorID() const override;

private:

    std::shared_ptr<OperatorContext<BlindRotator>> m_rotation_context;
    std::shared_ptr<LWEConversionContext> m_lwe_to_lwe_converter;

};

struct KSFunctionEvaluator : FunctionEvaluator {

    /** Function evaluation with a provided function mapping Z_{2N} to Z_Q
     * where N is the ring dimension and Q is the ring modulus
     *
     * @param result Buffer for output LWE
     * @param function Function to be evaluated
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(uint64_t* result, std::function<uint64_t(uint64_t)>& function, OutputKeyType output) override;

    /** Function evaluation with a provided function mapping Z_{2N} to Z_Q
     * where N is the ring dimension and Q is the ring modulus
     *
     * @param result Buffer for output LWE
     * @param function Function to be evaluated
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(std::vector<uint64_t>& result, std::function<uint64_t(uint64_t)>& function, OutputKeyType output) override;

    /** Function evaluation with a (possibly encrypted) LUT in lut_rlwe
     * where N is the ring dimension and Q is the ring modulus
     *
     * @param result Buffer for output LWE
     * @param function LUT to perform the evaluation on
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(uint64_t* RESTRICTED result, uint64_t* RESTRICTED lut_rlwe, OutputKeyType output) override;

    /** Function evaluation with a (possibly encrypted) LUT in lut_rlwe
     * where N is the ring dimension and Q is the ring modulus
     *
     * @param result Buffer for output LWE
     * @param function LUT to perform the evaluation on
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(std::vector<uint64_t>& result, std::vector<uint64_t>& lut_rlwe, OutputKeyType output) override;

    ~KSFunctionEvaluator() = default;

private:

    friend struct KSFunctionEvaluationContext;

    KSFunctionEvaluator(std::shared_ptr<const KSFunctionEvaluationContext> ctx,
                        std::unique_ptr<BlindRotator> rotator,
                        std::unique_ptr<LWEtoLWEConverter> conv);

    /* pointer to context that constructed the operator */
    std::shared_ptr<const KSFunctionEvaluationContext> m_context;
    /* internal blind-rotator */
    std::unique_ptr<BlindRotator> m_rotator;
    /* converter from LWE with ring dimension/modulus to target output LWE */
    std::unique_ptr<LWEtoLWEConverter> m_converter;
    /* temporary buffer for blind-rotation */
    AlignedVector m_acc_buffer;

};


#endif //TOOTHPASTE_KS_FUNCTIONAL_BOOTSTRAP_H
