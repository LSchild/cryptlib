//
// Created by leonard on 6/3/26.
//

#ifndef LARGE_FUNCTIONS_SAP_DECOMPOSITION_H
#define LARGE_FUNCTIONS_SAP_DECOMPOSITION_H

#include "interfaces/decomposition_operator.h"
#include "interfaces/blindrotation_operator.h"
#include "operators/lwe_to_rlwe_packing.h"

struct SAPDecomposer;
struct SAPDecompositionContext : OperatorContext<SAPDecomposer>,
public std::enable_shared_from_this<SAPDecompositionContext> {

    static const uint64_t IMPLICIT_SK_L0 = 32;

    SAPDecompositionContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
    std::shared_ptr<LWEtoRLWEPackingContext> packing_context,
            std::shared_ptr<LWEConversionParameters> lwe_conversion_context,
            uint64_t default_radix);

    /**
     * Computes output variance as a function of the input variance
     * In this case, the input is effectively ignored as the method is independent of the input variance
     * @param input_variance (ignored)
     * @return Output variance
     */
    long double ComputeOutputVariance(long double input_variance = 0) const override;

    /** Gets a pointer to the Rotation context
     *
     * @return BlindRotator
     */
    [[nodiscard]] std::shared_ptr<OperatorContext<BlindRotator>> GetBlindRotationContext() const;

    /** Gets a pointer to the Packing context
     *
     * @return Packer
     */
    [[nodiscard]] std::shared_ptr<LWEtoRLWEPackingContext> GetPackingContext() const;

    /** Gets a pointer to the Rotation context
     *
     * @return LWE to LWE conversion context
     */
    [[nodiscard]] std::shared_ptr<LWEConversionParameters> GetLWEConversionContext() const;

    /** Returns default radix for decomposition
     *
     * @return value of default radix
     */
    [[nodiscard]] const uint64_t GetDefaultRadix() const;

    /**
     * Returns the number of iterations after which we need to restart the accumulator
     * @return Restart iteration number
     */
    [[nodiscard]] const uint64_t GetRestartIteration() const;

    /**
     * Returns the number of iterations after which we need to reset the accumulator
     * as a function of a radix
     * @param radix
     * @param sk_hamming hamming weight of secret key
     * @return restart iteration number
     */
    [[nodiscard]] const uint64_t ComputeRestartIterationForRadix(uint64_t radix, uint64_t sk_hamming) const;

    /**
     * Sets m_rotation_context to new rotation context
     * @param new_rot_context new rotation context
     */
    void SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context);

    /**
     * Sets m_packing_context to new packing context
     * @param new_packing_context new packing context
     */
    void SetPackingContext(std::shared_ptr<LWEtoRLWEPackingContext> new_packing_context);

    /**
     * Sets m_conversion_context to new conversion context
     * @param new_packing_context new conversion context
     */
    void SetLWEConversionContext(std::shared_ptr<LWEConversionParameters> new_lwe_context);

    /**
     * Sets new default radix
     * @param new_radix new radix value
     */
    void SetDefaultRadix(uint64_t new_radix);

    /**
     * Returns pointer to container describing the input format
     * @return pointer to expected input container
     */
    [[nodiscard]] Container GetInputContainer() const override;

    /** Returns pointer to container describing the output as a function of the input
     *
     * @param container input container (ignored here)
     * @return output container
     */
    [[nodiscard]] Container GetOutputContainer(Container container) const override;

    /**
     * Constructs a new SAP decomposition context
     * @param keys Key bundle including keys[0] = key compatible with the given blind-rotation context
     * @return unique_ptr to SAPDecomposer instance
     */
    std::unique_ptr<SAPDecomposer> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    /**
     *
     * @return Operator id (DECOMP_SAP)
     */
    [[nodiscard]] OperatorID GetOperatorID() const override;

private:

    /* blind-rotation interface context */
    std::shared_ptr<OperatorContext<BlindRotator>> m_rotation_context;
    /* LWE to RLWE packing context */
    std::shared_ptr<LWEtoRLWEPackingContext> m_packing_context;
    /* LWE to LWE conversion context, used only for resetting */
    std::shared_ptr<LWEConversionParameters> m_conversion_context;

    uint64_t m_default_radix;
    uint64_t m_restart_iter;
};


struct SAPDecomposer : RadixDecomposer<SAPDecompositionContext> {


    /** Performs the decomposition
     *
     * Note that the radix parameter usually needs to be specified for
     * both the associated context, and the method.
     * This is due to the fact that often, for any fixed basis B, the operator will also work
     * for bases B' either s.t. B' | B or B' < B
     *
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer
     * @param radix Radix
     */
    void Decompose(uint64_t *output, const uint64_t *const input, uint64_t radix) override;

    /** Performs the decomposition
    *
    * Note that the radix parameter usually needs to be specified for
    * both the associated context, and the method.
    * This is due to the fact that often, for any fixed basis B, the operator will also work
    * for bases B' either s.t. B' | B or B' < B
    *
    * @param output Vector for output
    * @param input Vector for input
    * @param radix Radix
    */
    void Decompose(std::vector<uint64_t> &output, const std::vector<uint64_t> &input, uint64_t radix) override;

    /** Returns pointer to context that created the current operator.
    *
    * @return pointer to context
    */
    [[nodiscard]] const std::shared_ptr<SAPDecompositionContext> GetContext() const override;

    ~SAPDecomposer() override = default;

private:

    /**
     * Constructor for decomposition struct
     * @param ctx pointer to context / builder
     * @param rotator pointer to rotator, will be std::move-d !
     * @param packer pointer to lwe to rlwe packer, will be std::move-d !
     * @param max_radix maximum radix allowed during decomposition
     * @param lwe_sk_hamming_weight hamming weight of the LWE key
     * @param reset_period iteration number after which we need to refresh the accumulator
     */
    SAPDecomposer(std::shared_ptr<SAPDecompositionContext> ctx,
                  std::unique_ptr<BlindRotator> rotator,
                  std::unique_ptr<LWEtoRLWEPacker> packer,
                  std::unique_ptr<LWEtoLWEConverter> conv, uint64_t max_radix, uint64_t lwe_sk_hamming_weight,
                  uint64_t reset_period);

    /**
     * Given an input RLWE(X^c), returns RLWE(X^{floor(c / radix)})
     * @param output pointer to output buffer
     * @param input pointer to input buffer, assumed to be different from output
     * @param radix current radix, if smaller than m_max_radix may speed up this step
     * @param block_limit largest value of c/radix to consider assuming a bound on c is known
     */
    void HomTrunc(uint64_t* __restrict output, uint64_t* __restrict input, uint64_t radix, uint64_t block_limit);

    /**
     * After m_restart_iteration, the accumulator must be refreshed to allow the processing
     * of the next digits. This step, given a RLWE(X^c) with c < m_max_radix, returns a
     * new RLWE(X^c) with an error magnitude smaller than the input.
     * @param acc pointer to accumulator buffer, assumed to be given in COEF format
     */
    void ResetAccumulator(uint64_t* acc);

    /* context storing constants etc. */
    std::shared_ptr<SAPDecompositionContext> m_context;
    /* pointer to interface instance for blind-rotation */
    std::unique_ptr<BlindRotator> m_rotator;
    /* pointer to packing engine for truncation */
    std::unique_ptr<LWEtoRLWEPacker> m_packer;
    /* pointer to lwe-to-lwe keyswitch for reset/restart */
    std::unique_ptr<LWEtoLWEConverter> m_lwe_converter;

    /* hamming weight of a binary/ternary key */
    uint64_t m_beta;

    /* highest valid radix for decomposition */
    uint64_t m_max_radix;

    /* digit index after which we reset / restart the accumulator */
    uint64_t m_restart_iteration;

    /* temporary buffer required for packing */
    AlignedVector m_packing_buffer;
};

#endif //LARGE_FUNCTIONS_SAP_DECOMPOSITION_H
