//
// Created by leonard on 5/23/26.
//

#ifndef LARGE_FUNCTIONS_RADIX_DECOMPOSITION_H
#define LARGE_FUNCTIONS_RADIX_DECOMPOSITION_H

#include "common_types.h"
#include "blindrotation_common.h"
#include "automorphism_key.h"

template<typename BR>
struct SAPRadixDecomposer;

template<typename BR>
struct SAPRadixDecompositionContext : OperatorContext<SAPRadixDecomposer<BR>, BlindRotationKeys>,
 public std::enable_shared_from_this<SAPRadixDecompositionContext<BR>> {

    SAPRadixDecompositionContext(std::shared_ptr<OperatorContext<BR, BlindRotationKeys>>& rot_context,
                                 AutomorphismParameters params, uint64_t radix, KeyDistribution sk_distr, uint64_t sk_hamming)
                                 : m_br_context(rot_context), m_auto_params(params), m_lwe_sk_distr(sk_distr), m_radix(radix) {
        static_assert(std::is_base_of<BlindRotator, BR>::value, "Type BR is not a blind-rotation operator");

    }


    [[nodiscard]] const std::shared_ptr<OperatorContext<BR,BlindRotationKeys>>& GetBlindRotationContext() const {
        return m_br_context;
    }

    [[nodiscard]] const AutomorphismParameters& GetAutomorphismParameters() const {
        return m_auto_params;
    }

    [[nodiscard]] uint64_t GetRadix() const {
        return m_radix;
    }

    [[nodiscard]] KeyDistribution GetKeyDistribution() const {
        return m_lwe_sk_distr;
    }

    [[nodiscard]] uint64_t GetRestartIteration() const {
        return m_restart_iteration;
    }

    [[nodiscard]] uint64_t GetHammingWeight() const {
        return m_hamming;
    }

    void SetBlindRotationContext(std::shared_ptr<OperatorContext<BR,BlindRotationKeys>> new_rot_params) {
        m_br_context = new_rot_params;
    }

    void SetAutomorphismParameters(const AutomorphismParameters& new_auto_params) {
        m_auto_params = new_auto_params;
    }

    void SetRadix(uint64_t new_radix) {
        m_radix = new_radix;
    }

    void SetHamming(uint64_t new_hamming) {
        m_hamming = new_hamming;
    }

    void SetKeyDistribution(KeyDistribution new_distribution) {
        m_lwe_sk_distr = new_distribution;
    }

    void SetRestartIteration(uint64_t new_iter) {
        m_restart_iteration = new_iter;
    }

    [[nodiscard]] uint64_t DetermineRestartThreshold() const {
        // TODO
        return 0;
    }

    [[nodiscard]] Container GetInputContainer() const override {
        auto params_br = std::dynamic_pointer_cast<TupleContainerImpl>(m_br_context->GetInputContainer());
        return params_br->GetElem(0);
    }

    [[nodiscard]]  Container GetOutputContainer() const override {
        auto br_output = std::dynamic_pointer_cast<RLWEContainerImpl>(m_br_context->GetOutputContainer());

        // TODO: Fix
        RLWEContainer params_rlwe = std::dynamic_pointer_cast<RLWEContainerImpl>(br_output);

        return params_rlwe;
    }

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override {
        // TODO
        return 0.0;
    }

    [[nodiscard]] OperatorID GetOperatorID() const override {
        return DECOMP_SAP;
    }

    virtual std::shared_ptr<SAPRadixDecomposer<BR>> ConstructOperator(BlindRotationKeys& keys) const {

        // TODO: Fixme
        std::shared_ptr<SAPRadixDecomposer<BR>> op = std::make_shared<SAPRadixDecomposer<BR>>();

        return op;
    }


    std::shared_ptr<OperatorContext<BR, BlindRotationKeys>> m_br_context;
    AutomorphismParameters m_auto_params;

    KeyDistribution m_lwe_sk_distr;

    uint64_t m_radix;
    uint64_t m_hamming;
    uint64_t m_restart_iteration;

};

template<typename BR>
struct SAPRadixDecomposer : RadixDecomposer<SAPRadixDecompositionContext<BR>> {

    void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) {

    }

    void KeyGen(const std::vector<uint64_t>& source_key, const std::vector<uint64_t>&  target_key) {
        KeyGen(source_key.data(), target_key.data());
    }

    void Decompose(uint64_t* output, const uint64_t*const input, uint64_t radix) {

    }

    void Decompose(std::vector<uint64_t>& output, const std::vector<uint64_t>& input, uint64_t radix) {
        Decompose(output.data(), input.data(), radix);
    }

    [[nodiscard]] virtual Container GetInputContainer() const {
        return m_context->GetInputContainer();
    }

    [[nodiscard]] virtual Container GetOutputContainer() const {
        return m_context->GetOutputContainer();
    }

    void PhaseSplit(uint64_t* output_buffer, uint64_t* input_lwe);

    void HomTrunc(uint64_t* output_rlwe, uint64_t* input_rlwe);



    [[nodiscard]] const SAPRadixDecompositionContext<BR>& GetParams() const {
        return m_context;
    }

    bool m_params_set = false;
    bool m_keys_generated = false;

    std::shared_ptr<const SAPRadixDecompositionContext<BR>> m_context;

    std::shared_ptr<BR> m_rotator;
    std::vector<std::shared_ptr<AutomorphismEvaluator>> m_auto_converters;

    AlignedVector m_acc;
    AlignedVector m_phase_buffer;

};

#endif //LARGE_FUNCTIONS_RADIX_DECOMPOSITION_H
