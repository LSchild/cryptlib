//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_AUTOMORPHISM_EVALUATION_H
#define LARGE_FUNCTIONS_AUTOMORPHISM_EVALUATION_H

#include <cassert>
#include "operators/endo_glwe_conversion.h"

struct AutomorphismEvaluator;
struct AutomorphismContext : public OperatorContext<AutomorphismEvaluator>,
                             public std::enable_shared_from_this<AutomorphismContext> {

    /**
     * Constructor for automorphism X -> X^t evaluator over Z_Q[X]/(X^N + 1)
     * @param source_key_distribution Distribution of RLWE key (Binary/Ternary/Gaussian)
     * @param modulus Q
     * @param N N
     * @param basis Gadget Basis
     * @param digits Gadget Digits
     * @param std Standard Deviation for automorphism key
     * @param automorphism_index t with gcd(t, N) = 1
     */
    AutomorphismContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N, uint64_t basis, uint64_t digits, double std, uint32_t automorphism_index);

    /**
     * Constructor for automorphism X -> X^t evaluator over Z_Q[X]/(X^N + 1)
     *
     * Modulus and ring dimension are obtained from the ntt engine
     *
     * @param source_key_distribution Distribution of the RLWE key (Binary/Ternary/Gaussian)
     * @param ntt NTT Engine
     * @param basis Gadget Basis
     * @param digits Gadget Digits
     * @param std Standard Deviation for automorphism key
     * @param automorphism_index t with gcd(t, N) = 1
     */
    AutomorphismContext(KeyDistribution source_key_distribution, std::shared_ptr<MathWorker> ntt, uint64_t basis, uint64_t digits, double std, uint32_t automorphism_index);

    /**
     * Constructor for automorphism X -> X^t evaluator over Z_Q[X]/(X^N + 1)
     * @param other Other context
     */
    AutomorphismContext(const AutomorphismContext &other);

    /**
     * Modifies the index of the automorphism X -> X^t
     * @param idx New value of t
     */
    void SetAutomorphismIndex(uint32_t idx);

    /**
     *
     * @return Index t of the automorphism X-> X^t
     */
    [[nodiscard]] uint32_t GetAutomorphismIndex() const;

    /**
     * Computes variance induced by automorphism evaluation
     * @param input_variance Variance on the input container
     * @return Output variance as a function of the output variance
     */
    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    /** On input a given KeyBundle, constructs the operator
     *
     * @param keys The KeyBundle as required by the operator
     * @return Pointer to instance of the operator
     */
    [[nodiscard]] std::unique_ptr<AutomorphismEvaluator> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @return A Container describing the expected input
     */
    [[nodiscard]] Container GetInputContainer() const override;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @param input Input container description. Necessary as some operators are not additive in variance.
     *        if equal to nullptr, is ignored.
     * @return A Container describing the calculated output
     */
    [[nodiscard]] Container GetOutputContainer(Container input) const override;

    /** Returns an ID describing the Operator
     *
     * @return The ID
     */
    [[nodiscard]] OperatorID GetOperatorID() const override;

    /**
     *
     * @return Key distribution
     */
    [[nodiscard]] KeyDistribution GetSourceKeyDistribution() const;

    /**
     *
     * @return Ring Q of the ring R=Z_Q[X]/(X^N + 1)
     */
    [[nodiscard]] uint64_t GetModulus() const;

    /**
     *
     * @return Ring dimension of the ring R=Z_Q[X]/(X^N + 1)
     */
    [[nodiscard]] uint64_t GetDimension() const;

    /**
     *
     * @return Gadget Basis L
     */
    [[nodiscard]] uint64_t GetGadgetBasis() const;

    /**
     * Returns the bits/log2 of the gadget basis for efficient gadget computation
     * @return Base 2 logarithm of gadget basis
     */
    [[nodiscard]] uint64_t GetGadgetBasisLog2() const;

    /**
     *
     * @return Gadget Digits
     */
    [[nodiscard]] uint64_t GetGadgetDigits() const;

    /**
     *
     * @return Standard Deviation of the key
     */
    [[nodiscard]] double GetStd() const;

    /**
     *
     * @return Pointer to NTT engine
     */
    [[nodiscard]] std::shared_ptr<MathWorker> GetNTT() const;

    /**
     * Update declared source key distribution
     * @param distribution New distribution
     */
    void SetSourceKeyDistribution(KeyDistribution distribution);

    /**
     * Update declared Q of the ring R=Z_Q[X]/(X^N + 1)
     * @param modulus New Q Q
     */
    void SetModulus(uint64_t modulus);

    /**
     * Update declared ring dimension of the ring R=Z_Q[X]/(X^N + 1)
     * @param input_dimension New dimension N
     */
    void SetDimension(uint64_t input_dimension) ;

    /**
     * Update Gadget Basis
     * @param basis New basis
     */
    void SetGadgetBasis(uint64_t basis) ;

    /**
     * Update Gadget digits
     * @param digits New digits
     */
    void SetGadgetDigits(uint64_t digits) ;

    /**
     * Update key standard deviation
     * @param std new value
     */
    void SetStd(double std) ;

    /**
     * Update ntt/math engine
     * @param ntt engine
     */
    void SetNTT(std::shared_ptr<MathWorker> ntt);

private:

    // internally, an automorphism evaluator is near equivalent to a rlwe to rlwe switch to we keep a context here
    std::shared_ptr<RLWEConversionContext> m_rlwe_conversion;
    // index t of the automorphism X -> X^t
    uint32_t m_automorphism_index;
};

struct AutomorphismEvaluator {
    friend AutomorphismContext;

    // input assumed in [a=COEF|b=COEF] form
    /**
     * Performs the automorphism evaluation
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer, input assumed to be in coefficient form and different from output
     */
    void Eval(uint64_t *output, const uint64_t *const input);

    /**
    * Performs *only* the key-switching subprocedure of the automorphism evaluation
    * @param output Pointer to output buffer
    * @param input Pointer to input buffer, input=[A,B] assumed to be in [coefficient,ntt] form and different from output
    */
    void EvalSwitchOnly(uint64_t* output, const uint64_t* const input);

    /**
    * Performs the automorphism evaluation
    * @param output Reference to output vector
    * @param input Reference to input vector, input assumed to be in coefficient form
    */
    void Eval(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);


private:

    /**
     * Private constructor for the automorphism evaluator. Note that the RLWEtoRLWEConverter pointer will be std::move-d !
     * @param params Automorphism context
     * @param conv Pointer to RLWE converter, switching from a key s(X^t) to s(X) for an automorphism X -> X^t
     *
     */
    AutomorphismEvaluator(std::shared_ptr<const AutomorphismContext> params, std::unique_ptr<RLWEtoRLWEConverter> conv);

    /**
     * Evaluates the automorphism / negacyclic permutation
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer, must be different than input
     */
    void ApplyAutomorphism(uint64_t* output, const uint64_t *const input);

    /**
    * Evaluates the automorphism / negacyclic permutation
    * @param output Output vector
    * @param input Input Vector
    */
    void ApplyAutomorphism(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);

    // debug flags
    bool m_params_set = false;
    bool m_keys_generated = false;

    // pointer to parent context
    std::shared_ptr<const AutomorphismContext> m_automorphism_params;
    // pointer for switch from s(X^t) to s(X)
    std::unique_ptr<RLWEtoRLWEConverter> m_converter;
    // scratch space
    AlignedBuffer m_auto_buffer;
};

#endif //LARGE_FUNCTIONS_AUTOMORPHISM_EVALUATION_H
