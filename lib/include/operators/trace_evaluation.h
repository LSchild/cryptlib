//
// Created by leonard on 5/30/26.
//

#ifndef LARGE_FUNCTIONS_TRACE_EVALUATION_H
#define LARGE_FUNCTIONS_TRACE_EVALUATION_H

#include "backend/backend.h"

#include "interfaces/operator_context.h"
#include "operators/automorphism_evaluation.h"

struct TraceEvaluator;

struct TraceEvaluationContext : public OperatorContext<TraceEvaluator>,
public std::enable_shared_from_this<TraceEvaluationContext> {

    TraceEvaluationContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N, uint64_t basis, uint64_t digits, double std);

    TraceEvaluationContext(KeyDistribution source_key_distribution, std::shared_ptr<MathWorker> ntt, uint64_t basis, uint64_t digits, double std);

    TraceEvaluationContext(const TraceEvaluationContext& other);

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    /** On input a given KeyBundle, constructs the operator
     *
     * @param keys The KeyBundle as required by the operator
     * @return Pointer to instance of the operator
     */
    [[nodiscard]] std::unique_ptr<TraceEvaluator> ConstructOperator(const std::vector<GenericKey>& keys) const override;

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

    [[nodiscard]] KeyDistribution GetSourceKeyDistribution() const;

    [[nodiscard]] uint64_t GetModulus() const;

    [[nodiscard]] uint64_t GetDimension() const;

    [[nodiscard]] uint64_t GetGadgetBasis() const;

    [[nodiscard]] uint64_t GetGadgetBasisLog2() const;

    [[nodiscard]] uint64_t GetGadgetDigits() const;

    [[nodiscard]] double GetStd() const;

    [[nodiscard]] std::shared_ptr<MathWorker> GetNTT() const;

    void SetSourceKeyDistribution(KeyDistribution distribution);

    void SetModulus(uint64_t modulus);

    void SetDimension(uint64_t input_dimension) ;

    void SetGadgetBasis(uint64_t basis) ;

    void SetGadgetDigits(uint64_t digits) ;

    void SetStd(double std) ;

    void SetNTT(std::shared_ptr<MathWorker>);

private:

    std::shared_ptr<AutomorphismContext> m_auto_context;


};

struct TraceEvaluator {
    friend struct TraceEvaluationContext;

    // input assumed in [a=COEF|b=COEF] form
    /**
     * Evaluates the homomorphic trace, that is, given a RLWE(p) = [a,b] in [COEF, COEF] format,
     * returns a RLWE(p_0)
     * @param output RLWE(p)
     * @param input RLWE(p_0)
     */
    void Eval(uint64_t* output, const uint64_t* input);

    // input assumed in [a=COEF|b=COEF] form

    /**
     * Evaluates the homomorphic trace, that is, given a RLWE(p) = [a,b] in [COEF, COEF] format,
     * returns a RLWE(p_0)
     * @param output RLWE(p)
     * @param input RLWE(p_0)
     */
    void Eval(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);

    [[nodiscard]] std::shared_ptr<const TraceEvaluationContext> GetContext() const;

    void EvalAuto(uint64_t* output, const uint64_t* input, uint64_t auto_idx);

    void EvalAuto(std::vector<uint64_t >& output, const std::vector<uint64_t>& input, uint64_t auto_idx);
private:

    TraceEvaluator(std::shared_ptr<const TraceEvaluationContext> context, std::vector<std::unique_ptr<AutomorphismEvaluator>> evaluators);

    std::shared_ptr<const TraceEvaluationContext> m_context;
    std::vector<std::unique_ptr<AutomorphismEvaluator>> m_trace_evaluators;

};


#endif //LARGE_FUNCTIONS_TRACE_EVALUATION_H
