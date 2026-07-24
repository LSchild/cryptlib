//
// Created by leonard on 6/2/26.
//

#ifndef LARGE_FUNCTIONS_LWE_TO_RLWE_PACKING_H
#define LARGE_FUNCTIONS_LWE_TO_RLWE_PACKING_H

#include "interfaces/operator_context.h"
#include "operators/automorphism_evaluation.h"
#include "hexl/ntt/ntt.hpp"
#include <memory>

struct LWEtoRLWEPacker;

struct LWEtoRLWEPackingContext : public OperatorContext<LWEtoRLWEPacker>,
                                public std::enable_shared_from_this<LWEtoRLWEPackingContext> {

    LWEtoRLWEPackingContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N, uint64_t basis, uint64_t digits, double std);

    LWEtoRLWEPackingContext(KeyDistribution source_key_distribution, std::shared_ptr<MathWorker> ntt, uint64_t basis, uint64_t digits, double std);

    LWEtoRLWEPackingContext(const LWEtoRLWEPackingContext& other);

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    /** On input a given KeyBundle, constructs the operator
     *
     * @param keys The KeyBundle as required by the operator
     * @return Pointer to instance of the operator
     */
    [[nodiscard]] std::unique_ptr<LWEtoRLWEPacker> ConstructOperator(const std::vector<GenericKey>& keys) const override;

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

struct LWEtoRLWEPacker {
    friend struct LWEtoRLWEPackingContext;

    // input assumed in [a=COEF|b=COEF] form
    void Pack(uint64_t* output, const uint64_t* input, uint64_t len);

    // input assumed in [a=COEF|b=COEF] form
    void Pack(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);

    void PackConsecutively(uint64_t* output, const uint64_t* input, uint64_t len);

    [[nodiscard]] const std::shared_ptr<const LWEtoRLWEPackingContext> GetContext() const;
private:

    LWEtoRLWEPacker(std::shared_ptr<const LWEtoRLWEPackingContext> context, std::vector<std::unique_ptr<AutomorphismEvaluator>> evaluators);

    void L0Setup(uint64_t* __restrict output, const uint64_t* __restrict input);

    std::shared_ptr<const LWEtoRLWEPackingContext> m_context;
    std::vector<std::unique_ptr<AutomorphismEvaluator>> m_auto_evaluators;

};


#endif //LARGE_FUNCTIONS_LWE_TO_RLWE_PACKING_H
