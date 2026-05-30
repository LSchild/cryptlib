//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_AUTOMORPHISM_EVALUATION_H
#define LARGE_FUNCTIONS_AUTOMORPHISM_EVALUATION_H

#include <cassert>
#include "operators/endo_glwe_conversion.h"
#include "rlwe_key_switching.h"
#include "setup.h"

struct AutomorphismEvaluator;
struct AutomorphismContext : public OperatorContext<AutomorphismEvaluator>,
                             public std::enable_shared_from_this<AutomorphismContext> {

    AutomorphismContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N, uint64_t basis, uint64_t digits, double std, uint32_t automorphism_index);

    AutomorphismContext(KeyDistribution source_key_distribution, std::shared_ptr<intel::hexl::NTT>, uint64_t basis, uint64_t digits, double std, uint32_t automorphism_index);

    AutomorphismContext(const AutomorphismContext &other);

    void SetAutomorphismIndex(uint32_t idx);

    [[nodiscard]] uint32_t GetAutomorphismIndex() const;

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

    [[nodiscard]] KeyDistribution GetSourceKeyDistribution() const;

    [[nodiscard]] uint64_t GetModulus() const;

    [[nodiscard]] uint64_t GetDimension() const;

    [[nodiscard]] uint64_t GetGadgetBasis() const;

    [[nodiscard]] uint64_t GetGadgetBasisLog2() const;

    [[nodiscard]] uint64_t GetGadgetDigits() const;

    [[nodiscard]] double GetStd() const;

    [[nodiscard]] std::shared_ptr<intel::hexl::NTT> GetNTT() const;

    void SetSourceKeyDistribution(KeyDistribution distribution);

    void SetModulus(uint64_t modulus);

    void SetDimension(uint64_t input_dimension) ;

    void SetGadgetBasis(uint64_t basis) ;

    void SetGadgetDigits(uint64_t digits) ;

    void SetStd(double std) ;

    void SetNTT(std::shared_ptr<intel::hexl::NTT>);

private:

    std::shared_ptr<RLWEConversionParameters> m_rlwe_conversion;
    uint32_t m_automorphism_index;
};

struct AutomorphismEvaluator {
    friend AutomorphismContext;

    // input assumed in [a=COEF|b=COEF] form
    void Eval(uint64_t *output, const uint64_t *const input);

    // input assumed in [a=COEF|b=NTT] form
    void EvalSwitchOnly(uint64_t* output, const uint64_t* const input);

    // input assumed in [a=COEF|b=COEF] form
    void Eval(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);


private:




    AutomorphismEvaluator(std::shared_ptr<const AutomorphismContext> params, std::unique_ptr<RLWEtoRLWEConverter> conv);

    void ApplyAutomorphism(uint64_t* output, const uint64_t *const input);

    void ApplyAutomorphism(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);

    bool m_params_set = false;
    bool m_keys_generated = false;

    std::shared_ptr<const AutomorphismContext> m_automorphism_params;
    std::unique_ptr<RLWEtoRLWEConverter> m_converter;
    AlignedVector m_auto_buffer;
};

// OLD VERSION
class AutomorphismKey : public RLWEKeySwitchingKey {

public:

    AutomorphismKey(std::shared_ptr<RingGSWCryptoParams>& params, NativePoly sk, uint32_t automorphism_idx);

    RLWECiphertext AutomorphismTransform(const RLWECiphertext& source) const;

    RLWECiphertext AutomorphismTransform(const NativePoly& A, const NativePoly& B);

    static NativePoly ApplyAutomorphism(const NativePoly& A, uint32_t automorphism_idx);

private:

    uint32_t m_automorphism_idx;

};

class FastAutomorphismKey  {

public:

    FastAutomorphismKey(std::shared_ptr<intel::hexl::NTT> &ntt_engine, uint32_t L, NativePoly sk, uint32_t automorphism_idx);


    RLWECiphertext AutomorphismTransform(const NativePoly& A, const NativePoly& B);


private:

    uint32_t m_automorphism_idx;
    RGSWSample m_key;
    std::shared_ptr<intel::hexl::NTT> ntt_engine;


};


#endif //LARGE_FUNCTIONS_AUTOMORPHISM_EVALUATION_H
