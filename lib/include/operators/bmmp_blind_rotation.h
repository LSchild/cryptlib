//
// Created by leonard on 3/23/26.
//

#ifndef BMMP_BLIND_ROTATOR_H
#define BMMP_BLIND_ROTATOR_H

#include "backend/backend.h"
#include "static/common_types.h"
#include "interfaces/operator_context.h"
#include "interfaces/blindrotation_operator.h"
#include "base_crypto.h"
#include "mux_operator.h"

struct BMMPBlindRotator;


struct BMMPBlindRotationContext : public OperatorContext<BMMPBlindRotator>
        , public std::enable_shared_from_this<BMMPBlindRotationContext> {

    friend struct CGGIBlindRotator;

    explicit BMMPBlindRotationContext(KeyDistribution distr, uint64_t modulus, uint64_t ring_dim, uint64_t lwe_dim, uint64_t m_basis, uint64_t m_digits, double std, uint64_t step_size);

    explicit BMMPBlindRotationContext(KeyDistribution distr, std::shared_ptr<MathWorker> ntt, uint64_t lwe_dim, uint64_t m_basis, uint64_t m_digits, double std, uint64_t step_size);

    BMMPBlindRotationContext(const BMMPBlindRotationContext &other);

    void SetKeyDistribution(KeyDistribution dist);

    void SetModulus(uint64_t modulus);

    void SetRingDimension(uint64_t ring_dim);

    void SetStepSize(uint64_t step_size);

    void SetNTT(std::shared_ptr<MathWorker> ntt);

    void SetLWEDimension(uint64_t lwe_dim);

    void SetBlindRotationBasis(uint64_t L);

    void SetBlindRotationRGSWDigits(uint64_t digits);

    void SetStd(double std);

    [[nodiscard]] KeyDistribution GetKeyDistribution() const;

    [[nodiscard]] uint64_t GetModulus() const;

    [[nodiscard]] uint64_t GetRingDimension() const;

    [[nodiscard]] uint64_t GetStepSize() const;

    [[nodiscard]] std::shared_ptr<MathWorker> GetNTT() const;

    [[nodiscard]] uint64_t GetLWEDimension() const;

    [[nodiscard]] uint64_t GetBlindRotationBasis() const;

    [[nodiscard]] uint64_t GetBlindRotationBasisLog2() const;

    [[nodiscard]] uint64_t GetBlindRotationRGSWDigits() const;

    [[nodiscard]] double GetStd() const;

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    [[nodiscard]] std::unique_ptr<BMMPBlindRotator> ConstructOperator(const std::vector<GenericKey>& bundle) const override;

    [[nodiscard]] Container GetInputContainer() const override;

    [[nodiscard]] Container GetOutputContainer(Container input) const override;

    [[nodiscard]] OperatorID GetOperatorID() const override;

protected:

    KeyDistribution m_distribution;

    uint64_t m_modulus;
    uint64_t m_N;
    uint64_t m_n;

    uint64_t m_basis;
    uint64_t m_basis_log2;
    uint64_t m_digits;

    uint64_t m_step_size;

    double m_std;

    std::shared_ptr<MathWorker> m_ntt;

};



struct BMMPBlindRotator : public BlindRotator {

    friend struct BMMPBlindRotationContext;

    /** Builds a BMMP blind-rotation object for binary or ternary LWE keys.
     *
     * @param params Struct containing the parameters (modulus, ring dimension, ...)
     */
    explicit BMMPBlindRotator(const std::shared_ptr<const BMMPBlindRotationContext> &params);

    /** Performs BMMP blind-rotation
     *
     * @param lwe_vec vector containing the LWE sample [a_0, a_1, ..., a_{n - 1}, b]
     * @param rlwe_acc_vec RLWE accumulator in COEF format
     */
    void BlindRotate(std::vector<uint64_t> &result, const std::vector<uint64_t> &lwe_vec,
                     std::vector<uint64_t> &rlwe_acc_vec,
                     bool output_as_coefficients = false) override;

    /** Performs BMMP blind-rotation
     *
     * @param lwe_vec vector containing the LWE sample [a_0, a_1, ..., a_{n - 1}, b]
     * @param rlwe_acc_vec RLWE accumulator in COEF format
     */
    void
    BlindRotate(uint64_t *result, const uint64_t *lwe_vec, uint64_t *rlwe_acc_vec, bool output_as_coefficients = false) override;

    [[nodiscard]] std::shared_ptr<RLWEEncryptor> GetEncryptor() const;

    [[nodiscard]] const std::shared_ptr<const OperatorContext<BlindRotator>> GetContext() const override;

    ~BMMPBlindRotator() override = default;

private:

    /** Performs CGGI key generation
     *
     * @param lwe_key (Binary) vector representing the secret LWE key
     * @param rlwe_key RLWE key assumed to be in COEFFICIENT representation
     */
    void KeyGen(const std::vector<uint64_t>& lwe_key, const std::vector<uint64_t>& rlwe_key);

    /** Performs CGGI key generation
     *
     * @param lwe_key (Binary) vector representing the secret LWE key
     * @param rlwe_key RLWE key assumed to be in COEFFICIENT representation
    */
    void KeyGen(const uint64_t* __restrict lwe_key, const uint64_t* __restrict rlwe_key);

    /** Builds the weight polynomials for blind-rotation for the case step_size = 2.
     * these polynomials are given by X^a_0 - 1, X^a_1 - 1, X^(a_0
     *
     * @param poly_buffer
     * @param a_0
     * @param a_1
     */
    void BuildStepPolynomials(uint64_t* __restrict poly_buffer, uint64_t a_0, uint64_t a_1);

    void KeyGenBinary(const uint64_t* __restrict lwe_key, const uint64_t* __restrict rlwe_key);

    void BlindRotateBinary(const uint64_t* __restrict lwe_vec, uint64_t* __restrict rlwe_vec);

    void SetupMonomials();

    std::shared_ptr<const BMMPBlindRotationContext> m_params;

    std::shared_ptr<MathWorker> m_engine;
    std::unique_ptr<MuxOperator> m_mux;
    std::shared_ptr<RLWEEncryptor> m_encryptor;

    AlignedVector m_brk;
    AlignedVector m_monomials;
    AlignedVector m_accumulator;

    bool m_params_set = false;
    bool m_keys_generated = false;

};


#endif //BMMP_BLIND_ROTATOR_H
