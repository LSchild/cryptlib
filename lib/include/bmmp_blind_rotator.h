//
// Created by leonard on 3/23/26.
//

#ifndef BMMP_BLIND_ROTATOR_H
#define BMMP_BLIND_ROTATOR_H
#include "cggi_blind_rotator.h"

struct BMMPBlindRotatorParams : CGGIBlindRotatorParams {

    BMMPBlindRotatorParams(KeyDistribution distr, uint64_t modulus, uint64_t ring_dim, uint64_t lwe_dim, uint64_t m_basis, uint64_t m_digits, double std, uint64_t step_size);

    BMMPBlindRotatorParams(KeyDistribution distr, std::shared_ptr<intel::hexl::NTT>, uint64_t lwe_dim, uint64_t m_basis, uint64_t m_digits, double std, uint64_t step_size);

    BMMPBlindRotatorParams(const BMMPBlindRotatorParams &other);

    long double ComputeOutputVariance(long double input_variance = 0) const override;

    uint64_t GetStepSize() const;

protected:

    uint64_t m_step_size;

};

struct BMMPBlindRotator : public BlindRotator<BMMPBlindRotatorParams> {
    /** Builds a CGGI blind-rotation object for binary or ternary LWE keys.
     *
     * @param params Struct containing the parameters (modulus, ring dimension, ...)
     */
    explicit BMMPBlindRotator(const BMMPBlindRotatorParams &params);

    void SetParams(BMMPBlindRotatorParams& params) override;

    /** Performs CGGI key generation
     *
     * @param lwe_key (Binary) vector representing the secret LWE key
     * @param rlwe_key RLWE key assumed to be in COEFFICIENT representation
     */
    void KeyGen(const std::vector<uint64_t>& lwe_key, const std::vector<uint64_t>& rlwe_key) override;

    /** Performs CGGI key generation
     *
     * @param lwe_key (Binary) vector representing the secret LWE key
     * @param rlwe_key RLWE key assumed to be in COEFFICIENT representation
    */
    void KeyGen(const uint64_t* __restrict lwe_key, const uint64_t* __restrict rlwe_key) override;

    /** Performs BMMP blind-rotation
     *
     * @param lwe_vec vector containing the LWE sample [a_0, a_1, ..., a_{n - 1}, b]
     * @param rlwe_acc_vec RLWE accumulator in NTT format
     */
    void BlindRotate(const std::vector<uint64_t>& lwe_vec, std::vector<uint64_t>& rlwe_acc_vec) override;

    /** Performs BMMP blind-rotation
     *
     * @param lwe_vec vector containing the LWE sample [a_0, a_1, ..., a_{n - 1}, b]
     * @param rlwe_acc_vec RLWE accumulator in NTT format
     */
    void BlindRotate(const uint64_t* __restrict lwe_vec, uint64_t* __restrict rlwe_acc_vec) override;

    BlindRotationMethod GetMethod() override;

    [[nodiscard]] const BMMPBlindRotatorParams& GetParams() const override;

    [[nodiscard]] std::shared_ptr<RLWEEncryptor> GetEncryptor() const;

    [[nodiscard]] Container GetInputContainer() const override;

    [[nodiscard]] Container GetOutputContainer() const override;

private:

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

    BMMPBlindRotatorParams m_params;

    std::shared_ptr<intel::hexl::NTT> m_engine;
    std::unique_ptr<MuxOperator> m_mux;
    std::shared_ptr<RLWEEncryptor> m_encryptor;

    AlignedVector m_brk;
    AlignedVector m_monomials;
    AlignedVector m_accumulator;

    bool m_params_set = false;
    bool m_keys_generated = false;

};


#endif //BMMP_BLIND_ROTATOR_H
