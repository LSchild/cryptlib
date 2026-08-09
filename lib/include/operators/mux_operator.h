//
// Created by leonard on 3/23/26.
//

#ifndef MUX_OPERATOR_H
#define MUX_OPERATOR_H

#include "backend/backend.h"
#include "backend/aligned_vector.h"

struct MuxOperator {


    MuxOperator(std::shared_ptr<MathWorker> ntt, uint64_t basis_log2, uint64_t digits);

    /** Computes the usual external product result (+)= rgsw * rhs
     *
     *
     * @param result pointer to RLWE sample coefficients of result
     * @param rgsw (NTT) pointer to RGSW sample that will be multiplied
     * @param rhs (NTT) pointer to input RLWE sample
     * @param add_to_result if true compute result += rgsw * rhs, else result = rgsw * rhs
     */
    void ExternalProduct(uint64_t * __restrict result, const uint64_t* __restrict rgsw, const uint64_t* __restrict rhs, bool add_to_result = false);

    /** Computes the usual MUX operation: result = case_0 if rgsw_control == 0 else case_1
     *
     *
     * @param result (NTT) pointer to RLWE sample coefficients of result
     * @param rgsw_control (NTT) pointer to RGSW sample of a control bit B
     * @param case_0 (NTT) pointer to RLWE with value to return if B == 0
     * @param case_1 (NTT) pointer to RLWE with value to return if B == 1
     */
    void BinaryMux(uint64_t* __restrict result, const uint64_t* __restrict rgsw_control, uint64_t *__restrict case_0, uint64_t *__restrict case_1);

    /** Computes the Multi-CMUX operation rlwe_sample += sum_i Extprod(rgsw_control_bits[i], rhs) * weights[i]
     * If weights == nullptr, skip the final product
     *
     *
     * @param rlwe_sample pointer to RLWE sample coefficients of input and rlwe_sample
     * @param k number of RGSW control bits
     * @param rgsw_control_bits (NTT) start pointer to beginning of RGSW control bit array
     * @param weights (NTT) pointer to vector of polynomials used for post-external-prod multiplication
     * @param digit_buffer pointer to buffer to store digits, if not provided will be allocated
     */
    void MultiMux(uint64_t* __restrict rlwe_sample, uint64_t k, const uint64_t* __restrict rgsw_control_bits, const uint64_t*__restrict weights);

    /** Computes the usual CMUX operation: acc += rgsw_control * (X^a - 1) * acc
     *
     *
     * @param acc (NTT) pointer to input/output RLWE sample
     * @param rgsw_control (NTT)  pointer to rgsw control bit
     * @param monomial (NTT) pointer to NTTed monomial-minus-one X^a - 1
     */
    void BinaryCMux(uint64_t* __restrict acc, const uint64_t* __restrict rgsw_control, const uint64_t *__restrict monomial);

    /** Computes the ternary CMUX operation: acc = X^{a_i * s_i} * acc, where s_i in {-1, 0, 1}
     * s_i is encoded as a pair of RGSW sample hi, lo with
     * - hi,lo = RGSW(0), RGSW(0) if s_i = 0
     * - hi,lo = RGSW(1), RGSW(0) if s_i = 1
     * - hi,lo = RGSW(0), RGSW(1) if s_i = -1
     *
     *
     *
     * @param acc (NTT) pointer to input/output RLWE sample
     * @param rgsw_controls (NTT) pointer to pair of RGSW samples
     * @param monomial (NTT) pointer to pair of monomials required
     */
    void TernaryCMux(uint64_t* __restrict acc, const uint64_t* __restrict rgsw_controls, const uint64_t *__restrict monomial);

    /** Computes the RLWE-Prime x polynomial product: acc (+)= rlwe_prime * polynomial
     *
     *
     * @param acc pointer to output RLWE sample
     * @param rlwe_prime (NTT) pointer to RLWE-prime sample
     * @param polynomial (COEF) pointer to  polynomial
     * @param add_to_acc bool to decide whether to update acc or overwrite it
     */
    void RLWEPrimeProduct(uint64_t * __restrict acc, const uint64_t* __restrict rlwe_prime, const uint64_t* __restrict polynomial, bool add_to_acc);

    [[nodiscard]] std::shared_ptr<MathWorker> GetNTT() const;

    [[nodiscard]] uint64_t GetModulusBits() const;

    [[nodiscard]] uint64_t GetRGSWBasis() const;

    [[nodiscard]] uint64_t GetRGSWBasisLog2() const;

    [[nodiscard]] uint64_t GetRGSWDigits() const;

    [[nodiscard]] uint64_t GetMask() const;

    [[nodiscard]] uint64_t GetFirstShift() const;

private:

    AlignedBuffer m_scratch_space;
    std::shared_ptr<MathWorker> m_ntt;

    uint64_t m_modulus_bits;
    uint64_t m_basis;
    uint64_t m_basis_log2;
    uint64_t m_digits;

    uint64_t m_mask;
    uint64_t m_first_shift;
};

#endif //MUX_OPERATOR_H
