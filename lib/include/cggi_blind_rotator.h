//
// Created by leonard on 3/23/26.
//

#ifndef CGGI_BLIND_ROTATOR_H
#define CGGI_BLIND_ROTATOR_H

#include <cstdint>
#include <memory>
#include <vector>

#include "hexl/hexl.hpp"

#include "blindrotation_common.h"
#include "interfaces.h"
#include "mux_operator.h"
#include "base_crypto.h"

// forward decl.
struct CGGIBlindRotator;

struct CGGIBlindRotationContext : OperatorContext<CGGIBlindRotator, BlindRotationKeys>
, public std::enable_shared_from_this<CGGIBlindRotationContext> {

    friend struct CGGIBlindRotator;

    explicit CGGIBlindRotationContext(KeyDistribution distr, uint64_t modulus, uint64_t ring_dim, uint64_t lwe_dim, uint64_t m_basis, uint64_t m_digits, double std);

    explicit CGGIBlindRotationContext(KeyDistribution distr, std::shared_ptr<intel::hexl::NTT> ntt, uint64_t lwe_dim, uint64_t m_basis, uint64_t m_digits, double std);

    CGGIBlindRotationContext(const CGGIBlindRotationContext &other);

    void SetKeyDistribution(KeyDistribution dist);

    void SetModulus(uint64_t modulus);

    void SetRingDimension(uint64_t ring_dim);

    void SetNTT(std::shared_ptr<intel::hexl::NTT> ntt);

    void SetLWEDimension(uint64_t lwe_dim);

    void SetBlindRotationBasis(uint64_t L);

    void SetBlindRotationRGSWDigits(uint64_t digits);

    void SetStd(double std);

    [[nodiscard]] KeyDistribution GetKeyDistribution() const;

    [[nodiscard]] uint64_t GetModulus() const;

    [[nodiscard]] uint64_t GetRingDimension() const;

    [[nodiscard]] std::shared_ptr<intel::hexl::NTT> GetNTT() const;

    [[nodiscard]] uint64_t GetLWEDimension() const;

    [[nodiscard]] uint64_t GetBlindRotationBasis() const;

    [[nodiscard]] uint64_t GetBlindRotationBasisLog2() const;

    [[nodiscard]] uint64_t GetBlindRotationRGSWDigits() const;

    [[nodiscard]] double GetStd() const;

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    [[nodiscard]] std::shared_ptr<CGGIBlindRotator> ConstructOperator(BlindRotationKeys& bundle) const override;

    [[nodiscard]] Container GetInputContainer() const override;

    [[nodiscard]]  Container GetOutputContainer() const override;

    [[nodiscard]] OperatorID GetOperatorID() const override;

protected:

    KeyDistribution m_distribution;

    uint64_t m_modulus;
    uint64_t m_N;
    uint64_t m_n;

    uint64_t m_basis;
    uint64_t m_basis_log2;
    uint64_t m_digits;

    double m_std;

    std::shared_ptr<intel::hexl::NTT> m_ntt;

};

struct CGGIBlindRotator : public BlindRotator {

    friend struct CGGIBlindRotationContext;


    CGGIBlindRotator(const std::shared_ptr<const CGGIBlindRotationContext>& params);

    /** Performs CGGI blind-rotation
     *
     * @param lwe_vec vector containing the LWE sample [a_0, a_1, ..., a_{n - 1}, b]
     * @param rlwe_acc_vec RLWE accumulator in NTT format
     */
    void BlindRotate(const std::vector<uint64_t>& lwe_vec, std::vector<uint64_t>& rlwe_acc_vec) override;

    /** Performs CGGI blind-rotation
     *
     * @param lwe_vec vector containing the LWE sample [a_0, a_1, ..., a_{n - 1}, b]
     * @param rlwe_acc_vec RLWE accumulator in NTT format
     */
    void BlindRotate(const uint64_t* __restrict lwe_vec, uint64_t* __restrict rlwe_acc_vec) override;

    [[nodiscard]] std::shared_ptr<RLWEEncryptor> GetEncryptor() const;

private:




    /** Performs CGGI key generation
     *
     * @param lwe_key (Binary/Ternary) vector representing the secret LWE key
     * @param rlwe_key RLWE key assumed to be in COEFFICIENT representation
     */
    void KeyGen(const std::vector<uint64_t>& lwe_key, const std::vector<uint64_t>& rlwe_key);

    /** Performs CGGI key generation
     *
     * @param lwe_key (Binary/Ternary) vector representing the secret LWE key
     * @param rlwe_key RLWE key assumed to be in COEFFICIENT representation
    */
    void KeyGen(const uint64_t* __restrict lwe_key, const uint64_t* __restrict rlwe_key);

    void KeyGenBinary(const uint64_t* __restrict lwe_key, const uint64_t* __restrict rlwe_key);

    void KeyGenTernary(const uint64_t* __restrict lwe_key, const uint64_t* __restrict rlwe_key);

    void BlindRotateBinary(const uint64_t* __restrict lwe_vec, uint64_t* __restrict rlwe_vec);

    void BlindRotateTernary(const uint64_t* __restrict lwe_vec, uint64_t* __restrict rlwe_vec);

    void SetupMonomials();

    std::shared_ptr<const CGGIBlindRotationContext> m_params;

    std::shared_ptr<intel::hexl::NTT> m_engine;
    std::unique_ptr<MuxOperator> m_mux;
    std::shared_ptr<RLWEEncryptor> m_encryptor;

    AlignedVector m_brk;
    AlignedVector m_monomials;
    AlignedVector m_accumulator;

    bool m_params_set = false;
    bool m_keys_generated = false;

};

#endif //CGGI_BLIND_ROTATOR_H
