//
// Created by leonard on 3/24/26.
//

#ifndef BASE_CRYPTO_H
#define BASE_CRYPTO_H
#include <cstdint>

#include "backend/backend.h"

struct RLWEEncryptor {
    /** Creates and Encryption for RLWE over the ring Z_modulus/(X^ring_dimension + 1)
     *
     * @param modulus ring Q
     * @param ring_dimension ring dimension N
     * @param std error standard deviation
     */
    RLWEEncryptor(uint64_t modulus, uint32_t ring_dimension, double std);

    /** Creates and Encryption for RLWE over the ring Z_modulus/(X^ring_dimension + 1) with specific
     * root of unity.
     *
     * @param modulus ring Q
     * @param ring_dimension ring dimension N
     * @param root_of_unity ring root of unity to enforce compat
     * @param std error standard deviation
     */
    RLWEEncryptor(uint64_t modulus, uint32_t ring_dimension, uint64_t root_of_unity, double std);

    /** Creates and Encryption for RLWE over the ring Z_modulus/(X^ring_dimension + 1)
     *
     * @param ntt hexl ntt engine
     * @param std standard devation
     */
    RLWEEncryptor(std::shared_ptr<MathWorker> ntt, double std);

    /** Creates RLWE(msg)
     *
     * @param result pointer to array of size 2 * N
     * @param msg pointer to polynomial of size N
     * @param secret_ntt pointer to secret polynomial in NTT form, size = N
     * @param msg_is_ntt toggle whether messages is in NTT form or not
     */
    void MakeRLWE(uint64_t* result, uint64_t* msg, uint64_t* secret_ntt, bool msg_is_ntt = false);

    /** Creates RGSW(msg)
     *
     * @param result pointer to array of size rgsw_digits * 4 * N
     * @param msg pointer to message polynomial
     * @param secret_ntt pointer to secret polynomial in NTT form, size = N
     * @param rgsw_basis rgsw basis L for rows of the RGSW sample
     * @param rgsw_digits number of digits for rgsw
     * @param msg_is_ntt toggle whether messages is in NTT form or not
     */
    void MakeRGSW(uint64_t* result, uint64_t* msg, uint64_t* secret_ntt, uint64_t rgsw_basis, uint64_t rgsw_digits, bool msg_is_ntt = false);

    /** Returns msg in COEFFICIENT format, from RLWE(msg) and the secret
     *
     * @param result pointer to output buffer of size N
     * @param rlwe pointer to input RLWE of size 2 * N
     * @param secret_ntt pointer to secret in polynomial of size N
     */
    void PhaseRLWE(uint64_t* result, uint64_t* rlwe, uint64_t* secret_ntt);

    /** Returns msg in COEFFICIENT format, from RGSW(msg) and the secret
    *
    * @param result pointer to output buffer of size N
    * @param rgsw pointer to input RLWE of size 2 * N
    * @param secret_ntt pointer to secret in polynomial of size N
    */
    void PhaseRGSW(uint64_t* result, uint64_t* rgsw, uint64_t* secret_ntt, uint64_t rgsw_digits);

    uint64_t GetModulus();

    uint64_t GetDimension();

    double GetStd();

    std::shared_ptr<MathWorker> GetNTT();

private:

    double m_std;

    std::shared_ptr<MathWorker> m_ntt;
};


struct LWEEncryptor {

    LWEEncryptor(uint64_t modulus, uint64_t lwe_dimension, double std);

    void MakeLWE(uint64_t* result, uint64_t msg, const uint64_t *const secret) const;

    void PhaseLWE(uint64_t* result, uint64_t* lwe_vec, const uint64_t *const secret) const;

    uint64_t PhaseLWE(uint64_t* lwe_vec, const uint64_t *const secret) const;

    uint64_t GetModulus();

    uint64_t GetDimension();

    double GetStd();

private:

    uint64_t m_modulus;
    uint64_t m_dimension;
    double m_std;

};

#endif //BASE_CRYPTO_H
