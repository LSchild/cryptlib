//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_SCHEME_SWITCH_H
#define LARGE_FUNCTIONS_SCHEME_SWITCH_H

#include "setup.h"
#include "rlwe_key_switching.h"
#include "automorphism_key.h"
#include "blind_rotator.h"

/**
 * The scheme switch key (SSK) is a special case of a rlwe switching key
 * Given a RLWE(m), it enables us to compute RLWE(-s * m) where s is the secret key
 *
 * More specifically, the SSK is ksk from s^2 to s. Then, given
 * RLWE(m) = (a, b = a * s + m + e), we start by performing a rlwe key switch on (a, 0)
 * obtaining (a', a' * s + a * s^2 + e') and add b to the left hand side, getting
 * (a' + b, a' * s + a * s^2 + e') = (a' + a * s + m + e, a' * s - a * s^2 + e')
 * We note by computing the phase of the last equation
 * phase = a' * s + a * s^2 + e' - a' * s - a * s * s - m * s - e * s
 *       = -s * m - e * s - e'
 * i.e. it is a valid RLWE sample of -s * m
 *
 * By using sparse secrets, we will have that || e * s||_\inf will be rather small
 */
class SchemeSwitchKey : public RLWEKeySwitchingKey {


public:

    explicit SchemeSwitchKey() = default;

    SchemeSwitchKey(std::shared_ptr<RingGSWCryptoParams>& params, NativePoly& sk);

    RLWECiphertext SchemeSwitch(const RLWECiphertext& source);

    RLWECiphertext SchemeSwitch(const NativePoly& A, const NativePoly& B);

};


class SchemeSwitchEngine2 {

public:

    SchemeSwitchEngine2(std::shared_ptr<RingGSWCryptoParams> &ek_params, std::shared_ptr<RingGSWCryptoParams> &
    switch_params, lbcrypto::NativePoly &sk, ConstLWEPrivateKey &sk_lwe);

    RingGSWSample Switch(LWECiphertext &ct);

    RLWECiphertext DecimateAndGetSlot(RLWECiphertext& in, uint32_t slot_idx);

    BlindRotationKey m_ek;

    /* parameters for blind-rotation */
    std::shared_ptr<RingGSWCryptoParams> m_params;
    /* parameters for scheme switch e.g. RGSW basis, or basis for automorphism keys */
    std::shared_ptr<RingGSWCryptoParams> m_switch_params;

    /* squaring key with which we can do RLWE(m) -> RLWE(-sk * m) */
    SchemeSwitchKey m_squaring_key;
    /* automorphism key with which we can do RLWE(m_0 + X m_1 + ...) to RLWE(m_i)_{i in ...} */
    std::vector<AutomorphismKey> m_automorphism_key;

    /* debugging */
    NativePoly m_sk;
    NativePoly m_sk_NTT;
};

struct SchemeSwitchParameters {
    /* LWE dimension */
    uint32_t m_n;
    /* RLWE dimension */
    uint32_t m_N;
    /* LWE modulus */
    uint64_t m_q;
    /* RLWE modulus */
    uint64_t m_Q;

    /* standard deviation */
    double m_std;

    /* Number of least significant digits to be skipped during scheme switching */
    uint32_t skip_first_digits;

    /** bases **/
    /* for blind-rotation key */
    uint32_t m_br_basis;
    /* for squaring key */
    uint32_t m_squaring_key_basis;
    /* for automorphism key */
    uint32_t m_automorphism_key_basis;
    /* number of digits in the output RGSW sample */
    uint32_t m_rgsw_basis;
};

class SchemeSwitchEngine {

public:

    SchemeSwitchEngine() {};

    explicit SchemeSwitchEngine(SchemeSwitchParameters& params);

    void SetKeys(NativePoly& sk_rlwe, LWEPrivateKey& sk_lwe);

    /**
     * Given an LWE sample encoding a bit i.e. ct = LWE(q / 2 * bit), returns a RGSW(bit)
     * @param ct = LWE(q / 2 * bit)
     * @return RGSW(bit)
     */
    RingGSWSample SwitchToRGSW(ConstLWECiphertext& ct);



//private:

    /**
    * Given a LWE ciphertext LWE(q/2 * bit) returns RLWE(bit * X^e * pad_poly * \sum_i L_i) where
    * pad_poly corresponds to a sum of subsequent monomials = 1 + X + X^2 ...and the highest degree depends on the number of output digits
    * @param ct the LWE ciphertext
    * @return the described RLWE sample
    */
    RLWECiphertext CreateDigits(ConstLWECiphertext& ct);

    /**
     * Given the output of CreateDigits, outputs a vector of RLWE ciphertexts ct_i, where ct_i = RLWE(bit * L^(i + skip_digits_val))
     * @param ct_digits the ciphertext as described in CreateDigits
     * @return The vector
     */
    std::vector<RLWECiphertext> ExtractDigits(RLWECiphertext& ct_digits);

    /**
     * Given a RLWE sample of a polynomial p, outputs a RLWE sample containing p_0, i.e. the constant term of p
     * @param ct the RLWE ciphertext
     * @return a rlwe ciphertext of the constant term
     */
    RLWECiphertext GetConstantTerm(const RLWECiphertext& ct);


    /* Debugging */
    NativePoly m_sk_rlwe;
    NativePoly m_sk_rlwe_ntt;
    NativeVector m_sk_lwe;

    /* Parameters */
    SchemeSwitchParameters m_params;

    /* internal parameter set wrappers */
    std::shared_ptr<RingGSWCryptoParams> m_br_params;
    std::shared_ptr<RingGSWCryptoParams> m_squaring_params;
    std::shared_ptr<RingGSWCryptoParams> m_automorphism_params;
    std::shared_ptr<RingGSWCryptoParams> m_rgsw_params;


    /* Required keys */
    SchemeSwitchKey m_squaring_key;
    FastBlindRotationKey m_br_key;
    //BlindRotationKey m_br_key;
    std::vector<AutomorphismKey> m_automorphism_keys;

};


#endif //LARGE_FUNCTIONS_SCHEME_SWITCH_H
