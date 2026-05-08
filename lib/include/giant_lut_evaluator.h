//
// Created by leonard on 12/10/24.
//

#ifndef LARGE_FUNCTIONS_GIANT_LUT_EVALUATOR_H
#define LARGE_FUNCTIONS_GIANT_LUT_EVALUATOR_H

#include <vector>
#include "setup.h"
#include "rgsw.h"
#include "digit_decomposer.h"
#include "scheme_switch.h"
#include "lwe_key_switching.h"

struct HugeLUTConfig {
    // to perform the initial bit decomposition
    BitDecomposerParams m_decomp_params;
    // needs to following keys:

    LWEPrivateKey m_huge_lwe_key;

    // to create the selectors
    SchemeSwitchParameters m_switch_params;
    // needs the following keys:
    // first the lwe key for blind-rotation inputs
    LWEPrivateKey m_scheme_switch_lwe_key;
    // then the output key
    NativePoly m_scheme_switch_rlwe_key;

    // we also need to switch from m_huge_lwe_key to m_scheme_switch_lwe_key once we have the bits
    // members used: {Q_ks, baseKs, std}, everything else can set to 0
    std::shared_ptr<LWECryptoParams> m_ksk_params;
};

class GiantLutEvaluator {

    // Pipeline for params
    // - Digit decomp outputs digits with (N, q, t)
    // - Bit Decomposition: in_params (N,q,t) intermediate_params (N, Q, 2) -> key & modulus switch to (n, q, 2) , n << N
    // - RGSW transform: in_params (n, q, 2) out_params (N, Q, 2, L)
    // need to make sure that N | LUT_N ideally

    // Option 1: for constructor, give the different engines etc as params
    // Option 2: for constructor, give all parameters and let it build

public:

    GiantLutEvaluator(HugeLUTConfig& cfg);

    /**
     *
     * @param ct_idx LWE sample encoding the record index
     * @param records Database encoded into polynomials
     * @param basebits L, L = 2^{basebits}
     * @param digits digits = ceil( log_L(ring modulus) )
     * @return RLWE sample containing the record
     */
    RLWECiphertext EvalLUTPolyMulRaw(lbcrypto::ConstLWECiphertext &ct_idx, std::vector<uint64_t>& records, uint32_t basebits, uint32_t digits);

    RLWECiphertext QueryCRT(lbcrypto::ConstLWECiphertext& ct_idx, std::vector<uint32_t>& records, uint32_t basebits, uint32_t digits, uint32_t P, uint32_t Q);
    /**
     *
     * @param database The database stored as 2d tensor / matrix
     * @param basebits log2(basis), for the basis that will be used for the RLWE' construction
     * @param digits number of digits induced by the basis
     * @return a vector / continuous chunk of memory containing the NTTed polynomials that encode the database
     */
    std::vector<uint64_t> SetupRecordsRaw(std::vector<std::vector<uint32_t>>& database, uint32_t basebits, uint32_t digits);

    std::vector<uint32_t> SetupRecordsCRT(std::vector<std::vector<uint32_t>>& database, uint32_t basebits, uint32_t digits, uint32_t P, uint32_t Q);
    /**
     *
     * @param cts The RLWE' samples
     * @return the same rlwe' samples but stored in continuous memory
     */
    std::vector<uint64_t> SetupRLWEPrimeRaw(std::vector<std::vector<RLWECiphertext>>& cts);

    std::vector<uint32_t> SetupRLWEPrimeCRT(std::vector<std::vector<RLWECiphertext>>& cts, uint64_t P, uint64_t Q);

    /**
     *
     * @param rlwe_prime The RLWE' samples stored in the memory allocated by \SetupRLWEPrimeRaw
     * @param records The database as given through \SetupRecordsRaw
     * @param digits The number of digits for the RLWE'
     * @param columns The database is stored in a block of memory of size (2 * N * digits) * columns, i.e. columns
     * determines the number of RLWE sample we'll have after phase 1
     * @return The resulting RLWE samples output by phase 1
     */
    std::vector<RLWECiphertext>
    PerformPhase1Raw(std::vector<uint64_t> &rlwe_prime, std::vector<uint64_t> &records, uint32_t digits,
                     uint32_t columns);

    /**
     *
    * @param rlwe_prime The RLWE' samples stored in the memory allocated by \SetupRLWEPrimeRaw
     * @param records The database as given through \SetupRecordsRaw
     * @param digits The number of digits for the RLWE'
     * @param columns The database is stored in a block of memory of size (2 * N * digits) * columns, i.e. columns
     * determines the number of RLWE sample we'll have after phase 1
     * @param P 1st prime number forming the CRT modulus
     * @param Q 2nd prime number forming the CRT modulus
     * @return The resulting RLWE samples output by phase 1
     */
    std::vector<RLWECiphertext>
    PerformPhase1CRT(std::vector<uint32_t>& rlwe_prime, std::vector<uint32_t>& records, uint32_t digits, uint32_t columns, uint64_t P, uint64_t Q);

    /**
     * Performs phase 2 from respire
     * @param current_choices list of initial RLWE ciphertexts to be selected from
     * @param selectors RGSW sample that will be used to compute the mux gates via. RGSW(bit) * (choice_1 - choice_0) + choice_0
     * @param first_idx To avoid copying, we only look at a slice of the vector of RGSW samples, induced by \first_idx and \last_idx
     * @param last_idx
     * @return The final RLWE sample to be used in Phase 3
     */
    static RLWECiphertext PerformPhase2(std::vector<RLWECiphertext>& current_choices, std::vector<RingGSWSample>& selectors, uint32_t first_idx, uint32_t last_idx);

    /**
     *
     * @param selectors The first log2(N) RGSW samples we construct and will use for phase 1
     * @param log2N log2(N)
     * @param target_polys the basis polynomials used for RLWE' prime construction. In our case it's a set (p_i) 0 < i < digits s.t. p_i = 1 * basis^i
     * @return returns the set of RLWE' selectors
     */
    std::vector<std::vector<RLWECiphertext>> BuildRLWEPrimeSelectors(std::vector<RingGSWSample>& selectors, uint32_t log2N, std::vector<NativePoly>& target_polys);

    /**
     * Transforms LWE samples of bits into RGSW samples of bits
     * @param bits the vector of LWE samples
     * @return RGSW(b_i) for LWE(b_i) in \bits
     */
    std::vector<RingGSWSample> BuildSelectors(std::vector<LWECiphertext> &bits);

    BitDecomposer m_bitdecomp_engine;
    SchemeSwitchEngine m_switch_engine;

    LWEKeySwitchingKey m_lwe_ksk;
    std::vector<NativePoly> m_rotation_polys;

};

#endif //LARGE_FUNCTIONS_GIANT_LUT_EVALUATOR_H
