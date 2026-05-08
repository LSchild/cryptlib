//
// Created by leonard on 3/10/25.
//

#ifndef PIROUETTE_H
#define PIROUETTE_H

#include <cstdint>
#include <vector>

#include "digit_decomposer.h"
#include "scheme_switch.h"
#include "setup.h"

enum PIROUETTE_MODE {
    CRT,
    PRIME
};

struct PirouetteParams {
    // to perform the initial bit decomposition
    DigitDecompositionParameters m_decomp_params;
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

    uint32_t m_CRT_P, m_CRT_Q;
    uint32_t m_rlwe_prime_basis_bits, m_rlwe_prime_digits;

    uint64_t m_scaling_params;
};

class PirouettePIR {

public:

    explicit PirouettePIR(PIROUETTE_MODE& mode);

    void SetCryptoParams(PirouetteParams& params);

    void SetMetaParams(uint32_t mu_1, uint32_t mu_2, uint32_t mu_3, uint32_t total_bitwidth);

    void PrepareDB(uint32_t entries_bits, uint32_t record_size_bits, std::vector<std::vector<uint32_t>>& database);

    RLWECiphertext Query(LWECiphertext &record_index);

    RLWECiphertext Query(std::pair<uint64_t, uint64_t>& compressed_record);

    std::pair<uint64_t, uint64_t> CreateQuery(uint64_t idx, uint64_t seed);

    static std::unique_ptr<PirouettePIR> GenerateDefaultConfig(PIROUETTE_MODE);

    NativePoly GetOutputKey() const {
        return m_params.m_scheme_switch_rlwe_key;
    }

    RLWECiphertext QueryCRT(LWECiphertext &record_index);
    RLWECiphertext QueryPrime(LWECiphertext &record_index);

    std::vector<std::vector<RLWECiphertext>> ConstructKBitSelector(std::vector<RingGSWSample>& selectors, std::vector<NativePoly>& target_polys);
    /**
       *
       * @param cts The RLWE' samples
       * @return the same rlwe' samples but stored in continuous memory
       */
    void LayoutRLWEPrimeRaw(std::vector<std::vector<RLWECiphertext>>& cts);
    void LayoutRLWEPrimeCRT(std::vector<std::vector<RLWECiphertext>>& cts);

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
    PerformPhase1Prime(std::vector<uint64_t, intel::hexl::AlignedAllocator<uint64_t, 32>> &rlwe_prime);

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
    PerformPhase1CRT(std::vector<uint32_t, intel::hexl::AlignedAllocator<uint32_t, 32>>& rlwe_prime);

    /**
     * Performs phase 2 from respire
     * @param current_choices list of initial RLWE ciphertexts to be selected from
     * @param selectors RGSW sample that will be used to compute the mux gates via. RGSW(bit) * (choice_1 - choice_0) + choice_0
     * @param first_idx To avoid copying, we only look at a slice of the vector of RGSW samples, induced by \first_idx and \last_idx
     * @param last_idx
     * @return The final RLWE sample to be used in Phase 3
     */
    static RLWECiphertext PerformPhase2(std::vector<RLWECiphertext>& current_choices, std::vector<RingGSWSample>& selectors, uint32_t first_idx, uint32_t last_idx);

    void PrepareDBCRT(uint32_t entries_bits, uint32_t record_size_bits, std::vector<std::vector<uint32_t>>& database);
    void PrepareDBPrime(uint32_t entries_bits, uint32_t record_size_bits, std::vector<std::vector<uint32_t>>& database);


    /**
     * Transforms LWE samples of bits into RGSW samples of bits
     * @param bits the vector of LWE samples
     * @return RGSW(b_i) for LWE(b_i) in \bits
     */
    std::vector<RingGSWSample> BuildSelectors(std::vector<LWECiphertext> &bits);


    PIROUETTE_MODE m_mode;
    PirouetteParams m_params;

    bool m_crypto_params_set = false, m_meta_params_set = false, m_db_set = false;

    uint32_t m_mu_1, m_mu_2, m_mu_3;

    DigitDecomposer m_digit_decomposer;
    SchemeSwitchEngine m_switch_engine;

    LWEKeySwitchingKey m_ksk;

    uint64_t m_db_entries_bits, m_db_record_bits;
    uint32_t m_total_bitwidth;


    std::vector<uint64_t, intel::hexl::AlignedAllocator<uint64_t, 32>> m_current_db_prime;
    std::vector<uint32_t, intel::hexl::AlignedAllocator<uint32_t, 32>> m_current_db_crt;

    std::vector<uint64_t, intel::hexl::AlignedAllocator<uint64_t, 32>> m_rlwe_prime_buffer;
    std::vector<uint32_t, intel::hexl::AlignedAllocator<uint32_t, 32>> m_rlwe_prime_buffer_crt;
};

#endif //PIROUETTE_H
