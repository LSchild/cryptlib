//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_AUTOMORPHISM_KEY_H
#define LARGE_FUNCTIONS_AUTOMORPHISM_KEY_H

#include "glwe_conversion.h"
#include "rlwe_key_switching.h"
#include "setup.h"

struct AutomorphismParameters : public RLWEConversionParameters {

    AutomorphismParameters(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N, uint64_t basis, uint64_t digits, double std, uint32_t automorphism_index);

    AutomorphismParameters(KeyDistribution source_key_distribution, std::shared_ptr<intel::hexl::NTT>, uint64_t basis, uint64_t digits, double std, uint32_t automorphism_index);

    AutomorphismParameters(const AutomorphismParameters &other);

    void SetAutomorphismIndex(uint32_t idx);

    [[nodiscard]] uint32_t GetAutomorphismIndex() const;

private:

    uint32_t m_automorphism_index;
};

struct AutomorphismEvaluator : protected RLWEtoRLWEConverter {

    AutomorphismEvaluator(AutomorphismParameters& params);

    void SetParams(AutomorphismParameters params);

    void KeyGen(const std::vector<uint64_t>& key);

    void KeyGen(const uint64_t* const key);

    // input assumed in [a=COEF|b=COEF] form
    void Eval(uint64_t *output, const uint64_t *const input);

    // input assumed in [a=COEF|b=COEF] form
    void Eval(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);

    [[nodiscard]] const AutomorphismParameters& GetParams() const override;

    void ApplyAutomorphism(uint64_t* output, const uint64_t *const input);

    void ApplyAutomorphism(std::vector<uint64_t>& output, const std::vector<uint64_t>& input);

    AutomorphismParameters m_automorphism_params;
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


#endif //LARGE_FUNCTIONS_AUTOMORPHISM_KEY_H
