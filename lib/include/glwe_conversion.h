//
// Created by leonard on 4/22/26.
//

#ifndef LARGE_FUNCTIONS_GLWE_CONVERSION_H
#define LARGE_FUNCTIONS_GLWE_CONVERSION_H

#include "common_types.h"

struct LWEConversionParameters : OperationParameters {

    LWEConversionParameters(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t source_dimension, uint64_t target_dimension, uint64_t basis, uint64_t digits, double std);

    LWEConversionParameters(LWEConversionParameters&);

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance =0.0) const override;

    [[nodiscard]] KeyDistribution GetSourceKeyDistribution() const;

    [[nodiscard]] uint64_t GetModulus() const;

    [[nodiscard]] uint64_t GetSourceDimension() const;

    [[nodiscard]] uint64_t GetTargetDimension() const;

    [[nodiscard]] uint64_t GetGadgetBasis() const;

    [[nodiscard]] uint64_t GetGadgetBasisLog2() const;

    [[nodiscard]] uint64_t GetGadgetDigits() const;

    [[nodiscard]] double GetStd() const;

    void SetSourceKeyDistribution(KeyDistribution distribution);

    void SetModulus(uint64_t modulus);

    void SetSourceDimension(uint64_t input_dimension) ;

    void SetTargetDimension(uint64_t output_dimension) ;

    void SetGadgetBasis(uint64_t basis) ;

    void SetGadgetDigits(uint64_t digits) ;

    void SetStd(double std) ;

    KeyDistribution m_source_distribution;
    uint64_t m_modulus;
    uint64_t m_source_dimension;
    uint64_t m_target_dimension;
    uint64_t m_basis;
    uint64_t m_basis_log2;
    uint64_t m_digits;
    double m_std;

};

struct LWEtoLWEConverter : SchemeConverter<LWEConversionParameters> {

    LWEtoLWEConverter(LWEConversionParameters& params);

    void SetParams(LWEConversionParameters &params) override;

    void KeyGen(const std::vector<uint64_t> &source_key, const std::vector<uint64_t> &target_key) override;

    void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) override;

    void Convert(uint64_t *output, const uint64_t *const input) override;

    void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) override;

    [[nodiscard]] Container GetInputContainer() const override;

    [[nodiscard]] Container GetOutputContainer() const override;

    [[nodiscard]] const LWEConversionParameters& GetParams() const override;

    AlignedVector m_ksk;
    AlignedVector m_acc;
    LWEConversionParameters m_params;

    bool m_params_set = false;
    bool m_keys_generated = false;

};


struct RLWEConversionParameters : OperationParameters {

    RLWEConversionParameters(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N, uint64_t basis, uint64_t digits, double std);

    RLWEConversionParameters(KeyDistribution source_key_distribution, std::shared_ptr<intel::hexl::NTT>, uint64_t basis, uint64_t digits, double std);

    RLWEConversionParameters(const RLWEConversionParameters &other);

    long double ComputeOutputVariance(long double input_variance) const override;

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

    KeyDistribution m_source_distribution;

    uint64_t m_modulus;
    uint64_t m_N;
    uint64_t m_basis;
    uint64_t m_basis_log2;
    uint64_t m_digits;
    double m_std;

    std::shared_ptr<intel::hexl::NTT> m_ntt;

};

struct RLWEtoRLWEConverter : SchemeConverter<RLWEConversionParameters> {

    RLWEtoRLWEConverter(const RLWEConversionParameters &params);

    void SetParams(RLWEConversionParameters &params) override;

    void KeyGen(const std::vector<uint64_t> &source_key, const std::vector<uint64_t> &target_key) override;

    void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) override;

    // input assumed in [a=COEF|b=NTT] form
    void Convert(uint64_t *output, const uint64_t *const input) override;

    // input assumed in [a=COEF|b=NTT] form
    void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) override;

    [[nodiscard]] Container GetInputContainer() const override;

    [[nodiscard]] Container GetOutputContainer() const override;

    [[nodiscard]] const RLWEConversionParameters& GetParams() const override;

    AlignedVector m_ksk;
    AlignedVector m_acc;
    RLWEConversionParameters m_params;

    bool m_params_set = false;
    bool m_keys_generated = false;

};



#endif //LARGE_FUNCTIONS_GLWE_CONVERSION_H
