//
// Created by leonard on 4/22/26.
//

#ifndef LARGE_FUNCTIONS_ENDO_GLWE_CONVERSION_H
#define LARGE_FUNCTIONS_ENDO_GLWE_CONVERSION_H

#include "common_types.h"
#include "interfaces/enum_ids.h"
#include "interfaces/conversion_operator.h"
#include "backend/backend.h"

struct LWEtoLWEConverter;

struct LWEConversionContext : public OperatorContext<LWEtoLWEConverter>,
                              public std::enable_shared_from_this<LWEConversionContext> {

    /**
     * LWEConversionContext constructor for switch from Z_Q^{n_0 + 1} to Z_Q^{n_1 + 1}
     * @param source_key_distribution Source key distribution as it influences the error growth
     * @param modulus Modulus Q
     * @param source_dimension Dimension of the input LWE n_0
     * @param target_dimension Dimension of the output LWE n_1
     * @param basis Gadget basis
     * @param digits Gadget digits
     * @param std Standard deviation for key
     */
    LWEConversionContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t source_dimension, uint64_t target_dimension, uint64_t basis, uint64_t digits, double std);

    /**
     * Constructor from different instance
     * @param other Different context
     */
    LWEConversionContext(LWEConversionContext& other);

    /**
     * Computes variance of LWE to LWE switching / conversion
     * @param input_variance Variance of the input LWE
     * @return Variance of the output LWE
     */
    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    /**
     * Construct a LWE to LWE converter from given key(s)
     * @param keys Vector containing two keys: the source key S and target key T
     * @return Unique ptr to LWEtoLWEConverter from S to T
     */
    [[nodiscard]] std::unique_ptr<LWEtoLWEConverter> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @return A Container describing the expected input
     */
    [[nodiscard]] Container GetInputContainer() const override;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @param input Input container description. Necessary as some operators are not additive in variance.
     *        if equal to nullptr, is ignored.
     * @return A Container describing the calculated output
     */
    [[nodiscard]] Container GetOutputContainer(Container input) const override;

    /** Returns an ID describing the Operator
     *
     * @return The ID
     */
    [[nodiscard]] OperatorID GetOperatorID() const override;

    /**
     *
     * @return Key distribution
     */
    [[nodiscard]] KeyDistribution GetSourceKeyDistribution() const;
    /**
     *
     * @return Ring modulus of the ring R=Z_Q[X]/(X^N + 1)
     */
    [[nodiscard]] uint64_t GetModulus() const;

    /**
     *
     * @return Dimension of the input LWE
     */
    [[nodiscard]] uint64_t GetSourceDimension() const;

    /**
     *
     * @return Dimension of the output LWE
     */
    [[nodiscard]] uint64_t GetTargetDimension() const;

    /**
     *
     * @return Gadget Basis
     */
    [[nodiscard]] uint64_t GetGadgetBasis() const;

    /**
     *
     * @return Bits / Log2 of the gadget basis
     */
    [[nodiscard]] uint64_t GetGadgetBasisLog2() const;

    /**
     *
     * @return Gadget Digits
     */
    [[nodiscard]] uint64_t GetGadgetDigits() const;

    /**
     *
     * @return Key standard deviation
     */
    [[nodiscard]] double GetStd() const;

    /**
     * Update the source key distribution
     * @param distribution New distribution
     */
    void SetSourceKeyDistribution(KeyDistribution distribution);

    /**
     * Update input and output modulus
     * @param modulus New modulus
     */
    void SetModulus(uint64_t modulus);

    /**
     * Update input LWE dimension
     * @param input_dimension New dimension
     */
    void SetSourceDimension(uint64_t input_dimension) ;

    /**
     * Update target LWE dimension
     * @param output_dimension New dimension
     */
    void SetTargetDimension(uint64_t output_dimension) ;

    /**
     * Update Gadget basis
     * @param basis New basis
     */
    void SetGadgetBasis(uint64_t basis) ;

    /**
     * Update Gadget digits
     * @param digits New digits
     */
    void SetGadgetDigits(uint64_t digits) ;

    /**
     * Update key standard deviation
     * @param std
     */
    void SetStd(double std) ;

private:
    // input key distribution
    KeyDistribution m_source_distribution;
    // modulus Q
    uint64_t m_modulus;
    // input LWE dimension
    uint64_t m_source_dimension;
    // output LWE dimension
    uint64_t m_target_dimension;
    // Gadget basis
    uint64_t m_basis;
    // Log2 of gadget basis
    uint64_t m_basis_log2;
    // gadget digits
    uint64_t m_digits;
    // Key Standard deviation
    double m_std;

};

struct LWEtoLWEConverter : public SchemeConverter<LWEConversionContext> {
    friend class LWEConversionContext;

    /**
     * Performs the LWE to LWE conversion / switching
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer, must differ from input
     */
    void Convert(uint64_t *output, const uint64_t *const input) override;

    /**
     * Performs the LWE to LWE conversion / switching
     * @param output Reference to output vector
     * @param input Reference to input vector
     */
    void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) override;

    /**
     * Getter for parent context
     * @return Parent context
     */
    [[nodiscard]] const std::shared_ptr<const LWEConversionContext> GetContext() const override;

private:

    /**
     * Primary constructor
     * @param params LWEConversion context
     */
    LWEtoLWEConverter(std::shared_ptr<const LWEConversionContext> params);

    /**
     * Function to generate the conversion/switching key
     * @param source_key Key to switch FROM
     * @param target_key Key to switch TO
     */
    void KeyGen(const std::vector<uint64_t> &source_key, const std::vector<uint64_t> &target_key);

    /**
     * Function to generate the conversion/switching key
     * @param source_key Key to switch FROM
     * @param target_key Key to switch TO
     */
    void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key);

    // buffer that stores the key
    AlignedVector m_ksk;
    // additional scratch space
    AlignedVector m_acc;
    // pointer to parent context
    std::shared_ptr<const LWEConversionContext> m_params;

    // debug flags
    bool m_params_set = false;
    bool m_keys_generated = false;

};


struct RLWEtoRLWEConverter;

struct RLWEConversionContext : public OperatorContext<RLWEtoRLWEConverter>,
                               public std::enable_shared_from_this<RLWEConversionContext> {

    /**
     * Constructor for RLWE conversion over the ring Z_Q[X]/(X^N + 1)
     * @param source_key_distribution Distribution of the source RLWE key (Binary/Ternary/Gaussian)
     * @param modulus Modulus Q
     * @param N Dimension N
     * @param basis Gadget basis
     * @param digits Gadget digits
     * @param std Standard deviation for the conversion/switching key
     */
    RLWEConversionContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N, uint64_t basis, uint64_t digits, double std);

    /**
     * Constructor for RLWE conversion over the ring Z_Q[X]/(X^N + 1)
     *
     * Modulus and dimension are obtained from the ntt engine
     *
     * @param source_key_distribution Distribution of the source RLWE key (Binary/Ternary/Gaussian)
     * @param ntt NTT engine
     * @param basis Gadget basis
     * @param digits Gadget digits
     * @param std Standard deviation for the conversion/switching key
     */
    RLWEConversionContext(KeyDistribution source_key_distribution, std::shared_ptr<MathWorker> ntt, uint64_t basis, uint64_t digits, double std);

    /**
     * Constructor from different instance
     * @param other Different context
     */
    RLWEConversionContext(const RLWEConversionContext &other);

    /**
     * Computes the variance of the RLWE to RLWE key conversion
     * @param input_variance Variance of input RLWE sample
     * @return Output variance
     */
    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    /**
     * Constructs the RLWEtoRLWEConverter from RLWE under key S to RLWE under key T
     * @param keys Pointer to vector of keys = {S, T}
     * @return Unique pointer to RLWEtoRLWEConverter
     */
    [[nodiscard]] std::unique_ptr<RLWEtoRLWEConverter> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @return A Container describing the expected input
     */
    [[nodiscard]] Container GetInputContainer() const override;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @param input Input container description. Necessary as some operators are not additive in variance.
     *        if equal to nullptr, is ignored.
     * @return A Container describing the calculated output
     */
    [[nodiscard]] Container GetOutputContainer(Container input) const override;

    /** Returns an ID describing the Operator
     *
     * @return The ID
     */
    [[nodiscard]] OperatorID GetOperatorID() const override;

    /**
     *
     * @return Source key distribution
     */
    [[nodiscard]] KeyDistribution GetSourceKeyDistribution() const;

    /**
     *
     * @return Ring modulus of the ring R=Z_Q[X]/(X^N + 1)
     */
    [[nodiscard]] uint64_t GetModulus() const;

    /**
     *
     * @return Ring dimension of the ring R=Z_Q[X]/(X^N + 1)
     */
    [[nodiscard]] uint64_t GetDimension() const;

    /**
     *
     * @return Gadget Basis L
     */
    [[nodiscard]] uint64_t GetGadgetBasis() const;

    /**
     * Returns the bits/log2 of the gadget basis for efficient gadget computation
     * @return Base 2 logarithm of gadget basis
     */
    [[nodiscard]] uint64_t GetGadgetBasisLog2() const;

    /**
     *
     * @return Gadget Digits
     */
    [[nodiscard]] uint64_t GetGadgetDigits() const;

    /**
     *
     * @return Standard Deviation of the key
     */
    [[nodiscard]] double GetStd() const;

    /**
     * Getter for NTT engine
     * @return Pointer to NTT Engine
     */
    [[nodiscard]] std::shared_ptr<MathWorker> GetNTT() const;

    /**
     * Update declared source key distribution
     * @param distribution New distribution
     */
    void SetSourceKeyDistribution(KeyDistribution distribution);

    /**
     * Update declared modulus of the ring R=Z_Q[X]/(X^N + 1)
     * @param modulus New modulus Q
     */
    void SetModulus(uint64_t modulus);

    /**
     * Update declared ring dimension of the ring R=Z_Q[X]/(X^N + 1)
     * @param input_dimension New dimension N
     */
    void SetDimension(uint64_t input_dimension) ;

    /**
     * Update Gadget Basis
     * @param basis New basis
     */
    void SetGadgetBasis(uint64_t basis) ;

    /**
     * Update Gadget digits
     * @param digits New digits
     */
    void SetGadgetDigits(uint64_t digits) ;

    /**
    * Update key standard deviation
     * @param std new value
     */
    void SetStd(double std) ;

    /**
     * Update ntt/math engine
     * @param ntt engine
     */
    void SetNTT(std::shared_ptr<MathWorker> ntt);

private:
    // distribution of source key, relevant for noise growth
    KeyDistribution m_source_distribution;
    // ring modulus Q
    uint64_t m_modulus;
    // ring dimensio N
    uint64_t m_N;
    // gadget basis
    uint64_t m_basis;
    // log2(gadget basis)
    uint64_t m_basis_log2;
    // gadget digits
    uint64_t m_digits;
    // RLWE standard deviation
    double m_std;

    // pointer to NTT engine
    std::shared_ptr<MathWorker> m_ntt;

};

struct RLWEtoRLWEConverter : public SchemeConverter<RLWEConversionContext> {

    friend struct RLWEConversionContext;


    /**
     * Performs the conversion of a RLWE sample
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer, input = [A, B] assumed to be in [coefficient,ntt] format
     */
    void Convert(uint64_t *output, const uint64_t *const input) override;

    /**
    * Performs the conversion of a RLWE sample
    * @param output Output vector
    * @param input Input vector, input = [A, B] assumed to be in [coefficient,ntt] format
    */
    void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) override;

    /**
     *
     * @return Parent context
     */
    const std::shared_ptr<const RLWEConversionContext> GetContext() const override;

private:

    /**
     * Constructor from params / context
     * @param params Pointer to parent context
     */
    RLWEtoRLWEConverter(std::shared_ptr<const RLWEConversionContext> params);

    /**
     * Key generation functions ALL INPUTS ASSUMED TO BE IN COEFFICIENT FORM
     * @param source_key Source RLWE key
     * @param target_key Target RLWE key
     */
    void KeyGen(const std::vector<uint64_t> &source_key, const std::vector<uint64_t> &target_key);

    /**
    * Key generation functions ALL INPUTS ASSUMED TO BE IN COEFFICIENT FORM
    * @param source_key Source RLWE key
    * @param target_key Target RLWE key
    */
    void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key);

private:
    // key buffer
    AlignedVector m_ksk;
    // scratch space
    AlignedVector m_acc;
    // pointer to parent context
    std::shared_ptr<const RLWEConversionContext> m_params;
    // debug flags
    bool m_params_set = false;
    bool m_keys_generated = false;

};



#endif //LARGE_FUNCTIONS_ENDO_GLWE_CONVERSION_H
