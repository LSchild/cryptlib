//
// Created by leonard on 4/28/26.
//

#ifndef LARGE_FUNCTIONS_CONTAINER_TYPES_H
#define LARGE_FUNCTIONS_CONTAINER_TYPES_H

#include <cstdint>
#include <variant>
#include <vector>
#include <memory>

enum ContainerLabel {
    RLWE,
    RGSW,
    LWE,
    MLWE,
    GLWE,
    NTRU,
    NTRUVec,
    Tuple,
    Vector,
    Flag,
    Empty
};

/**
 * Base class for our containers. The purpose of such containers are to determine
 * output variances and can be used for analysis of computational graphs
 */
struct ContainerImpl {
    /**
     *
     * @return Container Label
     */
    virtual ContainerLabel GetLabel() {
        return Empty;
    }
};

using Container = std::shared_ptr<ContainerImpl>;

/**
 * Container for Ring Learning With Error (RLWE) inputs & outputs
 */
struct RLWEContainerImpl : public ContainerImpl {

    /**
     * Constructs RLWE Container over a cyclotomic ring with quotient X^N + 1
     * @param N ring dimension
     * @param Q ring modulus
     * @param var error variance
     */
    RLWEContainerImpl(uint64_t N, uint64_t Q, long double var) : m_Q(Q), m_N(N), m_var(var) {}

    /**
     *
     * @return Container Label (here RLWE)
     */
    ContainerLabel GetLabel() override {
        return RLWE;
    }

    /**
     *
     * @return Ring Modulus
     */
    [[nodiscard]] uint64_t GetQ() const {
        return m_Q;
    }

    /**
     *
     * @return Ring Dimension
     */
    [[nodiscard]] uint64_t GetN() const {
        return m_N;
    }

    /**
     *
     * @return Error Variance
     */
    [[nodiscard]] long double GetVariance() const {
        return m_var;
    }

private:

    // ring modulus
    uint64_t m_Q;
    // ring dimension
    uint64_t m_N;
    // error variance
    long double m_var;

};

using RLWEContainer = std::shared_ptr<RLWEContainerImpl>;

/**
 * Container for Learning With Error (LWE) inputs & outputs
 */
struct LWEContainerImpl : public ContainerImpl {

    /**
     * Constructs LWE Container over a rank-1 module Z_Q^(n + 1)
     * @param n module dimension
     * @param Q module modulus
     * @param var error variance
     */
    LWEContainerImpl(uint64_t n, uint64_t Q, long double var) : m_q(Q), m_n(n), m_var(var) {}

    /**
     *
     * @return Container Label (here LWE)
     */
    ContainerLabel GetLabel() override {
        return LWE;
    }

    /**
     *
     * @return Module Modulus
     */
    [[nodiscard]] uint64_t GetQ() const {
        return m_q;
    }

    /**
     *
     * @return Module Dimension
     */
    [[nodiscard]] uint64_t GetN() const {
        return m_n;
    }

    /**
     *
     * @return Error Variance
     */
    [[nodiscard]] long double GetVariance() const {
        return m_var;
    }

private:

    // module modulus
    uint64_t m_q;
    // module dimension
    uint64_t m_n;

    long double m_var;

};

using LWEContainer = std::shared_ptr<LWEContainerImpl>;

/**
 * Container for (ring) GSW inputs & outputs
 * For further details see [GSW13]
 */
struct RGSWContainerImpl : public ContainerImpl {

    /**
     * Constructs RGSW Container over a cyclotomic ring with quotient X^N + 1
     * @param N Ring dimension
     * @param Q Ring modulus
     * @param digits Number of digits for the gadget in the Lev (RLWE') sub-structure
     * @param basis Gadget basis in the Lev (RLWE') sub-structure
     * @param var error variance
     */
    RGSWContainerImpl(uint64_t N, uint64_t Q, uint64_t digits, uint64_t basis, long double var) : m_Q(Q), m_N(N), m_digits(digits), m_basis(basis), m_var(var) {}

    /**
     *
     * @return Container Label (here RGSW)
     */
    ContainerLabel GetLabel() override {
        return RGSW;
    }

    /**
     *
     * @return Ring Modulus
     */
    [[nodiscard]] uint64_t GetQ() const {
        return m_Q;
    }

    /**
     *
     * @return Ring Dimension
     */
    [[nodiscard]] uint64_t GetN() const {
        return m_N;
    }

    /**
     *
     * @return Gadget Digits
     */
    [[nodiscard]] uint64_t GetDigits() const {
        return m_digits;
    }

    /**
     *
     * @return Gadget Basis
     */
    [[nodiscard]] uint64_t GetBasis() const {
        return m_basis;
    }

    /**
     *
     * @return Error Variance
     */
    [[nodiscard]] long double GetVariance() const {
        return m_var;
    }

private:

    // Ring modulus
    uint64_t m_Q;
    // Ring dimension
    uint64_t m_N;
    // Gadget Digits
    uint64_t m_digits;
    // Gadget Basis
    uint64_t m_basis;

    // Error variance
    long double m_var;

};

using RGSWContainer = std::shared_ptr<RGSWContainerImpl>;

/**
 * Container class for tuples of containers
 */
struct TupleContainerImpl : public ContainerImpl {

    /**
     * Constructor
     * @param container_list List/Vector of containers that will form the tuple
     */
    TupleContainerImpl(const std::vector<Container>& container_list) {
        for(auto& elem : container_list) {
            m_elems.push_back(elem);
        }
    }

    /**
     *
     * @return Container Label (here Tuple)
     */
    ContainerLabel GetLabel() override {
        return Tuple;
    }

    /**
     * Lookup function
     * @param idx index of tuple element
     * @return Container at given index
     */
    Container GetElem(uint32_t idx) {
        return m_elems.at(idx);
    }

private:

    // vector of containers in the tuple
    std::vector<Container> m_elems;

};

using TupleContainer = std::shared_ptr<TupleContainerImpl>;

/**
 * Container class for vector of Containers of the same type
 */
struct VectorContainerImpl : public ContainerImpl {

    /**
     * Constructor
     * @param elem Representative Container for the Vector
     * @param max_len Maximum capacity
     */
    VectorContainerImpl(const Container elem, uint64_t max_len) : m_inner_type(elem), m_max_len(max_len) {

    }

    /**
     *
     * @return Maximum Length
     */
    [[nodiscard]] uint64_t GetMaxLen() const {
        return m_max_len;
    }

    /**
     *
     * @return Representative Container type
     */
    [[nodiscard]] Container GetElemType() const {
        return m_inner_type;
    };

    /**
     *
     * @return Container label (here Vector)
     */
    ContainerLabel GetLabel() override {
        return Vector;
    }

private:
    // max length
    uint64_t m_max_len;
    // internal, representative Container
    Container m_inner_type;

};

template<typename FlagValueType>
struct FlagContainerImpl : public ContainerImpl {

    FlagContainerImpl(std::string& flag_name, FlagValueType flag_value) : m_name(flag_name), m_value(flag_value) {

    }

    ContainerLabel GetLabel() override {
        return Flag;
    }

    [[nodiscard]] FlagValueType GetValue() const {
        return m_value;
    }

    [[nodiscard]] const std::string& GetName() const {
        return m_name;
    }

private:

    FlagValueType m_name;
    FlagValueType m_value;

};

template<typename FlagContainerType>
using FlagContainer = std::shared_ptr<FlagContainerImpl<FlagContainerType>>;

#endif //LARGE_FUNCTIONS_CONTAINER_TYPES_H
