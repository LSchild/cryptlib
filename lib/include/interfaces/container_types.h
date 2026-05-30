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
    Empty
};

struct ContainerImpl {
    virtual ContainerLabel GetLabel() {
        return Empty;
    }
};

using Container = std::shared_ptr<ContainerImpl>;


struct RLWEContainerImpl : public ContainerImpl {

    RLWEContainerImpl(uint64_t N, uint64_t Q, long double var) : m_Q(Q), m_N(N), m_var(var) {}

    ContainerLabel GetLabel() override {
        return RLWE;
    }

    [[nodiscard]] uint64_t GetQ() const {
        return m_Q;
    }

    [[nodiscard]] uint64_t GetN() const {
        return m_N;
    }

    [[nodiscard]] long double GetVariance() const {
        return m_var;
    }

private:

    uint64_t m_Q, m_N;
    long double m_var;

};

using RLWEContainer = std::shared_ptr<RLWEContainerImpl>;

struct LWEContainerImpl : public ContainerImpl {

    LWEContainerImpl(uint64_t n, uint64_t Q, long double var) : m_q(Q), m_n(n), m_var(var) {}

    ContainerLabel GetLabel() override {
        return LWE;
    }

    [[nodiscard]] uint64_t GetQ() const {
        return m_q;
    }

    [[nodiscard]] uint64_t GetN() const {
        return m_n;
    }

    [[nodiscard]] long double GetVariance() const {
        return m_var;
    }

private:

    uint64_t m_q, m_n;
    long double m_var;

};

using LWEContainer = std::shared_ptr<LWEContainerImpl>;

struct RGSWContainerImpl : public ContainerImpl {

    RGSWContainerImpl(uint64_t N, uint64_t Q, uint64_t digits, uint64_t basis, long double var) : m_Q(Q), m_N(N), m_digits(digits), m_basis(basis), m_var(var) {}

    ContainerLabel GetLabel() override {
        return RGSW;
    }

    [[nodiscard]] uint64_t getQ() const {
        return m_Q;
    }

    [[nodiscard]] uint64_t getN() const {
        return m_N;
    }

    [[nodiscard]] uint64_t GetDigits() const {
        return m_digits;
    }

    [[nodiscard]] uint64_t GetBasis() const {
        return m_basis;
    }

    [[nodiscard]] long double getVar() const {
        return m_var;
    }

private:

    uint64_t m_Q, m_N, m_digits, m_basis;
    long double m_var;

};

using RGSWContainer = std::shared_ptr<RGSWContainerImpl>;

struct TupleContainerImpl : public ContainerImpl {

    TupleContainerImpl(const std::vector<Container>& container_list) {
        for(auto& elem : container_list) {
            m_elems.push_back(elem);
        }
    }

    ContainerLabel GetLabel() override {
        return Tuple;
    }

    Container GetElem(uint32_t idx) {
        return m_elems.at(idx);
    }

private:

    std::vector<Container> m_elems;

};

using TupleContainer = std::shared_ptr<TupleContainerImpl>;



#endif //LARGE_FUNCTIONS_CONTAINER_TYPES_H
