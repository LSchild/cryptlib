//
// Created by leonard on 6/18/26.
//

#ifndef LARGE_FUNCTIONS_BACKEND_H
#define LARGE_FUNCTIONS_BACKEND_H

#include <cstdint>
#include <concepts>
#include <vector>

struct MathWorker {

    virtual void ForwardNTT(uint64_t* result, uint64_t* input) = 0;

    virtual void ForwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) = 0;

    virtual void BackwardNTT(uint64_t* result, uint64_t* input) = 0;

    virtual void BackwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) = 0;

    virtual uint64_t AddMod(uint64_t left, uint64_t right) = 0;

    virtual uint64_t SubMod(uint64_t left, uint64_t right) = 0;

    virtual uint64_t MulMod(uint64_t left, uint64_t right) = 0;

    virtual uint64_t FMAMod(uint64_t left, uint64_t right, uint64_t summand) = 0;

    virtual uint64_t InvMod(uint64_t left) = 0;

    virtual uint64_t PowMod(uint64_t base, uint64_t exponent) = 0;

    virtual void AddModElw(uint64_t* __restrict result, uint64_t* __restrict left, uint64_t* __restrict right, uint64_t len) = 0;

    virtual void AddEqModElw(uint64_t* __restrict left_out, uint64_t* __restrict right, uint64_t len) = 0;

    virtual void SubModElw(uint64_t* __restrict result, uint64_t* __restrict left, uint64_t* __restrict right, uint64_t len) = 0;

    virtual void SubEqModElw(uint64_t* __restrict left_out, uint64_t* __restrict right, uint64_t len) = 0;

    virtual void MulModElw(uint64_t* __restrict result, uint64_t* __restrict left, uint64_t* __restrict right, uint64_t len) = 0;

    virtual void MulEqModElw(uint64_t* __restrict left_out, uint64_t* __restrict right, uint64_t len) = 0;

    virtual void FMAModElw(uint64_t* __restrict result, uint64_t* __restrict left, uint64_t* __restrict right, uint64_t* __restrict summand, uint64_t len) = 0;

    virtual void FMAEqModElw(uint64_t* __restrict left_out, uint64_t* __restrict right, uint64_t* __restrict summand, uint64_t len) = 0;

    virtual void MulSModElw(uint64_t* __restrict result, uint64_t* __restrict left, uint64_t right, uint64_t len) = 0;

    virtual void MulSEqModElw(uint64_t* __restrict left_out, uint64_t right, uint64_t len) = 0;

    virtual uint64_t GetModulus() const = 0;

    virtual uint64_t GetDimension() const = 0;

    virtual ~MathWorker() = default;


};


template<typename T>
concept Backend = requires(T c, uint64_t x, uint64_t y) {
    {T::modulus } -> std::convertible_to<uint64_t>;
    {T::mod_add(x, y)} -> std::convertible_to<uint64_t>;
    {T::mod_sub(x, y)} -> std::convertible_to<uint64_t>;
    {T::mod_mul(x, y)} -> std::convertible_to<uint64_t>;
    {T::mod_exp(x, y)} -> std::convertible_to<uint64_t>;
};



#endif //LARGE_FUNCTIONS_BACKEND_H
