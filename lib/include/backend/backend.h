//
// Created by leonard on 6/18/26.
//

#ifndef LARGE_FUNCTIONS_BACKEND_H
#define LARGE_FUNCTIONS_BACKEND_H

#include <cstdint>
#include <concepts>

template<typename T>
concept Backend = requires(T c, uint64_t x, uint64_t y) {
    {T::modulus } -> std::convertible_to<uint64_t>;
    {T::mod_add(x, y)} -> std::convertible_to<uint64_t>;
    {T::mod_sub(x, y)} -> std::convertible_to<uint64_t>;
    {T::mod_mul(x, y)} -> std::convertible_to<uint64_t>;
    {T::mod_exp(x, y)} -> std::convertible_to<uint64_t>;
};



#endif //LARGE_FUNCTIONS_BACKEND_H
