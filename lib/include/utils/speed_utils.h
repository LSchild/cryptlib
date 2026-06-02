//
// Created by leonard on 4/16/26.
//

#ifndef LARGE_FUNCTIONS_SPEED_UTILS_H
#define LARGE_FUNCTIONS_SPEED_UTILS_H

#include <cstdint>
#include <cstring>

#define ZERO_UINT64_ARR(start, N) std::memset(start, 0, N * sizeof(uint64_t))

template <class T>
T reverse_bits(T n) {
    short bits = sizeof(n) * 8;
    T mask = ~T(0); // equivalent to uint32_t mask = 0b11111111111111111111111111111111;

    while (bits >>= 1) {
        mask ^= mask << (bits); // will convert mask to 0b00000000000000001111111111111111;
        n = (n & ~mask) >> bits | (n & mask) << bits; // divide and conquer
    }

    return n;
}


#endif //LARGE_FUNCTIONS_SPEED_UTILS_H
