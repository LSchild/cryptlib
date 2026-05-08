//
// Created by leonard on 3/13/25.
//

#include <bitset>
#include <cstdint>
#include <string>
#include <vector>
#include <chrono>
#include "aesctr.h"
#include <iostream>

std::vector<uint64_t> respire_query_pack(uint32_t dim2_bits, uint32_t dim3_bits, uint32_t N, uint32_t dim_1_index, uint32_t dim_2_index, uint32_t dim_3_index, uint32_t basis_bits, uint32_t digits) {
    // note that we cannot statically allocate a buffer as the size of the query strongly depends
    // on N
    // we pack the respire query into bits as follows
    // for the 1st dimension, the result is a Nbit string K with K_{dim_1_index} = 1 and 0 elsewhere
    // for the second and third dimension we store the gadgets required for each rgsw as bits

    // calculate size necessary in bits
    auto dim23_bit_size = dim2_bits + dim3_bits;
    auto total_size = (N >> 6) + (dim23_bit_size) * digits * basis_bits;
    std::cerr << N + ((dim23_bit_size) * digits * basis_bits * 64) / 56 << " ";
    std::vector<uint64_t> buffer(total_size, 0);

    // dim1
    buffer[(dim_1_index >> 6) % (N >> 6)] = dim_1_index & (64 - 1);
    auto current_offset = (dim_1_index >> 6) % (N >> 6);
    for (uint32_t i = 0; i < dim2_bits; i++) {
        auto bvi = (dim_2_index >> i) & 1;
        for (uint32_t j = 0; j < digits; j++) {
            buffer[current_offset] = bvi << (j * basis_bits);
            current_offset++;
        }
    }
    // dim1
    for (uint32_t i = 0; i < dim3_bits; i++) {
        auto bvi = (dim_3_index >> i) & 1;
        for (uint32_t j = 0; j < digits; j++) {
            buffer[current_offset] = bvi << (j * basis_bits);
            current_offset++;
        }
    }

    return buffer;
}

int main(int argc, char** argv) {

    uint64_t mu1 = std::stoll(argv[1]);
    uint64_t mu2 = std::stoll(argv[2]);
    uint64_t mu3 = std::stoll(argv[3]);

    uint64_t index_ = std::stoll(argv[4]);
    uint64_t N = std::stoll(argv[5]);
    uint64_t L = std::stoll(argv[6]);
    uint64_t l = std::stoll(argv[7]);

    uint64_t rounds = std::stoll(argv[8]);
    uint64_t seed_lo = std::atoll(argv[9]);
    uint64_t seed_hi = std::atoll(argv[10]);

    __uint128_t full_seed = (__uint128_t(seed_hi) << 64) + seed_lo;
    aes::AesCtr<> prng;
    prng.seed(full_seed);

    std::vector<uint64_t> results(rounds, 0);

    auto start = std::chrono::high_resolution_clock::now();

    for (uint32_t j = 0; j < rounds; j++) {

        auto index = index_ |= (1 << j);

        auto input = respire_query_pack(mu2, mu3, N, index >> (mu3 + mu2), index >> (mu2), index, L, l);
        for (auto&v : input) {
            v &= prng();
        }

        results[j] = input.back();
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

    for (auto& v: results) {
        std::cout << v << std::endl;
    }

    std::cout << "Took elasped = " << elapsed << "us. " << std::endl;


}
