//
// Created by leonard on 3/13/25.
//

#include <chrono>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include "aesctr.h"

int main(int argc, char** argv) {

    // argv[1] n
    // argv[2] q_bits
    // argv[3] msg
    // argv[4] runs
    // stddev = 3.19 always

    uint32_t n = std::stoi(argv[1]);
    uint64_t q_bits = std::stoi(argv[2]);
    uint64_t q = 1ull << q_bits;
    uint64_t q_mask = (1ull << q_bits) - 1;
    uint64_t msg = std::atoll(argv[3]);
    uint64_t runs = std::atoll(argv[4]);
    uint64_t seed_lo = std::atoll(argv[5]);
    uint64_t seed_hi = std::atoll(argv[6]);

    __uint128_t full_seed = (__uint128_t(seed_hi) << 64) + seed_lo;

    // generate some key
    // in an implementation that would be an input

    std::vector<uint64_t> sk(n, 0);
    for (uint32_t i = 0; i < n; i++)
        sk[i] = std::rand() & 1;

    std::default_random_engine gen;
    std::normal_distribution<float> distr(q >> 1, 3.19);

    aes::AesCtr<> prng;
    prng.seed(full_seed);

    std::vector<uint64_t> results(runs, 0);

    auto start = std::chrono::high_resolution_clock::now();

    for (uint32_t j = 0; j < runs; j++) {

        uint64_t running_sum = msg + (uint64_t(distr(gen)));
        for (uint32_t n_i = 0; n_i < n; n_i++) {
            auto a_i = prng();
            running_sum = running_sum + sk[n_i] & a_i;
        }
        results[j] = running_sum;
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(stop-start).count();

    std::cout << "Took elasped = " << elapsed << "us. " << std::endl;

    for (auto& v : results) {
        std::cout << v << std::endl;
    }

}
