//
// Created by leonard on 11/12/24.
//

#include <gtest/gtest.h>
#include "hexl/hexl.hpp"

TEST(DUMMY,DUMMY) {

    volatile uint32_t N = 2048;
    volatile uint64_t Q = 1125899906949121;

    std::vector<uint64_t> buffer(N * 1000, 0);

    std::srand(time(nullptr));

    for (long i = 0; i < N * 1000; i++) {
        buffer[i] = rand() % Q;
    }

    auto ntt = intel::hexl::NTT(N, Q);

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < 1000; i++) {
        ntt.ComputeForward(buffer.data() + i * N, buffer.data() + i * N, 1, 1);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
    std::cerr << "Took " << elapsed / 1000.0 << std::endl;

    EXPECT_NE(buffer[0], buffer[buffer.size() - 1]);
    EXPECT_NE(buffer[123], buffer[432]);
}