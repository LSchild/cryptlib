//
// Created by leonard on 5/23/26.
//

#ifndef LARGE_FUNCTIONS_BLINDROTATION_COMMON_H
#define LARGE_FUNCTIONS_BLINDROTATION_COMMON_H

#include <cstdint>

struct BlindRotationKeys {
    const uint64_t* lwe_sk;
    const uint64_t* rlwe_sk;
};

#endif //LARGE_FUNCTIONS_BLINDROTATION_COMMON_H
