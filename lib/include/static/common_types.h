//
// Created by leonard on 3/23/26.
//

#ifndef INTERFACES_H
#define INTERFACES_H

#include <vector>
#include <cstdint>
#include "hexl/hexl.hpp"
#include "static/container_types.h"

using AlignedVector = std::vector<uint64_t, intel::hexl::AlignedAllocator<uint64_t, 32>>;

#endif //INTERFACES_H
