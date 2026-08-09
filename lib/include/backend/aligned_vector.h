//
// Created by leonard on 3/23/26.
//

#ifndef INTERFACES_H
#define INTERFACES_H

#include <vector>
#include <cstdint>
#include "hexl/hexl.hpp"
#include "static/container_types.h"

static const size_t ALIGN_TO = 32;

using AlignedBuffer = std::vector<uint64_t, intel::hexl::AlignedAllocator<uint64_t, ALIGN_TO>>;


#endif //INTERFACES_H
