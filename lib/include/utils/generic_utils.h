//
// Created by leonard on 4/16/26.
//

#ifndef LARGE_FUNCTIONS_SPEED_UTILS_H
#define LARGE_FUNCTIONS_SPEED_UTILS_H

#include <cstdint>
#include <cstring>
#include <memory>

#define ZERO_UINT64_ARR(start, N) std::memset(start, 0, N * sizeof(uint64_t))

template<typename Base, typename Derived> requires std::is_base_of_v<Base, Derived>
bool HasUnderlyingType(std::shared_ptr<Base> ptr) {
    return dynamic_cast<Derived*>(ptr.get()) != nullptr;
}

#endif //LARGE_FUNCTIONS_SPEED_UTILS_H
