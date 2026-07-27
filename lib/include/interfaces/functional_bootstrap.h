//
// Created by leonard on 7/27/26.
//

#ifndef TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H
#define TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H

#include <cstdint>
#include <functional>

#include "backend/definitions.h"

enum struct OutputKeyType {
    BlindRotationKey,
    LWEKey
};

struct FunctionEvaluator {

    virtual void EvalFunc(uint64_t* result, std::function<uint64_t(uint64_t)>& function, OutputKeyType output) = 0;

    virtual void EvalFunc(std::vector<uint64_t>& result, std::function<uint64_t(uint64_t)>& function, OutputKeyType output) = 0;

    virtual void EvalFunc(uint64_t* RESTRICTED result, uint64_t* RESTRICTED lut_rlwe, OutputKeyType output) = 0;

    virtual void EvalFunc(std::vector<uint64_t>& result, std::vector<uint64_t>& lut_rlwe, OutputKeyType output) = 0;

    ~FunctionEvaluator() = default;
};

#endif //TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H
