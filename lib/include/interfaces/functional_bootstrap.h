//
// Created by leonard on 7/27/26.
//

#ifndef TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H
#define TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H

#include <cstdint>
#include <functional>
#include <string>

#include "backend/definitions.h"

enum OutputKeyType {
    BlindRotationKey = 0,
    LWEKey = 1
};

static const std::string OutputKeyTypeLabels[2] = {"BlindRotationKey", "LWEKey"};

struct FunctionEvaluator {

    virtual void EvalFunc(uint64_t* RESTRICTED result, const uint64_t* RESTRICTED input_lwe, std::function<uint64_t(uint64_t)>& function, OutputKeyType output) = 0;

    virtual void EvalFunc(std::vector<uint64_t>& result, const std::vector<uint64_t>& input_lwe, std::function<uint64_t(uint64_t)>& function, OutputKeyType output) = 0;

    virtual void EvalFunc(uint64_t* RESTRICTED result, const uint64_t* RESTRICTED input_lwe, uint64_t* RESTRICTED lut_rlwe, OutputKeyType output) = 0;

    virtual void EvalFunc(std::vector<uint64_t>& result, const std::vector<uint64_t>& input_lwe, std::vector<uint64_t>& lut_rlwe, OutputKeyType output) = 0;

    virtual void Finalize(uint64_t* RESTRICTED result, const uint64_t* RESTRICTED input) = 0;

    virtual void Finalize(std::vector<uint64_t>& result, const std::vector<uint64_t>& input) = 0;

    ~FunctionEvaluator() = default;
};

#endif //TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H
