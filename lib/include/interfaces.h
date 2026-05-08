//
// Created by leonard on 3/23/26.
//

#ifndef INTERFACES_H
#define INTERFACES_H

#include <vector>
#include <cstdint>
#include "hexl/hexl.hpp"
#include "container_types.h"

using AlignedVector = std::vector<uint64_t, intel::hexl::AlignedAllocator<uint64_t, 32>>;

enum KeyDistribution {
    BINARY,
    TERNARY,
    GAUSSIAN
};

enum BlindRotationMethod {
    CGGI,
    BMMP
};

struct OperationParameters {
    virtual long double ComputeOutputVariance(long double input_variance = 0.0) const {
        return 0.0;
    }
};


template<typename T>
struct BlindRotator {

    typedef T ParameterType;

    virtual void SetParams(ParameterType& params) = 0;

    virtual void KeyGen(const std::vector<uint64_t>& lwe_key, const std::vector<uint64_t>& rlwe_key) = 0;

    virtual void KeyGen(const uint64_t* __restrict lwe_key, const uint64_t* __restrict rlwe_key) = 0;

    virtual void BlindRotate(const std::vector<uint64_t>& lwe_vec, std::vector<uint64_t>& rlwe_acc_vec) = 0;

    virtual void BlindRotate(const uint64_t *__restrict lwe_vec, uint64_t* __restrict rlwe_acc_vec) = 0;

    virtual BlindRotationMethod GetMethod() = 0;

    virtual const ParameterType& GetParams() const = 0;

    virtual Container GetInputContainer() const = 0;

    virtual Container GetOutputContainer() const = 0;

    //virtual ~BlindRotator();

};

template<typename SchemeConversionParams>
struct SchemeConverter {

    virtual void SetParams(SchemeConversionParams& params) = 0;

    virtual void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) = 0;

    virtual void KeyGen(const std::vector<uint64_t>& source_key, const std::vector<uint64_t>&  target_key) = 0;

    virtual void Convert(uint64_t* output, const uint64_t*const input) = 0;

    virtual void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) = 0;

    virtual Container GetInputContainer() const = 0;

    virtual Container GetOutputContainer() const = 0;

    virtual const SchemeConversionParams& GetParams() const = 0;
};

#endif //INTERFACES_H
