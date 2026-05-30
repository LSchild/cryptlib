//
// Created by leonard on 5/27/26.
//

#ifndef LARGE_FUNCTIONS_BLINDROTATION_OPERATOR_H
#define LARGE_FUNCTIONS_BLINDROTATION_OPERATOR_H

#include <cstdint>
#include <vector>

#include "interfaces/operator_context.h"

/**
 * In general, BlindRotator rotate the coefficients of a polynomial P
 * as a function of a LWE ciphertext c.
 * The interfaces makes no difference between accumulators of different types.
 * The most common blind-rotation computes P * X^{-Phase(x)}. Over the cyclic ring, this would correspond to moving
 * P[Phase(x)] into the constant coefficient.
 */
struct BlindRotator {

    /** Performs the operation and writes the result to the first argument
     *
     * @param result Buffer for the result.
     * @param lwe_ciphertext Vector containing the coefficients of the LWE sample. We expect for LWE(c) = [vec{a}, b], i.e.
     * b is in the coefficient with index n.
     * @param accumulator Vector containing the accumulator
     */
    virtual void BlindRotate(std::vector<uint64_t> &result, const std::vector<uint64_t> &lwe_ciphertext,
                             std::vector<uint64_t> &accumulator, bool output_as_coefficients = false) = 0;

    /** Performs the operation and writes the result to the first argument if non-null, overwise overwrite the accumulator.
     *
     * @param result Pointer to the result buffer. Can be null.
     * @param lwe_ciphertext Pointer to the coefficients of the LWE sample. We expect for LWE(c) = [vec{a}, b], i.e.
     * b is in the coefficient with index n.
     * @param accumulator Pointer the accumulator
     */
    virtual void
    BlindRotate(uint64_t *result, const uint64_t *lwe_vec, uint64_t *accumulator, bool output_as_coefficients = false) = 0;

    /** Returns pointer to context that created the current operator.
    *
    * @return pointer to context
    */
    virtual const std::shared_ptr<const OperatorContext<BlindRotator>> GetContext() const = 0;

    virtual ~BlindRotator() = default;
};

#endif //LARGE_FUNCTIONS_BLINDROTATION_OPERATOR_H
