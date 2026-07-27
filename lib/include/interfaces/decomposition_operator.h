//
// Created by leonard on 5/27/26.
//

#ifndef LARGE_FUNCTIONS_DECOMPOSITION_OPERATOR_H
#define LARGE_FUNCTIONS_DECOMPOSITION_OPERATOR_H

#include <cstdint>
#include <vector>

#include "static/container_types.h"
#include "interfaces/operator_context.h"
/** A radix decomposer takes as input an encryption of a datum
 * and returns encryptions of its message or phase as digits in a given radix
 *
 */
 template<typename ContextType>
struct RadixDecomposer {

    /** Performs the decomposition
     *
     * Note that the radix parameter usually needs to be specified for
     * both the associated context, and the method.
     * This is due to the fact that often, for any fixed basis B, the operator will also work
     * for bases B' either s.t. B' | B or B' < B
     *
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer
     * @param radix Radix
     */
    virtual void Decompose(uint64_t* output, const uint64_t *const input, uint64_t radix) = 0;

    /** Performs the decomposition
     *
     * Note that the radix parameter usually needs to be specified for
     * both the associated context, and the method.
     * This is due to the fact that often, for any fixed basis B, the operator will also work
     * for bases B' either s.t. B' | B or B' < B
     *
     * @param output Vector for output
     * @param input Vector for input
     * @param radix Radix
     */
    virtual void Decompose(std::vector<uint64_t>& output, const std::vector<uint64_t>& input, uint64_t radix) = 0;

    /** Returns pointer to context that created the current operator.
    *
    * @return pointer to context
    */
    virtual const std::shared_ptr<const ContextType> GetContext() const = 0;

    /**
     * Destructor, marked as virtual as we often treat them as black-boxes
     */
    virtual ~RadixDecomposer() = default;
};


#endif //LARGE_FUNCTIONS_DECOMPOSITION_OPERATOR_H
