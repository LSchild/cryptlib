//
// Created by leonard on 5/27/26.
//

#ifndef LARGE_FUNCTIONS_CONVERSION_OPERATOR_H
#define LARGE_FUNCTIONS_CONVERSION_OPERATOR_H

#include <cstdint>
#include <vector>

#include "interfaces/operator_context.h"

/** The purpose of scheme converters are to convert instance of
 * one scheme into another.
 * Examples include, but are not limited to
 * - LWE to RLWE switching
 * - LWE to RGSW switching (circuit bootstrapping)
 *
 */
struct SchemeConverter {

    /** Performs the conversion.
     *
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer
     */
    virtual void Convert(uint64_t* output, const uint64_t* const input) = 0;

    /** Performs the conversion.
    *
    * @param output Vector for output
    * @param input Vector for input
    */
    virtual void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) = 0;

    /** Returns pointer to context that created the current operator.
     *
     * @return pointer to context
     */
    virtual const std::shared_ptr<const OperatorContext<SchemeConverter>> GetContext() const = 0;
};


#endif //LARGE_FUNCTIONS_CONVERSION_OPERATOR_H
