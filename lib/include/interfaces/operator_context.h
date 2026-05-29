//
// Created by leonard on 5/27/26.
//

#ifndef LARGE_FUNCTIONS_OPERATOR_CONTEXT_H
#define LARGE_FUNCTIONS_OPERATOR_CONTEXT_H

#include "interfaces/container_types.h"
#include "interfaces/enum_ids.h"

#include "keys.h"

/**
 *  Every Operator is associated with a context as we are using a builder-type pattern.
 *  The point here is that every type of operation will most likely pop up with different keys, slightly different params
 *  which we unify in this way. Note, this is not strictly better than a class & params approach.
 *
 *
 * @tparam OperatorType Type of Operator the context will construct.
 * @tparam KeyBundle Type of keys the operator works with.
 */
template<typename OperatorType>
struct OperatorContext {

    /**
     * Almost every FHE operator induces an ouput variance.
     * This function returns said variance as a function of the input variance
     * @param input_variance Variance of input (LWE/RLWE/...)
     * @return Variance induced by operator
     */
    [[nodiscard]] virtual long double ComputeOutputVariance(long double input_variance = 0.0) const = 0;

    /** On input a given KeyBundle, constructs the operator
     *
     * @param keys The KeyBundle as required by the operator
     * @return Pointer to instance of the operator
     */
    [[nodiscard]] virtual std::shared_ptr<OperatorType> ConstructOperator(const std::vector<GenericKey>& keys) const = 0;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @return A Container describing the expected input
     */
    [[nodiscard]] virtual Container GetInputContainer() const = 0;

    /** Every Operator expects a certain input format (e.g. LWE)
     * and possible constraints. The container returned describes
     * the constraints in detail.
     *
     * @param input Input container description. Necessary as some operators are not additive in variance.
     *        if equal to nullptr, is ignored.
     * @return A Container describing the calculated output
     */
    [[nodiscard]] virtual Container GetOutputContainer(Container input) const = 0;

    /** Returns an ID describing the Operator
     *
     * @return The ID
     */
    [[nodiscard]] virtual OperatorID GetOperatorID() const = 0;
};


#endif //LARGE_FUNCTIONS_OPERATOR_CONTEXT_H
