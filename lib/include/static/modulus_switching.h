//
// Created by leonard on 7/6/26.
//

#ifndef LARGE_FUNCTIONS_MODULUS_SWITCHING_H
#define LARGE_FUNCTIONS_MODULUS_SWITCHING_H

#include "enum_ids.h"

/**
 * Computes the output variance after switching moduli
 * @param input_variance Variance of input before the switch
 * @param source_modulus Initial Q of the input
 * @param target_modulus Output Q of the operation
 * @param kd Secret Key Distribution
 * @param expected_l0 Mean/Expected value of the hamming weight
 * @return An upper bound on the Modulus Switching variance
 */
long double EstimateModulusSwitchingVariance(long double input_variance, uint64_t source_modulus, uint64_t target_modulus, KeyDistribution kd, uint64_t expected_l0);

/**
 * Class describing different Q switching approaches
 */
enum struct ModulusSwitchType {
    ROUND,
    FLOOR,
    CEIL,
    RANDOM // coin toss to determine whether to round up or down
};

/**
 * Modulus switch function, for each element of vec computes
 *
 * vec[i] = ROUND_FUNCTION(target_modulus * vec[i] / source_modulus)
 * where ROUND_function depends on the given Q switch type
 *
 * @param vec
 * @param n
 * @param source_modulus
 * @param target_modulus
 * @param type
 */
void ModulusSwitch(uint64_t* vec, uint64_t n, uint64_t source_modulus, uint64_t target_modulus, ModulusSwitchType type);

#endif //LARGE_FUNCTIONS_MODULUS_SWITCHING_H
