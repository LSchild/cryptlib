//
// Created by leonard on 7/6/26.
//

#ifndef LARGE_FUNCTIONS_MODULUS_SWITCHING_H
#define LARGE_FUNCTIONS_MODULUS_SWITCHING_H

#include "interfaces/enum_ids.h"

/**
 * Computes the output variance after switching moduli
 * @param input_variance Variance of input before the switch
 * @param source_modulus Initial modulus of the input
 * @param target_modulus Output modulus of the operation
 * @param kd Secret Key Distribution
 * @param expected_l0 Mean/Expected value of the hamming weight
 * @return An upper bound on the Modulus Switching variance
 */
long double EstimateModulusSwitchingVariance(long double input_variance, uint64_t source_modulus, uint64_t target_modulus, KeyDistribution kd, uint64_t expected_l0);

#endif //LARGE_FUNCTIONS_MODULUS_SWITCHING_H
