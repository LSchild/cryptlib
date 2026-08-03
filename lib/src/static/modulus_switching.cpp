//
// Created by leonard on 7/6/26.
//

#include <cstdint>
#include <iostream>
#include <functional>
#include <cmath>
#include <random>
#include "static/modulus_switching.h"

long double MSErrorBT(long double input_variance, uint64_t source_modulus, uint64_t target_modulus, uint64_t expected_l0) {
    long double Qs = source_modulus;
    auto Qs2 = Qs * Qs;
    long double Qt = target_modulus;
    auto Qt2 = Qt * Qt;

    long double principal = (Qs2 * input_variance) / Qt2;

    long double hamming_mean = expected_l0;

    return principal + hamming_mean / 12.0;
}

long double EstimateModulusSwitchingVariance(long double input_variance, uint64_t source_modulus, uint64_t target_modulus, KeyDistribution kd, uint64_t expected_l0) {

    if (kd == BINARY or kd == TERNARY) {
        return MSErrorBT(input_variance, source_modulus, target_modulus, expected_l0);
    } else {
        std::cerr << "Modulus Switching Variance computation implemented for Binary/Ternary keys only for now" << std::endl;
        std::exit(1);
    }

}

void ModulusSwitch(uint64_t* vec, uint64_t n, uint64_t source_modulus, uint64_t target_modulus, ModulusSwitchType type) {

    std::function<uint64_t(uint64_t)> round_f;

    long double sM = source_modulus;
    long double tM = target_modulus;

    switch (type) {
        case ModulusSwitchType::FLOOR: {
            round_f = [sM,tM] (uint64_t v) -> uint64_t {
                return std::floorl((tM * (long double)(v)) / sM);
            };
            break;
        }
        case ModulusSwitchType::CEIL: {
            round_f = [sM,tM] (uint64_t v) -> uint64_t {
                return std::ceill((tM * (long double)(v)) / sM);
            };
            break;
        }
        case ModulusSwitchType::ROUND : {
            round_f = [sM,tM] (uint64_t v) ->  uint64_t  {
                auto hi = tM * static_cast<long double>(v);
                auto res = hi / sM;
                return std::roundl(res);
            };
            break;
        }

        case ModulusSwitchType::RANDOM : {

            std::random_device random_device;
            std::mt19937 engine(random_device());
            std::uniform_real_distribution<float> unif_1_0(0.0, 1.0);
            unif_1_0(engine);

            round_f = [sM, tM, unif_1_0, engine] (uint64_t v) mutable -> uint64_t {
                auto vv = (tM * (long double)(v)) / sM;
                if (unif_1_0(engine) >= 0.5) {
                    return std::floorl(vv);
                } else {
                    return std::ceill(vv);
                }
            };
            break;
        }

    }

    for(uint64_t i = 0; i < n; i++) {
        vec[i] = round_f(vec[i]);
        if (vec[i] >= target_modulus) {
            vec[i] -= target_modulus;
        }
    }


}