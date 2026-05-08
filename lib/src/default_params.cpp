//
// Created by leonard on 3/10/25.
//

#include "pirouette.h"

std::unique_ptr<PirouettePIR> PirouettePIR::GenerateDefaultConfig(PIROUETTE_MODE mode) {
    std::srand(0);
    // set up parameters

    // variance
    double std = 0;//3.19; // for testing
    double ksk_var_after_bitdecomp = 0;//3.19;

    // moduli
    uint64_t ring_modulus = 72057594037641217;
    uint64_t small_lwe_modulus = 1 << 12;
    uint32_t large_lwe_modulus_bits = 32;
    uint64_t large_lwe_modulus = 1ull << large_lwe_modulus_bits;
    uint64_t ksk_modulus_during_dd = 1ull << 42;
    uint64_t ksk_after_bit_decomp_modulus = 1 << 30;

    uint32_t large_lwe_n = 1800;
    uint32_t lwe_n = 600;

    // (r)lwe dimensions
    uint32_t N_bits = std::round(std::log2((long double) small_lwe_modulus)) - 1;
    uint32_t N = 1u << N_bits;

    // bound values
    uint32_t error_width = 7;
    uint32_t alpha = 1 << error_width;
    uint32_t beta = 1 << (error_width - 2);
    uint32_t clear_bits = (N_bits + 1) - error_width;

    // bases
    uint32_t decomp_br_basis = 1 << 12;
    uint32_t decomp_auto_basis = 1 << 4;

    uint32_t sswitch_br_basis = 1 << 10;
    uint32_t sswitch_auto_basis = 1 << 4;
    uint32_t sswitch_squaring_basis = 1 << 2;
    // Basis for output RGSW
    uint32_t sswitch_rgsw_basis = 1 << 4;
    uint32_t ksk_basis = 1 << 3;

    // Number of Digits for output RGSW
    uint32_t skip_n_first_digits = 6;

    // crt
    uint32_t P = 268496897;
    uint32_t Q = 268460033;

    // set up parameter sets

    DigitDecompositionParameters decomposerParams{large_lwe_n, N, small_lwe_modulus,
                                         large_lwe_modulus, large_lwe_modulus_bits, ring_modulus,
        std, decomp_br_basis, clear_bits, error_width, ksk_modulus_during_dd, ksk_basis};

    SchemeSwitchParameters switchParameters{lwe_n, N, small_lwe_modulus, ring_modulus, std, skip_n_first_digits, sswitch_br_basis, sswitch_squaring_basis, sswitch_auto_basis, sswitch_rgsw_basis};

    auto ksk_params = std::make_shared<lbcrypto::LWECryptoParams>(lwe_n, N, small_lwe_modulus, ring_modulus, ksk_after_bit_decomp_modulus, ksk_var_after_bitdecomp, ksk_basis);
    auto genT = BinaryUniformGeneratorImpl<NativeVector>();
    NativeVector huge_key_vec = genT.GenerateVector(large_lwe_n, large_lwe_modulus);

    // there's no nicer way to set the hamming weight...
    auto one_ctr = 0;
    for (uint32_t i = 0; i < large_lwe_n; i++) {
        if (huge_key_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                huge_key_vec[i] = 0;
            }

        }
    }
    NativeVector ring_sk_vec = genT.GenerateVector(N, ring_modulus);
    one_ctr = 0;
    for (uint32_t i = 0; i < large_lwe_n; i++) {
        if (ring_sk_vec[i] == 1) {
            if (one_ctr < beta) {
                one_ctr++;
            } else {
                ring_sk_vec[i] = 0;
            }

        }
    }

    NativeInteger rmod = ring_modulus;
    auto pp = std::make_shared<ILNativeParams>(2 * N, rmod);
    NativePoly sk_out(pp, COEFFICIENT, true);
    sk_out.SetValues(ring_sk_vec, COEFFICIENT);

    auto huge_lwe_key = std::make_shared<LWEPrivateKeyImpl>(huge_key_vec);
    auto small_lwe_key_vec = genT.GenerateVector(lwe_n, small_lwe_modulus);
    auto small_lwe_key = std::make_shared<LWEPrivateKeyImpl>(small_lwe_key_vec);


    PirouetteParams pirouette_params {
        decomposerParams,
        huge_lwe_key,
        switchParameters,
        small_lwe_key,
        sk_out,
        ksk_params,
        P,
        Q,
        4,
        2,
        alpha
    };

    auto pir = std::make_unique<PirouettePIR>(mode);
    pir->SetCryptoParams(pirouette_params);

    return pir;
}
