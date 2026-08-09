//
// Created by leonard on 8/9/26.
//

#include <cstdint>

#ifndef TOOTHPASTE_SOLINAS_2_32_1_M1_AVX_H
#define TOOTHPASTE_SOLINAS_2_32_1_M1_AVX_H

#ifdef AAAAA

#include <immintrin.h>

#define U64x4 __m256i

#define AddU64x4(a,b) _mm256_add_epi64(a,b)
#define SubU64x4(a,b) _mm256_sub_epi64(a,b)

// Source: Peter Cordes @ stackoverflow
U64x4 mul64_avx2 (U64x4 a, U64x4 b)
{
    // There is no vpmullq until AVX-512. Split into 32-bit multiplies
    // Given a and b composed of high<<32 | low  32-bit halves
    // a*b = a_low*(u64)b_low  + (u64)(a_high*b_low + a_low*b_high)<<32;  // same for signed or unsigned a,b since we aren't widening to 128
    // the a_high * b_high product isn't needed for non-widening; its place value is entirely outside the low 64 bits.

    U64x4 b_swap  = _mm256_shuffle_epi32(b, _MM_SHUFFLE(2,3, 0,1));   // swap H<->L
    U64x4 crossprod  = _mm256_mullo_epi32(a, b_swap);                 // 32-bit L*H and H*L cross-products

    U64x4 prodlh = _mm256_slli_epi64(crossprod, 32);          // bring the low half up to the top of each 64-bit chunk
    U64x4 prodhl = _mm256_and_si256(crossprod, _mm256_set1_epi64x(0xFFFFFFFF00000000)); // isolate the other, also into the high half were it needs to eventually be
    U64x4 sumcross = _mm256_add_epi32(prodlh, prodhl);       // the sum of the cross products, with the low half of each u64 being 0.

    U64x4 prodll  = _mm256_mul_epu32(a,b);                  // widening 32x32 => 64-bit  low x low products
    U64x4 prod    = _mm256_add_epi32(prodll, sumcross);     // add the cross products into the high half of the result
    return  prod;
}

#else



#endif

#endif //TOOTHPASTE_SOLINAS_2_32_1_M1_AVX_H
