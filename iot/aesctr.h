#ifndef AESCTR_H
#define AESCTR_H

// Taken from https://github.com/lemire/testingRNG
// Added C++ interface compatible with std::shuffle, &c.

// contributed by Samuel Neves
/*
                          Apache License
                           Version 2.0, January 2004
                        http://www.apache.org/licenses/

   TERMS AND CONDITIONS FOR USE, REPRODUCTION, AND DISTRIBUTION

   1. Definitions.

      "License" shall mean the terms and conditions for use, reproduction,
      and distribution as defined by Sections 1 through 9 of this document.

      "Licensor" shall mean the copyright owner or entity authorized by
      the copyright owner that is granting the License.

      "Legal Entity" shall mean the union of the acting entity and all
      other entities that control, are controlled by, or are under common
      control with that entity. For the purposes of this definition,
      "control" means (i) the power, direct or indirect, to cause the
      direction or management of such entity, whether by contract or
      otherwise, or (ii) ownership of fifty percent (50%) or more of the
      outstanding shares, or (iii) beneficial ownership of such entity.

      "You" (or "Your") shall mean an individual or Legal Entity
      exercising permissions granted by this License.

      "Source" form shall mean the preferred form for making modifications,
      including but not limited to software source code, documentation
      source, and configuration files.

      "Object" form shall mean any form resulting from mechanical
      transformation or translation of a Source form, including but
      not limited to compiled object code, generated documentation,
      and conversions to other media types.

      "Work" shall mean the work of authorship, whether in Source or
      Object form, made available under the License, as indicated by a
      copyright notice that is included in or attached to the work
      (an example is provided in the Appendix below).

      "Derivative Works" shall mean any work, whether in Source or Object
      form, that is based on (or derived from) the Work and for which the
      editorial revisions, annotations, elaborations, or other modifications
      represent, as a whole, an original work of authorship. For the purposes
      of this License, Derivative Works shall not include works that remain
      separable from, or merely link (or bind by name) to the interfaces of,
      the Work and Derivative Works thereof.

      "Contribution" shall mean any work of authorship, including
      the original version of the Work and any modifications or additions
      to that Work or Derivative Works thereof, that is intentionally
      submitted to Licensor for inclusion in the Work by the copyright owner
      or by an individual or Legal Entity authorized to submit on behalf of
      the copyright owner. For the purposes of this definition, "submitted"
      means any form of electronic, verbal, or written communication sent
      to the Licensor or its representatives, including but not limited to
      communication on electronic mailing lists, source code control systems,
      and issue tracking systems that are managed by, or on behalf of, the
      Licensor for the purpose of discussing and improving the Work, but
      excluding communication that is conspicuously marked or otherwise
      designated in writing by the copyright owner as "Not a Contribution."

      "Contributor" shall mean Licensor and any individual or Legal Entity
      on behalf of whom a Contribution has been received by Licensor and
      subsequently incorporated within the Work.

   2. Grant of Copyright License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      copyright license to reproduce, prepare Derivative Works of,
      publicly display, publicly perform, sublicense, and distribute the
      Work and such Derivative Works in Source or Object form.

   3. Grant of Patent License. Subject to the terms and conditions of
      this License, each Contributor hereby grants to You a perpetual,
      worldwide, non-exclusive, no-charge, royalty-free, irrevocable
      (except as stated in this section) patent license to make, have made,
      use, offer to sell, sell, import, and otherwise transfer the Work,
      where such license applies only to those patent claims licensable
      by such Contributor that are necessarily infringed by their
      Contribution(s) alone or by combination of their Contribution(s)
      with the Work to which such Contribution(s) was submitted. If You
      institute patent litigation against any entity (including a
      cross-claim or counterclaim in a lawsuit) alleging that the Work
      or a Contribution incorporated within the Work constitutes direct
      or contributory patent infringement, then any patent licenses
      granted to You under this License for that Work shall terminate
      as of the date such litigation is filed.

   4. Redistribution. You may reproduce and distribute copies of the
      Work or Derivative Works thereof in any medium, with or without
      modifications, and in Source or Object form, provided that You
      meet the following conditions:

      (a) You must give any other recipients of the Work or
          Derivative Works a copy of this License; and

      (b) You must cause any modified files to carry prominent notices
          stating that You changed the files; and

      (c) You must retain, in the Source form of any Derivative Works
          that You distribute, all copyright, patent, trademark, and
          attribution notices from the Source form of the Work,
          excluding those notices that do not pertain to any part of
          the Derivative Works; and

      (d) If the Work includes a "NOTICE" text file as part of its
          distribution, then any Derivative Works that You distribute must
          include a readable copy of the attribution notices contained
          within such NOTICE file, excluding those notices that do not
          pertain to any part of the Derivative Works, in at least one
          of the following places: within a NOTICE text file distributed
          as part of the Derivative Works; within the Source form or
          documentation, if provided along with the Derivative Works; or,
          within a display generated by the Derivative Works, if and
          wherever such third-party notices normally appear. The contents
          of the NOTICE file are for informational purposes only and
          do not modify the License. You may add Your own attribution
          notices within Derivative Works that You distribute, alongside
          or as an addendum to the NOTICE text from the Work, provided
          that such additional attribution notices cannot be construed
          as modifying the License.

      You may add Your own copyright statement to Your modifications and
      may provide additional or different license terms and conditions
      for use, reproduction, or distribution of Your modifications, or
      for any such Derivative Works as a whole, provided Your use,
      reproduction, and distribution of the Work otherwise complies with
      the conditions stated in this License.

   5. Submission of Contributions. Unless You explicitly state otherwise,
      any Contribution intentionally submitted for inclusion in the Work
      by You to the Licensor shall be under the terms and conditions of
      this License, without any additional terms or conditions.
      Notwithstanding the above, nothing herein shall supersede or modify
      the terms of any separate license agreement you may have executed
      with Licensor regarding such Contributions.

   6. Trademarks. This License does not grant permission to use the trade
      names, trademarks, service marks, or product names of the Licensor,
      except as required for reasonable and customary use in describing the
      origin of the Work and reproducing the content of the NOTICE file.

   7. Disclaimer of Warranty. Unless required by applicable law or
      agreed to in writing, Licensor provides the Work (and each
      Contributor provides its Contributions) on an "AS IS" BASIS,
      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
      implied, including, without limitation, any warranties or conditions
      of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A
      PARTICULAR PURPOSE. You are solely responsible for determining the
      appropriateness of using or redistributing the Work and assume any
      risks associated with Your exercise of permissions under this License.

   8. Limitation of Liability. In no event and under no legal theory,
      whether in tort (including negligence), contract, or otherwise,
      unless required by applicable law (such as deliberate and grossly
      negligent acts) or agreed to in writing, shall any Contributor be
      liable to You for damages, including any direct, indirect, special,
      incidental, or consequential damages of any character arising as a
      result of this License or out of the use or inability to use the
      Work (including but not limited to damages for loss of goodwill,
      work stoppage, computer failure or malfunction, or any and all
      other commercial damages or losses), even if such Contributor
      has been advised of the possibility of such damages.

   9. Accepting Warranty or Additional Liability. While redistributing
      the Work or Derivative Works thereof, You may choose to offer,
      and charge a fee for, acceptance of support, warranty, indemnity,
      or other liability obligations and/or rights consistent with this
      License. However, in accepting such obligations, You may act only
      on Your own behalf and on Your sole responsibility, not on behalf
      of any other Contributor, and only if You agree to indemnify,
      defend, and hold each Contributor harmless for any liability
      incurred by, or claims asserted against, such Contributor by reason
      of your accepting any such warranty or additional liability.

   END OF TERMS AND CONDITIONS

   APPENDIX: How to apply the Apache License to your work.

      To apply the Apache License to your work, attach the following
      boilerplate notice, with the fields enclosed by brackets "[]"
      replaced with your own identifying information. (Don't include
      the brackets!)  The text should be enclosed in the appropriate
      comment syntax for the file format. We also recommend that a
      file or class name and description of purpose be included on the
      same "printed page" as the copyright notice for easier
      identification within third-party archives.

   Copyright [yyyy] [name of copyright owner]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.


 */
#include <cassert>
#include <cstddef>
#include <limits>
#include <cstdint>
#include <cstring>
#include <array>
#include <type_traits>
#include <immintrin.h>
#include <tmmintrin.h>


#if __cplusplus >= 201703L
#define AES_MAYBE_UNUSED [[maybe_unused]]
#else
#define AES_MAYBE_UNUSED
#endif

#ifndef TYPES_TEMPLATES
#define TYPES_TEMPLATES
namespace types {
    template<typename T>
    struct is_integral: std::false_type {};
    template<>struct is_integral<unsigned char>: std::true_type {};
    template<>struct is_integral<signed char>: std::true_type {};
    template<>struct is_integral<unsigned short>: std::true_type {};
    template<>struct is_integral<signed short>: std::true_type {};
    template<>struct is_integral<unsigned int>: std::true_type {};
    template<>struct is_integral<signed int>: std::true_type {};
    template<>struct is_integral<unsigned long>: std::true_type {};
    template<>struct is_integral<signed long>: std::true_type {};
    template<>struct is_integral<unsigned long long>: std::true_type {};
    template<>struct is_integral<signed long long>: std::true_type {};
#if __cplusplus >= 201703L
    template<class T> inline constexpr bool is_integral_v = is_integral<T>::value;
#endif

    template<typename T> struct is_simd: std::false_type {};
    template<typename T> struct is_simd_int: std::false_type {};
    template<typename T> struct is_simd_float: std::false_type {};

#if __SSE2__
    template<>struct is_simd<__m128i>: std::true_type {};
    template<>struct is_simd<__m128>:  std::true_type {};
    template<>struct is_simd_int<__m128i>: std::true_type {};
    template<>struct is_simd_float<__m128>: std::true_type {};
#endif
#if __AVX2__
    template<>struct is_simd<__m256i>: std::true_type {};
    template<>struct is_simd<__m256>:  std::true_type {};
    template<>struct is_simd_int<__m256i>: std::true_type {};
    template<>struct is_simd_float<__m256>: std::true_type {};
#endif
#if __AVX512__
    template<>struct is_simd<__m512i>: std::true_type {};
    template<>struct is_simd<__m512>:  std::true_type {};
    template<>struct is_simd_int<__m512i>: std::true_type {};
    template<>struct is_simd_float<__m512>: std::true_type {};
#endif
#if __cplusplus >= 201703L
    template<class T> inline constexpr bool is_simd_v = is_simd<T>::value;
    template<class T> inline constexpr bool is_simd_int_v = is_simd_int<T>::value;
    template<class T> inline constexpr bool is_simd_float_v = is_simd_float<T>::value;
#endif
} // namespace types
#endif



namespace aes {

using std::uint64_t;
using std::uint8_t;
using std::size_t;


#define AES_ROUND(rcon, index)                                                 \
  do {                                                                         \
    __m128i k2 = _mm_aeskeygenassist_si128(k, rcon);                           \
    k = _mm_xor_si128(k, _mm_slli_si128(k, 4));                                \
    k = _mm_xor_si128(k, _mm_slli_si128(k, 4));                                \
    k = _mm_xor_si128(k, _mm_slli_si128(k, 4));                                \
    k = _mm_xor_si128(k, _mm_shuffle_epi32(k2, _MM_SHUFFLE(3, 3, 3, 3)));      \
    seed_[index] = k;                                                    \
  } while (0)

#ifndef HAS_AVX_512
#  define HAS_AVX_512 (_FEATURE_AVX512F || _FEATURE_AVX512ER || _FEATURE_AVX512PF || _FEATURE_AVX512CD || __AVX512BW__ || __AVX512CD__ || __AVX512F__)
#endif
#if HAS_AVX_512
#define VEC_ALIGNMENT_FOR_BUFFER 64
#elif __AVX2__
#define VEC_ALIGNMENT_FOR_BUFFER 32
#else
#define VEC_ALIGNMENT_FOR_BUFFER 16
#endif



template<typename GeneratedType=uint64_t, size_t UNROLL_COUNT=4,
         typename=typename std::enable_if<
            types::is_integral<GeneratedType>::value || types::is_simd_int<GeneratedType>::value
            >::type
        >
class AesCtr {
    static const size_t AESCTR_ROUNDS = 10;
    uint8_t state_[sizeof(__m128i) * UNROLL_COUNT] __attribute__ ((aligned (VEC_ALIGNMENT_FOR_BUFFER)));
    __m128i ctr_[UNROLL_COUNT];
    __m128i seed_[AESCTR_ROUNDS + 1];
    __m128i work[UNROLL_COUNT];
    unsigned offset_;

    // Unrollers
    template<size_t ind, size_t todo>
    struct aes_unroll_impl {
        void operator()(__m128i *ret, AesCtr &state) const {
            ret[ind] = _mm_xor_si128(state.ctr_[ind], state.seed_[0]);
            aes_unroll_impl<ind + 1, todo - 1>()(ret, state);
        }
        void aesenc(__m128i *ret, __m128i subkey) const {
            ret[ind] = _mm_aesenc_si128(ret[ind], subkey);
            aes_unroll_impl<ind + 1, todo - 1>().aesenc(ret, subkey);
        }
        template<size_t NUMROLL>
        void round_and_enc(__m128i *ret, AesCtr &state) const {
            const __m128i subkey = state.seed_[ind];
            aes_unroll_impl<0, NUMROLL>().aesenc(ret, subkey);
            aes_unroll_impl<ind + 1, todo - 1>().template round_and_enc<NUMROLL>(ret, state);
        }
        void add_store(__m128i *work, AesCtr &state) const {
          state.ctr_[ind] =
              _mm_add_epi64(state.ctr_[ind], _mm_set_epi64x(0, UNROLL_COUNT));
              _mm_store_si128(
                  reinterpret_cast<__m128i *>(&state.state_[16 * ind]),
                  _mm_aesenclast_si128(work[ind], state.seed_[AESCTR_ROUNDS]));
          aes_unroll_impl<ind + 1, todo - 1>().add_store(work, state);
        }
    };
    // Termination conditions
    template<size_t ind>
    struct aes_unroll_impl<ind, 0> {
        void operator()(AES_MAYBE_UNUSED __m128i *ret, AES_MAYBE_UNUSED AesCtr &state) const {}
        void aesenc(AES_MAYBE_UNUSED __m128i *ret, AES_MAYBE_UNUSED __m128i subkey) const {}
        template<size_t NUMROLL>
        void round_and_enc(AES_MAYBE_UNUSED __m128i *ret, AES_MAYBE_UNUSED AesCtr &state) const {}
        void add_store(AES_MAYBE_UNUSED __m128i *work, AES_MAYBE_UNUSED AesCtr &state) const {}
    };

public:
    using result_type = GeneratedType;
    constexpr AesCtr(uint64_t seedval=0) {
        seed(seedval);
    }
    void generate_new_values() {
        aes_unroll_impl<0, UNROLL_COUNT>()(work, *this);
        aes_unroll_impl<1, AESCTR_ROUNDS - 1>().template round_and_enc<UNROLL_COUNT>(work, *this);
        aes_unroll_impl<0, UNROLL_COUNT>().add_store(work, *this);
        offset_ = 0;
    }
    result_type operator()() {
        if (__builtin_expect(offset_ >= sizeof(__m128i) * UNROLL_COUNT, 0))
            generate_new_values(); // sets offset_ to 0.
        result_type ret;
        std::memcpy(&ret, &state_[offset_], sizeof(ret));
        offset_ += sizeof(result_type);
        return ret;
    }
    static constexpr result_type max() {return std::numeric_limits<result_type>::max();}
    static constexpr result_type min() {return std::numeric_limits<result_type>::min();}
    void seed(uint64_t k) {
        seed(_mm_set_epi64x(0, k));
    }
    void seed(__m128i k) {
      seed_[0] = k;
      // D. Lemire manually unrolled following loop since _mm_aeskeygenassist_si128
      // requires immediates

      AES_ROUND(0x01, 1);
      AES_ROUND(0x02, 2);
      AES_ROUND(0x04, 3);
      AES_ROUND(0x08, 4);
      AES_ROUND(0x10, 5);
      AES_ROUND(0x20, 6);
      AES_ROUND(0x40, 7);
      AES_ROUND(0x80, 8);
      AES_ROUND(0x1b, 9);
      AES_ROUND(0x36, 10);

      for (unsigned i = 0; i < UNROLL_COUNT; ++i) ctr_[i] = _mm_set_epi64x(0, i);
      offset_ = sizeof(__m128i) * UNROLL_COUNT;
    }
    result_type operator[](size_t count) const {
        static constexpr unsigned DIV   = sizeof(__m128i) / sizeof(result_type);
        static constexpr unsigned BMASK = DIV - 1;
        const unsigned offset_(count & BMASK);
        result_type ret[DIV];
        count /= DIV;
        __m128i tmp(_mm_xor_si128(_mm_set_epi64x(0, count), seed_[0]));
        for (unsigned r = 1; r <= AESCTR_ROUNDS - 1; tmp = _mm_aesenc_si128(tmp, seed_[r++]));
        _mm_store_si128(reinterpret_cast<__m128i *>(ret), _mm_aesenclast_si128(tmp, seed_[AESCTR_ROUNDS]));
        return ret[offset_];
    }
    static constexpr size_t BUFSIZE = sizeof(state_);
    const uint8_t *buf() const {return &state_[0];}
    using ThisType = AesCtr<GeneratedType, UNROLL_COUNT>;

    template<typename T, bool manual_override=false,
             typename=typename std::enable_if<
                manual_override || types::is_integral<T>::value || types::is_simd_int<T>::value
                >::type
             >
    class buffer_view {
        ThisType &ref;
    public:
        buffer_view(ThisType &ctr): ref{ctr} {}
        using const_pointer = const T *;
        using pointer       = T *;
        const_pointer cbegin() const {
            return reinterpret_cast<const_pointer>(&ref.state_[0]);
        }
        const_pointer cend() const {
            return reinterpret_cast<const_pointer>(&ref.state_[BUFSIZE]);
        }
        pointer begin() {
            return reinterpret_cast<pointer>(&ref.state_[0]);
        }
        pointer end() {
            return reinterpret_cast<pointer>(&ref.state_[BUFSIZE]);
        }
    };
    template<typename T, bool manual_override=false>
    buffer_view<T, manual_override> view() {return buffer_view<T, manual_override>(*this);}
};
#undef AES_ROUND


template<typename size_type, size_t arrsize>
constexpr std::array<size_type, arrsize> seed_to_array(size_type seedseed) {
    std::array<size_type, arrsize> ret{};
    aes::AesCtr<size_type> gen(seedseed);
    for(auto &el: ret) el = gen();
    return ret;
}

template<typename T>
struct is_aes: std::false_type {};

template<typename T, size_t n>
struct is_aes<AesCtr<T, n>>: std::true_type {};

} // namespace aes

#undef AESCTR_UNROLL
#undef AESCTR_ROUNDS
#undef AES_MAYBE_UNUSED

#endif
