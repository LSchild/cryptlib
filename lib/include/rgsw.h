//
// Created by leonard on 11/12/24.
//

#ifndef LARGE_FUNCTIONS_RGSW_H
#define LARGE_FUNCTIONS_RGSW_H

#include "binfhecontext.h"
#include "hexl/hexl.hpp"

using namespace lbcrypto;

class RingGSWSample {
public:

    RingGSWSample() = default;

    RingGSWSample(const std::shared_ptr<RingGSWCryptoParams>& params, std::vector<RLWECiphertext>& msm, std::vector<RLWECiphertext>& m);

    RingGSWSample(const RingGSWSample&);

    [[nodiscard]] RLWECiphertext mul(const RLWECiphertext & rhs) const;

    void mul_inplace(RLWECiphertext& rhs) const;

    RLWECiphertext poly_mul_with_precomp(std::vector<NativePoly>& precomp_digs);

    void cmux(RLWECiphertext& in, uint32_t idx) const;

    void cmux(RLWECiphertext& in, NativePoly& mon) const;

    RingGSWSample create_one();

    RingGSWSample flip_bit(RingGSWSample& one);

    std::vector<RLWECiphertext> MultiCMUX(std::vector<RLWECiphertext> &choice_0);

//private:
    std::shared_ptr<RingGSWCryptoParams> m_params;
    std::vector<RLWECiphertext> m_msm;
    std::vector<RLWECiphertext> m_m;
};

class RGSWSample {
public:
    RGSWSample() = default;

    RGSWSample(std::shared_ptr<intel::hexl::NTT>& ntt_engine, uint32_t L, uint32_t digits, const std::vector<RLWECiphertext> &msm, const std::vector<RLWECiphertext> &m);
    RGSWSample(std::shared_ptr<intel::hexl::NTT>& ntt_engine, uint32_t L, uint32_t digits, const uint64_t* msm, const uint64_t* m);

    static RGSWSample convert(const RingGSWSample& old_version);

    RLWECiphertext mul(const RLWECiphertext & rhs) ;
    RLWECiphertext mul_poly(const NativePoly& poly) ;

    /**
     * Important: we assume that the INTT has been applied to the RLWE operand
     *
     * Performs a RGSW-RLWE product (external product)
     * @param result output buffer of size 2 * N
     * @param rhs RLWE sample encoded in buffer of size 2 * N
     * @param add_to_result if true, we compute result += this * rhs, if false result = this * rhs
     */
    void mul(uint64_t * __restrict result, const uint64_t* __restrict rhs, bool add_to_result = false) ;

    void mul_poly(uint64_t* __restrict result, const uint64_t* __restrict poly, bool add_to_result = false) ;

    void cmux(RLWECiphertext& in, uint32_t idx);

    void cmux(uint64_t* in_out, uint32_t idx);

    void cmux(RLWECiphertext& in, NativePoly& mon);
    /** Given the (NTTed) RLWE sample RLWE(m) and polynomial X^a - 1 computes RLWE(m x X^{a * b}) where b is the bit encoded
     * in the current RGSW sample. The result is written back to the input.
     *
     * NOTE: This function is only correct if the RGSW sample contains a BIT
     *
     * @param in_out pointer to memory block of size 2 * N containing the RLWE sample
     * @param mon pointer to the polynomial mentioned above
     */
    void cmux(uint64_t* in_out, uint64_t* mon);

    void ternary_mux(uint64_t* in_out, uint64_t* mon_pos, uint64_t* mon_neg, RGSWSample& high_bit);

    void mul_dir(uint64_t* __restrict res, const uint64_t* __restrict poly, bool mul_left = false, bool add_to_result = false);

    std::shared_ptr<intel::hexl::NTT> m_engine;
    uint32_t m_basis;
    uint32_t m_basis_bits;
    uint32_t m_mask;
    uint32_t m_first_shift;


    uint64_t m_digits;

    std::vector<uint64_t> m_key;
    std::vector<uint64_t> m_scratch_space;

};

#endif //LARGE_FUNCTIONS_RGSW_H
