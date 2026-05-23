//
// Created by leonard on 4/28/26.
//

#ifndef LARGE_FUNCTIONS_RGSW_CONVERSION_H
#define LARGE_FUNCTIONS_RGSW_CONVERSION_H

#include "container_types.h"
#include "interfaces.h"
#include "glwe_conversion.h"
#include "automorphism_key.h"
#include "math_utils.h"
#include "speed_utils.h"

template<typename BlindRotationParameterType>
struct RGSWConversionParams : OperationParameters {

    RGSWConversionParams(BlindRotationParameterType br_params,
                         AutomorphismParameters auto_params,
                         RLWEConversionParameters squaring_params, uint64_t digits, uint64_t basis) : m_rot_params(br_params),
                                                                                                      m_auto_params(auto_params),
                                                                                                      m_squaring_params(squaring_params),
                                                                                                      m_output_digits(digits), m_output_basis(basis),
                                                                                                      m_output_basis_log2(IntLog2(basis)){
    }


    long double ComputeOutputVariance(long double input_variance = 0) const override {
        // todo
        auto test = static_cast<OperationParameters>(m_rot_params).ComputeOutputVariance(0.0);
        return test;
    }

    [[nodiscard]] const BlindRotationParameterType& GetBlindRotationParameters() const {
        return m_rot_params;
    }

    [[nodiscard]] const AutomorphismParameters& GetAutomorphismParameters() const {
        return m_auto_params;
    }

    [[nodiscard]] const RLWEConversionParameters& GetSquaringParameters() const {
        return m_squaring_params;
    }

    [[nodiscard]] uint64_t GetOutputDigits() const {
        return m_output_digits;
    }

    [[nodiscard]] uint64_t GetOutputBasis() const {
        return m_output_basis;
    }

    [[nodiscard]] uint64_t GetOutputBasisLog2() const {
        return m_output_basis_log2;
    }

    void SetBlindRotationParameters(BlindRotationParameterType new_rot_params) {
        m_rot_params = new_rot_params;
    }

    void SetAutomorphismParameters(const AutomorphismParameters& new_auto_params) {
        m_auto_params = new_auto_params;
    }

    void SetSquaringParameters(const RLWEConversionParameters& new_squaring_params) {
        m_squaring_params = new_squaring_params;
    }

    void SetOutputDigits(uint64_t new_digits) {
        m_output_digits = new_digits;
    }

    void SetOutputBasis(uint64_t new_basis) {
        m_output_basis = new_basis;
        m_output_basis_log2 = IntLog2(new_basis);
    }


    BlindRotationParameterType m_rot_params;
    AutomorphismParameters m_auto_params;
    RLWEConversionParameters m_squaring_params;

    uint64_t m_output_digits;
    uint64_t m_output_basis;
    uint64_t m_output_basis_log2;

};

template<typename BRParamType, typename BRType>
struct LWEtoRGSWConverter : SchemeConverter<RGSWConversionParams<BRParamType>> {

    LWEtoRGSWConverter(RGSWConversionParams<BRParamType> params) : m_params(params) {

        // TODO: should we get rid of this assertion
        static_assert(std::is_base_of<BlindRotator<BRParamType>, BRType>::value, "Incompatible blind-rotation parameter and blind-rotation type");

        m_rotator = std::make_shared<BRType>(params.GetBlindRotationParameters());
        m_squaring_converter = std::make_shared<RLWEtoRLWEConverter>(params.GetSquaringParameters());

        auto output = m_rotator->GetOutputContainer();
        auto output_rlwe = std::dynamic_pointer_cast<RLWEContainerImpl>(output);

        for(uint32_t i = 2; i <= output_rlwe->getN(); i *= 2) {
            AutomorphismParameters auto_params = m_params.GetAutomorphismParameters();
            auto_params.SetAutomorphismIndex(i + 1);
            m_auto_converters.push_back(std::make_shared<AutomorphismEvaluator>(auto_params));
        }
        m_params_set = true;

    }

    void SetParams(RGSWConversionParams<BRParamType>& params) override {
        m_params = params;
        m_rotator = std::make_shared<BRType>(params.GetBlindRotationParameters());
        m_squaring_converter = std::make_shared<RLWEtoRLWEConverter>(params.GetSquaringParameters());

        auto output = m_rotator->GetOutputContainer();
        auto output_rlwe = std::dynamic_pointer_cast<RLWEContainerImpl>(output);

        AutomorphismParameters auto_params = m_params.GetAutomorphismParameters();
        for(uint32_t i = 2; i <= output_rlwe->getN(); i *= 2) {
            auto_params.SetAutomorphismIndex(i + 1);
            m_auto_converters.push_back(std::make_shared<AutomorphismEvaluator>(auto_params));
        }
        m_params_set = true;

    }

    void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) override {
        assert(m_params_set);
        m_rotator->KeyGen(source_key, target_key);

        // need to square the key
        auto N = m_params.GetSquaringParameters().GetDimension();
        auto Q = m_params.GetSquaringParameters().GetModulus();

        m_acc.resize(2 * N);
        m_extraction_buffer.resize(4 * N);

        AlignedVector tmp_keys(2 * N, 0);
        auto ntt = m_params.GetSquaringParameters().GetNTT();
        ntt->ComputeForward(tmp_keys.data() + N, target_key, 1, 1);
        intel::hexl::EltwiseMultMod(tmp_keys.data(), tmp_keys.data() + N, tmp_keys.data() + N, N,Q ,1);
        ntt->ComputeInverse(tmp_keys.data(), tmp_keys.data(), 1, 1);

        m_squaring_converter->KeyGen(tmp_keys.data(), target_key);
        /* (a,b ) -> (a, 0) -> Switch-> Negate -> (-a', -b'= -a' * s + a*s^2) -> (-a' + b, -b') -> - s * m
         *
         *
         */
        for(auto& conv : m_auto_converters) {
            conv->KeyGen(target_key);
        }


        // assumption is input is a LWE of q/2 * b + e
        // max |error| for q = 2N / 4block_size = N/2block_size
        const auto block_size = N / (m_params.GetOutputDigits());
        const auto digits = m_params.GetOutputDigits();
        const auto basebits = m_params.GetOutputBasisLog2();
        uint64_t max_digits = std::ceil(std::log2((long double)Q)/(long double)basebits);

        uint64_t scale = 1 << (m_params.GetOutputBasisLog2() * (max_digits - digits));

        // we shift by -N/(2 * block_size) so that it is centered properly
        // set first digit

        // NOTE: Q - 2 instead of two since we're fusing the multiplication by -1 here
        // NOTE: mult. by -1 is necessary to map 1->L^i/2 instead of 1-> -L^i/2
        uint64_t two_inverse = intel::hexl::InverseMod(2, Q);
        for(uint32_t i = 0; i < block_size/2; i++) {
            m_acc[N + i] = scale ;
            m_acc[N + N - i - 1] = (Q - scale);
        }

        for(uint32_t block_idx = 1; block_idx < digits; block_idx++) {
            scale <<= basebits;

            for(uint32_t j = 0; j < block_size; j++) {
                m_acc[N + block_idx * block_size + j - (block_size / 2)] = scale;
            }
        }

        intel::hexl::EltwiseFMAMod(m_acc.data() + N, m_acc.data() + N, two_inverse, nullptr, N, Q,1);
        intel::hexl::EltwiseFMAMod(m_acc.data(), m_acc.data() + N, Q - 1, nullptr, N, Q,1);


        ntt->ComputeForward(m_acc.data() + N, m_acc.data() + N, 1, 1);
        ntt->ComputeForward(m_acc.data(), m_acc.data(), 1, 1);

        m_keys_generated = true;
    }

    void KeyGen(const std::vector<uint64_t>& source_key, const std::vector<uint64_t>&  target_key) override {
        KeyGen(source_key.data(), target_key.data());
    }

    virtual void Convert(uint64_t* output, const uint64_t*const input) override {
        assert(m_keys_generated);
        // need to square the key
        auto square_params = m_params.GetSquaringParameters();
        auto N = square_params.GetDimension();
        auto Q = square_params.GetModulus();
        auto out_digits = m_params.GetOutputDigits();
        auto block_size = N / out_digits;
        auto ntt = square_params.GetNTT();

        // initial accumulator is perturbed already
        uint64_t * const lev_left = output;
        uint64_t * const lev_right = output + 2 * N * out_digits;

        std::copy(m_acc.data(), m_acc.data() + N, lev_right + N);
        m_rotator->BlindRotate(input, lev_right);

        intel::hexl::EltwiseAddMod(lev_right + N, lev_right + N, m_acc.data() + N, N, Q);

        // Lev extraction
        AlignedVector shift_poly(N, 0);
        shift_poly[N - block_size] = Q - 1;

        ntt->ComputeForward(shift_poly.data(),shift_poly.data(), 1, 1);

        uint64_t* buffer = m_extraction_buffer.data();

        uint64_t Ninv = intel::hexl::InverseMod(N, Q);
        intel::hexl::EltwiseFMAMod(lev_right,lev_right,Ninv, nullptr,2*N,Q,1);

        for(uint32_t d_i = 1; d_i < out_digits; d_i++) {
            auto row_d_im1 = lev_right + (d_i - 1) * 2 * N;
            auto row_d_i = lev_right + d_i * 2 * N;
            intel::hexl::EltwiseMultMod(row_d_i, row_d_im1, shift_poly.data(), N, Q, 1);
            intel::hexl::EltwiseMultMod(row_d_i + N, row_d_im1 + N, shift_poly.data(), N, Q, 1);
        }

        for(uint32_t d_i = 0; d_i < out_digits; d_i++) {
            auto row_d_i = output + 2 * N * out_digits + d_i * 2 * N;
            for (auto key = m_auto_converters.rbegin(); key != m_auto_converters.rend(); key++) {
                ntt->ComputeInverse(buffer, row_d_i, 1, 1);
                ntt->ComputeInverse(buffer + N, row_d_i + N, 1,1);
                ZERO_UINT64_ARR(buffer + 2 * N, 2 * N);
                (*key)->Eval(buffer + 2 * N, buffer);
                intel::hexl::EltwiseAddMod(row_d_i, buffer + 2 * N, row_d_i, 2 * N, Q);
            }
        }


        // now we just need to compute the RLWE(-s * m) from RLWE(m)
        ZERO_UINT64_ARR(buffer, 4 * N);

        for(uint32_t i = 0; i < out_digits; i++) {
            uint64_t* rlwe_m = lev_right + 2 * N * i;
            uint64_t* rlwe_sm = lev_left + 2 * N * i;

            // rlwe key switching expects [COEF, NTT] format
            // std::copy(rlwe_m, rlwe_m + N, output_buffer);
            ntt->ComputeInverse(buffer + N, rlwe_m, 1,1);
            // rlwe_sm = RLWE(-a * s^2)
            m_squaring_converter->Convert(rlwe_sm, buffer);
            // rlwe_sm = RLWE( a* s^2) = [a', b']
            intel::hexl::EltwiseSubMod(rlwe_sm, buffer + 2 * N, rlwe_sm, 2 * N, Q);
            // rlwe_sm = [a' + b, b'] = RLWE(b' - (a' + b) * s) = RLWE(a*s^2 - (a*s + m) * s) = RLWE(-s * m)
            intel::hexl::EltwiseAddMod(rlwe_sm, rlwe_m + N, rlwe_sm, N, Q);
        }
    };

    void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) override {
        Convert(output.data(), input.data());
    }

    [[nodiscard]] Container GetInputContainer() const override {
        auto params_br = std::dynamic_pointer_cast<TupleContainerImpl>(m_rotator->GetInputContainer());
        return params_br->GetElem(0);
    }

    [[nodiscard]]  virtual Container GetOutputContainer() const override {
        auto params_br = std::dynamic_pointer_cast<TupleContainerImpl>(m_rotator->GetInputContainer());

        auto tmp = params_br->GetElem(1);

        RLWEContainer params_rlwe = std::dynamic_pointer_cast<RLWEContainerImpl>(tmp);
        auto var = m_params.ComputeOutputVariance();
        auto dig = m_params.GetOutputDigits();
        auto basis = m_params.GetOutputBasis();
        return std::make_shared<RGSWContainerImpl>(params_rlwe->getQ(), params_rlwe->getN(), dig, basis, var);
    }

    const RGSWConversionParams<BRParamType>& GetParams() const override {
        return m_params;
    }


    bool m_params_set = false;
    bool m_keys_generated = false;

    RGSWConversionParams<BRParamType> m_params;

    std::shared_ptr<BlindRotator<BRParamType>> m_rotator;
    std::shared_ptr<RLWEtoRLWEConverter> m_squaring_converter;
    std::vector<std::shared_ptr<AutomorphismEvaluator>> m_auto_converters;

    AlignedVector m_acc;
    AlignedVector m_extraction_buffer;
};

#endif //LARGE_FUNCTIONS_RGSW_CONVERSION_H
