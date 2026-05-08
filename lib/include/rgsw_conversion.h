//
// Created by leonard on 4/28/26.
//

#ifndef LARGE_FUNCTIONS_RGSW_CONVERSION_H
#define LARGE_FUNCTIONS_RGSW_CONVERSION_H

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

    LWEtoRGSWConverter(RGSWConversionParams<BRParamType> params) : m_rotator(params.GetBlindRotationParameters()),
                                                                   m_squaring_converter(params.GetSquaringParameters()) {

        static_assert(std::is_base_of<BlindRotator<BRParamType>, BRType>::value, "Incompatible blind-rotation parameter and blind-rotation type");

        auto output = m_rotator->GetOutputContainer();
        auto output_rlwe = dynamic_cast<RLWEContainer>(output);

        for(uint32_t i = 2; i <= output_rlwe->getN(); i *= 2) {
            AutomorphismParameters auto_params = m_params.GetAutomorphismParameters();
            auto_params.SetAutomorphismIndex(i + 1);
            m_auto_converters.push_back(std::make_shared<AutomorphismEvaluator>(auto_params));
        }

    }

    void SetParams(RGSWConversionParams<BRParamType>& params) override {
        m_params = params;
        m_rotator = std::make_shared<BRType>(params.GetBlindRotationParameters());
        m_squaring_converter = std::make_shared<RLWEtoRLWEConverter>(params.GetSquaringParameters());

        auto output = m_rotator->GetOutputContainer();
        auto output_rlwe = dynamic_cast<RLWEContainer>(output);

        for(uint32_t i = 2; i <= output_rlwe->getN(); i *= 2) {
            AutomorphismParameters auto_params = m_params.GetAutomorphismParameters();
            auto_params.SetAutomorphismIndex(i + 1);
            m_auto_converters.push_back(std::make_shared<AutomorphismEvaluator>(auto_params));
        }
    }

    void KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) override {
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
        for(uint32_t i = 2; i <= N; i *= 2) {
            m_auto_converters[i]->KeyGen(target_key);
        }

        // assumption is input is a LWE of q/2 * b + e
        // max |error| for q = 2N / 4block_size = N/2block_size
        const auto block_size = N / (m_params.GetOutputDigits());
        const auto digits = m_params.GetOutputDigits();
        const auto basebits = m_params.GetOutputBasisLog2();
        uint64_t max_digits = std::ceil(std::log2((long double)Q)/(long double)basebits);

        auto scale = 1 << (m_params.GetOutputBasisLog2() * (max_digits - digits));

        // we shift by -N/(2 * block_size) so that it is centered properly
        // set first digit
        for(uint32_t i = 0; i < block_size/2; i++) {
            m_acc[i] = scale;
            m_acc[N - i - 1] = (Q - scale);
        }

        for(uint32_t block_idx = 1; block_idx < digits; block_idx+=block_size) {
            for(uint32_t j = 0; j < block_size; j++) {
                m_acc[N + block_idx + j - (block_idx / 2)] = scale;
            }
            scale <<= basebits;
        }
    }

    void KeyGen(const std::vector<uint64_t>& source_key, const std::vector<uint64_t>&  target_key) override {
        KeyGen(source_key.data(), target_key.data());
    }

    virtual void Convert(uint64_t* output, const uint64_t*const input) override {
        // need to square the key
        auto square_params = m_params.GetSquaringParameters();
        auto N = square_params.GetDimension();
        auto Q = square_params.GetModulus();
        auto out_digits = m_params.GetOutputDigits();
        auto block_size = N / out_digits;
        auto ntt = square_params.GetNTT();

        // initial accumulator is perturbed already
        std::copy(m_acc.data(), m_acc.data() + 2 * N, output);
        m_rotator->BlindRotate(input, output);

        // Lev extraction
        AlignedVector shift_poly(N, 0);
        shift_poly[N - block_size] = Q - 1;
        ntt->ComputeForward(shift_poly.data(),shift_poly.data(), 1, 1);

        uint64_t* output_buffer = m_extraction_buffer.data();

        for(uint32_t i = 0; i < out_digits; i++) {

            uint64_t* input_buffer = output + 4 * N * out_digits - (i + 1) * 2 * N;
            std::copy(output, output + 2 * N, input_buffer);

            // TODO: premul by N^{-1}
            for(auto key : m_auto_converters) {
                // auto eval input expected in [COEF, COEF] format
                ntt->ComputeInverse(input_buffer, input_buffer, 1,1);
                ntt->ComputeInverse(input_buffer + N, input_buffer + N, 1, 1);
                key->Eval(output_buffer, input_buffer);
                intel::hexl::EltwiseAddMod(input_buffer, input_buffer, output_buffer, 2 * N, Q);
            }

            intel::hexl::EltwiseMultMod(output, output, shift_poly.data(), N, Q, 1);
            intel::hexl::EltwiseMultMod(output + N, output + N, shift_poly.data(), N, Q, 1);
        }

        // now we just need to compute the RLWE(-s * m) from RLWE(m)
        for(uint32_t i = 0; i < out_digits; i++) {
            uint64_t* rlwe_m = output + 2 * N * out_digits + 2 * N * i;
            uint64_t* rlwe_sm = output + 2 * N *i;

            // rlwe key switching expects [COEF, NTT] format
            std::copy(rlwe_m, rlwe_m + N, output_buffer);
            ntt->ComputeInverse(output_buffer, output_buffer, 1,1);
            // rlwe_sm = RLWE(-a * s^2)
            m_squaring_converter->Convert(rlwe_sm, output_buffer);
            ZERO_UINT64_ARR(output_buffer, 2 * N);
            // rlwe_sm = RLWE( a* s^2) = [a', b']
            intel::hexl::EltwiseSubMod(rlwe_sm, output_buffer, rlwe_sm, 2 * N, Q);
            // rlwe_sm = [a' + b, b'] = RLWE(b' - (a' + b') * s) = RLWE(a*s^2 - (a*s + m) * s) = RLWE(-s * m)
            intel::hexl::EltwiseAddMod(rlwe_sm, rlwe_m + N, N, Q);
        }
    };

    void Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) override {
        Convert(output.data(), input.data());
    }

    [[nodiscard]] Container GetInputContainer() const override {
        auto params_br = dynamic_cast<TupleContainer>(m_rotator->GetInputContainer());
        return params_br->GetElem(0);
    }

    [[nodiscard]]  virtual Container GetOutputContainer() const override {
        auto params_br = dynamic_cast<TupleContainer>(m_rotator->GetInputContainer());
        RLWEContainer params_rlwe = dynamic_cast<RLWEContainer>(params_br->GetElem(1));
        auto var = m_params.ComputeOutputVariance();
        auto dig = m_params.GetOutputDigits();
        auto basis = m_params.GetOutputBasis();
        return std::make_shared<RGSWContainerImpl>(params_rlwe->getQ(), params_rlwe->getN(), dig, basis, var);
    }

    const RGSWConversionParams<BRParamType>& GetParams() const {
        return m_params;
    }


    RGSWConversionParams<BRParamType> m_params;

    std::shared_ptr<BlindRotator<BRParamType>> m_rotator;
    std::shared_ptr<RLWEtoRLWEConverter> m_squaring_converter;
    std::vector<std::shared_ptr<AutomorphismEvaluator>> m_auto_converters;

    AlignedVector m_acc;
    AlignedVector m_extraction_buffer;
};

#endif //LARGE_FUNCTIONS_RGSW_CONVERSION_H
