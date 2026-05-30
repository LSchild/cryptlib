//
// Created by leonard on 5/30/26.
//

#include "operators/lwe_to_rgsw_conversion.h"


[[nodiscard]] uint64_t SchemeSwitchingContext::GetOutputDigits() const {
    return m_output_digits;
}

[[nodiscard]] uint64_t SchemeSwitchingContext::GetOutputBasis() const {
    return m_output_basis;
}

[[nodiscard]] uint64_t SchemeSwitchingContext::GetOutputBasisLog2() const {
    return m_output_basis_log2;
}

void SchemeSwitchingContext::SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context) {
    m_rotation_context = new_rot_context;
}

void SchemeSwitchingContext::SetTraceContext(std::shared_ptr<TraceEvaluationContext> new_trace_context) {
    m_trace_params = new_trace_context;
}

void SchemeSwitchingContext::SetSquaringParameters(std::shared_ptr<RLWEConversionParameters> new_squaring_params) {
    m_squaring_params = new_squaring_params;
}

void SchemeSwitchingContext::SetOutputDigits(uint64_t new_digits) {
    m_output_digits = new_digits;
}

void SchemeSwitchingContext::SetOutputBasis(uint64_t new_basis) {
    m_output_basis = new_basis;
    m_output_basis_log2 = IntLog2(new_basis);
}

long double SchemeSwitchingContext::ComputeOutputVariance(long double input_variance) const {
    auto o_var = m_rotation_context->ComputeOutputVariance(input_variance);
    auto trace_var = m_trace_params->ComputeOutputVariance(o_var);
    auto square_var = m_squaring_params->ComputeOutputVariance(trace_var);

    return square_var;
}

[[nodiscard]] Container SchemeSwitchingContext::GetInputContainer() const  {
auto params_br = std::dynamic_pointer_cast<TupleContainerImpl>(m_rotation_context->GetInputContainer());
return params_br->GetElem(0);
}

[[nodiscard]]  Container SchemeSwitchingContext::GetOutputContainer(Container container) const  {

auto N = m_squaring_params->GetDimension();
auto Q = m_squaring_params->GetModulus();
auto var = ComputeOutputVariance();
auto dig = GetOutputDigits();
auto basis = GetOutputBasis();
return std::make_shared<RGSWContainerImpl>(N, Q, dig, basis, var);
}

const std::shared_ptr<OperatorContext<BlindRotator>> SchemeSwitchingContext::GetBlindRotationContext() const {
    return m_rotation_context;
}

const std::shared_ptr<TraceEvaluationContext> SchemeSwitchingContext::GetTraceContext() const {
    return m_trace_params;
}

const std::shared_ptr<RLWEConversionParameters> SchemeSwitchingContext::GetSquaringContext() const {

    return m_squaring_params;
}

std::unique_ptr<LWEtoRGSWConverter> SchemeSwitchingContext::ConstructOperator(const std::vector<GenericKey>& keys) const {

std::vector<GenericKey> trace_key(keys.begin() + 1, keys.end());

auto rotator = m_rotation_context->ConstructOperator(keys);
auto tracer = m_trace_params->ConstructOperator(trace_key);
auto squarer = m_squaring_params->ConstructOperator(trace_key);

std::unique_ptr<LWEtoRGSWConverter> op = std::unique_ptr<LWEtoRGSWConverter>(new LWEtoRGSWConverter(shared_from_this(), std::move(rotator), std::move(tracer), std::move(squarer)));

return op;
}

[[nodiscard]] OperatorID SchemeSwitchingContext::GetOperatorID() const  {
    return CONV_LWE_RGSW;
}



LWEtoRGSWConverter::LWEtoRGSWConverter(std::shared_ptr<const SchemeSwitchingContext> context, std::unique_ptr<BlindRotator> rotator,
        std::unique_ptr<TraceEvaluator> trace_eval, std::unique_ptr<RLWEtoRLWEConverter> squaring_eval) :
        m_params(context),
        m_rotator(std::move(rotator)),
        m_trace_eval(std::move(trace_eval)),
        m_square_eval(std::move(squaring_eval))
        {

    auto trace_context = m_params->GetTraceContext();
    const auto N = trace_context->GetDimension();
    const auto Q = trace_context->GetModulus();

    m_acc.resize(2 * N);
    m_extraction_buffer.resize(4 * N);

    const auto digits = m_params->GetOutputDigits();
    const auto basebits = m_params->GetOutputBasisLog2();

    const auto block_size = N / digits;
    uint64_t max_digits = std::ceil(std::log2((long double)Q)/(long double)basebits);

    uint64_t scale = 1 << (basebits * (max_digits - digits));

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

    m_params_set = true;
    m_keys_generated = true;

}

void LWEtoRGSWConverter::Convert(uint64_t* output, const uint64_t*const input)  {
    assert(m_keys_generated);
    // need to square the key
    auto square_params = m_square_eval->GetContext();
    auto N = square_params->GetDimension();
    auto Q = square_params->GetModulus();
    auto out_digits = m_params->GetOutputDigits();
    auto block_size = N / out_digits;
    auto ntt = square_params->GetNTT();

    // initial accumulator is perturbed already
    uint64_t * const lev_left = output;
    uint64_t * const lev_right = output + 2 * N * out_digits;

    std::copy(m_acc.data(), m_acc.data() + N, lev_right + N);
    m_rotator->BlindRotate(lev_right, input, m_acc.data(), true);

    intel::hexl::EltwiseAddMod(lev_right + N, lev_right + N, m_acc.data() + N, N, Q);



    // Lev extraction
    AlignedVector shift_poly(N, 0);
    shift_poly[N - block_size] = Q - 1;

    ntt->ComputeForward(shift_poly.data(),shift_poly.data(), 1, 1);

    uint64_t* buffer = m_extraction_buffer.data();

    uint64_t Ninv = intel::hexl::InverseMod(N, Q);
    intel::hexl::EltwiseFMAMod(lev_right,lev_right,Ninv, nullptr,2*N,Q,1);

    // TODO: the rotation can be done without ntt...
    // TODO: High prio
    for(uint32_t d_i = 1; d_i < out_digits; d_i++) {
        auto row_d_im1 = lev_right + (d_i - 1) * 2 * N;
        auto row_d_i = lev_right + d_i * 2 * N;
        intel::hexl::EltwiseMultMod(row_d_i, row_d_im1, shift_poly.data(), N, Q, 1);
        intel::hexl::EltwiseMultMod(row_d_i + N, row_d_im1 + N, shift_poly.data(), N, Q, 1);
    }

    for(uint32_t d_i = 0; d_i < out_digits; d_i++) {
        auto row_d_i = output + 2 * N * out_digits + d_i * 2 * N;
        ntt->ComputeInverse(row_d_i, row_d_i, 1, 1);
        ntt->ComputeInverse(row_d_i, row_d_i, 1, 1);
        m_trace_eval->Eval(row_d_i, row_d_i);
    }


    // now we just need to compute the RLWE(-s * m) from RLWE(m)
    ZERO_UINT64_ARR(buffer, 4 * N);

    for(uint32_t i = 0; i < out_digits; i++) {
        uint64_t* rlwe_m = lev_right + 2 * N * i;
        uint64_t* rlwe_sm = lev_left + 2 * N * i;

        // rlwe key switching expects [COEF, NTT] format
        // std::copy(rlwe_m, rlwe_m + N, output_buffer);
        ntt->ComputeInverse(buffer, rlwe_m, 1,1);
        // rlwe_sm = RLWE(-a * s^2)
        m_square_eval->Convert(rlwe_sm, buffer);
        // rlwe_sm = RLWE( a* s^2) = [a', b']
        intel::hexl::EltwiseSubMod(rlwe_sm, buffer + 2 * N, rlwe_sm, 2 * N, Q);
        // rlwe_sm = [a' + b, b'] = RLWE(b' - (a' + b) * s) = RLWE(a*s^2 - (a*s + m) * s) = RLWE(-s * m)
        intel::hexl::EltwiseAddMod(rlwe_sm, rlwe_m + N, rlwe_sm, N, Q);
    }
};

void LWEtoRGSWConverter::Convert(std::vector<uint64_t>& output, const std::vector<uint64_t>& input) {

    Convert(output.data(), input.data());
}



[[nodiscard]] const std::shared_ptr<const SchemeSwitchingContext> LWEtoRGSWConverter::GetContext() const  {
    return m_params;
}