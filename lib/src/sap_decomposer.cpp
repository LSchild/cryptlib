//
// Created by leonard on 6/3/26.
//

#include "operators/sap_decomposition.h"
#include "utils/speed_utils.h"
#include "utils/math_utils.h"

SAPDecompositionContext::SAPDecompositionContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
std::shared_ptr<LWEtoRLWEPackingContext> packing_context,
        std::shared_ptr<LWEConversionParameters> lwe_conversion_context,
uint64_t default_radix) :
m_rotation_context(std::move(blind_rotation_context)),
m_packing_context(std::move(packing_context)),
m_conversion_context(std::move(lwe_conversion_context)), m_default_radix(default_radix) {

    m_restart_iter = ComputeRestartIterationForRadix(default_radix, IMPLICIT_SK_L0);
}

std::shared_ptr<LWEConversionParameters> SAPDecompositionContext::GetLWEConversionContext() const {
    return m_conversion_context;
}

std::shared_ptr<OperatorContext<BlindRotator>> SAPDecompositionContext::GetBlindRotationContext() const {
    return m_rotation_context;
}

std::shared_ptr<LWEtoRLWEPackingContext> SAPDecompositionContext::GetPackingContext() const {
    return m_packing_context;
}

Container SAPDecompositionContext::GetInputContainer() const {
    auto dummy = m_conversion_context->GetInputContainer();
    auto out = m_conversion_context->GetOutputContainer(dummy);
    auto out_lwe = std::dynamic_pointer_cast<LWEContainerImpl>(out);

    auto in_n = out_lwe->GetN();
    auto in_q = UINT64_MAX;

    return std::make_shared<LWEContainerImpl>(in_n, in_q, 0.0);
}

long double SAPDecompositionContext::ComputeOutputVariance(long double input_variance) const {

    auto N = m_packing_context->GetDimension();

    long double out_var = 0.0;

    long double cluster_var = 0.0;
    for(uint32_t i = 0; i < m_restart_iter; i++) {
        auto tmp_var = m_rotation_context->ComputeOutputVariance(out_var);
        cluster_var = tmp_var * N * std::pow(m_default_radix, 2.0) / 12.0;
        out_var = m_conversion_context->ComputeOutputVariance(cluster_var);
    }

    return cluster_var;
}

Container SAPDecompositionContext::GetOutputContainer(Container container) const {
    auto dummy_rlwe_in = m_rotation_context->GetInputContainer();
    auto dummy_rlwe_out = std::dynamic_pointer_cast<RLWEContainerImpl>(m_rotation_context->GetOutputContainer(dummy_rlwe_in));

    auto o_var = ComputeOutputVariance(0.0);

    return std::make_shared<LWEContainerImpl>(dummy_rlwe_out->GetN(), dummy_rlwe_out->GetQ(), o_var);
}

const uint64_t SAPDecompositionContext::GetDefaultRadix() const {
    return m_default_radix;
}

const uint64_t SAPDecompositionContext::GetRestartIteration() const {
    return m_restart_iter;
}

const uint64_t SAPDecompositionContext::ComputeRestartIterationForRadix(uint64_t radix, uint64_t sk_hamming) const {

    long double log_max_p_fail = -40;
    long double current_var = 0;

    auto Q = m_packing_context->GetModulus();
    auto N = m_packing_context->GetDimension();

    auto iter = 0;
    while (true) {

        // compute iteration variance
        current_var = m_rotation_context->ComputeOutputVariance(current_var);
        current_var = m_packing_context->ComputeOutputVariance(current_var);

        auto extraction_var = (current_var * N * radix*radix) / 12.0;
        auto p_fail = EstimateFailureProbability(extraction_var, Q, radix);

        // we apply log next so we have to make sure we're within range
        // if we're at or below 0, we're fine
        if (p_fail <= 0.0) {
            iter++;
            continue;
        }

        auto log_p_fail= std::log2l(p_fail);
        if (log_p_fail <= log_max_p_fail) {
            iter++;
            continue;
        } else {
            break;
        }
    }

    auto reset_var = m_conversion_context->ComputeOutputVariance(current_var);
    auto lwe_var = std::pow((long double)Q, 2.0) * reset_var / std::pow((long double)(2 * N), 2.0);
    lwe_var += (long double)(std::pow(sk_hamming, 2.0) + 1) / 3.0;

    auto lwe_p_fail = EstimateFailureProbability(lwe_var, 2 * N, radix);
    if (lwe_p_fail <= 0.0)
        return iter;

    auto log_lwe_p_fail = std::log2(lwe_p_fail);
    if (log_lwe_p_fail <= log_max_p_fail) {
        return iter;
    } else {
        std::cerr << __FUNCTION__ << " Warning: LWE KeySwitch for restart induces large variance. Check your hamming weight & other params." << std::endl;
        std::cerr << __FUNCTION__ << " Warning: This is probably fine though." << std::endl;
        return iter - 1;
     }
}

OperatorID SAPDecompositionContext::GetOperatorID() const {
    return DECOMP_SAP;
}

void SAPDecompositionContext::SetDefaultRadix(uint64_t new_radix) {
    m_default_radix = new_radix;
    // TODO update restart iter
}

void SAPDecompositionContext::SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context) {
    m_rotation_context = new_rot_context;
}

void SAPDecompositionContext::SetPackingContext(std::shared_ptr<LWEtoRLWEPackingContext> new_trace_context) {
    m_packing_context = new_trace_context;
}

void SAPDecompositionContext::SetLWEConversionContext(std::shared_ptr<LWEConversionParameters> new_lwe_context) {
    m_conversion_context = new_lwe_context;
}

std::unique_ptr<SAPDecomposer> SAPDecompositionContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    // TODO
}

void SAPDecomposer::HomTrunc(uint64_t *__restrict output, uint64_t *__restrict input, uint64_t radix, uint64_t block_limit) {
    // Assume input is given in COEF form
    auto packing_cont = m_context->GetPackingContext();
    auto N = packing_cont->GetDimension();
    auto Q = packing_cont->GetModulus();
    auto packing_p = m_packing_buffer.data();
    ZERO_UINT64_ARR(packing_p, m_packing_buffer.size());

    AlignedVector a_rev(2 * N, 0);
    auto a_rev_p = a_rev.data();
    // evaluate reversal automorphism first
    a_rev[N] = input[N];
    for(uint64_t idx = 1; idx < N; idx++) {
        a_rev[N + idx] = intel::hexl::SubUIntMod(0, input[N - idx], Q);
    }
    // mirror so now a_rev = [-Auto(input; -1), Auto(input; -1)]
    // more specifically for any k < N, a_rev[N-k:2*N-k] = Auto(input * X^{-k}, -1)
    intel::hexl::EltwiseSubMod(a_rev.data(), a_rev.data(), a_rev.data() + N, N, Q);

    // todo: bench
    for(uint32_t block_idx = 0; block_idx < block_limit; block_idx++) {
        auto p_i = packing_p + block_idx * (N + 1);
        std::copy(input + block_idx * radix, input + (block_idx + 1) * radix, p_i);
        p_i[N] = input[block_idx];
        for(uint64_t b_i = 1; b_i < radix; b_i++) {
            auto shift_idx = block_idx * radix + b_i;
            intel::hexl::EltwiseAddMod(p_i, p_i, a_rev_p + shift_idx, N, Q);
            p_i[N] = intel::hexl::AddUIntMod(p_i[N], input[shift_idx], Q);
        }
    }

    // packing buffer should now be an array of block-wise added coefs
    // we now re-pack consecutively
    m_packer->PackConsecutively(output, packing_p, block_limit * radix);
}

void SAPDecomposer::ResetAccumulator(uint64_t *acc) {
    // TODO
}

void SAPDecomposer::Decompose(uint64_t *output, const uint64_t *const input, uint64_t radix) {

    // We extract expected input shapes
    auto br_input_container = m_rotator->GetContext()->GetInputContainer();
    auto br_input_tuple = std::dynamic_pointer_cast<TupleContainerImpl>(br_input_container);
    // LWE input
    auto lwe_container = std::dynamic_pointer_cast<LWEContainerImpl>(br_input_tuple->GetElem(0));
    auto lwe_dim_in = lwe_container->GetN();
    auto lwe_mod_in = lwe_container->GetQ();
    // RLWE input
    auto rlwe_container = std::dynamic_pointer_cast<RLWEContainerImpl>(br_input_tuple->GetElem(1));
    auto rlwe_mod = rlwe_container->GetQ();
    auto rlwe_dim_in = rlwe_container->GetN();
    // NTT
    auto ntt = m_context->GetPackingContext()->GetNTT();

    // set up constants
    auto perturb_lo = radix * m_beta;
    auto perturb_hi = m_beta;
    auto radix_log2 = IntLog2(radix);
    // assumes radix is power of two
    auto radix_mask = radix - 1;

    // buffer for accumulator and extraction poly
    AlignedVector rlwe_scratch(0, 5 * rlwe_dim_in);
    auto acc_p = rlwe_scratch.data();
    // acc = Q / B * 1
    acc_p[0] = rlwe_mod >> radix_log2;
    // extraction poly encodes map x -> x mod radix, encoded in reverse
    auto extract_poly = rlwe_scratch.data() + 2 * rlwe_dim_in;
    extract_poly[0] = 0;
    for(uint64_t i = 1; i < rlwe_dim_in; i++) {
        extract_poly[rlwe_dim_in - i] = i & radix_mask;
    }
    ntt->ComputeForward(extract_poly, extract_poly, 1, 1);
    // extraction output buffer
    auto ext_buffer = rlwe_scratch.data() + 3 * rlwe_dim_in;

    // buffers for current phase digit and updated input
    std::vector<uint64_t> scratch(0, 2 * (lwe_dim_in + 1));
    auto input_copy = scratch.data();
    auto current_sub_phase = scratch.data() + lwe_dim_in + 1;
    std::copy(input, input + lwe_dim_in + 1, input_copy);

    // dynamically determine input modulus
    auto current_modulus = 1ull << IntLog2(input[0]);
    for(uint64_t i = 1; i <= lwe_dim_in; i++) {
        if (current_modulus <= input[i]) [[unlikely]] {
            current_modulus <<= 1;
        }
    }

    auto restart_iter = radix == m_max_radix ? m_restart_iteration : m_context->ComputeRestartIterationForRadix(radix, m_beta);
    auto current_output_digit = 0;
    while (current_modulus > 0) {
        /* start by obtaining current phase digit */
        for(uint64_t i = 0; i < lwe_dim_in + 1; i++) {
            current_sub_phase[i] = input_copy[i] & radix_mask;
            input_copy[i] >>= radix_log2;
        }
        current_modulus >>= radix_log2;

        // perturb as required
        current_sub_phase[lwe_dim_in] = intel::hexl::AddUIntMod(current_sub_phase[lwe_dim_in], perturb_lo, lwe_mod_in);
        input_copy[lwe_dim_in] = intel::hexl::SubUIntMod(input_copy[lwe_dim_in], perturb_hi, current_modulus);

        m_rotator->BlindRotate(nullptr, current_sub_phase, acc_p);

        // apply extraction
        intel::hexl::EltwiseMultMod(ext_buffer, acc_p, extract_poly, rlwe_dim_in, rlwe_mod, 1);
        intel::hexl::EltwiseMultMod(ext_buffer + rlwe_dim_in, acc_p + rlwe_dim_in,extract_poly, rlwe_dim_in, rlwe_mod, 1);
        ntt->ComputeInverse(ext_buffer, ext_buffer, 1, 1);
        ntt->ComputeInverse(ext_buffer + rlwe_dim_in, ext_buffer + rlwe_dim_in, 1, 1);

        // actually extract
        auto current_digit_buf = output + current_output_digit * (rlwe_dim_in + 1);
        current_digit_buf[0] = ext_buffer[0];
        current_digit_buf[rlwe_dim_in] = ext_buffer[rlwe_dim_in];
        std::reverse_copy(ext_buffer + 1, ext_buffer + rlwe_dim_in, current_digit_buf + 1);
        current_output_digit++;

        // truncate
        ntt->ComputeInverse(ext_buffer, acc_p, 1, 1);
        ntt->ComputeInverse(ext_buffer + rlwe_dim_in, acc_p + rlwe_dim_in, 1, 1);

        if ((current_output_digit + 1) % restart_iter == 0 and current_modulus != 0) {
            ResetAccumulator(ext_buffer);
        }
        // todo: is m_beta ok ?
        HomTrunc(acc_p, ext_buffer, radix, m_beta);
    }
}

void SAPDecomposer::Decompose(std::vector<uint64_t> &output, const std::vector<uint64_t> &input, uint64_t radix) {
    Decompose(output.data(), input.data(), radix);
}

const std::shared_ptr<SAPDecompositionContext> SAPDecomposer::GetContext() const {
    return m_context;
}

SAPDecomposer::SAPDecomposer(std::shared_ptr<SAPDecompositionContext> ctx, std::unique_ptr<BlindRotator> rotator,
                             std::unique_ptr<LWEtoRLWEPacker> packer, std::unique_ptr<LWEtoLWEConverter> conv, uint64_t max_radix,
                             uint64_t lwe_sk_hamming_weight, uint64_t reset_period) : m_context(std::move(ctx)),
                             m_rotator(std::move(rotator)), m_packer(std::move(packer)), m_lwe_converter(std::move(conv)),
                             m_max_radix(max_radix), m_beta(lwe_sk_hamming_weight), m_restart_iteration(reset_period) {
    auto N = ctx->GetPackingContext()->GetDimension();
    auto max_block = N / max_radix;
    m_packing_buffer.resize(max_block * (N + 1));
}