//
// Created by leonard on 6/2/26.
//

#include "operators/lwe_to_rlwe_packing.h"

#include "utils/speed_utils.h"
#include "utils/math_utils.h"

LWEtoRLWEPackingContext::LWEtoRLWEPackingContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N,
                                               uint64_t basis, uint64_t digits, double std) {
    m_auto_context = std::make_shared<AutomorphismContext>(source_key_distribution, modulus, N, basis, digits, std, 3);

}

LWEtoRLWEPackingContext::LWEtoRLWEPackingContext(KeyDistribution source_key_distribution,
                                               std::shared_ptr<intel::hexl::NTT> ntt, uint64_t basis, uint64_t digits,
                                               double std) {
    m_auto_context = std::make_shared<AutomorphismContext>(source_key_distribution, ntt, basis, digits, std, 3);
}

LWEtoRLWEPackingContext::LWEtoRLWEPackingContext(const LWEtoRLWEPackingContext &other) {
    m_auto_context = std::make_shared<AutomorphismContext>(other.GetSourceKeyDistribution(), other.GetNTT(), other.GetGadgetBasis(), other.GetGadgetDigits(), other.GetStd(), 3);
}

Container LWEtoRLWEPackingContext::GetInputContainer() const {
    return std::make_shared<RLWEContainerImpl>(m_auto_context->GetDimension(), m_auto_context->GetModulus(), 0.0);
}

long double LWEtoRLWEPackingContext::ComputeOutputVariance(long double input_variance) const {
    long double N = m_auto_context->GetDimension();
    auto auto_result_rlwe = std::dynamic_pointer_cast<RLWEContainerImpl>(m_auto_context->GetInputContainer());
    return auto_result_rlwe->GetVariance() * (std::pow(N, 2) / 3);

}

Container LWEtoRLWEPackingContext::GetOutputContainer(Container input) const {
    long double N = m_auto_context->GetDimension();
    auto input_rlwe =  std::dynamic_pointer_cast<RLWEContainerImpl>(m_auto_context->GetOutputContainer(input));
    // todo: implement the better one
    auto o_var = ComputeOutputVariance(input_rlwe->GetVariance());
    return std::make_shared<RLWEContainerImpl>(N, m_auto_context->GetModulus(), o_var);
}

std::unique_ptr<LWEtoRLWEPacker> LWEtoRLWEPackingContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    std::vector<std::unique_ptr<AutomorphismEvaluator>> evaluators;
    auto N = GetDimension();
    for(uint32_t i = 2; i <= N; i *= 2) {
        auto tmp_context = std::make_shared<AutomorphismContext>(*m_auto_context);
        tmp_context->SetAutomorphismIndex(i + 1);
        evaluators.push_back(tmp_context->ConstructOperator(keys));
    }

    auto op = std::unique_ptr<LWEtoRLWEPacker>(new LWEtoRLWEPacker(shared_from_this(), std::move(evaluators)));

    return op;
}

OperatorID LWEtoRLWEPackingContext::GetOperatorID() const {
    return PACK_LWE_TO_RLWE;
}

KeyDistribution LWEtoRLWEPackingContext::GetSourceKeyDistribution() const {
    return m_auto_context->GetSourceKeyDistribution();
}

uint64_t LWEtoRLWEPackingContext::GetModulus() const {
    return m_auto_context->GetModulus();
}

uint64_t LWEtoRLWEPackingContext::GetDimension() const {
    return m_auto_context->GetDimension();
}

uint64_t LWEtoRLWEPackingContext::GetGadgetBasis() const {
    return m_auto_context->GetGadgetBasis();
}

uint64_t LWEtoRLWEPackingContext::GetGadgetBasisLog2() const {
    return m_auto_context->GetGadgetBasisLog2();
}

uint64_t LWEtoRLWEPackingContext::GetGadgetDigits() const {
    return m_auto_context->GetGadgetDigits();
}

double LWEtoRLWEPackingContext::GetStd() const {
    return m_auto_context->GetStd();
}

std::shared_ptr<intel::hexl::NTT> LWEtoRLWEPackingContext::GetNTT() const {
    return m_auto_context->GetNTT();
}

void LWEtoRLWEPackingContext::SetSourceKeyDistribution(KeyDistribution distribution) {
    m_auto_context->SetSourceKeyDistribution(distribution);
}

void LWEtoRLWEPackingContext::SetModulus(uint64_t modulus) {
    m_auto_context->SetModulus(modulus);
}

void LWEtoRLWEPackingContext::SetDimension(uint64_t input_dimension) {
    m_auto_context->SetDimension(input_dimension);
}

void LWEtoRLWEPackingContext::SetGadgetBasis(uint64_t basis) {
    m_auto_context->SetGadgetBasis(basis);
}

void LWEtoRLWEPackingContext::SetGadgetDigits(uint64_t digits) {
    m_auto_context->SetGadgetDigits(digits);
}

void LWEtoRLWEPackingContext::SetStd(double std) {
    m_auto_context->SetStd(std);
}

void LWEtoRLWEPackingContext::SetNTT(std::shared_ptr<intel::hexl::NTT> ntt) {
    m_auto_context->SetNTT(std::move(ntt));
}

LWEtoRLWEPacker::LWEtoRLWEPacker(std::shared_ptr<const LWEtoRLWEPackingContext> context,
                                 std::vector<std::unique_ptr<AutomorphismEvaluator>> evaluators) : m_context(std::move(context)), m_auto_evaluators(std::move(evaluators)) {

}

void LWEtoRLWEPacker::Pack(uint64_t *output, const uint64_t *input, uint64_t len) {
    // note that the point of doing the packing like this is to save space
    // if we apply the usual, automorphism based packing
    // we'd need to first allocate len * RLWE buffers / or stacks of vector,
    // either allocating too much memory or losing memory continuity
    // In this version, elements that are processed together are always next to each other in memory.

    auto Q = m_context->GetModulus();
    auto N = m_context->GetDimension();
    auto log_N = IntLog2(N);
    auto log_l = IntLog2(len);
    auto ntt = m_context->GetNTT();

    // set up index reversal lambda
    auto brev_index = [len, log_l] (uint64_t idx) {return (len - 1) ^ (reverse_bits<uint64_t>(idx) >> (64 - log_l)); };

    assert(log_l >= 1);
    assert((1 << log_l) == len);
    // in NTT format
    AlignedVector stack_elem((log_N + 1) * N * 2, 0);
    AlignedVector rotated(2 * N, 0); auto rotated_ptr = rotated.data();
    const auto stack_ptr = stack_elem.data();

    std::vector<std::pair<uint64_t, uint64_t>> stack(log_l + 2);
    uint64_t stack_head = 0;
    stack[stack_head] = {0, 0};
    stack[++stack_head] = {0, 1};

    auto current_elem = 0;

    // Copy first 2 LWE samples and transpose

    L0Setup(stack_ptr, input + brev_index(current_elem) * (N + 1));
    current_elem++;

    L0Setup(stack_ptr + 2 * N, input + brev_index(current_elem) * (N + 1));
    current_elem++;

    uint64_t head_level, neck_level;

    while (stack[0].first != log_l) {

        if (stack_head > 0) {
            head_level = stack[stack_head].first;
            neck_level = stack[stack_head - 1].first;
        } else {
            head_level = UINT64_MAX;
            neck_level = 0;
        }

        if (head_level == neck_level) {
            auto output_level = head_level + 1;

            auto shift = N >> output_level;
            // even
            auto head_index = stack_head;
            auto head = stack_ptr + head_index * 2 * N;
            // odd
            auto neck_index = stack_head - 1;
            auto neck = stack_ptr + neck_index * 2 * N;

            // rotate without mult
            ZERO_UINT64_ARR(rotated_ptr, 2 * N);
            intel::hexl::EltwiseSubMod(rotated_ptr, rotated_ptr, neck + N - shift, shift, Q);
            std::copy(neck, neck + N - shift, rotated_ptr + shift);
            intel::hexl::EltwiseSubMod(rotated_ptr + N, rotated_ptr + N, neck + 2 * N - shift, shift, Q);
            std::copy(neck + N, neck + 2 * N - shift, rotated_ptr + N + shift);

            // left operand for packing
            // head + X^N/2^l * neck
            intel::hexl::EltwiseAddMod(neck, rotated_ptr, head, 2 * N, Q);

            // right operand
            // head - X^N/2^l * neck
            intel::hexl::EltwiseSubMod(head, head, rotated_ptr, 2 * N, Q);
            ZERO_UINT64_ARR(rotated_ptr, 2 * N);

            m_auto_evaluators[output_level - 1]->Eval(rotated_ptr, head);
            ntt->ComputeInverse(rotated_ptr, rotated_ptr, 1, 1);
            ntt->ComputeInverse(rotated_ptr + N, rotated_ptr + N, 1, 1);
            intel::hexl::EltwiseAddMod(neck, neck, rotated_ptr, 2 * N, Q);

            stack_head--;
            stack[stack_head] = {output_level, neck_index};

        } else {
            auto c_elem = stack_head + 1;// stack[stack_head].second + 1;
            auto stack_elem_ptr = stack_ptr + c_elem * 2 * N;
            stack[++stack_head] = {0, c_elem};

            ZERO_UINT64_ARR(stack_elem_ptr, 2 * N);
            L0Setup(stack_elem_ptr, input + brev_index(current_elem) * (N + 1));
            current_elem++;
        }
    }

    // now the partial solution is in stack at index 0, if l = N, we are done, otherwise we need to apply more traces

    for(uint32_t i = 0; i < (log_N - log_l); i++) {
        ZERO_UINT64_ARR(rotated_ptr, 2 * N);
        m_auto_evaluators[log_N - i - 1]->Eval(rotated_ptr, stack_ptr);
        ntt->ComputeInverse(rotated_ptr, rotated_ptr, 1,1 );
        ntt->ComputeInverse(rotated_ptr + N, rotated_ptr + N, 1, 1);
        intel::hexl::EltwiseAddMod(stack_ptr, rotated_ptr, stack_ptr, 2 * N, Q);
    }
    std::copy(stack_ptr, stack_ptr + 2 * N, output);
}

void LWEtoRLWEPacker::L0Setup(uint64_t *__restrict output, const uint64_t *__restrict input) {
    auto N = m_context->GetDimension();
    auto Q = m_context->GetModulus();
    auto N_inverse = intel::hexl::InverseMod(N, Q);

    output[0] = input[0];
    for(uint64_t i = 1; i < (N >> 1); i++) {
        output[i] = intel::hexl::SubUIntMod(0, input[N - i], Q);
        output[N - i] = intel::hexl::SubUIntMod(0, input[i], Q);
    }
    output[N >> 1] = intel::hexl::SubUIntMod(0, input[N >> 1], Q);
    output[N] = input[N];
    intel::hexl::EltwiseFMAMod(output, output, N_inverse, nullptr, N + 1, Q, 1);
}

void LWEtoRLWEPacker::Pack(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    auto N = m_context->GetDimension();
    auto n = input.size() / N;

    Pack(output.data(), input.data(), n);
}

void LWEtoRLWEPacker::PackConsecutively(uint64_t *output, const uint64_t *input, uint64_t len) {
    // the previous packing function, packs $len samples into N coefficients
    // such that the stride between elements is N/$len
    // This function takes care of the case where we desire a consecutive packing


    auto Q = m_context->GetModulus();
    auto N = m_context->GetDimension();
    auto log_N = IntLog2(N);
    auto log_l = IntLog2(len);
    auto ntt = m_context->GetNTT();

    // set up index reversal lambda
    auto brev_index = [N, log_N] (uint64_t idx) {return (N - 1) ^ (reverse_bits<uint64_t>(idx) >> (64 - log_N)); };

    assert(log_l >= 1);
    //assert((1 << log_l) == len);
    // in NTT format
    AlignedVector stack_elem((log_N + 1) * N * 2, 0);
    AlignedVector rotated(2 * N, 0); auto rotated_ptr = rotated.data();
    const auto stack_ptr = stack_elem.data();

    std::vector<std::pair<uint64_t, uint64_t>> stack(log_N + 2);
    uint64_t stack_head = 0;

    auto current_elem = 0;

    // Copy first 2 LWE samples and transpose
    auto current_elem_read_idx = brev_index(current_elem);

    if (current_elem_read_idx < len) {
        L0Setup(stack_ptr, input + brev_index(current_elem) * (N + 1));
        stack[stack_head] = {0, 0};
    } else {
        ZERO_UINT64_ARR(stack_ptr, 2 * N);
        stack[stack_head] = {0, 1};
    }

    ++current_elem;
    ++stack_head;

    current_elem_read_idx = brev_index(current_elem);
    if (current_elem_read_idx < len) {
        L0Setup(stack_ptr + 2 * N, input + brev_index(current_elem) * (N + 1));
        stack[stack_head] = {0, 0};
    } else {
        ZERO_UINT64_ARR(stack_ptr + 2 * N, 2 * N);
        stack[stack_head] = {0, 1};
    }
    current_elem++;

    uint64_t auto_counter = 0;
    uint64_t head_level, neck_level;

    while (stack[0].first != log_N) {

        if (stack_head > 0) {
            head_level = stack[stack_head].first;
            neck_level = stack[stack_head - 1].first;
        } else {
            head_level = UINT64_MAX;
            neck_level = 0;
        }

        if (head_level == neck_level) {

            if (stack[stack_head].second == 1 and stack[stack_head - 1].second == 1) {
                stack_head--;
                stack[stack_head] = {head_level + 1, 1};
                continue;
            }

            auto output_level = head_level + 1;

            auto shift = N >> output_level;
            // even
            auto head_index = stack_head;
            auto head = stack_ptr + head_index * 2 * N;
            // odd
            auto neck_index = stack_head - 1;
            auto neck = stack_ptr + neck_index * 2 * N;

            // rotate without mult
            ZERO_UINT64_ARR(rotated_ptr, 2 * N);
            intel::hexl::EltwiseSubMod(rotated_ptr, rotated_ptr, neck + N - shift, shift, Q);
            std::copy(neck, neck + N - shift, rotated_ptr + shift);
            intel::hexl::EltwiseSubMod(rotated_ptr + N, rotated_ptr + N, neck + 2 * N - shift, shift, Q);
            std::copy(neck + N, neck + 2 * N - shift, rotated_ptr + N + shift);

            // left operand for packing
            // head + X^N/2^l * neck
            intel::hexl::EltwiseAddMod(neck, rotated_ptr, head, 2 * N, Q);

            // right operand
            // head - X^N/2^l * neck
            intel::hexl::EltwiseSubMod(head, head, rotated_ptr, 2 * N, Q);
            ZERO_UINT64_ARR(rotated_ptr, 2 * N);

            m_auto_evaluators[output_level - 1]->Eval(rotated_ptr, head);
            auto_counter++;
            ntt->ComputeInverse(rotated_ptr, rotated_ptr, 1, 1);
            ntt->ComputeInverse(rotated_ptr + N, rotated_ptr + N, 1, 1);
            intel::hexl::EltwiseAddMod(neck, neck, rotated_ptr, 2 * N, Q);

            stack_head--;
            stack[stack_head] = {output_level, 0};

        } else {

            current_elem_read_idx = brev_index(current_elem);
            auto stack_elem_ptr = stack_ptr + (stack_head + 1) * 2 * N;
            ZERO_UINT64_ARR(stack_elem_ptr, 2 * N);

            if (current_elem_read_idx < len) {
                L0Setup(stack_elem_ptr, input + current_elem_read_idx * (N + 1));
                stack[++stack_head] = {0, 0};
            } else {
                stack[++stack_head] = {0, 1};
            }
            current_elem++;
        }
    }

    std::cerr << "Called " << auto_counter << std::endl;
    std::copy(stack_ptr, stack_ptr + 2 * N, output);
}