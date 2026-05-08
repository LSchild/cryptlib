//
// Created by leonard on 4/15/26.
//
#include <cassert>
#include <iostream>
#include <utility>
#include "bmmp_blind_rotator.h"
#include "speed_utils.h"

BMMPBlindRotatorParams::BMMPBlindRotatorParams(KeyDistribution distr, uint64_t modulus, uint64_t ring_dim, uint64_t lwe_dim, uint64_t m_basis, uint64_t m_digits, double std, uint64_t step_size)
: CGGIBlindRotatorParams(distr, modulus, ring_dim, lwe_dim, m_basis, m_digits, std), m_step_size(step_size) {
    if (GetLWEDimension() % m_step_size != 0) {
        std::cerr << "Unsupported parameter: LWE dimension not divisible by 2" << std::endl;
        assert(false);
    }

    if (m_step_size != 2) {
        std::cerr << "Only step size of 2 is allowed (for now)" << std::endl;
        assert(false);
    }
}

BMMPBlindRotatorParams::BMMPBlindRotatorParams(KeyDistribution distr, std::shared_ptr<intel::hexl::NTT> ntt,
                                               uint64_t lwe_dim, uint64_t basis, uint64_t digits, double std,
                                               uint64_t step_size) : CGGIBlindRotatorParams(distr, std::move(ntt), lwe_dim,basis,digits,std), m_step_size(step_size) {
    if (GetLWEDimension() % m_step_size != 0) {
        std::cerr << "Unsupported parameter: LWE dimension not divisible by 2" << std::endl;
        assert(false);
    }

    if (m_step_size != 2) {
        std::cerr << "Only step size of 2 is allowed (for now)" << std::endl;
        assert(false);
    }
}

BMMPBlindRotatorParams::BMMPBlindRotatorParams(BMMPBlindRotatorParams& other) : CGGIBlindRotatorParams(other), m_step_size(other.m_step_size) {
    if (GetLWEDimension() % m_step_size != 0) {
        std::cerr << "Unsupported parameter: LWE dimension not divisible by 2" << std::endl;
        assert(false);
    }

    if (m_step_size != 2) {
        std::cerr << "Only step size of 2 is allowed (for now)" << std::endl;
        assert(false);
    }
}


uint64_t BMMPBlindRotatorParams::GetStepSize() const {
    return m_step_size;
}

long double BMMPBlindRotatorParams::ComputeOutputVariance(long double input_variance) const {
    // TODO
    return 0.0;
}

BMMPBlindRotator::BMMPBlindRotator(BMMPBlindRotatorParams &params) : m_params(params), m_params_set(true) {
    m_engine = params.GetNTT();
    m_mux = std::make_unique<MuxOperator>(m_engine, m_params.GetBlindRotationBasisLog2(), m_params.GetBlindRotationRGSWDigits());
    m_encryptor = std::make_shared<RLWEEncryptor>(m_engine, m_params.GetStd());
    m_accumulator.resize(3 * m_params.GetRingDimension());
    SetupMonomials();
    m_params_set = true;
}

void BMMPBlindRotator::SetParams(BMMPBlindRotatorParams &params) {
    m_params = params;
    m_engine = m_params.GetNTT();
    m_mux = std::make_unique<MuxOperator>(m_engine, m_params.GetBlindRotationBasisLog2(), m_params.GetBlindRotationRGSWDigits());
    m_encryptor = std::make_shared<RLWEEncryptor>(m_engine, m_params.GetStd());
    m_accumulator.resize(2 * m_params.GetRingDimension());

    m_params_set = true;
    SetupMonomials();
}

BlindRotationMethod BMMPBlindRotator::GetMethod() {
    return BMMP;
}

const BMMPBlindRotatorParams& BMMPBlindRotator::GetParams() const {
    return m_params;
}

void BMMPBlindRotator::KeyGen(const std::vector<uint64_t> &lwe_key, const std::vector<uint64_t> &rlwe_key) {
    assert(m_params_set);
    KeyGen(lwe_key.data(), rlwe_key.data());
}

void BMMPBlindRotator::KeyGen(const uint64_t *__restrict lwe_key, const uint64_t *__restrict rlwe_key) {
    assert(m_params_set);
    if (m_params.GetKeyDistribution() == BINARY) {
        KeyGenBinary(lwe_key, rlwe_key);
    } else {
        std::cerr << "Unsupported key distribution (for now)" << std::endl;
        assert(false);
    }
    m_keys_generated = true;
}

void BMMPBlindRotator::KeyGenBinary(const uint64_t *__restrict lwe_key, const uint64_t *__restrict rlwe_key) {
    auto n = m_params.GetLWEDimension();
    auto N = m_params.GetRingDimension();
    auto L = m_params.GetBlindRotationBasis();
    auto l = m_params.GetBlindRotationRGSWDigits();

    auto sk_ntt = AlignedVector(N);
    m_engine->ComputeForward(sk_ntt.data(), rlwe_key, 1, 1);
    auto data = AlignedVector(N, 0);

    m_brk.resize(3 * (n >> 1) * (4 * N * l));
    for(uint64_t i = 0; i < n; i+=2) {
        auto key_chunk = m_brk.data() + 3 * (i >> 1) * (4 * N * l);
        auto s_i = lwe_key[i];
        auto s_ip1 = lwe_key[i + 1];

        // used for poly = X^{a_0} - 1
        data[0] = s_i;
        m_encryptor->MakeRGSW(key_chunk, data.data(), sk_ntt.data(), L, l);
        // used for poly = X^{a_1} - 1
        data[0] = s_ip1;
        key_chunk += 4 * N * l;
        m_encryptor->MakeRGSW(key_chunk, data.data(), sk_ntt.data(), L, l);
        // used for poly = (X^{a_0 + a_1} - X^a_0 - X^{a_1} + 1
        data[0] = s_ip1 * s_i;
        key_chunk += 4 * N * l;
        m_encryptor->MakeRGSW(key_chunk, data.data(), sk_ntt.data(), L, l);
        /**
         * s_0 s_1 * (X^{a_0 + a_1} - 1) * acc + s_0 (1 - s_1) * X^{a_0} * acc + (1 - s_0) s_1 * X^{a_1} * acc + (1 - s_0)(1 - s_1) * acc
         * "" + acc + s_0 s_1 * acc - s_0 acc -s_1 * acc
         */
    }
    m_keys_generated = true;
}

void BMMPBlindRotator::SetupMonomials() {
    auto N = m_params.GetRingDimension();

    m_monomials.resize(N * (N + 2));
    std::fill(m_monomials.begin(), m_monomials.end(), 0);
    for(uint64_t i = 0; i < N; i++) {
        auto idx = i * N;
        m_monomials[idx + i] = 1;
        m_engine->ComputeForward(m_monomials.data() + idx, m_monomials.data() + idx, 1, 1);
    }
}

std::shared_ptr<RLWEEncryptor> BMMPBlindRotator::GetEncryptor() const {
    return m_encryptor;
}

void BMMPBlindRotator::BuildStepPolynomials(uint64_t *__restrict poly_buffer, uint64_t a_0, uint64_t a_1) {
    const auto Q = m_params.GetModulus();
    const auto N = m_params.GetRingDimension();
    const auto N2 = (N << 1);
    const auto mask = N2 - 1;
    // We assume poly_buffer has size 3 * N
    uint64_t sum = (a_0 + a_1) & mask;

    const auto* const monomial_list = m_monomials.data();
    const auto* const one_ptr = monomial_list;

    ZERO_UINT64_ARR(poly_buffer, 3 * N);

    const uint64_t coefs [] = {a_0, a_1, sum};

    for(uint32_t i = 0; i < 3; i++) {
        if (coefs[i] >= N) {
            intel::hexl::EltwiseSubMod(poly_buffer + i * N, poly_buffer + i * N, monomial_list + (coefs[i] - N) * N, N, Q);
        } else {
            std::copy(monomial_list + coefs[i] * N, monomial_list + coefs[i] * N + N, poly_buffer + i * N);
            //intel::hexl::EltwiseAddMod(poly_buffer + i * N, poly_buffer + i * N, monomial_list + coefs[i] * N, N, Q);
        }
    }

    intel::hexl::EltwiseSubMod(poly_buffer + 2 * N, poly_buffer + 2 * N, poly_buffer, N, Q);
    intel::hexl::EltwiseSubMod(poly_buffer + 2 * N, poly_buffer + 2 * N, poly_buffer + N, N, Q);

    intel::hexl::EltwiseSubMod(poly_buffer, poly_buffer, one_ptr, N, Q);
    intel::hexl::EltwiseSubMod(poly_buffer + N, poly_buffer + N, one_ptr, N, Q);
    intel::hexl::EltwiseAddMod(poly_buffer + 2 * N, poly_buffer + 2 * N, one_ptr, N, Q);

    /*
    auto ntt = m_engine;
    std::cerr << "a0 = " << a_0 << " a1 = " << a_1 << " a0 + a1 = " << sum << std::endl;
    for(uint32_t j = 0; j < 3; j++) {
        ntt->ComputeInverse(poly_buffer + j * N, poly_buffer + j * N, 1, 1);

        for (uint32_t i = 0; i < N; i++) {
            std::cerr << poly_buffer[i + j * N] << ", ";
        }
        std::cerr << std::endl;
    } */

}

void BMMPBlindRotator::BlindRotate(const uint64_t *__restrict lwe_vec, uint64_t *__restrict rlwe_acc_vec) {
    BlindRotateBinary(lwe_vec, rlwe_acc_vec);
}

void BMMPBlindRotator::BlindRotate(const std::vector<uint64_t> &lwe_vec, std::vector<uint64_t> &rlwe_acc_vec) {
    BlindRotateBinary(lwe_vec.data(), rlwe_acc_vec.data());
}

void BMMPBlindRotator::BlindRotateBinary(const uint64_t *const __restrict lwe_vec, uint64_t *__restrict rlwe_vec) {
    // LWE dimension
    const auto n = m_params.GetLWEDimension();
    // RLWE dimension
    const auto N = m_params.GetRingDimension();
    const auto N2 = (N << 1);
    const auto N2mask = N2 - 1;
    // RLWE modulus
    const auto Q = m_params.GetModulus();
    // RGSW digits
    const auto l = m_params.GetBlindRotationRGSWDigits();
    // RGSW size in bytes
    const auto rgsw_size = 4 * N * l;

    // pointer to scratch space
    const auto m_acc_p = m_accumulator.data();
    const auto* const m_mon_p = m_monomials.data();
    const auto* const m_zero_p = m_mon_p + N * N;


    // Multiply by X^b
    uint64_t b = lwe_vec[n];
    if (b >= N) {
        // b >= N => X^b = -X^{b - N}
        b -= N;
        // m_acc = [-rlwe_A, -rlwe_B]
        intel::hexl::EltwiseSubMod(m_acc_p, m_zero_p, rlwe_vec, 2 * N, Q);
        // monomial ptr to X^{b - N}
        auto monomial = m_mon_p + b * N;
        // m_acc = [-rlwe_A * X^{b - N}, -rlwe_B * X^{b - N}]
        //       = [rlwe_A * X^b, rlwe_B * X^b]
        intel::hexl::EltwiseMultMod(m_acc_p, m_acc_p, monomial, N, Q, 1);
        intel::hexl::EltwiseMultMod(m_acc_p + N, m_acc_p + N, monomial, N, Q, 1);
    } else {
        // monomial ptr to X^b
        auto monomial = m_mon_p + b * N;
        // m_acc = [rlwe_A * X^b, rlwe_B * X^b]
        intel::hexl::EltwiseMultMod(m_acc_p, rlwe_vec, monomial, N, Q, 1);
        intel::hexl::EltwiseMultMod(m_acc_p + N, rlwe_vec + N, monomial, N, Q, 1);
    }

    // copy back to recycle scratch space
    std::copy(m_acc_p, m_acc_p + 2 * N, rlwe_vec);

    ZERO_UINT64_ARR(m_accumulator.data(), 3 * N);

    // pointer to blind-rotation key
    const auto* const m_brk_p = m_brk.data();
    for(uint32_t i = 0; i < n; i+=2) {

        const auto* rgsw_samples = m_brk_p + rgsw_size * 3 * (i >> 1);
        const auto a_i = (N2 - lwe_vec[i]) & N2mask;
        const auto a_i1 = (N2 - lwe_vec[i + 1]) & N2mask;

        BuildStepPolynomials(m_acc_p, a_i, a_i1);

        // TODO: do we gain anything by skipping for a_i == 0 or a_i1 == 0
        m_mux->MultiMux(rlwe_vec, 3, rgsw_samples, m_acc_p);
    }


}

Container BMMPBlindRotator::GetInputContainer() const {
    auto n = m_params.GetLWEDimension();
    auto N = m_params.GetRingDimension();
    auto q = 2 * N;
    auto Q = m_params.GetModulus();

    Container lwecont = std::make_shared<LWEContainerImpl>(q, n, 0.0);
    Container rlwecont = std::make_shared<RLWEContainerImpl>(Q, N, 0.0);
    std::vector<Container> args = {lwecont, rlwecont};
    return std::make_shared<TupleContainerImpl>(args);
}

Container BMMPBlindRotator::GetOutputContainer() const {
    auto N = m_params.GetRingDimension();
    auto Q = m_params.GetModulus();
    auto var = m_params.ComputeOutputVariance();

    return std::make_shared<RLWEContainerImpl>(Q, N, var);
}