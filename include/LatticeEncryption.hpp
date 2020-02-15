#pragma once

#include <tuple>
#include <algorithm>
#include <gmpxx.h>

#include "Common.hpp"
#include <nfl.hpp>


namespace ligero::lattice {

template <typename T, size_t Degree, size_t NbPrimes>
using pub = std::pair<nfl::poly_p<T, Degree, NbPrimes>, nfl::poly_p<T, Degree, NbPrimes>>;

template <typename T, size_t Degree, size_t NbPrimes>
using sec = nfl::poly_p<T, Degree, NbPrimes>;

template <typename T, size_t Degree, size_t NbPrimes>
using cipher = pub<T, Degree, NbPrimes>;

template <typename T, size_t Degree, size_t NbPrimes>
using key_pair = std::pair<pub<T, Degree, NbPrimes>, sec<T, Degree, NbPrimes>>;

template <typename T, size_t Degree, size_t NbPrimes>
cipher<T, Degree, NbPrimes> operator+(const cipher<T, Degree, NbPrimes>& a, const cipher<T, Degree, NbPrimes>& b) {
    return { a.first + b.first, a.second + b.second };
}

template <typename Q>
std::pair<Q, Q> pair_sub(const std::pair<Q, Q>& a, const std::pair<Q, Q>& b) {
    return { a.first - b.first, a.second - b.second };
}

template <typename Q>
std::pair<Q, Q> pair_add(const std::pair<Q, Q>& a, const std::pair<Q, Q>& b) {
    return { a.first + b.first, a.second + b.second };
}

template <typename Q>
std::pair<Q, Q> operator-(const std::pair<Q, Q>& a, const std::pair<Q, Q>& b) {
    return { a.first - b.first, a.second - b.second };
}

template <typename Q>
std::pair<Q, Q> operator+(const std::pair<Q, Q>& a, const std::pair<Q, Q>& b) {
    return { a.first + b.first, a.second + b.second };
}

template <typename Q>
std::pair<Q, Q> pair_mul_q(const std::pair<Q, Q>& a, const Q& b) {
    return { a.first * b, a.second * b };
}
template <typename T, size_t Degree, size_t NbPrimes>
cipher<T, Degree, NbPrimes> operator*(const cipher<T, Degree, NbPrimes>& a, const nfl::poly_p<T, Degree, NbPrimes>& b) {
    return { a.first * b, a.second * b };
}

/* Round input number @in to nearest quotient of @div
 * Set the @out number to the quotient (or +1) of @in / @div
 * WARNING: @out number must be initialized.
 */
void roundNearest(mpz_t out, mpz_t in, mpz_t div) {
    mpz_t q, r;
    mpz_init(q);
    mpz_init(r);
    mpz_set_ui(q, 0);
    mpz_set_ui(r, 0);

    mpz_tdiv_qr(q, r, in, div);
    mpz_mul_ui(r, r, 2);
    if (mpz_cmp(r, div) >= 0) {
        mpz_add_ui(q, q, 1);
    }
    mpz_set(out, q);
    mpz_clears(q, r, NULL);  // NULl terminating
}

/* out = pow(op1, op2)
 * WARNING: @out must be initialized
 */
void mpz_pow_d(mpz_t out, mpz_t op1, double op2) {
    mpfr_t x, exp;
    mpfr_init_set_z(x, op1, MPFR_RNDN);
    mpfr_init_set_d(exp, op2, MPFR_RNDN);

    mpfr_pow(x, x, exp, MPFR_RNDN);
    mpfr_get_z(out, x, MPFR_RNDN);

    mpfr_clear(x);
    mpfr_clear(exp);
}

template <typename T, size_t Degree, size_t NbPrimes>
std::array<mpz_class, Degree> generateNoise(const mpz_class& bound) {

    gmp_randstate_t state;
    gmp_randinit_mt(state);

    mpz_class bound2 = bound * 2;

    mpz_class tmp;
    std::array<mpz_class, Degree> arr;

    for (auto it = arr.begin(); it != arr.end(); it++) {
        mpz_urandomm(tmp.get_mpz_t(), state, bound2.get_mpz_t());

        if(tmp > bound) {
            tmp = mpz_class(nfl::poly_p<T, Degree, NbPrimes>::moduli_product()) - (tmp - bound);
        }
        *it = tmp;
    }
    gmp_randclear(state);
    return arr;
}

template <typename T, size_t Degree, size_t NbPrimes>
std::array<mpz_class, Degree> generateNoiseWithMultiple(const std::array<mpz_class, Degree>& bound, const std::array<mpz_class, Degree>& tau) {

    gmp_randstate_t state;
    gmp_randinit_mt(state);

    std::array<mpz_class, Degree> bound2;
    std::transform(bound.begin(), bound.end(), bound2.begin(), [](const auto& x) { return x * 2; });

    mpz_class tmp;
    std::array<mpz_class, Degree> arr;

    for (auto it = arr.begin(), index = 0; it != arr.end(); it++, index++) {
        if(bound2[index] != 0)
          mpz_urandomm(tmp.get_mpz_t(), state, bound2[index].get_mpz_t());
        else tmp = 0;


        if(tmp > bound[index]) {
            *it = mpz_class(nfl::poly_p<T, Degree, NbPrimes>::moduli_product()) - ((tmp - bound[index])*tau[index]);
        }
	else
            *it = tmp*tau[index];
    }
    gmp_randclear(state);
    return arr;
}

template <typename T, size_t Degree, size_t NbPrimes>
std::array<mpz_class, Degree> generateNoise(const std::array<mpz_class, Degree>& bound) {

    gmp_randstate_t state;
    gmp_randinit_mt(state);

    std::array<mpz_class, Degree> bound2;
    std::transform(bound.begin(), bound.end(), bound2.begin(), [](const auto& x) { return x * 2; });

    mpz_class tmp;
    std::array<mpz_class, Degree> arr;

    for (auto it = arr.begin(), index = 0; it != arr.end(); it++, index++) {
        if(bound2[index] != 0)
           mpz_urandomm(tmp.get_mpz_t(), state, bound2[index].get_mpz_t());
        else 
          tmp = 0;

        if(tmp > bound[index]) {
            tmp = mpz_class(nfl::poly_p<T, Degree, NbPrimes>::moduli_product()) - (tmp - bound[index]);
        }
        *it = tmp;
    }
    gmp_randclear(state);
    return arr;
}

/*
 * generate an array of uniformly distributed random noise
 * ranging from [±R], where R = (2 ^ λ) * (2 * σ * P * N^1.5 * n)
 */
template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::array<mpz_class, Degree> generateR(const ProtocolConfig<T>& config) {
    mpz_class pown = config.numParties();
    mpz_pow_d(pown.get_mpz_t(), pown.get_mpz_t(), 1.5);

    mpz_class bound = mpz_class(1) << config.lambda();
    bound *= mpz_class(2) * mpz_class(config.sigma()) * mpz_class(pown) * mpz_class(Degree) * mpz_class(nfl::poly_p<T, Degree, NbPrimesP>::moduli_product());

    return generateNoise<T, Degree, NbPrimesQ>(bound);
}

/*
 * generate an array of uniformly distributed random noise
 * ranging from [±Z], where Z = τ * N * (2 ^ λ)
 */
template <typename T, size_t Degree, size_t NbPrimes>
std::array<mpz_class, Degree> _generateZ(const ProtocolConfig<T>& config, mpz_class tau) {
    mpz_class bound = tau * mpz_class(config.numParties()) * (mpz_class(1) << config.lambda());
    return generateNoise<T, Degree, NbPrimes>(bound);
}

template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
nfl::poly_p<T, Degree, NbPrimesP> _generateZ(const ProtocolConfig<T>& config, std::array<mpz_class, Degree>& tau) {
    using P = nfl::poly_p<T, Degree, NbPrimesP>;

    mpz_class b = mpz_class(config.numParties()) * (mpz_class(1) << config.lambda());
    std::array<mpz_class, Degree> bnd;

    for(int i = 0; i < Degree; i++) bnd[i] = tau[i]*b;
    //std::transform(tau.begin(), tau.end(), tau.begin(), [&](const auto& t) { return t * b; });
    auto zi = generateNoiseWithMultiple<T, Degree, NbPrimesP>(bnd,tau);
    //std::transform(zi.begin(), zi.end(), tau.begin(), zi.begin(), [](const auto& z, const auto& t) { return z * t; });
    P zp;
    zp.set_mpz(zi);
    return zp;
}

template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::tuple<cipher<T, Degree, NbPrimesQ>, nfl::poly_p<T, Degree, NbPrimesQ>, nfl::poly_p<T, Degree, NbPrimesQ>, nfl::poly_p<T, Degree, NbPrimesQ>> _encrypt(
    nfl::poly_p<T, Degree, NbPrimesP>& message,
    const pub<T, Degree, NbPrimesQ>& pubkey,
    nfl::gaussian<uint16_t, T, 2>& chi,
    const mpz_class& q_over_p)
{
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

    Q ui(chi), vi(chi), wi(chi);
    // keep everything in FFT domain
    ui.ntt_pow_phi();
    vi.ntt_pow_phi();
    wi.ntt_pow_phi();

    // transform Q/P to FFT polynomial
    Q qdivp{q_over_p};
    qdivp.ntt_pow_phi();

    // transform message from mod P to mod Q
    message.invntt_pow_invphi();
    Q message_q;
    auto msg_mpz = message.poly2mpz();
    message_q.mpz2poly(msg_mpz);

    // avoid memory leak
    std::for_each(msg_mpz.begin(), msg_mpz.end(), mpz_clear);

    message_q.ntt_pow_phi();
    Q fst = ui * std::get<0>(pubkey) + vi;
    Q sec = ui * std::get<1>(pubkey) + wi + qdivp * message_q;

    return {{fst, sec}, ui, vi, wi};
}

template <typename T, size_t Degree, size_t NbPrimesQ>
std::tuple<cipher<T, Degree, NbPrimesQ>, nfl::poly_p<T, Degree, NbPrimesQ>, nfl::poly_p<T, Degree, NbPrimesQ>, nfl::poly_p<T, Degree, NbPrimesQ>> _encrypt(
    const std::vector<mpz_class>& message,
    const pub<T, Degree, NbPrimesQ>& pubkey,
    nfl::gaussian<uint16_t, T, 2>& chi,
    const mpz_class& q_over_p)
{
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

    assert(message.size() <= Degree && "encrypt: message size must be LE than poly degree");

    Q ui(chi), vi(chi), wi(chi);
    // keep everything in FFT domain
    ui.ntt_pow_phi();
    vi.ntt_pow_phi();
    wi.ntt_pow_phi();

    // transform Q/P to FFT polynomial
    Q qdivp{q_over_p};
    qdivp.ntt_pow_phi();

    // preparing the message polynomial
    auto *ptr = new std::array<mpz_class, Degree>{};
    std::copy(message.begin(), message.end(), ptr->begin());
    Q message_q;
    message_q.set_mpz(*ptr);
    delete(ptr);

    Q fst = ui * std::get<0>(pubkey) + vi;
    Q sec = ui * std::get<1>(pubkey) + wi + qdivp * message_q;

    return {{fst, sec}, ui, vi, wi};
}

template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::pair<nfl::poly_p<T, Degree, NbPrimesQ>,nfl::poly_p<T, Degree, NbPrimesQ>> _partial_decrypt(
    const cipher<T, Degree, NbPrimesQ>& message,
    const sec<T, Degree, NbPrimesQ>& seckey,
    const ProtocolConfig<T>& config)
{
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;
    Q r{generateR<T, Degree, NbPrimesP, NbPrimesQ>(config)};
    r.ntt_pow_phi();
    return {std::get<1>(message) - std::get<0>(message) * seckey + r, r};
}

/*
 * Evaluating the polynomial once get the sum of partial decrypt.
 * Input message must in FFT domain, return an array of original message.
 */
template <typename TP, typename TQ, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::array<mpz_t, Degree> _eval_poly(nfl::poly_p<TQ, Degree, NbPrimesQ>& partial, mpz_class& q_over_p, const std::vector<mpz_class>& tauVector) {
    partial.invntt_pow_invphi();
    auto coeff = partial.poly2mpz();

    for (auto it = coeff.begin(); it != coeff.end(); it++) {
        roundNearest(*it, *it, q_over_p.get_mpz_t());
    };

    nfl::poly_p<TP, Degree, NbPrimesP> msg;
    msg.mpz2poly(coeff);
    std::for_each(coeff.begin(), coeff.end(), mpz_clear);
    msg.ntt_pow_phi();
    coeff = msg.poly2mpz();

    size_t i = 0;
    for (auto it = coeff.begin(); it != coeff.end(); it++) {
        mpz_class Pby2 = mpz_class(nfl::poly_p<TP, Degree, NbPrimesP>::moduli_product()) / 2;
        mpz_class tauminusoneP = (tauVector[i] - 1) * mpz_class(nfl::poly_p<TP, Degree, NbPrimesP>::moduli_product());
        if(mpz_cmp(*it, Pby2.get_mpz_t()) > 0) mpz_add(*it, *it, tauminusoneP.get_mpz_t());
        mpz_tdiv_r(*it, *it, tauVector[i].get_mpz_t());
        ++i;
    };

    return coeff;
}
/*
 * Evaluating the polynomial once get the sum of partial decrypt.
 * Input message must in FFT domain, return an array of original message.
 */
template <typename TP, typename TQ, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::array<mpz_t, Degree> _eval_poly(nfl::poly_p<TQ, Degree, NbPrimesQ>& partial, mpz_class& q_over_p, const mpz_class& tau) {
    partial.invntt_pow_invphi();
    auto coeff = partial.poly2mpz();

    for (auto it = coeff.begin(); it != coeff.end(); it++) {
        roundNearest(*it, *it, q_over_p.get_mpz_t());
    };

    nfl::poly_p<TP, Degree, NbPrimesP> msg;
    msg.mpz2poly(coeff);
    std::for_each(coeff.begin(), coeff.end(), mpz_clear);
    msg.ntt_pow_phi();
    coeff = msg.poly2mpz();

    for (auto it = coeff.begin(); it != coeff.end(); it++) {
        mpz_class Pby2 = mpz_class(nfl::poly_p<TP, Degree, NbPrimesP>::moduli_product()) / 2;
        mpz_class tauminusoneP = (tau - 1) * mpz_class(nfl::poly_p<TP, Degree, NbPrimesP>::moduli_product());
        if(mpz_cmp(*it, Pby2.get_mpz_t()) > 0) mpz_add(*it, *it, tauminusoneP.get_mpz_t());
        mpz_tdiv_r(*it, *it, tau.get_mpz_t());
    };

    return coeff;
}

template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
std::tuple<cipher<T, Degree, NbPrimesQ>,
    nfl::poly_p<T, Degree, NbPrimesQ>> // obsolete: zi NOTE: zi is in P not Q
    // nfl::poly_p<T, Degree, NbPrimesQ>, // u
    // nfl::poly_p<T, Degree, NbPrimesQ>, // v
    // nfl::poly_p<T, Degree, NbPrimesQ>> // w
cipher_mul_scalar(
    const cipher<T, Degree, NbPrimesQ>& ciph,
    nfl::poly_p<T, Degree, NbPrimesP>& scalar,
    const ProtocolConfig<T>& config,
    std::array<mpz_class, Degree>& tau,
    const pub<T, Degree, NbPrimesQ>& key,
    nfl::gaussian<uint16_t, T, 2>& chi,
    const mpz_class& q_over_p
) {
    using P = nfl::poly_p<T, Degree, NbPrimesP>;
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;

    scalar.invntt_pow_invphi();
    auto s_coeff = scalar.poly2mpz();
    Q sq;
    sq.set_mpz(s_coeff);
    std::for_each(s_coeff.begin(), s_coeff.end(), mpz_clear);
    sq.ntt_pow_phi();

    P zi = _generateZ<T, Degree, NbPrimesP, NbPrimesQ>(config, tau);
    Q promotedZi;

    zi.invntt_pow_invphi();
    auto zi_coeff = zi.poly2mpz();
    promotedZi.set_mpz(zi_coeff);
    std::for_each(zi_coeff.begin(), zi_coeff.end(), mpz_clear);
    promotedZi.ntt_pow_phi();

    // auto [ezi, uz, vz, wz] = _encrypt<T, Degree, NbPrimesP, NbPrimesQ>(zi, key, chi, q_over_p);
    cipher<T, Degree, NbPrimesQ> temp;
    temp.second = promotedZi;

    return {ciph*sq + temp, promotedZi};
}

template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
class LatticeEncryption {

public:
    using P = nfl::poly_p<T, Degree, NbPrimesP>;
    using Q = nfl::poly_p<T, Degree, NbPrimesQ>;
    using pubKey = pub<T, Degree, NbPrimesQ>;
    using secKey = sec<T, Degree, NbPrimesQ>;
    using cipherText = cipher<T, Degree, NbPrimesQ>;

    using num_type = mpz_class;

    LatticeEncryption(ProtocolConfig<T>& config)
        : _config(config), _chi(_init_gaussian()), _q_over_p(_init_q_over_p())
        { }

    ~LatticeEncryption() { delete(_fg); }

public:

    std::tuple<cipher<T, Degree, NbPrimesQ>, Q, Q, Q> encrypt(P& message, const pubKey& key) {
        return _encrypt(message, key, _chi, _q_over_p);
    }

    std::tuple<cipher<T, Degree, NbPrimesQ>, Q, Q, Q> encrypt(const std::vector<num_type>& message, const pubKey& key) {
        return _encrypt(message, key, _chi, _q_over_p);
    }

    auto partial_decrypt(const cipherText& message, const secKey& key) {
        return _partial_decrypt<T, Degree, NbPrimesP, NbPrimesQ>(message, key, _config);
    }

    std::array<mpz_t, Degree> eval_poly(Q& poly, const num_type& tau) {
        return _eval_poly<T, T, Degree, NbPrimesP, NbPrimesQ>(poly, _q_over_p, tau);
    }

    std::array<mpz_t, Degree> eval_poly(Q& poly, const std::vector<num_type>& tauVector) {
        return _eval_poly<T, T, Degree, NbPrimesP, NbPrimesQ>(poly, _q_over_p, tauVector);
    }


public:
    P generateZ(const num_type& tau) {
        P noise{_generateZ<T, Degree, NbPrimesP>(_config, tau)};
        return noise;
    }

    nfl::gaussian<uint16_t, T, 2>& chi() {
        return _chi;
    }

    std::tuple<cipher<T, Degree, NbPrimesQ>, Q> product(const cipherText& cipher, P& scalar, std::array<num_type, Degree>& tau, const pubKey& key) {
        return cipher_mul_scalar<T, Degree, NbPrimesP, NbPrimesQ>(cipher, scalar, _config, tau, key, _chi, _q_over_p);
    }

private:
    mpz_class _init_q_over_p() {
        mpz_class acc = 1;
        for (auto i = NbPrimesP; i < NbPrimesQ; i++) {
            acc *= nfl::params<T>::P[i];
        }
        return acc;
    }

    nfl::gaussian<uint16_t, T, 2> _init_gaussian() {
        _fg = new nfl::FastGaussianNoise<uint16_t, T, 2>(_config.sigma(), _config.lambda(), Degree);
        return nfl::gaussian(_fg);
    }

private:
    ProtocolConfig<T>& _config;
    nfl::FastGaussianNoise<uint16_t, T, 2> *_fg;
    nfl::gaussian<uint16_t, T, 2> _chi;
    mpz_class _q_over_p;
};

} // namespace ligero::lattice
