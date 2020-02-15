/*
 *  This file contains helper functions for testing encryption/decryption
 *  methods locally. No transport layer needed.
 */

#pragma once

#include <nfl.hpp>
#include <mpfr.h>
#include <gmpxx.h>

#include <array>
#include <vector>
#include <tuple>
#include <algorithm>

namespace ligero::lwe {

    template <typename _T>
    struct _params {};

    template <>
    struct _params<uint16_t> {
        constexpr static auto p = 1;
        constexpr static auto q = 2;
        constexpr static auto degree = 16;
        using Chi = nfl::gaussian<uint16_t, uint16_t, 2>;
        using poly_p = nfl::poly_p<uint16_t, degree, p>;
        using poly_q = nfl::poly_p<uint16_t, degree, q>;
        using poly_p_ntt = poly_p;
        using poly_q_ntt = poly_q;
    };

    template <>
    struct _params<uint32_t> {
        constexpr static auto p = 6;
        constexpr static auto q = 15;
        constexpr static auto degree = 128;
        using Chi = nfl::gaussian<uint16_t, uint32_t, 2>;
        using poly_p = nfl::poly_p<uint32_t, degree, p>;
        using poly_q = nfl::poly_p<uint32_t, degree, q>;
        using poly_p_ntt = poly_p;
        using poly_q_ntt = poly_q;
    };

    template <>
    struct _params<uint64_t> {
        constexpr static auto p = 6;
        constexpr static auto q = 15;
        constexpr static auto degree = 1<<15;
        using Chi = nfl::gaussian<uint16_t, uint64_t, 2>;
        using poly_p = nfl::poly_p<uint64_t, degree, p>;
        using poly_q = nfl::poly_p<uint64_t, degree, q>;
        using poly_p_ntt = poly_p;
        using poly_q_ntt = poly_q;
    };

    template <typename _T>
    using pub_key = std::pair<typename _params<_T>::poly_q, typename _params<_T>::poly_q>;

    template <typename _T>
    using sec_key = typename _params<_T>::poly_q;

    template <typename _T>
    using keys = std::pair<pub_key<_T>, sec_key<_T>>;

    template <typename _T>
    using cipher = pub_key<_T>;

    template <typename _T>
    using many_keys = std::pair<pub_key<_T>, std::vector<sec_key<_T>>>;

    /* provide operator+ overloading for pairs. We require the first and the second
     * element has same type (i.e. poly_q)
     */
    template <typename _T>
    std::pair<_T, _T> operator+(const std::pair<_T, _T>& _p1, const std::pair<_T, _T>& _p2) {
        return { _p1.first + _p2.first, _p1.second + _p2.second };
    }

    /* provide operator* overloading for pairs.
     * (a, b) * k == (a * k, b * k)
     */
    template <typename _T>
    std::pair<_T, _T> operator*(const std::pair<_T, _T>& _p1, const _T& konst) {
        return { _p1.first * konst, _p1.second * konst };
    }

    /* provide operator* overloading for pairs.
     * k * (a, b) == (a * k, b * k)
     */
    template <typename _T>
    std::pair<_T, _T> operator*(const _T& konst, const std::pair<_T, _T>& _p1) {
        return { _p1.first * konst, _p1.second * konst };
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


    template <typename _T>
    class EncryptionEngine {
    public:
        using value_type = _T;
        using P = typename _params<_T>::poly_p;
        using P_fft = typename _params<_T>::poly_p_ntt;
        using Q = typename _params<_T>::poly_q;
        using Q_fft = typename _params<_T>::poly_q_ntt;

        EncryptionEngine(size_t num_parties): EncryptionEngine(num_parties, 8, 80, 1000) { }
        EncryptionEngine(size_t num_parties, size_t sigma, size_t lambda, size_t tau)
            : _num_parties(num_parties)
            , _sigma(sigma), _lambda(lambda), _tau(tau)
            , _q_over_p(_calc_q_over_p()), _chi(_init_gaussian())
        { }
        ~EncryptionEngine() { delete(_fg); }


    // random noise generators
    public:
    /*
     * generate an array of uniformly distributed random noise
     * ranging from [±Z], where Z = τ * N * (2 ^ λ)
		 * 0 to 2R
		 * if x in [R, 2R], subtract R, do Q - (x - R)
     */

    std::array<mpz_class, _params<_T>::degree> genZ() {
        mpz_class bound = mpz_class(_tau) * mpz_class(_num_parties) * (mpz_class(1) << _lambda);
        mpz_class bound2 = bound * 2;

        gmp_randstate_t state;
        gmp_randinit_mt(state);

        mpz_class tmp{1};
        std::array<mpz_class, _params<_T>::degree> arr;

        for (auto it = arr.begin(); it != arr.end(); it++) {
            mpz_urandomm(tmp.get_mpz_t(), state, bound2.get_mpz_t());
            // tmp = tmp - bound;
        	if (mpz_cmp(tmp.get_mpz_t(),bound.get_mpz_t())> 0 ){
        	//if (tmp > bound){
        	    tmp = mpz_class(P::moduli_product()) - ((tmp-bound) * _tau);
            } else {
                tmp = tmp * _tau;
            }
            //tmp = 0; //TODO: remove and debug
            *it = tmp;
            //*it = _tau;
        }
        gmp_randclear(state);
        return arr;
    }

    /*
     * generate an array of uniformly distributed random noise
     * ranging from [±R], where R = (2 ^ λ) * (2 * σ * P * N^1.5 * n)
     */
    std::array<mpz_class, _params<_T>::degree> genR() {
        mpz_class pown = _num_parties;
        mpz_pow_d(pown.get_mpz_t(), pown.get_mpz_t(), 1.5);

		mpz_class bound = mpz_class(1) << _lambda;
        bound *= 2 * mpz_class(_sigma) * mpz_class(P::moduli_product()) * pown * _params<_T>::degree;
		mpz_class bound2 = bound * 2;

        gmp_randstate_t state;
        gmp_randinit_mt(state);

        mpz_class tmp{1};
        std::array<mpz_class, _params<_T>::degree> arr;

        for (auto it = arr.begin(); it != arr.end(); it++) {
            mpz_urandomm(tmp.get_mpz_t(), state, bound2.get_mpz_t());
            if (mpz_cmp(tmp.get_mpz_t(),bound.get_mpz_t()) > 0){
            //if (tmp > bound){
                tmp = mpz_class(Q::moduli_product()) - tmp + bound;
            }
    		//tmp = 0; //TODO: remove and debug
            *it = tmp;
        }
        gmp_randclear(state);
        return arr;
    }

    // encryption and decryption
    public:

        /*
         * generate a public key (a, b) and a secret key s
         * return ((a, b), s)
         */
        keys<_T> gen_keypair() {
            // first uniformly choose a
            Q a(nfl::uniform{});

            // gaussian dist
            Q s(_chi);
            Q e(_chi);

            // keep everything in FFT domain
            a.ntt_pow_phi();
            s.ntt_pow_phi();
            e.ntt_pow_phi();

            Q_fft b = s * a + e;
            // return (public Cipher, private Cipher)
            return {{a, b}, s};
        }


        /*
         * generate one public key and n secret keys indicated by num_parties
         * return ((a, b), {s0, s1 ..})
         */
        many_keys<_T> gen_keys() {
            std::vector<Q_fft> vec(_num_parties);
            Q_fft a;
            a.ntt_pow_phi();
            for (size_t i = 0; i < _num_parties; i++) {
                Q u(nfl::uniform{});
                u.ntt_pow_phi();
                a = a + u;
            }

            Q_fft b;
            b.ntt_pow_phi();
            for (size_t i = 0; i < _num_parties; i++) {
                Q si(_chi);
                Q ei(_chi);
                si.ntt_pow_phi();
                ei.ntt_pow_phi();
                vec[i] = si;

                Q_fft bi = si * a + ei;
                b = b + bi;
            }
            return {{a, b}, vec};
        }

        /*
         * Encrypt a coefficient domain polynomial using pub key k
         * return FFT domain cipher text
         */
        cipher<_T> encrypt(P& message, const cipher<_T>& k) {
            Q u(_chi), v(_chi), w(_chi);

            // keep everything in FFT domain
            u.ntt_pow_phi();
            v.ntt_pow_phi();
            w.ntt_pow_phi();

            // transform Q/P to FFT polynomial
            Q qdivp{_q_over_p};
            qdivp.ntt_pow_phi();

            // transform message from mod P to mod Q
            // TODO: find a more efficient way to do this
            Q message_q;
            auto msg_mpz = message.poly2mpz();
            message_q.mpz2poly(msg_mpz);

            // avoid memory leak
            std::for_each(msg_mpz.begin(), msg_mpz.end(), mpz_clear);

            message_q.ntt_pow_phi();
            Q_fft fst = u * std::get<0>(k) + v;
            Q_fft sec = u * std::get<1>(k) + w + qdivp * message_q;

            return {fst, sec};
        }

        /*
         * Decrypt a ciphertext.
         * When decrypting, only one party needs to decrypt the sum of ciphertext
         * (c, d), and all other parties should decrypt (c, 0) in order to get
         * correct result
         */
        Q_fft partial_decrypt(const cipher<_T>& message, const sec_key<_T>& s) {
            Q r{genR()};
            r.ntt_pow_phi();

            Q_fft dec = message.second - message.first * s + r;
			return dec;
        }

        std::array<mpz_t, _params<_T>::degree> eval_poly(Q_fft& partial) {
            partial.invntt_pow_invphi();
            auto coeff = partial.poly2mpz();

            for (auto it = coeff.begin(); it != coeff.end(); it++) {
                roundNearest(*it, *it, _q_over_p.get_mpz_t());
            };

            P msg;
            msg.mpz2poly(coeff);
            std::for_each(coeff.begin(), coeff.end(), mpz_clear);
            msg.ntt_pow_phi();
            coeff = msg.poly2mpz();

            mpz_class mpz_tau(_tau);
            for (auto it = coeff.begin(); it != coeff.end(); it++) {
                mpz_class Pby2 = mpz_class(P::moduli_product()) / 2;
                mpz_class tauminusoneP = (mpz_tau - 1) * mpz_class(P::moduli_product());
                if(mpz_cmp(*it,Pby2.get_mpz_t()) > 0) mpz_add(*it,*it,tauminusoneP.get_mpz_t()); 
                mpz_tdiv_r(*it, *it, mpz_tau.get_mpz_t());
            };

            return coeff;
        }

    public:
        size_t getTau() {
            return _tau;
        }

    private:
        mpz_class _calc_q_over_p() {
            mpz_class acc = 1;
            for (auto i = _params<_T>::p; i < _params<_T>::q; i++) {
                acc *= nfl::params<_T>::P[i];
            }
            return acc;
        }

        typename _params<_T>::Chi _init_gaussian() {
            _fg = new nfl::FastGaussianNoise<uint16_t, _T, 2>(_sigma, _lambda, _params<_T>::degree);
            return nfl::gaussian(_fg);
        }

    private:
        size_t _num_parties;
        size_t _sigma;
        size_t _lambda;
        // TODO: we might want tau > 64 bits
        size_t _tau;
        mpz_class _q_over_p;
        nfl::FastGaussianNoise<uint16_t, _T, 2> *_fg;
        typename _params<_T>::Chi _chi;
    };
} // namespace ligero::lwe
