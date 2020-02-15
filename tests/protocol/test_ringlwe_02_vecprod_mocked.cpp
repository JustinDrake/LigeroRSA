#include "../encryption/TestRingLWE.hpp"
#include "easylogging++.h"
#include "gtest/gtest.h"

INITIALIZE_EASYLOGGINGPP

using namespace ligero::lwe;

TEST(RingLWEMultiplication, 2PartiesMockedMultiplication01) {

    using T = uint64_t;

    using message = typename _params<T>::poly_p;
    using cipher = typename _params<T>::poly_q;

    EncryptionEngine<T> e(2);  // init with 2 parties

    auto [pub, sec] = e.gen_keys();

    message a0{1, 2, 3, 4, 5, 6};
    message a1{0};

    a0.invntt_pow_invphi();
    a1.invntt_pow_invphi();

    cipher b0{9};
    cipher b1{1};

    b0.ntt_pow_phi();
    b1.ntt_pow_phi();

    auto ea0 = e.encrypt(a0, pub);
    auto ea1 = e.encrypt(a1, pub);

    auto result = ea0 + ea1;

    // ===== noise =====
    message z0{e.genZ()};
    message z1{e.genZ()};
    message tau{e.getTau()};
    tau.ntt_pow_phi();
    //z0 = z0 * tau;
    //z1 = z1 * tau;

    z0.invntt_pow_invphi();
    z1.invntt_pow_invphi();

    auto ez0 = e.encrypt(z0, pub);
    auto ez1 = e.encrypt(z1, pub);

    auto m0 = result * b0 + ez0;
    auto m1 = result * b1 + ez1;

    auto m = m0 + m1;

    auto d0 = e.partial_decrypt(m, sec[0]);
    auto d1 = e.partial_decrypt({m.first, 0}, sec[1]);  // only one party decrypt the sum

    cipher d = d0 + d1;

    auto plain = e.eval_poly(d);  // get decrypted message

    for (auto i = 0; i < 6; i++) {
        EXPECT_EQ(mpz_get_ui(plain[i]), size_t(i+1) * 10);
    }
}

TEST(RingLWEMultiplication, 100PartiesMockedMultiplication01) {

    using T = uint64_t;

    using message = typename _params<T>::poly_p;
    using cipher = typename _params<T>::poly_q;

	auto numParties = 100;
    EncryptionEngine<T> e(numParties);  // init with numParties parties

    auto [pub, sec] = e.gen_keys();

	// each party samples a_i and encrypts
	message a{1, 2, 3, 4, 5, 6}; // a_0
	a.invntt_pow_invphi();
	auto result = e.encrypt(a, pub); // result = Enc(a) = \sum Enc(a_i)

	for (auto i = 1; i < numParties; i++) { // a_{1-999}
	    a = {0};
		a.invntt_pow_invphi();
		auto ea = e.encrypt(a, pub);
		result = result + ea;
	}

	// noise
	// tau
    message tau{e.getTau()};
    tau.ntt_pow_phi();

	// party 0
    message zi{e.genZ()};
	// zi = zi * tau;
	zi.invntt_pow_invphi();
	auto ezi = e.encrypt(zi, pub); // Enc(z_i * tau)

	cipher bi{9};
	bi.ntt_pow_phi();
	auto m = result * bi + ezi; // sum of all parties' Enc(a)*b_i + Enc(z_i * tau)

	for (auto i = 1; i < numParties; i++) { // party i=[1,999]
	    message zi{e.genZ()};
		// zi = zi * tau;
		zi.invntt_pow_invphi();
		auto ezi = e.encrypt(zi, pub); // Enc(z_i * tau)

		// sample b_i
		if (i == 1) {
			bi = {1};
		} else {
			bi = {0};
		}

		bi.ntt_pow_phi();
		auto mi  = result * bi + ezi; // Enc(a)*b_i + Enc(z_i * tau)
		m = m + mi;
	}

	cipher d;
	for (auto i = 0; i < numParties; i++) {
        if (i == 0) {
			auto di = e.partial_decrypt(m, sec[i]); // only one party decrypts (c,d)
			d = d + di;
		} else {
			auto di = e.partial_decrypt({m.first, 0}, sec[i]); // all other parties decrypt (c,0)
			d = d + di;
		}
	}

    auto plain = e.eval_poly(d);  // get decrypted message

    for (auto i = 0; i < 6; i++) {
        EXPECT_EQ(mpz_get_ui(plain[i]), size_t(i+1) * 10);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
