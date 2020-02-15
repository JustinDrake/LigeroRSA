#include "TestRingLWE.hpp"
#include "easylogging++.h"
#include "gtest/gtest.h"

INITIALIZE_EASYLOGGINGPP

using namespace ligero::lwe;

TEST(RingLWEAddition, 2PartiesMockedAddition01) {

    using T = uint64_t;

    using message = typename _params<T>::poly_p;
    using cipher = typename _params<T>::poly_q;

    EncryptionEngine<T> e(2);  // init with 2 parties

    auto [pub, sec] = e.gen_keys();

    message number{1, 2, 3, 4, 5, 6};
    message zero{0};

    number.invntt_pow_invphi();
    zero.invntt_pow_invphi();

    auto ezeros = e.encrypt(number, pub);
    auto eone = e.encrypt(zero, pub);

    auto result = ezeros + eone;

    auto d0 = e.partial_decrypt(result, sec[0]);
    auto d1 = e.partial_decrypt({result.first, 0}, sec[1]);  // only one party decrypt the sum

    cipher d = d0 + d1;

    auto plain = e.eval_poly(d);  // get decrypted message

    for (auto i = 0; i < 6; i++) {
        EXPECT_EQ(mpz_get_ui(plain[i]), size_t(i+1));
    }
}

TEST(RingLWEAddition, 2PartiesMockedAddition02) {

    using T = uint64_t;

    using message = typename _params<T>::poly_p;
    using cipher = typename _params<T>::poly_q;

    EncryptionEngine<T> e(2);  // init with 2 parties

    auto [pub, sec] = e.gen_keys();

    message number{1, 2, 3, 4, 5, 6};
    message nine{9, 9, 9, 9, 9, 9};

    number.invntt_pow_invphi();
    nine.invntt_pow_invphi();

    auto ezeros = e.encrypt(number, pub);
    auto eone = e.encrypt(nine, pub);

    auto result = ezeros + eone;

    auto d0 = e.partial_decrypt(result, sec[0]);
    auto d1 = e.partial_decrypt({result.first, 0}, sec[1]);  // only one party decrypt the sum

    cipher d = d0 + d1;

    auto plain = e.eval_poly(d);  // get decrypted message

    for (auto i = 0; i < 6; i++) {
        EXPECT_EQ(mpz_get_ui(plain[i]), size_t(i+10));
    }
}

TEST(RingLWEAddition, 100PartiesMockedAddition01) {

    using T = uint64_t;

    using message = typename _params<T>::poly_p;
    using cipher = typename _params<T>::poly_q;

	auto numParties = 100;
    EncryptionEngine<T> e(numParties);  // init with numParties parties

    auto [pub, sec] = e.gen_keys();

	// party 0's message
	message num{1, 2, 3, 4, 5, 6};
	num.invntt_pow_invphi();
	auto result = e.encrypt(num, pub);

	for (auto i = 1; i < numParties; i++) { // add party 1-999's message
	    message num{0};
		num.invntt_pow_invphi();
		auto e_num = e.encrypt(num, pub);
		result = result + e_num;
	}

	cipher d;
	for (auto i = 0; i < numParties; i++) {
		if (i == 0) {
			auto di = e.partial_decrypt(result, sec[i]); // only one party decrypts (c,d)
			d = d + di;
		} else {
			auto di = e.partial_decrypt({result.first, 0}, sec[i]); // all other parties decrypt (c,0)
			d = d + di;
		}
	}

    auto plain = e.eval_poly(d);  // get decrypted message

    for (auto i = 0; i < 6; i++) {
        EXPECT_EQ(mpz_get_ui(plain[i]), size_t(i+1));
    }
}

TEST(RingLWEAddition, 100PartiesMockedAddition02) {

    using T = uint64_t;

    using message = typename _params<T>::poly_p;
    using cipher = typename _params<T>::poly_q;

	auto numParties = 100;
    EncryptionEngine<T> e(numParties);  // init with numParties parties

    auto [pub, sec] = e.gen_keys();

	// party 0's message
	message num{1, 2, 3, 4, 5, 6};
	num.invntt_pow_invphi();
	auto result = e.encrypt(num, pub);

	for (auto i = 1; i < numParties; i++) {
		if (i == 1) { // party i==1 adds {9, 9, 9, 9, 9, 9}
		    message num{9, 9, 9, 9, 9, 9};
			num.invntt_pow_invphi();
			auto e_num = e.encrypt(num, pub);
			result = result + e_num;
		} else { // party 1<i<1000 adds {0}
		    message num{0};
			num.invntt_pow_invphi();
			auto e_num = e.encrypt(num, pub);
			result = result + e_num;
		}
	}

	cipher d;
	for (auto i = 0; i < numParties; i++) {
		if (i == 0) {
			auto di = e.partial_decrypt(result, sec[i]); // only one party decrypts (c,d)
			d = d + di;
		} else {
			auto di = e.partial_decrypt({result.first, 0}, sec[i]); // all other parties decrypt (c,0)
			d = d + di;
		}
	}

    auto plain = e.eval_poly(d);  // get decrypted message

    for (auto i = 0; i < 6; i++) {
        EXPECT_EQ(mpz_get_ui(plain[i]), size_t(i+10));
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
