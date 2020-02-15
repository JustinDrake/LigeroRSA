#include <thread>
#include <future>
#include <boost/serialization/utility.hpp>

#include "LatticeEncryption.hpp"
#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"
#include "gtest/gtest.h"

using namespace ligero;
using ligero::lattice::operator*;
using ligero::lattice::operator+;

constexpr auto degree = 1 << 16;
constexpr auto p = 13;
constexpr auto q = 29;
constexpr auto parties = 2;
constexpr auto sigma = 8;
constexpr auto lambda = 80;
constexpr auto tau_limit_bit = 340;
constexpr auto pbs = 1000;
mpz_class tau = 1000;

constexpr auto ip = "tcp://127.0.0.1:5555";

using P = nfl::poly_p<uint64_t, degree, p>;
using Q = nfl::poly_p<uint64_t, degree, q>;

std::pair<Q, Q> pair_add(const std::pair<Q, Q>& a, const std::pair<Q, Q>& b) {
    return { a.first + b.first, a.second + b.second };
}

void serve_vecplus() {
    // Registration
    ProtocolConfig<uint64_t> config(parties, p, q, degree, sigma, lambda, tau_limit_bit, pbs);
    ZeroMQCoordinatorTransport trans(ip, parties);
    auto ids = hostRegistration(trans, config);
    // End registration

    // Host keygen
    hostGenerateKeyPair<uint64_t, degree, q>(trans, ids);
    // End host keygen

    std::pair<Q, Q> acc = {0, 0};
    auto c = trans.awaitAggregateVectorInput<std::pair<Q, Q>>(MessageType::MESSAGE, ids, pair_add, acc);
    trans.broadcast(MessageType::MESSAGE, ids, c);
    auto d = trans.awaitAggregateVectorInput<Q>(MessageType::MESSAGE, ids, std::plus<Q>(), 0);
    trans.broadcast<Q>(MessageType::MESSAGE, ids, d);
}

void join_vecplus(int index) {
    auto id = boost::uuids::random_generator()();
    ZeroMQClientTransport trans(ip, id);
    auto config = registerAndAwaitConfiguration<uint64_t>(trans, ip);
    auto e = lattice::LatticeEncryption<uint64_t, degree, p, q>(config);

    EncryptedClient<uint64_t, uint64_t, degree, p, q> ec;
    auto [pub, sec] = ec.generateKeyPair(trans, e.chi());

    P message{1, 2, 3, 4, 5, 6, 7, 8};

    auto c = e.encrypt(message, pub);

    trans.send(MessageType::MESSAGE, c);
    auto sum = trans.awaitReply<lattice::cipher<uint64_t, degree, q>>(MessageType::MESSAGE);

    Q dec;
    if (index == 0) {
        dec = e.partial_decrypt(sum, sec);
    }
    else {
        dec = e.partial_decrypt({sum.first, 0}, sec);
    }

    trans.send(MessageType::MESSAGE, dec);
    auto m = trans.awaitReply<Q>(MessageType::MESSAGE);

    auto mm = e.eval_poly(m, tau);

    for (auto i = 0; i < 8; i++) {
        ASSERT_EQ(mpz_get_ui(mm[i]), (parties * (i + 1)) % tau);
    }
}

void serve_vecprod() {
    // Registration
    ProtocolConfig<uint64_t> config(parties, p, q, degree, sigma, lambda, tau_limit_bit, pbs);
    ZeroMQCoordinatorTransport trans(ip, parties);
    auto ids = hostRegistration(trans, config);
    // End registration

    // Host keygen
    hostGenerateKeyPair<uint64_t, degree, q>(trans, ids);
    // End host keygen

    std::pair<Q, Q> acc = {0, 0};
    auto c = trans.awaitAggregateVectorInput<std::pair<Q, Q>>(MessageType::MESSAGE, ids, pair_add, acc);
    trans.broadcast(MessageType::MESSAGE, ids, c);
    auto c1 = trans.awaitAggregateVectorInput<std::pair<Q, Q>>(MessageType::MESSAGE, ids, pair_add, acc);
    trans.broadcast(MessageType::MESSAGE, ids, c1);
    auto d = trans.awaitAggregateVectorInput<Q>(MessageType::MESSAGE, ids, std::plus<Q>(), Q{0});
    trans.broadcast(MessageType::MESSAGE, ids, d);
}

void join_vecprod(int index) {
    auto id = boost::uuids::random_generator()();
    ZeroMQClientTransport trans(ip, id);
    auto config = registerAndAwaitConfiguration<uint64_t>(trans, ip);
    auto e = lattice::LatticeEncryption<uint64_t, degree, p, q>(config);

    EncryptedClient<uint64_t, uint64_t, degree, p, q> ec;
    auto [pub, sec] = ec.generateKeyPair(trans, e.chi());

    P message{1, 2, 3, 4, 5, 6, 7, 8};
    P _product{1};

    auto c = e.encrypt(message, pub);

    trans.send(MessageType::MESSAGE, c);
    auto sum = trans.awaitReply<lattice::cipher<uint64_t, degree, q>>(MessageType::MESSAGE);

    std::array<mpz_class, degree> taus;
    std::fill(taus.begin(), taus.end(), tau);
    sum = e.product(sum, _product, taus, pub);

    trans.send(MessageType::MESSAGE, sum);
    sum = trans.awaitReply<lattice::cipher<uint64_t, degree, q>>(MessageType::MESSAGE);

    Q dec;
    if (index == 0) {
        dec = e.partial_decrypt(sum, sec);
    }
    else {
        dec = e.partial_decrypt({sum.first, 0}, sec);
    }

    trans.send(MessageType::MESSAGE, dec);
    auto m = trans.awaitReply<Q>(MessageType::MESSAGE);

    auto mm = e.eval_poly(m, tau);

    ASSERT_EQ(mpz_get_ui(mm[0]), (parties * parties) % tau);
    for (auto i = 1; i < degree; i++) {
        ASSERT_EQ(mpz_get_ui(mm[i]), 0);
    }
}

void join_vecprod_vect(int index) {
    auto id = boost::uuids::random_generator()();
    ZeroMQClientTransport trans(ip, id);
    auto config = registerAndAwaitConfiguration<uint64_t>(trans, ip);
    auto e = lattice::LatticeEncryption<uint64_t, degree, p, q>(config);

    EncryptedClient<uint64_t, uint64_t, degree, p, q> ec;
    auto [pub, sec] = ec.generateKeyPair(trans, e.chi());

    std::vector<mpz_class> msg(degree, 1); // {1, 1, 1, 1, 1, 1, 1, 1}
    P _product;
    auto *p_arr = new std::array<mpz_class, degree>{};
    std::iota(p_arr->begin(), p_arr->end(), 1);
    _product.set_mpz(*p_arr);
    delete(p_arr);

    auto c = e.encrypt(msg, pub);

    trans.send(MessageType::MESSAGE, c);
    auto sum = trans.awaitReply<lattice::cipher<uint64_t, degree, q>>(MessageType::MESSAGE);

    std::array<mpz_class, degree> taus;
    std::fill(taus.begin(), taus.end(), tau);
    sum = e.product(sum, _product, taus, pub);

    trans.send(MessageType::MESSAGE, sum);
    sum = trans.awaitReply<lattice::cipher<uint64_t, degree, q>>(MessageType::MESSAGE);

    Q dec;
    if (index == 0) {
        dec = e.partial_decrypt(sum, sec);
    }
    else {
        dec = e.partial_decrypt({sum.first, 0}, sec);
    }

    trans.send(MessageType::MESSAGE, dec);
    auto m = trans.awaitReply<Q>(MessageType::MESSAGE);

    auto mm = e.eval_poly(m, tau);

    for (auto i = 1; i < degree; i++) {
        ASSERT_EQ(mpz_get_ui(mm[i]), parties * parties * (i + 1) % tau);
    }
}

TEST(DistributedEncryption, 2PartiesAddition) {

    std::vector<std::thread> threads;

    for (auto i = 0; i < parties; i++) {
        threads.emplace_back(join_vecplus, i);
    }

    threads.emplace_back(serve_vecplus);

    for (auto i = 0; i < parties + 1; i++) {
        threads[i].join();
    }
}

TEST(DistributedEncryption, 2PartiesProduct) {

    std::vector<std::thread> threads;

    for (auto i = 0; i < parties; i++) {
        threads.emplace_back(join_vecprod, i);
    }

    threads.emplace_back(serve_vecprod);

    for (auto i = 0; i < parties + 1; i++) {
        threads[i].join();
    }
}

TEST(DistributedEncryption, 2PartiesVectorProduct) {

    std::vector<std::thread> threads;

    for (auto i = 0; i < parties; i++) {
        threads.emplace_back(join_vecprod_vect, i);
    }

    threads.emplace_back(serve_vecprod);

    for (auto i = 0; i < parties + 1; i++) {
        threads[i].join();
    }
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
