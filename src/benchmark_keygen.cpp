#include <thread>
#include <chrono>

#include "Transport.hpp"
#include "EncryptedClient.hpp"
#include "EncryptedCoordinator.hpp"
#include "Ligero.hpp"

#define NDEBUG 1

/** Run the protocol over a given number type. */
template <typename NumberType>
void runTest (int numParties, int numBits, int numCandidates) {
    std::vector<std::thread> threads;
    static constexpr int n = 1024;
    static constexpr int m = 1024;
    static constexpr double chiStd = 3.2;
    const NumberType q = (NumberType(1) << 150) - 1;
    const int packingFactor = numCandidates;

    // We'll generate for this test a shared public key and a set of partial
    // secrete keys for our participants.
    ligero::Matrix<NumberType> Aprime(n, m);
    ligero::Matrix<NumberType> b(packingFactor, m);
    std::vector<ligero::SecretKey<NumberType>> secretKeys;
    std::vector<ligero::SocketId> ids;

    for (int i = 0; i < numParties; ++i) {
        auto Aprime_shares = ligero::encryption::generateAprime(n, m, q);
        Aprime += Aprime_shares;
    }

    auto e = ligero::math::sampleRandomMatrix<NumberType>(packingFactor, m, q, chiStd);
    auto sprime = ligero::encryption::generateSprime(n, packingFactor, q);

    auto timeBegin = std::chrono::high_resolution_clock::now();
    auto bi = ligero::encryption::computeB(q, e, sprime, Aprime);
    auto timeEnd = std::chrono::high_resolution_clock::now();

    long result = std::chrono::duration_cast<std::chrono::microseconds>( timeEnd - timeBegin ).count();

    std::cout << "EXECUTION TIME = " << result << std::endl;
}

//==============================================================================
int main(int argc, char** argv)
{
  runTest<ligero::MPInt>(1, 100, 5);
  return 0;
}
