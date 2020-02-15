#include <boost/python.hpp>
#include "EncryptionEngine.hpp"
#include <sstream>

//using Ligero = EncryptionEngine<2,3>;

class Ligero {
public: 
        Ligero () {
            // p = 2**128 - 1
            const BigInt p1 = BigInt("340282366920938463463374607431768211455"); 
            // q = 2**512 - 1
            const BigInt q1 = BigInt("1942668892225729070919461906823518906642406839052139521251812409738904285205208498175");
            ee = EncryptionEngine<2, 3>(p1, q1); 
            generateKeys();
        }

        void generateKeys() {
            keys = ee.generateKeys();
        }
        
        std::string getSecretKey() {
            auto s = keys.first;
            std::stringstream ss;
            LOG(DEBUG) << s;
            ss << s;
            return ss.str();
        }

        std::string getPublicKey() {
            auto s = keys.second;
            std::stringstream ss;
            LOG(DEBUG) << s;
            ss << s;
            return ss.str();
        }

        std::string encrypt(std::string msg) {
            encryptedMsg = ee.encrypt(BigInt(msg), keys.second);
            LOG(DEBUG) << encryptedMsg;
            std::stringstream ss;
            ss << encryptedMsg;
            return ss.str();
        }

        std::string decrypt() {
            decryptedMsg = ee.decrypt(encryptedMsg, keys.first);
            LOG(DEBUG) << encryptedMsg;
            std::stringstream ss;
            ss << decryptedMsg;
            return ss.str();
        }

private:
        Eigen::Matrix<BigInt, 1, Eigen::Dynamic> encryptedMsg;
        BigInt decryptedMsg;
        EncryptionEngine<2, 3> ee;
        KeyPair keys;
};

BOOST_PYTHON_MODULE(ligero)
{
    using namespace boost::python;
    class_<Ligero>("Ligero")
        .def("generateKeys", &Ligero::generateKeys)
        .def("getSecretKey", &Ligero::getSecretKey)
        .def("getPublicKeyKey", &Ligero::getPublicKey)
        .def("encrypt", &Ligero::encrypt)
        .def("decrypt", &Ligero::decrypt)
    ;
}
