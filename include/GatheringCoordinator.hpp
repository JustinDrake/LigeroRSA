#pragma once
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/utility.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

/* Ligero */

#include <tuple>

/* Ligero */
#include "Common.hpp"
#include "ZkArgument.hpp"
#include "Transport.hpp"

namespace ligero
{
    using pd = PublicData;
    namespace GatheringCoordinator
    {
        typedef std::vector<SocketId>::const_iterator SocketIterator;

        template<typename FieldT>
        void gatherAndSaveArgumentsOfKnowledge(ZeroMQCoordinatorTransport& transport, const std::string filename, size_t testingModulis = 0, bool testing = false)  
        {
            std::string stagPublicData("_PublicData_");
            std::string stagProof("_Proof_");
            std::vector<SocketId>& ids = transport.ids;
            size_t nbParties = ids.size();

            if (testing) {
                for (auto id : ids) {
                    system((std::string("./distributed_verifier Test test &").c_str()));
                }                
            } else {
                for (auto id : ids) {
                    std::stringstream ssId;
                    ssId << id;
                    std::string partyIdSs(ssId.str());
                    system((std::string("./distributed_verifier Prod ") + partyIdSs + std::string(" &")).c_str());
                }                
            }

            {
                boost::dynamic_bitset<> hasReceivedInput (nbParties, false);
                bool complete = false;

                // Tell all parties to send their public data
                transport.broadcast(MessageType::GATHER_PUBLIC_DATA, ids);
                do {
                    std::stringstream ssId, ss;
                    auto [socketId, publicData] = transport.awaitNextInputST<PublicData>(MessageType::GATHER_PUBLIC_DATA);
                    
                    ssId << socketId;
                    std::string ssocketId(ssId.str());

                    if (testing) {ssocketId.assign("test");}

                    sendToFile(publicData, filename + stagPublicData + ssocketId);

                    if (!testing) {
                        SocketIterator it = std::find(ids.begin(), ids.end(), socketId);
                        int socketIndex = std::distance<SocketIterator>(ids.begin(), it);

                        DBG("Received Public Data from Participant (" << socketId << ")");
                        assert (it != ids.end());
                        assert (!hasReceivedInput[socketIndex]);

                        hasReceivedInput[socketIndex] = true;
                        complete = hasReceivedInput.all();

                    } else {
                        ids.emplace_back(socketId);
                        complete = true;
                    }
                } while (!complete);
            }

            // Tell all parties to send their proofs now
            {
                transport.broadcast(MessageType::GATHER_PROOFS, ids);
            }

            // Receive all proofs
            {
                // std::vector<std::vector<ligero::zksnark::FullFSTranscript<FieldT>>> aggregatedProofs(ids.size());
                std::vector<size_t> proofsReceived(nbParties, 0);
                bool complete = false;
                size_t testCount = 0;

                if (testing) {
                    std::string command = std::string("rm -rf ") + filename + stagProof + std::string("test");
                    system(command.c_str());
                    }

                do {
                    std::stringstream ssId, ss;
                    auto [socketId, argumentOfKnowledge] = transport.awaitNextInputST<ligero::zksnark::FullFSTranscript<FieldT>>(MessageType::GATHER_PROOFS);

                    if (!testing) {
                        SocketIterator it = std::find(ids.begin(), ids.end(), socketId);
                        int socketIndex = std::distance<SocketIterator>(ids.begin(), it);

                        DBG("Received Argument of Knowledge from Participant (" << socketId << ")");
                        assert (it != ids.end());

                        // aggregatedProofs[socketIndex].emplace_back(argumentOfKnowledge);
                        proofsReceived[socketIndex]++;

                        complete = true;
                        for (auto id: ids) {
                            SocketIterator it = std::find(ids.begin(), ids.end(), id);
                            int socketIndex = std::distance<SocketIterator>(ids.begin(), it);
                            if (proofsReceived[socketIndex]<21) {
                                complete = false;
                                break;
                            }
                        }

                        {
                            std::stringstream ssId;
                            ssId << socketId;
                            std::string ssocketId(ssId.str());
                            std::string proofNb = std::to_string(proofsReceived[socketIndex]);

                            //std::cout << filename + stagProof + ssocketId + std::string("_") + proofNb << std::endl;
                            sendToFile(argumentOfKnowledge, filename + stagProof + ssocketId + std::string("_") + proofNb);
                        }

                    } else {
                      //  aggregatedProofs[0].emplace_back(argumentOfKnowledge);
                        testCount++;
                        if (testCount == testingModulis) complete = true;

                        {
                            std::string ssocketId("test");
                            std::string proofNb = std::to_string(testCount);

                            //std::cout << filename + stagProof + ssocketId + std::string("_") + proofNb << std::endl;
                            sendToFile(argumentOfKnowledge, filename + stagProof + ssocketId + std::string("_") + proofNb);
                        }
                    }

                } while (!complete);

                size_t size_modulis = ids.size();
                if (testing) size_modulis = 1;
            }

            return;
        }

    } // namespace coordinator

} // namespace ligero
