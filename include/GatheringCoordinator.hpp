#pragma once
/** Third Party */
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/unordered_set.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

/** STL */
#include <chrono>
#include <set>
#include <tuple>

/** Ligero */
#include "Common.hpp"
#include "Transport.hpp"
#include "ZkArgument.hpp"

using namespace std::chrono_literals;

/** Selecting a port for distributed verifying */
constexpr char* port = "5556";
constexpr int numberOfProofs = 21;

/** Various timeout */
const auto register_total_timeout = 10s;
const auto total_timeout = 60min;
constexpr uint64_t layeringDelayUs = 50000ull;

/** Exporting to S3? */
const bool s3export = false;

namespace ligero {
namespace GatheringCoordinator {
typedef std::vector<SocketId>::const_iterator SocketIterator;
std::unordered_map<SocketId, SocketId, boost::hash<SocketId>> taskAllocation;

/** This function serializes public data
 * @param socketId unique identifier attached to the data
 * @param pdvalues pointer to the public data content to be serialized
 * @param outMessage reference of a string that will contain the resulting
 * serialized message
 */
void serializePublicData(SocketId& socketId, PublicData* pdValues,
                         std::string& outMessage) {
  timers myTimers;
  myTimers.initialize_timer();

  myTimers.begin(1, "serialize_single", 0);

  std::stringstream ss;
  boost::archive::binary_oarchive oa(ss);

  oa << static_cast<int>(MessageType::GATHER_PUBLIC_DATA);
  oa << (*pdValues);

  outMessage = std::move(ss.str());

  myTimers.end(1, "serialize_single", 0, true);
}

/** This function gathers public data from the passive protocol
 * @param encryptedCoordinator the coordinator object handling the ceremony
 * (passive security)
 * @param id unique identifier for this party/verifier
 * @param outValue pointer to a structure that will host the public data
 * gathered
 */
template <typename T, size_t Degree, size_t NbPrimesP, size_t NbPrimesQ>
void generatePublicData(
    EncryptedCoordinator<T, Degree, NbPrimesP, NbPrimesQ>& encryptedCoordinator,
    SocketId& id, PublicData* outValue) {
  timers myTimers;
  myTimers.initialize_timer();

  myTimers.begin(1, "handle_single", 0);

  bool goodAlloc = false;

  while (!goodAlloc) {
    goodAlloc = true;
    try {
      encryptedCoordinator.updatePublicData(id, outValue);
    } catch (std::bad_alloc& ba) {
      std::cerr << "bad_alloc caught, prop detection: " << ba.what()
                << " client = " << id << '\n';
      goodAlloc = false;
    }
  }

  myTimers.end(1, "handle_single", 0, true);
}

/** Gathers ZKPs from parties and dispatches them to the designated distributed
 * verifiers
 * @param encryptedCoordinator the coordinator object handling the ceremony with
 * passive security
 * @param transport object handling communication with the parties
 * @param filename
 * @param bucketName
 * @param testingModulis
 * @param testing
 * @param ip
 * @param numVerifiers
 * @return Unit if successful
 */
template <typename FieldT, typename T, size_t Degree, size_t NbPrimesP,
          size_t NbPrimesQ>
expected<Unit> gatherAndVerifyArgumentsOfKnowledge(
    EncryptedCoordinator<T, Degree, NbPrimesP, NbPrimesQ> encryptedCoordinator,
    ZeroMQCoordinatorTransport& transport, const std::string filename,
    const std::string bucketName, size_t testingModulis = 0,
    bool testing = false, const char* ip = "127.0.0.1",
    size_t numVerifiers = 2) {
  taskAllocation.clear();

  std::chrono::milliseconds totalTimeRemaining = total_timeout;
  auto runIdentifier =
      boost::uuids::to_string(boost::uuids::random_generator()());

  std::vector<SocketId>& ids = transport.ids;
  size_t nbThreads = 90;
  size_t nbParties = ids.size();

  timers myTimers;
  myTimers.initialize_timer();

  SocketId socketId = boost::uuids::random_generator()();
  ZeroMQGatheringTransport verificationTransport(
      std::string("tcp://") + std::string(ip) + std::string(":") +
          std::string(port),
      numVerifiers);
  LOG(INFO) << "Gathering connecting "
            << std::string("tcp://") + std::string(ip) + std::string(":") +
                   std::string(port);

  /** Party's socket id that missing a verifier */
  std::set<SocketId> failedVerifier;

  if (testing) {
  } else {
    DBG("Parties: " << transport.ids.size());

    std::chrono::milliseconds timeRemaining = register_total_timeout;

    size_t count = 0;
    for (auto& partyId : transport.ids) {
      MessageType t;

      DBG("Party: " << partyId);

      auto startTime = std::chrono::steady_clock::now();
      if (timeRemaining < std::chrono::milliseconds(0))
        timeRemaining = std::chrono::milliseconds(0);
      auto maybeResponse = verificationTransport.awaitShortRegistration(
          MessageType::VERIFIER_READY, timeRemaining);

      // tick tock
      auto timeSpent = std::chrono::steady_clock::now() - startTime;
      timeRemaining -=
          std::chrono::duration_cast<std::chrono::milliseconds>(timeSpent);
      totalTimeRemaining -=
          std::chrono::duration_cast<std::chrono::milliseconds>(timeSpent);

      if (!hasError(maybeResponse)) {
        count++;
        auto [verifierId, verifierIPAddress] = getResult(maybeResponse);
        taskAllocation[partyId] = verifierId;

        LOG(INFO) << "Allocating Verifier " << verifierId << " with IP "
                  << verifierIPAddress << " to Party " << partyId;
      } else {
        failedVerifier.insert(partyId);
        LOG(INFO) << "Verifier for " << partyId << " failed to register";
      }
    }
    LOG(INFO) << count << " verifiers Registered.";
  }

  std::unordered_map<SocketId, PublicData*, boost::hash<SocketId>>
      publicDataValues;
  std::unordered_map<SocketId, std::string, boost::hash<SocketId>>
      publicDataMessages;
  std::unordered_map<SocketId, std::string, boost::hash<SocketId>>
      sigmaDataMessages;

  std::vector<PublicData> publicData(ids.size());
  std::vector<std::vector<zmq::message_t*>> proofs(ids.size());
  std::vector<SocketId> restartIds;
  std::set<SocketId> goodParties;

  std::vector<std::string> sigmaDataFileNames;

  /** Gather And Send Sigma Protocol Data */
  /** ======================================================================= */
  {
    for (auto& id : transport.ids) {
      auto startTime = std::chrono::steady_clock::now();
      if (totalTimeRemaining < std::chrono::milliseconds(0))
        totalTimeRemaining = std::chrono::milliseconds(0);
      auto responseMaybe = transport.awaitNextInput(
          totalTimeRemaining);  // MessageType::GATHER_SIGMA_DATA

      auto timeSpent = std::chrono::steady_clock::now() - startTime;
      totalTimeRemaining -=
          std::chrono::duration_cast<std::chrono::milliseconds>(timeSpent);

      if (hasError(responseMaybe)) {
        /** one or more parties is timed out, kick out and restart the protocol
         * again
         */
        std::vector<SocketId> goodPartyList(goodParties.begin(),
                                            goodParties.end());
        transport.broadcast(MessageType::PROTOCOL_RESTART, goodPartyList);
        transport.parties() = goodPartyList.size();

        /** also notify the verifier */
        std::vector<SocketId> goodVerifierList;
        for (auto& id : goodParties) {
          if (failedVerifier.find(id) == failedVerifier.end()) {
            goodVerifierList.emplace_back(taskAllocation[id]);
          }
        }
        verificationTransport.broadcast(MessageType::PROTOCOL_RESTART,
                                        goodVerifierList);
        verificationTransport.ids.clear();
        verificationTransport.verifiers() = goodVerifierList.size();
        return Error::RESTART;
      } else {
        auto response = getResult(responseMaybe);
        socketId = response.first;
        auto& sigmaProtocolPublicData = response.second;

        goodParties.insert(socketId);

        if (!testing) {
          SocketIterator it = std::find(ids.begin(), ids.end(), socketId);
          int socketIndex = std::distance<SocketIterator>(ids.begin(), it);

          DBG("Received Argument of Knowledge from Participant (" << socketId
                                                                  << ")");
          assert(it != ids.end());

          sigmaDataMessages[socketId].assign(
              sigmaProtocolPublicData->data<char>(),
              sigmaProtocolPublicData->size());

          /** Send Sigma data only if there is a verifier */
          if (std::find(failedVerifier.begin(), failedVerifier.end(),
                        socketId) == failedVerifier.end()) {
            LOG(INFO) << "Send SigmaData to Verifier ("
                      << taskAllocation[socketId] << " for Party " << socketId;
            verificationTransport.sendRaw(taskAllocation[socketId],
                                          sigmaProtocolPublicData);
          } else {
            LOG(INFO) << "Send " << boost::uuids::to_string(socketId)
                      << "'s sigma protocol data to local disk";
            std::string sigmaDataFileName =
                filename + "-" + runIdentifier + "-" +
                boost::uuids::to_string(socketId) + ".sigmadata";
            sendToFileRaw(sigmaProtocolPublicData, sigmaDataFileName);
            sigmaDataFileNames.emplace_back(std::move(sigmaDataFileName));
          }
        }
      }
    }
  }

  /** Tell all parties to send their proofs now */
  { transport.broadcast(MessageType::GATHER_PROOFS, ids); }

  /** Gather And Send Public Data */
  /** ==================================================================== */
  {
    if (testing) {
    } else {
      {
        {
          std::string timerVerificationPublicData =
              std::string("3.a. publicdata, public");
          myTimers.begin(1, timerVerificationPublicData.c_str(), 0);
          std::vector<std::thread> threads;
          size_t last = 0;

          for (size_t i = 0; i < ids.size(); i++) {
            try {
              publicDataValues[ids[i]] = new PublicData();
            } catch (...) {
              std::cout << "we failed to allocate public data" << std::endl;
              exit(1);
            }

            threads.emplace_back(
                generatePublicData<T, Degree, NbPrimesP, NbPrimesQ>,
                std::ref(encryptedCoordinator), std::ref(ids[i]),
                std::ref(publicDataValues[ids[i]]));

            usleep(layeringDelayUs);
          }

          for (size_t j = last; j < ids.size(); j++) {
            threads[j].join();
          }

          myTimers.end(1, timerVerificationPublicData.c_str(), 0, true);
        }

        {
          std::string timerVerificationSerialization =
              std::string("3.a. publicdata, serialization");
          myTimers.begin(1, timerVerificationSerialization.c_str(), 0);

          std::vector<std::thread> threads;
          size_t last = 0;

          for (size_t i = 0; i < ids.size(); i++) {
            threads.emplace_back(serializePublicData, std::ref(ids[i]),
                                 std::ref(publicDataValues[ids[i]]),
                                 std::ref(publicDataMessages[ids[i]]));

            if ((i + 1) % nbThreads == 0) {
              for (size_t j = last; j <= i; j++) {
                threads[j].join();
              }
              last = i + 1;
            }
          }

          for (size_t j = last; j < ids.size(); j++) {
            threads[j].join();
          }

          myTimers.end(1, timerVerificationSerialization.c_str(), 0, true);
        }

        {
          /** send them all */
          for (auto it = publicDataMessages.begin();
               it != publicDataMessages.end(); ++it) {
            /** if verifier failed in registration, send the public data to file
             */
            auto sid = it->first;
            if (failedVerifier.find(sid) != failedVerifier.end()) {
              LOG(INFO) << "Sending public data to file";
              std::string pdFileName = filename + "-" + runIdentifier + "-" +
                                       boost::uuids::to_string(sid) +
                                       ".publicdata";
              sendToFile(it->second, pdFileName);

              if (s3export) {
                std::string awsCommand = "aws s3 cp " + pdFileName + " s3://" +
                                         bucketName + "/" + pdFileName;
                system(awsCommand.c_str());
              }

              continue;
            } else {
              LOG(INFO) << "Sending serialized public data for party "
                        << it->first << " to verifier "
                        << taskAllocation[it->first];
              verificationTransport.sendMsg(taskAllocation[it->first],
                                            it->second);
            }
          }
        }
      }
    }
  }

  /** Receive all proofs */

  std::vector<size_t> proofsReceived(nbParties, 0);
  bool complete = false;
  size_t testCount = 0;

  goodParties.clear();

  do {
    std::stringstream ssId, ss;

    auto startTime = std::chrono::steady_clock::now();
    if (totalTimeRemaining < std::chrono::milliseconds(0))
      totalTimeRemaining = std::chrono::milliseconds(0);
    auto responseMaybe = transport.awaitNextInput(
        totalTimeRemaining);  // MessageType::GATHER_PROOFS

    auto timeSpent = std::chrono::steady_clock::now() - startTime;
    totalTimeRemaining -=
        std::chrono::duration_cast<std::chrono::milliseconds>(timeSpent);

    if (hasError(responseMaybe)) {
      LOG(INFO) << "Got error while receiving proofs";

      auto err = getError(responseMaybe);

      /** awaitNextInput() will only return timed out or unknown error */
      if (err == Error::TIMED_OUT) {
        /** one or more parties is timed out, kick out and restart the protocol
         * again
         */
        std::vector<SocketId> goodPartyList(goodParties.begin(),
                                            goodParties.end());
        transport.broadcast(MessageType::PROTOCOL_RESTART, goodPartyList);
        transport.parties() = goodPartyList.size();

        // also notify the verifier
        std::vector<SocketId> goodVerifierList;
        for (auto& id : goodParties) {
          if (failedVerifier.find(id) == failedVerifier.end()) {
            goodVerifierList.emplace_back(taskAllocation[id]);
          }
        }
        verificationTransport.broadcast(MessageType::PROTOCOL_RESTART,
                                        goodVerifierList);
        verificationTransport.ids.clear();
        verificationTransport.verifiers() = goodVerifierList.size();
        return Error::RESTART;
      }
      LOG(INFO) << "Caught unexpected error when receving proofs";
      return err;
    }

    /** No error happened */
    else {
      auto response = getResult(responseMaybe);
      socketId = response.first;
      auto& argumentOfKnowledge = response.second;

      goodParties.insert(socketId);

      if (!testing) {
        SocketIterator it = std::find(ids.begin(), ids.end(), socketId);

        /** if we are receiving from unknown party */
        if (it == ids.end()) {
          LOG(INFO) << "Receiving from unknown party "
                    << boost::uuids::to_string(socketId);
          continue;
        }

        int socketIndex = std::distance<SocketIterator>(ids.begin(), it);
        DBG("Received Argument of Knowledge from Participant (" << socketId
                                                                << ")");

        proofsReceived[socketIndex]++;
        proofs[socketIndex].push_back(argumentOfKnowledge);

        complete = std::all_of(proofsReceived.begin(), proofsReceived.end(),
                               [](size_t received) { return received == 21; });

        {
          /** Verifier is not avaliable */
          if (failedVerifier.find(socketId) == failedVerifier.end()) {
            LOG(INFO) << "Send Proof for Modulus "
                      << proofsReceived[socketIndex] - 1 << " to Verifier ("
                      << taskAllocation[socketId] << " for Party " << socketId;
            verificationTransport.sendRaw(taskAllocation[socketId],
                                          argumentOfKnowledge);
          }

          std::stringstream ssId;
          ssId << socketId;
          std::string ssocketId(ssId.str());
          std::string proofNb = std::to_string(proofsReceived[socketIndex]);

          printAtXY(1, 4,
                    std::string("verifying proof:") + proofNb +
                        std::string("_") + ssocketId + std::string("_") +
                        proofNb);
        }
      }
    }
  } while (!complete);

  /** everyone is good now */
  goodParties.clear();

  std::vector<size_t> proofsSuccess(numberOfProofs, 0);
  std::vector<bool> goodVerifierFailingToReport(
      nbParties, true); /** We are deemed to have failed unless and until we */
                        /** successfully verify */

  verificationTransport.broadcast(MessageType::GATHER_REPORTS,
                                  verificationTransport.ids);

  size_t verificationSize = ids.size() - failedVerifier.size();
  if (!testing) {
    /** first receive report from good verifiers */
    printAtXY(1, 1 + firstline,
              std::string("Verifiers Reports: (") +
                  std::to_string(verificationSize) + std::string(" parties)"));
    printAtXY(1, 2 + firstline, std::string("=================="));
    printAtXY(1, 3 + firstline, std::string("Moduli:"));
    for (size_t i = 0; i < numberOfProofs; i++) {
      printAtXY(i * displayStride + displayOffset, 3 + firstline,
                std::to_string(i));
    }

    printAtXY(1, 4 + firstline, std::string("Successful Verifications:"));

    for (size_t i = 0; i < numberOfProofs; i++) {
      printAtXY(i * displayStride + displayOffset, 4 + firstline,
                std::to_string(0));
    }

    LOG(INFO) << "Receiving Verifier Reports";

    for (size_t idx = 0; idx < verificationSize; idx++) {
      auto startTime = std::chrono::steady_clock::now();
      if (totalTimeRemaining < std::chrono::milliseconds(0))
        totalTimeRemaining = std::chrono::milliseconds(0);
      auto maybeResponse = verificationTransport.awaitNextInput(
          totalTimeRemaining);  // MessageType::VERIFIER_REPORT

      auto timeSpent = std::chrono::steady_clock::now() - startTime;
      totalTimeRemaining -=
          std::chrono::duration_cast<std::chrono::milliseconds>(timeSpent);

      if (hasError(maybeResponse)) {
        /** we will handle failed verifier later */
        LOG(INFO) << "Unresponsive Verifier Detected.";
        break;
      } else {
        auto response = getResult(maybeResponse);

        LOG(INFO) << "Received Valid Verifier Report From " << response.first;
        std::vector<int> report;

        auto it = std::find_if(ids.begin(), ids.end(), [&](const auto& id) {
          return taskAllocation[id] == response.first;
        });

        if (it == ids.end()) {
          LOG(INFO) << "Cannot Identify Verifier.";
          /** don't mark it as successful and handle it later */
          continue;
        }

        auto socketIndex = std::distance(ids.begin(), it);

        LOG(INFO) << "Report received from Verifier " << response.first
                  << " was for Party " << ids[idx];

        try {
          MessageType t;
          std::stringstream ss;
          ss.write(response.second->data<char>(), response.second->size());

          boost::archive::binary_iarchive ia(ss);
          ia >> t;

          if (t != MessageType::VERIFIER_REPORT) {
            LOG(INFO) << "Out of sync header received";
            continue;
          }

          ia >> report;
        } catch (...) {
          LOG(INFO) << "Error deserializing report";
          continue;
        }

        LOG(INFO) << "Deserialized Report From " << response.first;

        goodVerifierFailingToReport[socketIndex] = false;
        LOG(INFO) << "Verifier " << response.first << " reporting for party "
                  << ids[socketIndex];

        auto verifySuccess =
            std::all_of(report.begin(), report.end(), [](int x) { return x; });
        if (verifySuccess) {
          goodParties.insert(ids[idx]);
        } else {
          LOG(INFO) << "Verification for party " << ids[idx]
                    << " failed. Adding to kickout list";
        }

        for (size_t i = 0; i < numberOfProofs; i++) {
          proofsSuccess[i] += report[i];
        }

        std::string msg("");
        for (size_t i = 0; i < numberOfProofs; i++) {
          msg += std::to_string(i) + std::string(": ") +
                 std::to_string(proofsSuccess[i]);
          if (i < 20) msg += std::string(", ");
        }

        LOG(INFO) << msg;

        for (size_t i = 0; i < numberOfProofs; i++) {
          printAtXY(i * displayStride + displayOffset, 4 + firstline,
                    std::to_string(proofsSuccess[i]));
        }
      }
    }

    auto reportReceiptNb = std::count(goodVerifierFailingToReport.begin(),
                                      goodVerifierFailingToReport.end(), false);
    LOG(INFO) << "Successfully received report from " << reportReceiptNb
              << " verifiers";

    /** check if some report failed and need restart */
    if ((goodParties.size() != reportReceiptNb)) {
      LOG(INFO) << "Restart protocol with " << goodParties.size() << " parties";

      std::vector<SocketId> empty;

      std::vector<SocketId> restartVerifiers;
      std::transform(goodParties.begin(), goodParties.end(),
                     std::back_inserter(restartVerifiers),
                     [&](const auto& id) { return taskAllocation[id]; });
      verificationTransport.broadcast(MessageType::PROTOCOL_RESTART,
                                      restartVerifiers);
      verificationTransport.ids.clear();
      verificationTransport.verifiers() = restartVerifiers.size();

      std::vector<SocketId> goodPartyList(goodParties.begin(),
                                          goodParties.end());
      transport.broadcast(MessageType::PROTOCOL_RESTART, goodPartyList);
      transport.parties() = goodPartyList.size();
      return Error::RESTART;
    } else {
      transport.broadcast(MessageType::END, transport.ids);
      verificationTransport.broadcast(MessageType::END,
                                      verificationTransport.ids);
    }

    /** Save failed sigma data to local disk */
    for (const std::string& file : sigmaDataFileNames) {
      LOG(INFO) << "Pushing " << file << " AWS";
      std::string awsCommand = "aws s3 cp " + file + " s3://" + bucketName +
                               "/" + file +
                               std::string(" > s3_upload.log 2>&1");
      system(awsCommand.c_str());
    }

    /** Some verifiers failed, saving proofs to local disk */
    for (size_t i = 0; i < ids.size(); i++) {
      if (goodVerifierFailingToReport[i]) {
        auto it = sigmaDataMessages.find(ids[i]);

        if (it != sigmaDataMessages.end()) {
          LOG(INFO) << "Sending sigma data to file";
          std::string pdFileName = filename + "-" + runIdentifier + "-" +
                                   boost::uuids::to_string(ids[i]) +
                                   ".sigmadata";
          sendToFileString((const std::string&)it->second, pdFileName);

          if (s3export) {
            std::string awsCommand = "aws s3 cp " + pdFileName + " s3://" +
                                     bucketName + "/" + pdFileName +
                                     std::string(" > s3_upload.log 2>&1");
            system(awsCommand.c_str());
          }
        }

        if (it != publicDataMessages.end()) {
          LOG(INFO) << "Sending public data to file";
          std::string pdFileName = filename + "-" + runIdentifier + "-" +
                                   boost::uuids::to_string(ids[i]) +
                                   ".publicdata";
          sendToFileString((const std::string&)it->second, pdFileName);

          if (s3export) {
            std::string awsCommand = "aws s3 cp " + pdFileName + " s3://" +
                                     bucketName + "/" + pdFileName +
                                     std::string(" > s3_upload.log 2>&1");
            system(awsCommand.c_str());
          }
        }

        for (size_t j = 0; j < proofs[i].size(); j++) {
          std::string proofFileName = filename + "-" + runIdentifier + "-" +
                                      boost::uuids::to_string(ids[i]) +
                                      ".proof." + std::to_string(j);
          LOG(INFO) << "Saving proofs " << j << "to disk";
          sendToFileRaw(proofs[i][j], proofFileName);

          if (s3export) {
            std::string awsCommand = "aws s3 cp " + proofFileName + " s3://" +
                                     bucketName + "/" + proofFileName +
                                     std::string(" > s3_upload.log 2>&1");
            LOG(INFO) << "Pushing proofs " << j << "to s3 storage";
            system(awsCommand.c_str());
          }
        }
      }
    }
  }

  return Unit{};
}

}  // namespace GatheringCoordinator

}  // namespace ligero
