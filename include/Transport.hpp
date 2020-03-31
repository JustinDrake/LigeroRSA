#pragma once

#include "Common.hpp"

#include <algorithm>
#include <bitset>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

#include <boost/dynamic_bitset.hpp>
#include <boost/multiprecision/gmp.hpp>
#include <boost/program_options.hpp>
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>
#include <boost/asio.hpp>
#include <cassert>
#include <csignal>
#include <chrono>
#include <functional>
#include <string>
#include <set>
#include <thread>
#include <future>

#include <unordered_map>
#include <zmq.hpp>

#include "nfl.hpp"

#include "Math.hpp"

namespace ligero
{
typedef std::vector<SocketId>::const_iterator SocketIterator;

struct throughput_test_config {
    std::chrono::milliseconds duration;
    std::chrono::milliseconds wait_time;
    std::chrono::milliseconds cleanup_time;
    size_t message_size; /* in Bytes */
};


//==============================================================================
/** An interfacee for coordinating an MPC protocol over ZeroMQ. */
class ZeroMQCoordinatorTransport
{
public:
    timers myTimers;
    bool receiving;
    static std::vector<SocketId> ids;
    static std::unordered_map<SocketId, PartyCommitment, boost::hash<SocketId>> partiesCommitments;

    static boost::filesystem::ofstream ofs;
    static boost::archive::binary_oarchive scriptOArchive;
                
    //==============================================================================
    ZeroMQCoordinatorTransport(const std::string& uri, int parties
        , std::chrono::milliseconds receive_timeout = std::chrono::minutes(10))
        : context(kDefaultZeroMQIOThreads, 10),
          socket(context, ZMQ_ROUTER),
          numParties(parties),
          receive_timeout(receive_timeout)
    {
        socket.setsockopt(ZMQ_TCP_KEEPALIVE, 1);
        socket.setsockopt(ZMQ_TCP_KEEPALIVE_INTVL, 30);
        socket.setsockopt(ZMQ_LINGER, 0);

        socket.bind(uri);
        restartCount = 0;
    }

    ~ZeroMQCoordinatorTransport()  {
        ofs.close();
    }

    /** Waits for registration messages over the socket; clients registering for the protocol. */
    expected<Unit> awaitRegistration() {
        int registeredParties = 0;

        auto registrationEndTime = std::chrono::steady_clock::now() + kDefaultRegistrationWaitTime;

        auto timeNow = std::chrono::steady_clock::now();

        std::vector<zmq::pollitem_t> items = {
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };

        // clear ids before start
        ids.clear();

        while (timeNow < registrationEndTime
            && static_cast<int>(ids.size()) < numParties) {

            timeNow = std::chrono::steady_clock::now();
            zmq::message_t identity;
            zmq::message_t request;

            SocketId sid;
            std::stringstream ss;

            {

                auto rc = zmq::poll(items, kDefaultRegistrationWaitTime);
                if (rc && (items[0].revents & ZMQ_POLLIN)) {
                    socket.recv(&identity);
                } else {
                    LOG(INFO) << "Registration Timedout: no data transfer to coordinator";
                    continue;
                }

                // Registration IP Address
                // Do not count it toward communication
                socket.recv(&request);

                // Assign the socket id
                memcpy(&sid, identity.data<char>(), identity.size());

                // check the id is one of survivor
                auto sid_pos = std::find(survivor.begin(), survivor.end(), sid);
                if (sid_pos == survivor.end()) {
                    continue;
                }

                // Read the message data
                ss.write(request.data<char>(), request.size());
                
                MessageType header;
                boost::archive::binary_iarchive ia(ss);
                PartyCommitment input;

                try {
                    ia >> header;
                    ia >> input;
                }
                catch (...) {
                    LOG(INFO) << "Error deserializing message, kicking out party " << sid;
                    // since `sid` is in the survivor list, kick it out
                    survivor.erase(sid_pos);
                    continue;
                }
                
                if (header != MessageType::ID_PARTY) {
                    LOG(INFO) << "Party sent wrong header, adding to kickout list";
                    survivor.erase(sid_pos);
                    continue;
                }

                // Cross-referencing IDs and IPs
                //LOG(INFO) << "CROSS_REFERENCING,ID," << sid << "," << input.ipAddress();
                //LOG(INFO) << "REGISTRATION, PROCESS 1, PARTIES REMAINING = " << (numParties - (++registeredParties));

                partiesCommitments[sid] = std::move(input);
                received[sid] = true;
                ids.push_back(std::move(sid));
            }
        }

        if (ids.size() < kMinAmountOfPartiesRequired) {
            LOG(INFO) << "Registration failed: Too few parties left";
            return Error::TOO_FEW_PARTIES;
        }

        LOG(INFO) << "Registration completed for " << ids.size() << " out of " << numParties;
        numParties = ids.size();
        DBG("Writing numParties");
        scriptOArchive << numParties;

        return Unit{};
    }

    std::vector<SocketId> hostThroughputTest(
        throughput_test_config config,
        std::function<bool(size_t)> predicate  /* byte/s -> bool, test if client survive */
    ) {
        using count_map = std::unordered_map<SocketId, size_t, boost::hash<SocketId>>;

        count_map package_received;

        zmq::message_t identity;
        zmq::message_t request;
        SocketId sid;

        std::vector<zmq::pollitem_t> items = { 
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };

        auto test_end_time = std::chrono::steady_clock::now() + config.wait_time;
        int connected_clients = 0;
        while (std::chrono::steady_clock::now() < test_end_time 
            && package_received.size() < numParties
            && connected_clients < numParties) {

            auto nb_msg_queued = zmq::poll(items, 100);
            if (nb_msg_queued == 0) {
                continue;
            }

            socket.recv(&identity);
            socket.recv(&request);

            // Assign the socket id
            memcpy(&sid, identity.data<char>(), identity.size());

            // already receiving packages
            if (package_received.find(sid) != package_received.end()) {
                DBG("WARNING: duplicate registration detected");
            }
            else {
                // joining test, perform sanity check
                std::stringstream ss;
                ss.write(request.data<char>(), request.size());
                MessageType header;
                boost::archive::binary_iarchive ia(ss);

                try {
                    ia >> header;
                }
                catch(...) {
                    LOG(INFO) << "Deserialization fail detected";
                    continue;
                }

                if (header == MessageType::THROUGHPUT_TEST_START) {
                    connected_clients++;
                    // done handshaking, initialize test data
                    package_received[sid] = 1;

                    // wait untill all parties have joined
                    LOG(INFO) << "THROUGHPUT TEST: new client joined! (ID: " << sid << "), remaining: " << numParties - connected_clients;
                    continue;
                }
                else {
                    // simply discard unknown message
                    DBG("WARNING: coordinator discard unknown message");
                    LOG(INFO) << "WARNING: coordinator discard unknown message";
                    numParties--;
                    continue;
                }
            }
        }
        LOG(INFO) << connected_clients << " parties are joining throughput test";
        numParties = connected_clients;

        // get all participants
        std::vector<SocketId> participants(numParties);
        std::transform(package_received.begin(), package_received.end(), 
            participants.begin(), [](const std::pair<SocketId, size_t>& x) { return x.first; });

        // quit gracefully if no enough participants
        if (participants.size() < kMinAmountOfPartiesRequired) {
            LOG(INFO) << "Expect at least 2 parties, got " << participants.size();
            ids = std::move(participants);
            return ids;
        }

        // notify everyone to start throughput test
        broadcast(MessageType::THROUGHPUT_TEST_START, participants);
        
        // launch timer before start test
        auto end_time = std::chrono::steady_clock::now() + config.duration;

        DBG(">>> starting throughput test <<<");

        // All parties have joined, start test now
        while (std::chrono::steady_clock::now() < end_time) {
            // receive throughput test packages
            auto rc = zmq::poll(items, 100);
            if (rc == 0) continue;

            socket.recv(&identity);
            socket.recv(&request);

            // check where the package from
            memcpy(&sid, identity.data<char>(), identity.size());

            if (package_received.find(sid) != package_received.end()) {
                package_received[sid]++;
            }
            else {
                DBG("WARNING: receive message from unregisted client");
                continue;
            }
        }

        broadcast(MessageType::THROUGHPUT_TEST_STOP, participants);

        std::for_each(package_received.begin(), package_received.end(), [&](auto& count) { 
                count.second = count.second * config.message_size * 1024 / config.duration.count();
                DBG("throughput: " << count.second << "byte/s");
        });

        std::vector<SocketId> kickout;

        for (auto& KBps : package_received) {
            if (predicate(KBps.second)) {
                survivor.emplace_back(KBps.first);
            }
            else {
                LOG(INFO) << "party " << KBps.first << " was kicked out";
                kickout.emplace_back(KBps.first);
            }
        }

        // clear the receive buffer
        int rc = 1;
        while (rc > 0) {
            rc = zmq::poll(items, config.cleanup_time);
            if (items[0].revents & ZMQ_POLLIN) {
                socket.recv(&request);
            }
        }

        LOG(INFO) << "Parties kicked out: " << kickout.size() << ", parties remaining: " << survivor.size();
        broadcast(MessageType::THROUGHPUT_TEST_SURVIVE, survivor);
        broadcast(MessageType::THROUGHPUT_TEST_KICKOUT, kickout);
        numParties = survivor.size();

        return survivor;
    }

    expected<std::pair<SocketId, zmq::message_t*>> awaitNextInput() {
        return awaitNextInput(receive_timeout);
    }

    /** Waits for the next message over the socket.
     *
     *  Returns a deserialized SocketId and message pair from the
     *  sending party.
     */
    expected<std::pair<SocketId, zmq::message_t*>> awaitNextInput(std::chrono::milliseconds timeout) {
        DBG("In awaitNextInput");
        zmq::message_t identity;
        zmq::message_t *request = new zmq::message_t();

        SocketId sid;
        SocketIterator it;

        std::vector<zmq::pollitem_t> items = { 
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };

        try {
            do {
                // Block, waiting for the next incoming request. The ZeroMQ context
                // automatically handles reading from the TCP socket on a separate thread and
                // preparing/queueing incoming messages. Thus from the perspective
                // of our coordinator code we just wait synchronously until we have
                // the data necessary to proceed.
                DBG("Polling...");
                auto rc = zmq::poll(items, timeout);
                DBG("items[0].revents & ZMQ_POLLIN = " << static_cast<bool>(items[0].revents & ZMQ_POLLIN));
                if (rc && (items[0].revents & ZMQ_POLLIN)) {
                    socket.recv(&identity);
                } else {
                    LOG(INFO) << " Timedout: no data transfer to coordinator";
                    delete request;
                    return Error::TIMED_OUT;
                }

                if (!receiving) {
                    myTimers.end(4,"Idle (Waiting)",communicationCost);
                    myTimers.begin(4,"Actual transfer",communicationCost);
                    receiving = true;
                }

                // Receive remaining data
                socket.recv(request);

                // Communication Cost
                communicationCost += identity.size() + request->size();

                // Assign the socket id
                memcpy(&sid, identity.data<char>(), identity.size());

                it = std::find(ids.begin(), ids.end(), sid);

            } while (it == ids.end());
        } catch(...) {
            LOG(INFO) << "Caught exception";
            return Error::UNKNOWN_ERROR;
        }

        received[sid] = true;

        DBG("Finished awaitNextInput");
        return std::make_pair(sid,request);
    }

    template <typename T>
    static std::pair<std::vector<SocketId>, std::vector<SocketId>>
    processBatch1(
            MessageType t,
            std::function<T(T, T)> op,
            T& accumulator,
            std::vector< std::pair<SocketId,zmq::message_t*> > *accumulated,
            std::unordered_map<SocketId, T, boost::hash<SocketId>>& storage
            ) {
        bool initializing = true;

        DBG("In processBatch1");
        DBG("accumulated->size:" << accumulated->size());

        std::vector<SocketId> kickouts;
        std::vector<SocketId> restarts;

        size_t subIndex = 0;
        for (std::pair<SocketId,zmq::message_t*>& elt: (*accumulated)) {
            try {
                std::stringstream ss;
                MessageType header;
                T input;

                ss.write(elt.second->data<char>(), elt.second->size());

                boost::archive::binary_iarchive ia(ss);
                ia >> header;
                if (initializing) {
                    ia >> accumulator;
                    storage[elt.first] = accumulator;
                    subIndex++;
                } else {
                    ia >> input;
                    storage[elt.first] = input;
                    subIndex++;
                }

                if(t != header)
                {
                    DBG("Received shares from Participant (" << elt.first << ")");
                    DBG("t:" << (int)t);
                    DBG("t:" << msgs[1 + int(t) - int(MessageType::ID_PARTY)]);
                    DBG("header:" << msgs[1 + int(header) - int(MessageType::ID_PARTY)]);

                    kickouts.push_back(elt.first);
                    delete elt.second;
                    continue;
                }

                // Check parties commitments
                if (t == MessageType::GAMMA_RANDOM_SEED_SHARES 
                    || t == MessageType::PUBLIC_KEY_A_SHARES
                    || t == MessageType::GCD_AX_BY_SHARES) {

                    auto pcit = partiesCommitments.find(elt.first);
                    if (pcit != partiesCommitments.end()) {
                        std::string expectedHash;
                        if (t == MessageType::GAMMA_RANDOM_SEED_SHARES) {
                            expectedHash = pcit->second.seed1hash();
                        } else if (t == MessageType::PUBLIC_KEY_A_SHARES) {
                            expectedHash = pcit->second.aiHash();
                        } else if (t == MessageType::GCD_AX_BY_SHARES) {
                            expectedHash = pcit->second.seed2hash();
                        } 

                        std::string newHash;
                        if (initializing) {
                            newHash = math::computeHash(accumulator);
                        } else {
                            newHash = math::computeHash(input);
                        }

                        if (newHash != expectedHash) {
                            DBG("sid: " << elt.first);
                            DBG("Expected hash: " << expectedHash);
                            DBG("Got hash: " << newHash);
                            throw std::invalid_argument("hash mismatched");
                        }
                    } else {
                        LOG(ERROR) << "Something went wrong, party could not be found in parties commitments map";
                        throw std::invalid_argument("failed to find party in commitments");
                    }
                }

                if (initializing) {initializing = false;}
                else accumulator = op(accumulator, input);

                // Clean up the 0MQ message as it is no longer required
                delete elt.second;
                restarts.push_back(elt.first);
            } catch (...) {
                DBG("Participant (" << elt.first << "sent wrong data, adding to kickout list");
                kickouts.push_back(elt.first);
                delete elt.second;
            }
        }

        DBG("processBatch1 finished");

        delete accumulated;
        return {kickouts, restarts};
    }

    template<typename T>
    static void processBatch2( MessageType t, std::function<T(T, T)> op, T& accumulator, std::vector<T> *accumulated) {
        DBG("In procesBatch2");
        bool initializing = true;

        for (auto& elt: (*accumulated)) {
            if (initializing) {
                
                accumulator = elt;
                initializing = false;
                }
            else {
                accumulator = op(accumulator, elt);
            }
        }
        DBG("Finished procesBatch2");
        delete accumulated;
    }


    /** Waits for input from each of the given parties.
     *
     *  Aggregates the input by the given binary operand `op` and returns the
     *  result, as in a functional reduce or fold-left pattern.
     */
    template <typename T>
    expected<std::pair<T, std::unordered_map<SocketId, T, boost::hash<SocketId>>>> awaitAggregateVectorInput(
        MessageType t,
        const std::vector<SocketId>& ids,
        std::function<T(T, T)> op,
        T accumulator
        ) {

        // Record header of the phase into a file
        DBG("Writing to script file MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)]); 
        scriptOArchive << t;

        DBG("received.size() = " << received.size());

        boost::dynamic_bitset<> hasReceivedInput (numParties);
        receiving = false;
        myTimers.begin(3,"Receiving Data",communicationCost);
        myTimers.begin(4,"Idle (Waiting)",communicationCost);

        // Allocate tasks per thread
        //std::vector<std::thread> localThreads;
        std::vector<T> intermediateResultsPhaseIn;
        size_t entriesProcessed = 0;

        // Setting parameters
        this->nbThreads = std::min(NB_MAX_THREADS, static_cast<int>(ceil(static_cast<double>(numParties)/2.0)));
        this->entriesPerBatch = ceil(static_cast<double>(numParties)/this->nbThreads);
        std::vector<T> intermediateResultsPhaseOut(ceil(static_cast<double>(numParties)/this->entriesPerBatch));
        std::vector<std::unordered_map<SocketId, T, boost::hash<SocketId>>> storageVector(intermediateResultsPhaseOut.size());

        // Phase 1 Processing
        {
                    // clear all receive flags
            for (auto it = received.begin(); it != received.end(); it++) {
                it->second = false;
            }

            DBG("entriesPerBatch:" << this->entriesPerBatch);
            DBG("#threads:" << intermediateResultsPhaseOut.size());

            std::vector<std::future<std::pair<std::vector<SocketId>, std::vector<SocketId>>>> kickoutsAndRestarts;

            myTimers.begin(2,"3.a. Wall clock on time on deserializing " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost);

            for (size_t i = 0; i < intermediateResultsPhaseOut.size(); i++) {
                DBG("i = " << i);
                DBG(" intermediateResultsPhaseOut.size() = " << intermediateResultsPhaseOut.size());


                std::vector< std::pair<SocketId,zmq::message_t*> > *batch = new std::vector< std::pair<SocketId,zmq::message_t*> >;

                for (size_t j = 0; (j < this->entriesPerBatch) && (entriesProcessed < static_cast<unsigned int>(numParties));j++) {
                    DBG("Awaiting for input");

                    auto input = awaitNextInput();
                    if (!hasError(input)) {
                        DBG("Received input");
                        batch->push_back(getResult(input));
                        entriesProcessed++;
                    }
                    else {
                        DBG("Did not receive input");

                        std::set<SocketId> kickedout;
                        // Wait all threads end before return
                        if (i > 0) {
                            for ( auto& el : kickoutsAndRestarts ) {
                                auto [kick, success] = el.get();
                                if (!kick.empty()) {
                                    std::for_each(kick.begin(), kick.end(), [&](const SocketId& id) { kickedout.insert(id); });
                                }
                            }
                        }

                        // kick out timed out parties and broadcast restart
                        std::vector<SocketId> kickouts;
                        std::vector<SocketId> restarts;

                        for (auto it = received.begin(); it != received.end(); it++) {
                            if (it->second) {
                                if (!kickedout.empty() && (kickedout.find(it->first) != kickedout.end())) {
                                    continue;
                                }
                                LOG(INFO) << "Adding party to restart: " << it->first;
                                restarts.push_back(it->first);
                            } else {
                                LOG(INFO) << "Adding party to kickout: " << it->first;
                                kickouts.push_back(it->first);
                            }
                        }

                        LOG(INFO) << "kickouts.size() = " << kickouts.size();
                        LOG(INFO) << "restarts.size() = " << restarts.size();
                        DBG("Restarting protocol for " << restarts.size() << "parties");
                        // remove data and start from the beginning
                        restartCount++;
                        ofs.open(std::string("script") + std::to_string(restartCount) + std::string(".data"), std::ios::out | std::ios::binary);
                        DBG("Writing restart.size()");
                        scriptOArchive << (int)restarts.size();

                        // update remaining ids
                        auto success = update_ids(restarts);
                        if (hasError(success)) {
                            return getError(success);
                        }
                        else {
                            LOG(INFO) << "Update ids with " << getResult(success) << " parties";
                        }

                        delete batch;
                        return Error::RESTART;
                    }
                }

                DBG("Calling processBatch1...");

                auto f = std::async(&processBatch1<T>, t, op, std::ref(intermediateResultsPhaseOut[i]), batch, std::ref(storageVector[i]));
                kickoutsAndRestarts.push_back(std::move(f));
            }

            {
                DBG("Waiting for completion");

                std::vector<SocketId> restarts;
                bool restartProtocol = false;
                for ( auto& el : kickoutsAndRestarts ) {
                    auto ki = el.get();

                    DBG("Number of kickouts = " << ki.first.size());
                    if (ki.first.size() > 0) {
                        restartProtocol = true;
                    }

                    restarts.insert(restarts.end(), ki.second.begin(), ki.second.end());
                }

                if (restartProtocol) {
                    DBG("Restarting protocol for " << restarts.size() << " parties");
                    restartCount++;
                    ofs.open(std::string("script") + std::to_string(restartCount) + std::string(".data"), std::ios::out | std::ios::binary);
                    DBG("Writing restart.size()");
                    scriptOArchive << (int)restarts.size();

                    auto success = update_ids(restarts);
                    if (hasError(success)) {
                        return getError(success);
                    }
                    else {
                        LOG(INFO) << "Update ids with " << getResult(success) << " parties";
                    }
                    return Error::RESTART;
                }
            }
        }


        // Phase 2 Recursively Aggregating
        DBG("Phase 2");
        do {
            // Clean up
            //localThreads.clear();
            std::vector<std::thread> localThreads;
            intermediateResultsPhaseIn.swap(intermediateResultsPhaseOut);

            // Setting parameters for this round
            this->nbThreads = std::min(NB_MAX_THREADS, static_cast<int>(ceil(static_cast<double>(intermediateResultsPhaseIn.size())/2.0)));
            this->entriesPerBatch = ceil(static_cast<double>(intermediateResultsPhaseIn.size())/this->nbThreads);

            entriesProcessed = 0;
            intermediateResultsPhaseOut.clear();
            intermediateResultsPhaseOut.resize(ceil(static_cast<double>(intermediateResultsPhaseIn.size())/this->entriesPerBatch));

            DBG("entriesPerBatch:" << this->entriesPerBatch);
            DBG("#threads:" << intermediateResultsPhaseOut.size());

            for (size_t i = 0; i < intermediateResultsPhaseOut.size(); i++) {
                std::vector<T> *batch = new std::vector<T>;

                for (size_t j = 0; (j < this->entriesPerBatch) && (entriesProcessed < intermediateResultsPhaseIn.size());j++) {
                    batch->push_back(intermediateResultsPhaseIn[entriesProcessed]);
                    entriesProcessed++;
                }

                localThreads.emplace_back(processBatch2<T>, t, op, std::ref(intermediateResultsPhaseOut[i]), batch);
            }

            // Waiting for completion
            for( auto&& t : localThreads) {t.join();}


        } while (intermediateResultsPhaseOut.size() > MAX_SINGLE_THREADED_ITERATIONS);

        myTimers.end(2,"3.a. Wall clock on time on deserializing " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost, true);
        myTimers.end(4,"Actual transfer",communicationCost);

        // Phase 3 Aggregating                 
        {
            DBG("Phase 3");
            for (auto&& shares : intermediateResultsPhaseOut) accumulator = op(accumulator, shares);
            myTimers.end(3,"Receiving Data",communicationCost);
        }

        // combine all storage into one
        std::unordered_map<SocketId, T, boost::hash<SocketId>> storage;
        for (auto& storageVal : storageVector) {
            for (auto it = storageVal.begin(); it != storageVal.end(); ++it) {
                storage[it->first] = it->second;
            }
        }
        for (auto storageElement: storage) {
            DBG("Writing to scriptOA " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << "for party " << storageElement.first);
            scriptOArchive << std::pair{t, storageElement.second};
        }

        LOG(DEBUG) << "Writing to scriptOA " << msgs[1 + int(t) - int(MessageType::ID_PARTY)]; 
        scriptOArchive << std::pair{t, accumulator};
        return std::make_pair(accumulator, storage);
    }

    //==============================================================================
    /** Send an empty message of the given type to the given sockets. */
    void broadcast (MessageType t, const std::vector<SocketId>& sids) {
        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);

        oa << static_cast<int>(t);
        std::string msg = ss.str();

        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";

        for (const auto& id : sids) {
            bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
            if (!success) throw std::runtime_error("Failing to send message.");

            success = socket.send(msg.c_str(), msg.size());
            if (!success) throw std::runtime_error("Failing to send message.");

            // Communication Cost
            communicationCost += id.size() + msg.size();
        }
    }

    /** Send a set of values to each of the given sockets, prefaced with the given
     *  message type header.
     */
    template <typename T>
    void send (MessageType t, const SocketId& id, const T& values) {
        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);

        // Serilize column vector to string

        myTimers.begin(2,"3.a. Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost);
        oa << static_cast<int>(t);
        oa << values;
        myTimers.end(2,"3.a. Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost, true);
        std::string msg = ss.str();

        bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
        if (!success) throw std::runtime_error("Failing to send message.");

        success = socket.send(msg.c_str(), msg.size());
        if (!success) throw std::runtime_error("Failing to send message.");
        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";

    }

    /** Write data to the transcript file
     */
    template <typename T>
    void writeToTranscript (MessageType t, const T& values) {
        scriptOArchive << t;
        scriptOArchive << values;
    }

    /** Send a set of values to each of the given sockets, prefaced with the given
     *  message type header.
     */
    template <typename T>
    void broadcast (MessageType t, const std::vector<SocketId>& sids, const T& values) {
        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);

        myTimers.begin(3,"Broadcasting Data",communicationCost);

        // Serilize column vector to string
        myTimers.begin(2,"3.a Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost);
        oa << static_cast<int>(t);
        oa << values;
        myTimers.end(2,"3.a Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost, true);
        std::string msg = ss.str();

        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";

        for (const auto& id : sids) {
            bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
            if (!success) throw std::runtime_error("Failing to send message.");

            success = socket.send(msg.c_str(), msg.size());

            if (!success) throw std::runtime_error("Failing to send message.");

            // Communication Cost
            communicationCost += id.size() + msg.size();
        }

        myTimers.end(3,"Broadcasting Data",communicationCost);
    }

    std::string getCommunicationCost () {
        char buf[100];

        sprintf(buf,"Communication cost: %0.02f MB", (double)communicationCost/1e6);
        return std::string(buf);
    }

    expected<int> update_ids(const std::vector<SocketId>& newIds) {
        auto size = newIds.size();
        // no enough parties to participate
        DBG("Number of parties remaining: " << size);
        if (size < kMinAmountOfPartiesRequired) return Error::TOO_FEW_PARTIES;

        ids.clear();
        ids.reserve(size);

        numParties = size;

        ids = std::move(newIds);

        // clean received
        received.clear();
        for (auto it = ids.begin(); it != ids.end(); it++) {
            received[*it] = true;
        }
        return size;
    }

    expected<int> update_ids_from_received() {
        auto size = std::count_if(received.begin(), received.end(), [](const auto& p) { return p.second; });
        // no enough parties to participate
        DBG("Number of parties remaining: " << size);
        if (size < kMinAmountOfPartiesRequired) return Error::TOO_FEW_PARTIES;

        numParties = size;

        ids.clear();
        ids.reserve(size);
        for (auto it = received.begin(); it != received.end(); it++) {
            if (it->second) {
                ids.emplace_back(it->first);
            }
        }

        // clean received
        received.clear();
        for (auto it = ids.begin(); it != ids.end(); it++) {
            received[*it] = true;
        }
        return size;
    }

    int& parties() { return numParties; }
    const int& parties() const { return numParties; }

    //==============================================================================
    /** Run the protocol... */
    virtual void host() {};

    size_t communicationCost = 0;

protected:
    //==============================================================================
    zmq::context_t context;
    zmq::socket_t socket;
    int numParties;
    std::vector<SocketId> survivor;
    std::chrono::milliseconds receive_timeout;
    std::unordered_map<SocketId, bool, boost::hash<SocketId>> received;

    int restartCount;
    unsigned int nbThreads;
    unsigned int entriesPerBatch;
};

//==============================================================================
/** An interfacee for coordinating an MPC protocol over ZeroMQ. */
class ZeroMQGatheringTransport
{
public:
    timers myTimers;
    bool receiving;
    static std::vector<SocketId> ids;
    static std::unordered_map<SocketId, PartyCommitment, boost::hash<SocketId>> verifiersCommitments;
                
    //==============================================================================
    ZeroMQGatheringTransport(const std::string& uri, int verifiers)
        : context(kDefaultZeroMQIOThreads, 90),
          socket(context, ZMQ_ROUTER),
          numverifiers(verifiers)
    {
        socket.setsockopt(ZMQ_TCP_KEEPALIVE, 1);
        socket.setsockopt(ZMQ_TCP_KEEPALIVE_INTVL, 30);
        socket.setsockopt(ZMQ_LINGER, 0);

        socket.bind(uri);
    }

    //virtual ~ZeroMQCoordinatorTransport() = default;
    ~ZeroMQGatheringTransport()  {
    }

    /** Waits for the next message over the socket.
     *
     *  Returns a deserialized SocketId and message pair from the
     *  sending party.
     */
    expected<std::pair<SocketId, zmq::message_t*>> awaitNextInput(std::chrono::milliseconds timeout) {
        DBG("In awaitNextInput");
        zmq::message_t identity;
        zmq::message_t *request = new zmq::message_t();

        SocketId sid;
        SocketIterator it;

        std::vector<zmq::pollitem_t> items = { 
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };

        try {
            do {
                // Block, waiting for the next incoming request. The ZeroMQ context
                // automatically handles reading from the TCP socket on a separate thread and
                // preparing/queueing incoming messages. Thus from the perspective
                // of our coordinator code we just wait synchronously until we have
                // the data necessary to proceed.
                DBG("Polling...");
                auto rc = zmq::poll(items, timeout);
                DBG("items[0].revents & ZMQ_POLLIN = " << static_cast<bool>(items[0].revents & ZMQ_POLLIN));
                if (rc && (items[0].revents & ZMQ_POLLIN)) {
                    socket.recv(&identity);
                } else {
                    LOG(INFO) << " Timedout: no data transfer to coordinator";
                    delete request;
                    return Error::TIMED_OUT;
                }

                if (!receiving) {
                    myTimers.end(4,"Idle (Waiting)",communicationCost);
                    myTimers.begin(4,"Actual transfer",communicationCost);
                    receiving = true;
                }
                // Receive remaining data
                socket.recv(request);

                // Communication Cost
                communicationCost += identity.size() + request->size();

                // Assign the socket id
                memcpy(&sid, identity.data<char>(), identity.size());

                it = std::find(ids.begin(), ids.end(), sid);

                /*if (it == ids.end()) {
                  DBG("Msg from Participant (" << sid << ") Rejected");
                  } else {
                  DBG("Msg from Participant (" << sid << ") Accepted");
                  }*/
            } while (it == ids.end());
        } catch(...) {
            LOG(INFO) << "Caught exception";
            return Error::DESERIALIZE_FAIL;
        }

        DBG("Finished awaitNextInput");
        return std::make_pair(sid,request);
    }

    expected<std::pair<SocketId, std::string>> awaitShortRegistration(MessageType t, std::chrono::milliseconds register_timeout) {
        DBG("In awaitNextInput");
        SocketId sid;
        zmq::message_t identity;
        zmq::message_t *rec = new zmq::message_t();
        MessageType msgType;
        std::string verifierIPAddress = "UNKNOWN";

        std::vector<zmq::pollitem_t> items = { 
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };

        try {
            // Block, waiting for the next incoming request. The ZeroMQ context
            // automatically handles reading from the TCP socket on a separate thread and
            // preparing/queueing incoming messages. Thus from the perspective
            // of our coordinator code we just wait synchronously until we have
            // the data necessary to proceed.
            DBG("Polling...");
            auto rc = zmq::poll(items, register_timeout);
            DBG("items[0].revents & ZMQ_POLLIN = " << static_cast<bool>(items[0].revents & ZMQ_POLLIN));

            if (rc && (items[0].revents & ZMQ_POLLIN)) {
                socket.recv(&identity);

                // Assign the socket id
                memcpy(&sid, identity.data<char>(), identity.size());
                ids.push_back(sid);

                socket.recv(rec);

                // Check the msg
                std::stringstream ss;
                ss.write(rec->data<char>(), rec->size());
                delete rec;

                boost::archive::binary_iarchive ia(ss);
                ia >> msgType;
                ia >> verifierIPAddress;

                if (t != msgType) {
                    return Error::OUT_OF_SYNC;
                }

            } else {
                LOG(INFO) << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " Timedout: no data transfer to coordinator";
                delete rec;
                return Error::TIMED_OUT;
            }
        } catch(...) {
            LOG(INFO) << "Caught exception";
            delete rec;
            return Error::DESERIALIZE_FAIL;
        }

        DBG("Finished awaitShortRegistration");
        return std::make_pair(sid, verifierIPAddress);
    }

    //==============================================================================
    /** Send an empty message of the given type to the given sockets. */
    void broadcast (MessageType t, const std::vector<SocketId>& sids) {
        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);

        oa << static_cast<int>(t);
        std::string msg = ss.str();

        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";

        for (const auto& id : sids) {
            bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
            if (!success) throw std::runtime_error("Failing to send message.");

            success = socket.send(msg.c_str(), msg.size());
            if (!success) throw std::runtime_error("Failing to send message.");

            // Communication Cost
            communicationCost += id.size() + msg.size();
        }
    }

    /** Send a set of values to each of the given sockets, prefaced with the given
     *  message type header.
     */
    template <typename T>
    void send (MessageType t, const SocketId& id, const T& values) {
        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);

        // Serilize column vector to string
        myTimers.begin(2,"3.a. Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost);
        oa << static_cast<int>(t);
        oa << values;
        myTimers.end(2,"3.a. Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost, true);
        std::string msg = ss.str();

        bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
        if (!success) throw std::runtime_error("Failing to send message.");

        success = socket.send(msg.c_str(), msg.size());
        if (!success) throw std::runtime_error("Failing to send message.");
        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";

    }

    /** Send a set of values to each of the given sockets, prefaced with the given
     *  message type header.
     */
    void sendMsg (const SocketId& id, const std::string& msg) {
        LOG(INFO) << "party " << id << "msg.size()" << msg.size(); 

        bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
        if (!success) throw std::runtime_error("Failing to send message.");

        success = socket.send(msg.c_str(), msg.size());
        if (!success) throw std::runtime_error("Failing to send message.");
        LOG(INFO) << "MessageType: " << msgs[1 + int(MessageType::GATHER_PUBLIC_DATA) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";
    }

    /** Send a set of values to each of the given sockets, prefaced with the given
     *  message type header.
     */
    void sendRaw (const SocketId& id, zmq::message_t* values) {

        bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
        if (!success) throw std::runtime_error("Failing to send message.");

        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);

        success = socket.send(values->data(), values->size());
        if (!success) throw std::runtime_error("Failing to send message.");
        LOG(INFO) << "MessageType: " << msgs[1 + int(MessageType::GATHER_PROOFS) - int(MessageType::ID_PARTY)] << " message size: " << values->size() << " bytes";

    }

    /** Send a set of values to each of the given sockets, prefaced with the given
     *  message type header.
     */
    template <typename T>
    void broadcast (MessageType t, const std::vector<SocketId>& sids, const T& values) {
        std::stringstream ss;
        boost::archive::binary_oarchive oa(ss);

        myTimers.begin(3,"Broadcasting Data",communicationCost);

        // Serilize column vector to string
        myTimers.begin(3,"3.a. Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost);
        oa << static_cast<int>(t);
        oa << values;
        myTimers.end(3,"3.a. Serializing Data " + msgs[1 + int(t) - int(MessageType::ID_PARTY)], communicationCost, true);
        std::string msg = ss.str();

        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";
        for (const auto& id : sids) {
            bool success = socket.send(id.begin(), id.size(), ZMQ_SNDMORE);
            if (!success) throw std::runtime_error("Failing to send message.");

            success = socket.send(msg.c_str(), msg.size());
            if (!success) throw std::runtime_error("Failing to send message.");

            // Communication Cost
            communicationCost += id.size() + msg.size();
        }

        myTimers.end(3,"Broadcasting Data",communicationCost);
    }

    int& verifiers() { return numverifiers; }
    const int& verifiers() const { return numverifiers; }

    //==============================================================================
    /** Run the protocol... */
    virtual void host() {};

    size_t communicationCost = 0;

protected:
    //==============================================================================
    zmq::context_t context;
    zmq::socket_t socket;
    int numverifiers;
    std::vector<SocketId> survivor;
};

std::unordered_map<SocketId, PartyCommitment, boost::hash<SocketId>> ZeroMQCoordinatorTransport::partiesCommitments;
boost::filesystem::ofstream ZeroMQCoordinatorTransport::ofs("script.data", std::ios::out | std::ios::binary);
boost::archive::binary_oarchive ZeroMQCoordinatorTransport::scriptOArchive(ZeroMQCoordinatorTransport::ofs, boost::archive::no_header);


//==============================================================================
/** An interface for joining an MPC protocol over ZeroMQ. */
class ZeroMQClientTransport
{
public:
    timers myTimers;

    //==============================================================================
    // ZeroMQClientTransport(const std::string& uri)
    //     : context(kDefaultZeroMQIOThreads),
    //       socket(context, ZMQ_DEALER),
    //       socketId(boost::uuidFs::random_generator()())
    // {
    //     socket.setsockopt(ZMQ_IDENTITY, socketId);
    //     socket.connect(uri);
    // }

    ZeroMQClientTransport(const std::string& uri, const SocketId& _socketId, std::chrono::milliseconds receive_timeout)
        : context(kDefaultZeroMQIOThreads),
          socket(context, ZMQ_DEALER),
          socketId(_socketId),
          receive_timeout(receive_timeout)
    {
        socket.setsockopt(ZMQ_IDENTITY, socketId);
        socket.setsockopt(ZMQ_TCP_KEEPALIVE, 1);
        socket.setsockopt(ZMQ_TCP_KEEPALIVE_INTVL, 30);
        //socket.setsockopt(ZMQ_LINGER, 0);

        socket.connect(uri);
    }

    virtual ~ZeroMQClientTransport() = default;

    // joining throughput test. return a boolean indicating if this client can move on to registration
    // @wait_timeout: duration to wait when registration
    // @hard_timeout: duration of throughput test
    // @msg_size size of bytes of package to be sent
    expected<Unit> joinThroughputTest(size_t wait_timeout, size_t hard_timeout, size_t msg_size, size_t nb_max_send) {

        bool received_from_coordinator = false;
        std::vector<zmq::pollitem_t> items = { 
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };

        send(MessageType::THROUGHPUT_TEST_START);
        LOG(INFO) << "waiting to join throughput test...";

        auto end_time = std::chrono::steady_clock::now() + std::chrono::milliseconds(wait_timeout);
        while (std::chrono::steady_clock::now() < end_time) {
            zmq::poll(items, 100);
            if (items[0].revents & ZMQ_POLLIN) {
                auto header = awaitReply();

                if (hasError(header)) {
                    return getError(header);
                }

                received_from_coordinator = true;
                if (getResult(header) == MessageType::THROUGHPUT_TEST_START) {
                    break;
                }
                else {
                    LOG(ERROR) << "Cannot join throughput test: expect THROUGHPUT_TEST_START, got header index=" << static_cast<int>(getResult(header));
                    if (getResult(header) == MessageType::ABORT_PROTOCOL_VERSION_MISMATCH) {
                        return Error::KILLED_BY_COORDINATOR;
                    }
                    else if (getResult(header) == MessageType::PROTOCOL_KILLED) {
                        return Error::KILLED_BY_COORDINATOR;
                    }
                    else {
                        return Error::OUT_OF_SYNC;
                    }
                }
            }
        }

        if (!received_from_coordinator) {
            LOG(INFO) << "Not receiving message from coordinator after " << wait_timeout << " ms, terminating.";
            return Error::TIMED_OUT;
        }

        LOG(INFO) << "Joining throughput test";

        std::string message(msg_size, '9');

        end_time = std::chrono::steady_clock::now() + std::chrono::milliseconds(hard_timeout);
        received_from_coordinator = false;
        int message_sent = 0;
        while (std::chrono::steady_clock::now() < end_time) {
            auto rc = zmq::poll(items, 0);

            if (items[0].revents & ZMQ_POLLIN) {
                auto header = awaitReply();
                received_from_coordinator = true;
                if (getResult(header) == MessageType::THROUGHPUT_TEST_STOP) {
                    LOG(INFO) << "End message received, stop test";
                    break;
                }
                else if (getResult(header) == MessageType::PROTOCOL_KILLED) {
                    return Error::KILLED_BY_COORDINATOR;
                }
                else {
                    return Error::OUT_OF_SYNC;
                }
            }

            // don't send if already sent enought packages
            if (message_sent > nb_max_send) {
                std::this_thread::sleep_for(std::chrono::milliseconds(100));
                continue;
            }

            socket.send(message.c_str(), message.size(), int(zmq::send_flags::dontwait));
            message_sent++;
        }

        if (!received_from_coordinator) {
            LOG(INFO) << "Not receiving stop message from coordinator, terminating.";
            return Error::TIMED_OUT;
        }

        auto header = awaitReply();

        if (hasError(header)) { return getError(header); }

        if (getResult(header) == MessageType::THROUGHPUT_TEST_SURVIVE) {
            return Unit{};
        }

        LOG(ERROR) << "Client expect header SURVIVE or KICKOUT, got message index=" << (int)(getResult(header)) - (int)MessageType::ID_PARTY;
        return Error::OUT_OF_SYNC;
    }

    //==============================================================================
    /** Waits for an empty message over the socket and returns its type. */
    expected<MessageType> awaitReply () {
        zmq::message_t reply;

        std::vector<zmq::pollitem_t> items = {
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };
        auto rc = zmq::poll(items, receive_timeout);
        if (rc && (items[0].revents & ZMQ_POLLIN)) {
            socket.recv(&reply);
        }
        else {
            return Error::TIMED_OUT;
        }

        MessageType header;

        std::stringstream ss;
        ss.write(reply.data<char>(), reply.size());

        boost::archive::binary_iarchive ia(ss);

        try {
            ia >> header;
        } catch (...) {
            return Error::DESERIALIZE_FAIL;
        }

        if (header == MessageType::PROTOCOL_RESTART) {
            return Error::RESTART;
        }
        if (header == MessageType::PROTOCOL_KILLED) {
            return Error::KILLED_BY_COORDINATOR;
        }
        return header;
    }

    /** Waits for and deserializes a message of the given type over the socket. */
    template <typename T>
    expected<T> awaitReply(MessageType type, bool timers = false) {
        zmq::message_t reply;

        // Logging Performance data
        //TIMED_FUNC(timerObj);
        std::vector<zmq::pollitem_t> items = {
            {static_cast<void*>(socket), 0, ZMQ_POLLIN, 0}
        };

        auto rc = zmq::poll(items, receive_timeout);
        if (rc && (items[0].revents & ZMQ_POLLIN)) {
            socket.recv(&reply);
        }
        else {
            return Error::TIMED_OUT;
        }


        if (timers) LOG(INFO) << "Receiving Data " << msgs[1+(int)type-int(MessageType::ID_PARTY)];

        // Deserialize the message content to a column vector
        MessageType header;
        T body;

        std::stringstream ss;
        ss.write(reply.data<char>(), reply.size());

        boost::archive::binary_iarchive ia(ss);

        try {
            ia >> header;
        } catch (...) {
            return Error::DESERIALIZE_FAIL;
        }

        // check if kill message received
        if (header == MessageType::PROTOCOL_RESTART) {
            return Error::RESTART;
        }
        if (header == MessageType::PROTOCOL_KILLED) {
            return Error::KILLED_BY_COORDINATOR;
        }

        try  {
            ia >> body;
        } catch (...) {
            return Error::OUT_OF_SYNC;
        }

        if (timers) LOG(INFO) << "Deserialization Completed for " << msgs[1+(int)type-int(MessageType::ID_PARTY)];

        if (type != header) {
            LOG(INFO) << "Expected " << msgs[1+(int)type-int(MessageType::ID_PARTY)] << ", got " << msgs[1+(int)header-int(MessageType::ID_PARTY)];
            return Error::OUT_OF_SYNC;
        }
        return body;
    }

    // Sends an empty message of the given type over the socket. */
    void send (MessageType t) {
        std::stringstream ss;

        // Serialize the message header
        boost::archive::binary_oarchive oa(ss);
        oa << t;

        // Push it through the socket
        std::string msg = ss.str();
        bool success = socket.send(msg.c_str(), msg.size());

        if (!success) throw std::runtime_error("Failing to send message.");
        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";
    }

    void send (MessageType t, const PublicData& message) {
        std::stringstream ss;

        // Serilize column vector to string
        boost::archive::binary_oarchive oa(ss);
        oa << t;
        oa << message;

        // Push it through the socket
        std::string msg = ss.str();
        LOG(INFO) << "Pdata size = " << msg.size();
        bool success = socket.send(msg.c_str(), msg.size());
        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";
        if (!success) throw std::runtime_error("Failing to send message.");
    }


    // Sends a message of the given type over the socket. */
    template <typename T>
    void send (MessageType t, const T& message) {
        std::stringstream ss;

        // Serilize column vector to string
        boost::archive::binary_oarchive oa(ss);
        oa << t;
        oa << message;

        // Push it through the socket
        std::string msg = ss.str();
        bool success = socket.send(msg.c_str(), msg.size());
        LOG(INFO) << "MessageType: " << msgs[1 + int(t) - int(MessageType::ID_PARTY)] << " message size: " << msg.size() << " bytes";
        if (!success) throw std::runtime_error("Failing to send message.");
    }

    //==============================================================================
    /** Run the protocol... */
    virtual void start() {}

    SocketId getSocketId() { return socketId; }


protected:
    //==============================================================================
    zmq::context_t context;
    zmq::socket_t socket;
    boost::uuids::uuid socketId;
    std::chrono::milliseconds receive_timeout;
};

std::vector<SocketId> ZeroMQCoordinatorTransport::ids;
std::vector<SocketId> ZeroMQGatheringTransport::ids;
} // namespace ligero
