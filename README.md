# RSACeremony
Distributed RSA Modulus Generation

## Prerequisites
1. If you are on Mac OS X or on Windows please install Docker
2. If you are on Ubuntu 14.04 install dependencies and setup environment:
```bash
sudo apt-get install build-essential pkg-config libgmp-dev wget git libmpfr-dev libsodium-dev gcc-8 g++-8 libzmq3-dev
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 700 --slave /usr/bin/g++ g++ /usr/bin/g++-7
sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 800 --slave /usr/bin/g++ g++ /usr/bin/g++-8
wget --max-redirect 3 https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz
sudo mkdir -p /usr/include/boost && tar zxf boost_1_69_0.tar.gz
sudo bash -c 'cd boost_1_69_0 && ./bootstrap.sh --prefix=/usr/local && ./b2 --with=all install && echo "/usr/local/lib" >> /etc/ld.so.conf.d/local.conf && ldconfig'
sudo bash -c 'wget https://download.opensuse.org/repositories/network:/messaging:/zeromq:/release-stable/xUbuntu_18.04/Release.key -qO- | apt-key add'
sudo bash -c 'echo "deb https://download.opensuse.org/repositories/network:/messaging:/zeromq:/release-stable/xUbuntu_18.04 ./" >> /etc/apt/sources.list'
sudo bash -c 'wget -qO- https://github.com/zeromq/cppzmq/archive/v4.3.0.tar.gz | tar xvzf - -C /usr/local/include'
wget --max-redirect 3 https://github.com/Kitware/CMake/releases/download/v3.14.3/cmake-3.14.3-Linux-x86_64.tar.gz
tar -xzf cmake-3.14.3-Linux-x86_64.tar.gz
export PATH="/cmake-3.14.3-Linux-x86_64/bin:${PATH}"
```
3. Coordinator Node Settings
   1. The Protocol currently uses port 5555 (passive/active protocol) and 5556 (proof verification), both ports must be open for TCP traffic on the coordinator's instance.
   2. For a 1024-party RSA Ceremony, setting stack size to 100 MB (ulimit -s 100000) and number of open files to 10K (ulimit -n 10000) has proven sufficient in testing.
   3. Minimum hardware for 1024-party node: Memory 600 GB, HDD 500 GB.

## Fetch Dependencies
```bash
$ git submodule update --init --recursive
```

## Building

#### Build and Run Tests with CMake 

For Ubuntu 14.04 users:
```bash
$ cd aws_scripts && ./native_build.sh
```

For Mac OS X, Windows, or Linux users with Docker:
```bash
$ ./scripts/build-binary-docker.sh
```

## Running

### Ubuntu 14.04
Lets run an active computation for 10 parties: First start the coordinator.

```bash
$ ./build/src/coordinator_full_protocol --parties 10
```

Now we can start to connect parties and verifiers to this coordinator. We will need to 
provide the coordinator's IP address or by default our party will connect to 127.0.0.1:5555

```bash
$ ./build/src/party_full_protocol --ip 127.0.0.1
```

```bash
$ ./build/src/distributed_verifier --ip 127.0.0.1
```

### Passive Protocol with Docker
For the 10-party passive protocol in Docker, run the coordinator and parties as follow:

```bash
$ ./scripts/run-coordinator-in-docker.sh --parties 10
```

```bash
$ ./scripts/run-party-in-docker.sh
```

Once all ten parties connected or timeout for registration passed the protocol will start.

## Additional Scripts
Instructions here are for native Ubuntu runs.

### Validator
The coordinator exports a record of its public interactions with the parties in the _script.data_ file generated the directory where _coordinator_full_protocol_ is run. To verify the integrity of these interactions, run _./validator_ in the same directory as _script.data_.

### Replay mode
To record and replay an RSA ceremony with 10 parties using the same randomness: 

1. Run a ceremony adding command line option _--mode record_ and _--passive_:

```bash
$ ./coordinator_full_protocol --parties 10 --passive --mode record
```

2. Replay it using option _--mode replay_:

```bash
$./coordinator_full_protocol --parties 10 --passive --mode replay
```

Comments:
1. All party binaries must be run in both cases.
2. Intermediary data is stored in file _record.data_, and binaries must be run from the same directory when replaying.



## License 

Copyright 2020 Ligero, Inc.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

   http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.



This code was developed under a project supported by the VDF Alliance, 
and particularly the Ethereum Foundation and Protocol Labs, Inc.
