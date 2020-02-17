# RSACeremony
Distributed RSA Modulus Generation

## Prerequisites
1. If you are on Mac OS X or on Windows please install Docker
2. If you are on Ubuntu 14.04 install dependencies and setup environment:
```bash
$ sudo apt-get install build-essential pkg-config libgmp-dev wget git libmpfr-dev libsodium-dev gcc-8 g++-8 libzmq3-dev
$ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-7 700 --slave /usr/bin/g++ g++ /usr/bin/g++-7
$ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-8 800 --slave /usr/bin/g++ g++ /usr/bin/g++-8
$ wget --max-redirect 3 https://dl.bintray.com/boostorg/release/1.69.0/source/boost_1_69_0.tar.gz
$ sudo mkdir -p /usr/include/boost && tar zxf boost_1_69_0.tar.gz
$ sudo cd boost_1_69_0 && ./bootstrap.sh --prefix=/usr/local && ./b2 --with=all install && echo "/usr/local/lib" >> /etc/ld.so.conf.d/local.conf && ldconfig
$ sudo bash -c 'echo "deb https://download.opensuse.org/repositories/network:/messaging:/zeromq:/release-stable/xUbuntu_18.04 ./" >> /etc/apt/sources.list'
$ sudo wget https://download.opensuse.org/repositories/network:/messaging:/zeromq:/release-stable/xUbuntu_18.04/Release.key -qO- | apt-key add
$ sudo wget -qO- https://github.com/zeromq/cppzmq/archive/v4.3.0.tar.gz | tar xvzf - -C /usr/local/include
$ wget --max-redirect 3 https://github.com/Kitware/CMake/releases/download/v3.14.3/cmake-3.14.3-Linux-x86_64.tar.gz
$ tar -xzf cmake-3.14.3-Linux-x86_64.tar.gz
$ export PATH="/cmake-3.14.3-Linux-x86_64/bin:${PATH}"
```
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

Lets run a computation for 10 parties. First start coordinator.

For Ubuntu 14.04 users:
```bash
$ ./build/src/coordinator_full_protocol --parties 10
```

For the rest with Docker:
```bash
$ ./scripts/run-coordinator-in-docker.sh --parties 10
```

Now we can start connect parties to this coordinator. We will need to 
provide coordinators IP address or by default our party will connect to 127.0.0.1:5555

Again for Ubuntu 14.04 users:
```bash
$ ./build/src/party_full_protocol --ip 127.0.0.1
```

For the rest with Docker:
```bash
$ ./scripts/run-party-in-docker.sh
```

Once all ten parties connected or timeout for registration passed the protocol will start.


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
