#!/bin/bash

rm -rf out.log

 ./benchmark_coordinator_01_serialization "127.0.0.1" 4 out.log &
 ./benchmark_party_01_serialization "127.0.0.1" > log_party_1 &
 ./benchmark_party_01_serialization "127.0.0.1" > log_party_2 &
 ./benchmark_party_01_serialization "127.0.0.1" > log_party_3 &
 ./benchmark_party_01_serialization "127.0.0.1" > log_party_4 &
