project(benchmarks)

macro(add_benchmark_complete BENCHNAME)
    add_executable(${BENCHNAME} ${ARGN})
    add_dependencies(${BENCHNAME} nfl)
    target_sources(${BENCHNAME} PRIVATE ${EASYLOGGINGPP_CC_FILE})
    target_include_directories(${BENCHNAME} PRIVATE ../include ${CPPZMQ} ${DEPS_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${EASYLOGGINGPP_SRC_DIR} ${NFLLIB_INCLUDE_DIR})
    target_link_libraries(${BENCHNAME} gmp zmq pthread sodium ${Boost_LIBRARIES} ${NFLLIB_LIB_DIR}/libnfllib_static.a)
    add_test(NAME ${BENCHNAME} COMMAND ${BENCHNAME}  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro()

# Benchmark serialization on AWS
#add_benchmark_complete(benchmark_party_01_serialization benchmark_party_01_serialization.cpp)
#add_benchmark_complete(benchmark_coordinator_01_serialization benchmark_coordinator_01_serialization.cpp)
add_benchmark_complete(benchmark_fft_ntt_01 benchmark_fft_ntt_01.cpp)
add_benchmark_complete(benchmark_fft_ntt_02 benchmark_fft_ntt_02.cpp)
