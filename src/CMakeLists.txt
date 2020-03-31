project(src CXX)

# supress DEBUG for this folder, use INFO instead
add_compile_options(-DELPP_DISABLE_DEBUG_LOGS)

macro(add_binary TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    add_dependencies(${TESTNAME} nfl)
    target_sources(${TESTNAME} PRIVATE ${EASYLOGGINGPP_CC_FILE})
    target_include_directories(${TESTNAME} PRIVATE ../include ${CPPZMQ} ${DEPS_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${EASYLOGGINGPP_SRC_DIR} ${NFLLIB_INCLUDE_DIR})
    target_link_libraries(${TESTNAME} pthread gmp gmpxx sodium zmq ${Boost_LIBRARIES} mpfr ${NFLLIB_LIB_DIR}/libnfllib_static.a)
endmacro()

# Profiling of the protocol on AWS
#add_binary(benchmark_keygen benchmark_keygen.cpp)
add_binary(party_full_protocol party_full_protocol.cpp)
add_binary(coordinator_full_protocol coordinator_full_protocol.cpp)
add_binary(validator validator.cpp)
add_binary(distributed_verifier distributed_verifier.cpp)
add_binary(manual_verifier manual_verifier.cpp)
