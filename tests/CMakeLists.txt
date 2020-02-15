project(tests)

include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.8.0
)

FetchContent_GetProperties(googletest)
if(NOT googletest_POPULATED)
  FetchContent_Populate(googletest)
  add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})
endif()

macro(add_test_complete TESTNAME)
    add_executable(${TESTNAME} ${ARGN})
    add_dependencies(${TESTNAME} nfl)
    target_sources(${TESTNAME} PRIVATE ${EASYLOGGINGPP_CC_FILE})
    target_include_directories(${TESTNAME} PRIVATE ../include ${CPPZMQ} ${DEPS_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${EASYLOGGINGPP_SRC_DIR} ${NFLLIB_INCLUDE_DIR})
    target_link_libraries(${TESTNAME} gmp gmpxx zmq gtest gmock gtest_main sodium ${Boost_LIBRARIES} mpfr ${NFLLIB_LIB_DIR}/libnfllib_static.a)
    add_test(NAME ${TESTNAME} COMMAND ${TESTNAME}  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endmacro()

# Test Encrypt/Decrypt cycle 100 times
#add_test_complete(test_ringlwe_01_mocked encryption/test_ringlwe_01_mocked.cpp)

# Test Ligero Math Functions
#add_test_complete(test_math_01_bucketing math/test_math_01_bucketing.cpp)
add_test_complete(test_math_02_basic math/test_math_02_basic.cpp)

# Test Transport
#add_test_complete(test_transport_lwe_01_basic transport/test_transport_lwe_01_basic.cpp)
add_test_complete(test_transport_02_keygen_exception_handling transport/test_transport_02_keygen_exception_handling.cpp)
add_test_complete(test_transport_03_presieve_exception_handling transport/test_transport_03_presieve_exception_handling.cpp)
add_test_complete(test_transport_04_wrong_ai transport/test_transport_04_wrong_ai.cpp)
add_test_complete(test_transport_05_wrong_jacobi_seed transport/test_transport_05_wrong_jacobi_seed.cpp)
add_test_complete(test_transport_06_wrong_gcd_and_jacobi_seed transport/test_transport_06_wrong_gcd_and_jacobi_seed.cpp)
add_test_complete(test_transport_07_party_aborts_on_registration transport/test_transport_07_party_aborts_on_registration.cpp)

# Test the protocol over a local network
#add_test_complete(test_ringlwe_02_vecprod_mocked protocol/test_ringlwe_02_vecprod_mocked.cpp)
#add_test_complete(test_ringlwe_03_vecprod protocol/test_ringlwe_03_vecprod.cpp)
add_test_complete(test_ringlwe_04_presieve protocol/test_ringlwe_04_presieve.cpp)
add_test_complete(test_ringlwe_05_candidates protocol/test_ringlwe_05_candidates.cpp)
add_test_complete(test_ringlwe_06_jacobi protocol/test_ringlwe_06_jacobi.cpp)
add_test_complete(test_ringlwe_07_gcd protocol/test_ringlwe_07_gcd.cpp)

# Test functionality tied to the zero knowledge succint argument of knowledge (zkSNARK)
#add_test_complete(test_fft_01 zk/test_fft_01.cpp)
add_test_complete(test_mt_01 zk/test_mt_01.cpp)
add_test_complete(test_secret_sharing_01 zk/test_secret_sharing_01.cpp)
#add_test_complete(test_secret_sharing_02 zk/test_secret_sharing_02.cpp)
#add_test_complete(test_zk_01 zk/test_zk_01.cpp)
add_test_complete(test_zk_02 zk/test_zk_02.cpp)
add_test_complete(test_zk_03 zk/test_zk_03.cpp)
add_test_complete(test_zk_04 zk/test_zk_04.cpp)
add_test_complete(test_zk_05 zk/test_zk_05.cpp)
add_test_complete(test_zk_06 zk/test_zk_06.cpp)
