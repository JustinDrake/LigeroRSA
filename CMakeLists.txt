cmake_minimum_required (VERSION 3.11)
project (Encryption)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(GMP REQUIRED)
find_package(Threads REQUIRED)

find_package(PkgConfig)
pkg_check_modules(PC_ZeroMQ QUIET zmq)
find_path(ZMQ_INCLUDE_DIR
        NAMES zmq.hpp
        PATHS ${PC_ZeroMQ_INCLUDE_DIRS}
        )
find_library(ZMQ_LIBRARY
        NAMES zmq
        PATHS ${PC_ZeroMQ_LIBRARY_DIRS}
        )

set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
set(BOOST_ROOT "/usr/include/boost/")
find_package( Boost COMPONENTS program_options serialization random REQUIRED )

# grab git hash and pass it as a protocol version
execute_process(
  COMMAND bash version.sh
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)

include_directories( ${Boost_INCLUDE_DIR} )
add_definitions(-DELPP_FEATURE_PERFORMANCE_TRACKING  -DELPP_THREAD_SAFE -DELPP_FORCE_USE_STD_THREAD -DPROTOCOL_VERSION_NUMBER=0x${GIT_COMMIT_HASH} -O3 -w -std=c++1z)

enable_testing()

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
    add_definitions(-DELPP_NO_LOG_TO_FILE)
endif()

set(CPPZMQ "${PROJECT_SOURCE_DIR}/depends/cppzmq")
set(DEPS_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/depends/exrandom/include")
set(EASYLOGGINGPP_SRC_DIR "${PROJECT_SOURCE_DIR}/depends/easyloggingpp/src")
set(EASYLOGGINGPP_CC_FILE "${EASYLOGGINGPP_SRC_DIR}/easylogging++.cc")
set(NFLLIB_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/depends/NFLlib/include")
set(NFLLIB_LIB_DIR "${PROJECT_SOURCE_DIR}/lib/lib")

include(ExternalProject)
ExternalProject_Add(nfl 
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/depends/NFLlib
        CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DNFL_OPTIMIZED=ON -DCMAKE_INSTALL_PREFIX=${PROJECT_SOURCE_DIR}/lib)

include_directories(${PROJECT_SOURCE_DIR}/depends/NFLlib/include)
include_directories(${PROJECT_SOURCE_DIR}/depends/NFLlib/lib)

add_subdirectory (src)
