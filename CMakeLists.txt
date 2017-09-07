# Version number of the application
set (CTHYB_VERSION "1.5")

# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake ${TRIQS_PATH}/share/triqs/cmake)

# start configuration
cmake_minimum_required(VERSION 2.8)
project(cthyb CXX)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# Check that versions are compatible
if(NOT CTHYB_VERSION EQUAL TRIQS_VERSION)
  message(FATAL_ERROR "The application version is not compatible with the TRIQS library (TRIQS library version: ${TRIQS_VERSION} while this application version: ${CTHYB_VERSION})")
endif()

# Get hash
triqs_get_git_hash(${CMAKE_SOURCE_DIR} "CTHYB")
if(${GIT_RESULT} EQUAL 0)
  message(STATUS "Hash: ${CTHYB_GIT_HASH}")
endif(${GIT_RESULT} EQUAL 0)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Add the compiling options (-D... ) for C++
message(STATUS "TRIQS : Adding compilation flags detected by the library (C++11/14, libc++, etc...) ")
add_definitions(${TRIQS_CXX_DEFINITIONS})

# Options for the compilation ...
option(Tests "Enable Tests" ON)
option(HYBRIDISATION_IS_COMPLEX "If ON, the hybridization Delta(tau) is complex" OFF)
option(LOCAL_HAMILTONIAN_IS_COMPLEX "If ON, the H_loc is complex" OFF)
option(EXT_DEBUG "Enable extended debugging output [developers only]" OFF)
option(SAVE_CONFIGS "Save visited configurations to configs.h5 [developers only]" OFF)
option(LOOK_FOR_NFFT "Look for the NFFT library for the G2 frequency sampling" ON)

if(HYBRIDISATION_IS_COMPLEX AND NOT LOCAL_HAMILTONIAN_IS_COMPLEX)
 message(FATAL_ERROR "Combination of options HYBRIDISATION_IS_COMPLEX=ON "
                     "and LOCAL_HAMILTONIAN_IS_COMPLEX=OFF is not supported.")
endif(HYBRIDISATION_IS_COMPLEX AND NOT LOCAL_HAMILTONIAN_IS_COMPLEX)

# 2-particle GF measurement requires NFFT
if(LOOK_FOR_NFFT)
  find_package(NFFT)
else(LOOK_FOR_NFFT)
  message(STATUS "Not looking for NFFT")
endif(LOOK_FOR_NFFT)

if(NFFT_FOUND)
  add_definitions(-DCTHYB_G2_NFFT)
  include_directories(${NFFT_INCLUDE_DIR})
  link_libraries(${NFFT_LIBRARIES})
  if(NFFT_OLD_API)
    add_definitions(-DNFFT_OLD_API)
  endif(NFFT_OLD_API)
endif(NFFT_FOUND)

if(EXT_DEBUG)
 add_definitions(-DEXT_DEBUG)
endif(EXT_DEBUG)

if(SAVE_CONFIGS)
 set(NUM_CONFIGS_TO_SAVE 50000 CACHE STRING "Number of visited configurations to save [developers only]")
 add_definitions(-DSAVE_CONFIGS -DNUM_CONFIGS_TO_SAVE=${NUM_CONFIGS_TO_SAVE})
endif(SAVE_CONFIGS)

include_directories(${CMAKE_BINARY_DIR}/c++ c++)

# Compile C++ code
add_subdirectory(c++)

# Tests
 if (${Tests})
   add_subdirectory(test)
 endif()

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
 option(BUILD_DOC "Build documentation" OFF)
 if(${BUILD_DOC})
   if(NOT ${TRIQS_WITH_DOCUMENTATION})
     message(FATAL_ERROR "Error: TRIQS library has not been compiled with its documentation")
   endif()
   add_subdirectory(doc)
 endif(${BUILD_DOC})
endif()
