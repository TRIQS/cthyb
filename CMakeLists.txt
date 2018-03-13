# start configuration
cmake_minimum_required(VERSION 2.8.12)
project(cthyb CXX)

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} -------------")

# Options for the compilation ...
option(Hybridisation_is_complex "If ON, the hybridization Delta(tau) is complex" OFF)
option(Local_hamiltonian_is_complex "If ON, the H_loc is complex" OFF)
option(MeasureG2 "Measure the two particle object (requires the NFFT library)" ON)

# check that options are compatible
if(Hybridisation_is_complex AND NOT Local_hamiltonian_is_complex)
 message(FATAL_ERROR "Combination of options Hybridisation_is_complex=ON "
                     "and Local_hamiltonian_is_complex=OFF is not supported.")
endif()

# Use shared libraries
set(BUILD_SHARED_LIBS ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS 2.0 REQUIRED)
find_package(Cpp2Py 1.5 REQUIRED)

# We will need the FindNNFT of triqs, so we adjust CMake Module Path to find it
list(APPEND CMAKE_MODULE_PATH ${TRIQS_ROOT}/share/cmake)

# 2-particle GF measurement requires NFFT
if(MeasureG2)
 find_package(NFFT REQUIRED)
endif()

# Default Install directory to TRIQS_ROOT if not given. Checks an absolute name is given.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS " No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} -------------")

# Use git to detect Cthyb version and git hash
execute_process(COMMAND bash "-c" "${TRIQS_GIT_EXECUTABLE} describe | sed \"s/-\\([0-9]*\\)-.*/.\\1/g\""
 WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
 OUTPUT_VARIABLE   CTHYB_VERSION
 OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT CTHYB_VERSION)
  set(CTHYB_VERSION 1.5)
endif()
triqs_get_git_hash_of_source_dir(CTHYB_GIT_HASH)
MESSAGE(STATUS "Cthyb version : ${CTHYB_VERSION}")
message(STATUS "Git hash: ${CTHYB_GIT_HASH}")

option(Build_Tests "Enable Tests" ON)
option(Build_Documentation "Build documentation" OFF)

# All PRIVATE common options.
# The std for all targets
add_compile_options(-std=c++17 $<$<CONFIG:Release>:-march=native>)

# Common to all targets (ctint, python, tests).
include_directories(c++)

# To be fixed!!
triqs_set_rpath_for_all_targets()

add_subdirectory(c++/triqs_cthyb)

# Tests
if (${Build_Tests})
 enable_testing()
 add_subdirectory(test)
endif()

if (${TRIQS_WITH_PYTHON_SUPPORT})

 # Python interface
 add_subdirectory(python/triqs_cthyb)

 # Build the documentation
 if (${Build_Documentation})
   add_subdirectory(doc)
 endif()

endif()

option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  SET(CPACK_GENERATOR "DEB")
  SET(CPACK_PACKAGE_VERSION ${CTHYB_VERSION})
  SET(CPACK_PACKAGE_CONTACT "https://github.com/TRIQS/cthyb")
  EXECUTE_PROCESS(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(CPACK_DEBIAN_PACKAGE_DEPENDS "python, python-numpy, python-h5py, python-jinja2, python-mako, python-mpi4py, python-matplotlib, python-scipy, cpp2py (>= 1.5), triqs (>= 2.0)")
  SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  SET(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  INCLUDE(CPack)
endif()
