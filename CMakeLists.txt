# Start configuration
cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
project(triqs_cthyb CXX)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message( STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")

# Options for the compilation ...
option(Hybridisation_is_complex "If ON, the hybridization Delta(tau) is complex" OFF)
option(Local_hamiltonian_is_complex "If ON, the H_loc is complex" OFF)
option(MeasureG2 "Measure the two particle object (requires the NFFT library)" ON)

# check that options are compatible
if(Hybridisation_is_complex AND NOT Local_hamiltonian_is_complex)
 message(FATAL_ERROR "Combination of options Hybridisation_is_complex=ON "
                     "and Local_hamiltonian_is_complex=OFF is not supported.")
endif()

# Build static libraries
set(BUILD_SHARED_LIBS OFF)

# Export the list of compile-commands into compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Enable compiler warnings for the whole project
add_compile_options(-Wall
		     $<$<CONFIG:Debug>:-Og>
		     $<$<CONFIG:Debug>:-ggdb3>
		     )

# Load Dependencies
find_package(TRIQS 2.2 REQUIRED)
find_package(Cpp2Py 1.6 REQUIRED)
if(MeasureG2) # 2-particle GF measurement requires NFFT
 find_package(NFFT REQUIRED)
 install(TARGETS nfft EXPORT triqs_cthyb-targets)
endif()

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS " No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# Define the triqs_cthyb version numbers and get the git hash
set(CTHYB_VERSION_MAJOR 2)
set(CTHYB_VERSION_MINOR 2)
set(CTHYB_VERSION_PATCH 0)
set(CTHYB_VERSION ${CTHYB_VERSION_MAJOR}.${CTHYB_VERSION_MINOR}.${CTHYB_VERSION_PATCH})
triqs_get_git_hash_of_source_dir(CTHYB_GIT_HASH)
message(STATUS "triqs_cthyb version : ${CTHYB_VERSION}")
message(STATUS "Git hash: ${CTHYB_GIT_HASH}")

# Build and install the cthyb library
add_subdirectory(c++/triqs_cthyb)

# Tests
option(Build_Tests "Build tests" ON)
if(Build_Tests)
  enable_testing()
  add_subdirectory(test)
endif()

if(TRIQS_WITH_PYTHON_SUPPORT)
  # Python interface
  add_subdirectory(python/triqs_cthyb)

  # Build the documentation
  option(Build_Documentation "Build documentation" OFF)
  if(Build_Documentation)
    if(NOT TRIQS_WITH_DOCUMENTATION)
      message(WARNING "TRIQS library has not been compiled with its documentation! Cannot build documentation.")
    else()
      message(STATUS "Documentation Build enabled")
      add_subdirectory(doc)
    endif()
  endif()
else()
  message(WARNING "TRIQS library has been installed without Python support. Cannot build the Python Interface and Documentation.")
endif()

# Additional configuration files
add_subdirectory(share)

# Debian Package
option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  SET(CPACK_GENERATOR "DEB")
  SET(CPACK_PACKAGE_VERSION ${CTHYB_VERSION})
  SET(CPACK_PACKAGE_CONTACT "https://github.com/TRIQS/cthyb")
  EXECUTE_PROCESS(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
  SET(CPACK_DEBIAN_PACKAGE_DEPENDS "triqs (>= 2.2)")
  SET(CPACK_DEBIAN_PACKAGE_CONFLICTS "cthyb")
  SET(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  SET(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  INCLUDE(CPack)
endif()
