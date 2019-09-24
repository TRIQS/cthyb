# ##############################################################################
#
# Copyright (C) 2017-2018 N. Wentzell, O. Parcollet
# Copyright (C) 2019 The Simons Foundation
#  authors: N. Wentzell, P. Dumitrescu, H. Strand
#
# triqs_cthyb is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# triqs_cthyb is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# triqs_cthyb (in the file COPYING.txt in this directory). If not, see
# <http://www.gnu.org/licenses/>.
#
# ##############################################################################

cmake_minimum_required(VERSION 3.0.2 FATAL_ERROR)
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

# ############
# Define Project
project(CTHYB VERSION 2.2.0 LANGUAGES C CXX)

# ############
# Load TRIQS and CPP2PY
find_package(TRIQS 2.2 REQUIRED)
find_package(Cpp2Py 1.6 REQUIRED)

# Get the git hash & print status
triqs_get_git_hash_of_source_dir(PROJECT_GIT_HASH)
message(STATUS "${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "${PROJECT_NAME} Git hash: ${PROJECT_GIT_HASH}")

# Enforce Consistent Versioning
if(NOT ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR} VERSION_EQUAL ${TRIQS_VERSION_MAJOR}.${TRIQS_VERSION_MINOR})
  message(FATAL_ERROR "The ${PROJECT_NAME} version ${PROJECT_VERSION} is not compatible with TRIQS version ${TRIQS_VERSION}.")
endif()

# Default Install directory to TRIQS_ROOT if not given or invalid.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS "No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} --------")

# ############
# Options

# Default to Release build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Type of build" FORCE)
endif()
message(STATUS "-------- BUILD-TYPE: ${CMAKE_BUILD_TYPE} --------")

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

# Global compiler options
add_compile_options(
  $<$<CONFIG:Debug>:-Og>
  $<$<CONFIG:Debug>:-ggdb3>
)

# Create an Interface target for compiler warnings
add_library(project_warnings INTERFACE)
install(TARGETS project_warnings EXPORT triqs_cthyb-targets)
target_compile_options(project_warnings
  INTERFACE
    -Wall
    -Wextra
    -Wshadow
    -Wpedantic
    -Wno-sign-compare
    $<$<CXX_COMPILER_ID:Clang>:-Wno-gcc-compat>
)

# #############
# Build Project

# Install additional triqs headers provided by cthyb
add_subdirectory(c++/triqs)

# Build and install the triqs_cthyb library
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
      message(WARNING "TRIQS was installed without documentation! Cannot build documentation.")
    else()
      message(STATUS "Documentation Build enabled")
      add_subdirectory(doc)
    endif()
  endif()
else()
  message(WARNING "TRIQS was installed without Python support. Cannot build the Python Interface and Documentation.")
endif()

# Additional configuration files
add_subdirectory(share)

# #############
# Debian Package

option(BUILD_DEBIAN_PACKAGE "Build a deb package" OFF)
if(BUILD_DEBIAN_PACKAGE)
  if(NOT CMAKE_INSTALL_PREFIX STREQUAL "/usr")
    message(FATAL_ERROR "CMAKE_INSTALL_PREFIX must be /usr for packaging")
  endif()
  set(CPACK_PACKAGE_NAME triqs_cthyb)
  set(CPACK_GENERATOR "DEB")
  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
  set(CPACK_PACKAGE_CONTACT "https://github.com/TRIQS/triqs_cthyb")
  execute_process(COMMAND dpkg --print-architecture OUTPUT_VARIABLE CMAKE_DEBIAN_PACKAGE_ARCHITECTURE OUTPUT_STRIP_TRAILING_WHITESPACE)
  set(CPACK_DEBIAN_PACKAGE_DEPENDS "triqs (>= 2.2)")
  set(CPACK_DEBIAN_PACKAGE_CONFLICTS "cthyb")
  set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
  set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS ON)
  include(CPack)
endif()

