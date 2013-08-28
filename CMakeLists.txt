# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(ctkrylov CXX)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS REQUIRED)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Details of the krylov code
option(Extended_Debug "Enable extended debugging output" OFF)
option(Stats "Enable statistics collection" OFF)

if(Extended_Debug)
 add_definitions(-DEXT_DEBUG)
endif(Extended_Debug)

if(Stats)
 add_definitions(-DKRYLOV_STATS)
 set(Dims_Stats_File "dims.stats.dat" CACHE STRING "Statistics file: subspace dimensions")
 set(Boundary_Stats_File "boundary_states.stats.dat" CACHE STRING "Statistics file: boundary states")
 add_definitions(-DDIMS_STATS_FILE=\"${Dims_Stats_File}\")
 add_definitions(-DBOUNDARY_STATS_FILE=\"${Boundary_Stats_File}\")
endif(Stats)

# Compile C++ code
add_subdirectory(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
endif()

# Add tests
add_subdirectory(test)

# Link the pytriqs for convenience ?? 
execute_process(COMMAND ln -fs ${TRIQS_PATH}/bin .)


