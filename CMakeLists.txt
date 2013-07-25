# Append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

# start configuration 
cmake_minimum_required(VERSION 2.8)
project(ctseg CXX)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS, including all predefined variables from TRIQS installation
find_package(TRIQS)

# We want to be installed in the TRIQS tree
set(CMAKE_INSTALL_PREFIX ${TRIQS_PATH})

# Details of the krylov code
option(Krylov_Debug "Enable extended debugging output in Krylov solver" OFF)
option(Krylov_Stats "Enable statistics collection in Krylov solver" OFF)
mark_as_advanced(Krylov_Debug)
mark_as_advanced(Krylov_Stats)

if(Krylov_Debug)
add_definitions(-DKRYLOV_DEBUG)
endif(Krylov_Debug)

if(Krylov_Stats)
add_definitions(-DKRYLOV_STATS)
endif(Krylov_Stats)

# Compile C++ code
add_subdirectory(c++)

# Python interface
if (${TRIQS_WITH_PYTHON_SUPPORT})
 #add_subdirectory(python)
endif()

# Add tests
add_subdirectory(test)

# Link the pytriqs for convenience ?? 
execute_process(COMMAND ln -fs ${TRIQS_PATH}/bin .)


