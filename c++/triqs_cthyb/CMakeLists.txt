file(GLOB_RECURSE sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(REMOVE_ITEM sources impurity_trace.checks.cpp)
if(NOT MeasureG2)
  list(REMOVE_ITEM sources measures/G2_iw_acc.cpp measures/G2_iw.cpp measures/G2_iwll.cpp measures/G2_iw_nfft.cpp measures/G2_tau.cpp)
endif()

# Additional compilation options
option(EXT_DEBUG "Enable extended debugging output [developers only]" OFF)
option(TRACE_DEBUG_CHECKS "Enable various debug checks in the impurity trace computation [developers only]" OFF)
option(SAVE_CONFIGS "Save visited configurations to configs.h5 [developers only]" OFF)
set(NUM_CONFIGS_TO_SAVE 50000 CACHE STRING "Number of visited configurations to save [developers only]")

# ========= Additional Dependencies ==========

# 2-particle GF measurement requires NFFT
if(MeasureG2)
 find_package(NFFT REQUIRED)
endif()

# ========= Variant Library Function ==========

function(add_cthyb_library TARGET_NAME)
  cmake_parse_arguments(ARG "" "" "COMPILE_DEFS" ${ARGN})

  add_library(${TARGET_NAME} ${sources})
  add_library(${PROJECT_NAME}::${TARGET_NAME} ALIAS ${TARGET_NAME})

  # Link against triqs and enable warnings
  target_link_libraries(${TARGET_NAME} PUBLIC triqs PRIVATE $<BUILD_INTERFACE:${PROJECT_NAME}_warnings> $<BUILD_INTERFACE:${PROJECT_NAME}_performance_flags>)

  # Configure target and compilation
  set_target_properties(${TARGET_NAME} PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
  )
  target_include_directories(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
  target_include_directories(${TARGET_NAME} SYSTEM INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_PREFIX}/include>)
  target_link_directories(${TARGET_NAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
  target_link_directories(${TARGET_NAME} INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_LIBDIR}>)
  target_compile_definitions(${TARGET_NAME} PUBLIC
				CTHYB_GIT_HASH=${PROJECT_GIT_HASH}
				TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:CTHYB_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
				$<$<BOOL:${MeasureG2}>:CTHYB_G2_NFFT>
				$<$<BOOL:${EXT_DEBUG}>:EXT_DEBUG>
				$<$<BOOL:${TRACE_DEBUG_CHECKS}>:TRACE_DEBUG_CHECKS>
				$<$<BOOL:${SAVE_CONFIGS}>:SAVE_CONFIGS>
				$<$<BOOL:${SAVE_CONFIGS}>:NUM_CONFIGS_TO_SAVE=${NUM_CONFIGS_TO_SAVE}>
				${ARG_COMPILE_DEFS}
			  )

  # Link NFFT if needed
  if(MeasureG2)
    target_link_libraries(${TARGET_NAME} PUBLIC nfft)
  endif()

  install(TARGETS ${TARGET_NAME} EXPORT ${PROJECT_NAME}-targets DESTINATION ${CMAKE_INSTALL_LIBDIR})
endfunction()

# ========= Build Library Variants ==========

# Always build the real variant
add_cthyb_library(${PROJECT_NAME}_c)

# Optionally build complex hybridisation variant
if(Hybridisation_is_complex)
  add_cthyb_library(${PROJECT_NAME}_c_complex_hyb COMPILE_DEFS HYBRIDISATION_IS_COMPLEX)
endif()

# Optionally build complex hybridisation+hamiltonian variant
if(Local_hamiltonian_is_complex)
  if(NOT Hybridisation_is_complex)
    message(FATAL_ERROR "Local_hamiltonian_is_complex requires Hybridisation_is_complex=ON")
  endif()
  add_cthyb_library(${PROJECT_NAME}_c_complex_all COMPILE_DEFS HYBRIDISATION_IS_COMPLEX LOCAL_HAMILTONIAN_IS_COMPLEX)
endif()

# Install headers and NFFT target (once)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")
if(MeasureG2)
  install(TARGETS nfft EXPORT ${PROJECT_NAME}-targets)
endif()

# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(${PROJECT_NAME}_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
      TARGET ${PROJECT_NAME}_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++23
	  --template=gcc
	  --verbose
	  --force
	  --quiet
          ${sources}
	WORKING_DIRECTORY
	  ${CMAKE_CURRENT_SOURCE_DIR}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()

# ========= Dynamic Analyzer Checks ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
if(ASAN)
  if(NOT TARGET asan)
    find_package(sanitizer REQUIRED COMPONENTS asan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:asan>)
endif()

option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")
if(UBSAN)
  if(NOT TARGET ubsan)
    find_package(sanitizer REQUIRED COMPONENTS ubsan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:ubsan>)
endif()

option(MSAN OFF "Compile library and executables with LLVM Memory Sanitizer")
if(MSAN)
  if(NOT TARGET msan)
    find_package(sanitizer REQUIRED COMPONENTS msan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:msan>)
endif()

option(TSAN OFF "Compile library and executables with LLVM Thread Sanitizer")
if(TSAN)
  if(NOT TARGET tsan)
    find_package(sanitizer REQUIRED COMPONENTS tsan)
  endif()
  target_link_libraries(${PROJECT_NAME}_c PUBLIC $<BUILD_INTERFACE:tsan>)
endif()
