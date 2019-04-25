file(GLOB_RECURSE sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
list(REMOVE_ITEM sources impurity_trace.checks.cpp)
if(NOT MeasureG2)
  list(REMOVE_ITEM sources measures/G2_iw_acc.cpp measures/G2_iw.cpp measures/G2_iwll.cpp measures/G2_iw_nfft.cpp measures/G2_tau.cpp)
endif()
add_library(triqs_cthyb_c ${sources})

# Link against dependencies
target_link_libraries(triqs_cthyb_c PUBLIC triqs $<$<BOOL:${MeasureG2}>:nfft>)

# Additional compilation options
option(EXT_DEBUG "Enable extended debugging output [developers only]" OFF)
option(SAVE_CONFIGS "Save visited configurations to configs.h5 [developers only]" OFF)
set(NUM_CONFIGS_TO_SAVE 50000 CACHE STRING "Number of visited configurations to save [developers only]")

# Configure compilation
target_compile_options(triqs_cthyb_c PUBLIC -std=c++17 -fPIC)
target_include_directories(triqs_cthyb_c PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/c++>)
target_compile_definitions(triqs_cthyb_c PUBLIC
  				CTHYB_GIT_HASH=${CTHYB_GIT_HASH}
				TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:CTHYB_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
				$<$<BOOL:${Hybridisation_is_complex}>:HYBRIDISATION_IS_COMPLEX>
				$<$<BOOL:${Local_hamiltonian_is_complex}>:LOCAL_HAMILTONIAN_IS_COMPLEX>
				$<$<BOOL:${MeasureG2}>:CTHYB_G2_NFFT>
				$<$<BOOL:${EXT_DEBUG}>:EXT_DEBUG>
				$<$<BOOL:${SAVE_CONFIGS}>:SAVE_CONFIGS>
				$<$<BOOL:${SAVE_CONFIGS}>:NUM_CONFIGS_TO_SAVE=${NUM_CONFIGS_TO_SAVE}>
			  )

# Install library and headers
install(TARGETS triqs_cthyb_c EXPORT triqs_cthyb-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")


# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(triqs_cthyb_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET triqs_cthyb_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++14
	  --template=gcc
	  --verbose
	  --quiet
          ${sources}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()


# ========= Dynamic Analyzer Checks ==========

set(sanitizers "")
# Address Sanitizer
option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
if(ASAN)
  list(APPEND sanitizers asan)
  target_compile_options(triqs_cthyb_c PRIVATE -fsanitize=address -fno-omit-frame-pointer -ggdb3)
  target_link_libraries(triqs_cthyb_c INTERFACE "-fsanitize=address -fno-omit-frame-pointer")
  if(NOT DEFINED ENV{ASAN_OPTIONS})
    message(WARNING "ASAN_OPTIONS is not set. Consider setting ASAN_OPTIONS=symbolize=1:detect_leaks=0 when running tests")
  endif()
endif()

# Undefined Behavior Sanitizer
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")
if(UBSAN)
  list(APPEND sanitizers ubsan)
  target_compile_options(triqs_cthyb_c PUBLIC -fsanitize=undefined -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-omit-frame-pointer -ggdb3)
  target_link_libraries(triqs_cthyb_c INTERFACE "-fsanitize=undefined -fsanitize=float-divide-by-zero -fsanitize=float-cast-overflow -fno-omit-frame-pointer")
  if(NOT DEFINED ENV{UBSAN_OPTIONS})
    message(WARNING "UBSAN_OPTIONS is not set. Consider setting UBSAN_OPTIONS=symbolize=1:print_stacktrace=1:halt_on_error=1 when running tests")
  endif()
endif()

if(sanitizers)
  find_package(sanitizer REQUIRED ${sanitizers})
  if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
    target_link_libraries(triqs_cthyb_c INTERFACE "-fuse-ld=gold")
  endif()
endif()
