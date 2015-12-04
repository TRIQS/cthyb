# load triqs helper to set up tests
find_package(TriqsTest)

# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

link_libraries(cthyb_c ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
set(SIMPLE_TESTS sorted_test rbt)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    triqs_set_rpath_for_target(${t})
    triqs_add_cpp_test(${t})
endforeach(t)

# HDF5 tests
macro(add_test_defs name)
    set(full_name "${name}")
    if(${ARGC} GREATER 1)
        set(full_name "${name}_${ARGV1}")
    endif()
    add_executable(${full_name} "${name}.cpp")
    triqs_set_rpath_for_target(${full_name})
    if(${ARGC} GREATER 2)
        set(defs ${ARGV2})
        separate_arguments(defs)
        set_target_properties(${full_name} PROPERTIES COMPILE_DEFINITIONS "${defs}")
    endif()
    triqs_add_cpp_test(${full_name})
endmacro()

add_test_defs(anderson)
add_test_defs(anderson block "BLOCK")
add_test_defs(anderson qn "QN")
add_test_defs(anderson block_qn "BLOCK;QN")

add_test_defs(spinless)
add_test_defs(spinless "qn" "QN")

add_test_defs(kanamori)
add_test_defs(kanamori "qn" "QN")

add_test_defs(kanamori_offdiag)
add_test_defs(kanamori_offdiag "qn" "QN")

add_test_defs(legendre)
add_test_defs(atomic_gf)
