# load triqs helper to set up tests
find_package(TriqsTest)

link_libraries(cthyb_c ${TRIQS_LIBRARY_ALL})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})

# Simple tests
set(SIMPLE_TESTS hilbert_test sorted_test rbt space_partition)
foreach(t ${SIMPLE_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    triqs_set_rpath_for_target(${t})
    triqs_add_test_C_simple(${t})
endforeach(t)

# HDF5 tests
set(HDF5_TESTS
    anderson anderson_block anderson_qn anderson_block_qn
    anderson_new
    spinless spinless_qn
    kanamori kanamori_qn
    kanamori_offdiag kanamori_offdiag_qn)

foreach(t ${HDF5_TESTS})
    add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
    triqs_set_rpath_for_target(${t})
    triqs_add_test_C_hdf(${t} " -d 1.e-7")
endforeach(t)

#add_subdirectory(long_tests)
