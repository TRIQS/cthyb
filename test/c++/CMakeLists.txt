# load triqs helper to set up tests
find_package(TriqsTest)

if (${TRIQS_WITH_PYTHON_SUPPORT})
 #add_subdirectory(python)
endif()

link_libraries (cthyb_krylov ${TRIQS_LIBRARY_ALL} )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/../../c++ ${TRIQS_INCLUDE_ALL})
 
# Simple tests
add_executable( hilbert_test  ${CMAKE_CURRENT_SOURCE_DIR}/hilbert_test.cpp )
add_executable( operator_test  ${CMAKE_CURRENT_SOURCE_DIR}/operator_test.cpp )
add_executable( sorted_test  ${CMAKE_CURRENT_SOURCE_DIR}/sorted_test.cpp )
add_executable( krylov_test  ${CMAKE_CURRENT_SOURCE_DIR}/krylov_test.cpp )
triqs_add_test_C_simple( hilbert_test )
triqs_add_test_C_simple( operator_test )
triqs_add_test_C_simple( sorted_test )
triqs_add_test_C_simple( krylov_test )

# HDF5 tests
add_executable( anderson_block ${CMAKE_CURRENT_SOURCE_DIR}/anderson_block.cpp )
add_executable( anderson_block_qn ${CMAKE_CURRENT_SOURCE_DIR}/anderson_block_qn.cpp )
add_executable( anderson ${CMAKE_CURRENT_SOURCE_DIR}/anderson.cpp )
add_executable( anderson_qn ${CMAKE_CURRENT_SOURCE_DIR}/anderson_qn.cpp )
add_executable( kanamori ${CMAKE_CURRENT_SOURCE_DIR}/kanamori.cpp )
add_executable( kanamori_offdiag ${CMAKE_CURRENT_SOURCE_DIR}/kanamori_offdiag.cpp )
add_executable( kanamori_offdiag_qn ${CMAKE_CURRENT_SOURCE_DIR}/kanamori_offdiag_qn.cpp )
add_executable( kanamori_qn ${CMAKE_CURRENT_SOURCE_DIR}/kanamori_qn.cpp )
add_executable( spinless ${CMAKE_CURRENT_SOURCE_DIR}/spinless.cpp )
add_executable( spinless_qn ${CMAKE_CURRENT_SOURCE_DIR}/spinless_qn.cpp )
triqs_add_test_C_hdf( anderson_block " -d 1.e-7" )
triqs_add_test_C_hdf( anderson_block_qn " -d 1.e-7" )
triqs_add_test_C_hdf( anderson " -d 1.e-7" )
triqs_add_test_C_hdf( anderson_qn " -d 1.e-7" )
triqs_add_test_C_hdf( kanamori " -d 1.e-7" )
triqs_add_test_C_hdf( kanamori_offdiag " -d 1.e-7" )
triqs_add_test_C_hdf( kanamori_offdiag_qn " -d 1.e-7" )
triqs_add_test_C_hdf( kanamori_qn " -d 1.e-7" )
triqs_add_test_C_hdf( spinless " -d 1.e-7" )
triqs_add_test_C_hdf( spinless_qn " -d 1.e-7" )

#add_subdirectory(long_tests)
