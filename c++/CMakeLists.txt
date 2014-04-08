add_library(cthyb_c ctqmc_matrix.cpp sorted_spaces.cpp atomic_correlators_worker.cpp gf_block_structure.cpp)
target_link_libraries(cthyb_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(cthyb_c)

install(TARGETS cthyb_c DESTINATION lib)
