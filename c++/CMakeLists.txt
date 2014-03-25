add_library(cthyb_matrix2_c ctqmc_matrix.cpp sorted_spaces.cpp atomic_correlators_worker.cpp)
target_link_libraries(cthyb_matrix2_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL})
triqs_set_rpath_for_target(cthyb_matrix2_c)

install(TARGETS cthyb_matrix2_c DESTINATION lib)
