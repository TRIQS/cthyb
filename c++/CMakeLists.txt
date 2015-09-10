# The solver
add_library(cthyb_c solver_core.cpp sorted_spaces.cpp impurity_trace.cpp measure_density_matrix.cpp) 
target_link_libraries(cthyb_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(cthyb_c)

install(TARGETS cthyb_c DESTINATION lib)
