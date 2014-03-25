# where will the python end up in triqs?
set(python_destination applications/impurity_solvers/cthyb_matrix2)

# make a local pytriqs copy
triqs_prepare_local_pytriqs_merged_with_my_python(${python_destination})

# Build the cython module
triqs_cython_module(cthyb_matrix2 cthyb_matrix2 ${python_destination})
target_link_libraries(cthyb_matrix2 ${TRIQS_LIBRARY_ALL} cthyb_matrix2_c)
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_SOURCE_DIR})
triqs_set_rpath_for_target(cthyb_matrix2)

# to be able to run from toplevel
execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/cthyb_matrix2.so ${CMAKE_BINARY_DIR}/pytriqs/${python_destination})

