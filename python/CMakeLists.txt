# where will the python end up in triqs?
set(python_destination applications/impurity_solvers/cthyb_krylov)

# make a local pytriqs copy
triqs_prepare_local_pytriqs_merged_with_my_python(${python_destination})

# Build the cython module for many_body_operator
triqs_cython_module(many_body_operator many_body_operator ${python_destination})
target_link_libraries(many_body_operator ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_SOURCE_DIR})
triqs_set_rpath_for_target(many_body_operator)

execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/many_body_operator.so ${CMAKE_BINARY_DIR}/pytriqs/${python_destination})