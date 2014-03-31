# where will the python end up in triqs?
set(python_destination applications/impurity_solvers/cthyb)

# make a local pytriqs copy
triqs_prepare_local_pytriqs_merged_with_my_python(${python_destination})

# Build the cython module
triqs_cython_module(cthyb cthyb ${python_destination})
target_link_libraries(cthyb ${TRIQS_LIBRARY_ALL} cthyb_c)
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_SOURCE_DIR})
triqs_set_rpath_for_target(cthyb)

# to be able to run from toplevel
execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/cthyb.so ${CMAKE_BINARY_DIR}/pytriqs/${python_destination})

