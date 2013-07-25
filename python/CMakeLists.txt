# Build the cython module
triqs_add_rpath_for_cxx_code() # necessary for python module to load the C++ so file
triqs_cython_module(ctqmc_seg ctseg_solver applications/impurity_solvers/ctqmc_seg)
target_link_libraries(ctseg_solver triqs_solver_ctqmc_seg)
include_directories( ${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR}/../ )

# to be able to run from toplevel
execute_process(COMMAND ln -fs ${CMAKE_CURRENT_BINARY_DIR}/ctseg_solver.so ${CMAKE_BINARY_DIR})

# Install 
SET(PYTHON_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py )
install (FILES ${PYTHON_SOURCES} DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/impurity_solvers/ctqmc_seg)

#make sure that the __init__.py empty files are present ...
file (WRITE __init__.py.all_empty  "__all__=[]")
install (FILES __init__.py.all_empty DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/impurity_solvers/ RENAME __init__.py)
install (FILES __init__.py.all_empty DESTINATION ${TRIQS_PYTHON_LIB_DEST}/applications/                  RENAME __init__.py)


