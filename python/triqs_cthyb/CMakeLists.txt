# Build the python module
add_cpp2py_module(solver_core)

target_link_libraries(solver_core cthyb_c)

# We need to include the convertes.hxx files
target_include_directories(solver_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Configure the version
configure_file(version.py.in version.py)

# All Python file. Copy them in the build dir to have a complete package for the tests.
set(PYTHON_SOURCES __init__.py util.py solver.py tail_fit.py)

foreach(f ${PYTHON_SOURCES})
 configure_file(${f} ${f} COPYONLY)
endforeach()

# Install python module to proper location
set(PYTHON_LIB_DEST ${CPP2PY_PYTHON_LIB_DEST_ROOT}/triqs_cthyb)
install(TARGETS solver_core DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${PYTHON_SOURCES} ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})


