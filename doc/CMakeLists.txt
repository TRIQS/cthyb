if (NOT ${TRIQS_WITH_DOCUMENTATION})
     message("Error: TRIQS library has not been compiled with its documentation")
endif()

# generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py @ONLY)

# ---------------------------------
# Top Sphinx target
# ---------------------------------
# Sources
file(GLOB_RECURSE sources *.rst)

set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/contents.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${sources}
	COMMAND PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -j8 -b html ${CMAKE_CURRENT_BINARY_DIR} html)
add_custom_target(docs_sphinx ALL DEPENDS ${sphinx_top} ${CMAKE_CURRENT_BINARY_DIR})

# Application must be compiled first
add_dependencies(docs_sphinx solver_core) 

# ------------------------------------------------------------------------------------------------
# Copy the directory substructure and link all files with relevant extensions
# cp_rs is a script in cpp2py/bin, it mimics cp -rs on Linux
# and filters the relevant extension
# ------------------------------------------------------------------------------------------------

set(EXT_FOR_DOC "rst png txt css_t conf css js gif jpg py html bib h5 md ipynb")
execute_process(COMMAND cp_rs ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${EXT_FOR_DOC})

# ---------------------------------
# Install
# ---------------------------------
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/triqs_cthyb
  FILES_MATCHING
  REGEX "\\.(html|pdf|png|gif|jpg|svg|js|xsl|css|py|txt|inv|bib)$"
  PATTERN "_*"
)
