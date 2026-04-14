# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py")

# Add Python modules
c2py_add_module(configuration LINK_LIBRARIES ${PROJECT_NAME}_c ${PROJECT_NAME}_warnings)
c2py_add_module(solver_core LINK_LIBRARIES ${PROJECT_NAME}_c ${PROJECT_NAME}_warnings DEPENDS_ON_BINDINGS configuration)
install(TARGETS configuration solver_core DESTINATION ${PYTHON_LIB_DEST})

# ========= Complex Variant Modules ==========

# The complex variants reuse the bindings generated for the real variant and only differ
# in the cthyb library they link against. They are placed in their own subpackage so that
# the module names remain unchanged.
function(add_cthyb_variant VARIANT_SUFFIX VARIANT_LIB)
  set(output_dir ${CMAKE_CURRENT_BINARY_DIR}/_${VARIANT_SUFFIX})

  foreach(module_name configuration solver_core)
    set(variant_target ${module_name}_${VARIANT_SUFFIX})

    # Compile through a shim in the build tree. The binding generator replays every
    # compile command recorded for ${module_name}.cpp, so the source path of the
    # variants has to differ from the one of the real variant.
    set(variant_src ${output_dir}/${module_name}.cpp)
    file(WRITE ${variant_src} "#include \"${CMAKE_CURRENT_SOURCE_DIR}/${module_name}.cpp\"\n")

    Python_add_library(${variant_target} MODULE WITH_SOABI ${variant_src})
    set_target_properties(${variant_target} PROPERTIES
      OUTPUT_NAME ${module_name}
      LIBRARY_OUTPUT_DIRECTORY ${output_dir}
      CXX_VISIBILITY_PRESET hidden
      VISIBILITY_INLINES_HIDDEN ON
    )
    target_link_libraries(${variant_target} PRIVATE c2py::c2py ${VARIANT_LIB} ${PROJECT_NAME}_warnings)
    if(TARGET ${module_name}_bindings_generation)
      add_dependencies(${variant_target} ${module_name}_bindings_generation)
    endif()
    install(TARGETS ${variant_target} DESTINATION ${PYTHON_LIB_DEST}/_${VARIANT_SUFFIX})
  endforeach()

  # Generate and install __init__.py for the subpackage
  set(init_file ${output_dir}/__init__.py)
  file(WRITE ${init_file} "")
  install(FILES ${init_file} DESTINATION ${PYTHON_LIB_DEST}/_${VARIANT_SUFFIX})
endfunction()

if(Hybridisation_is_complex)
  add_cthyb_variant(complex_hyb ${PROJECT_NAME}_c_complex_hyb)
endif()

if(Local_hamiltonian_is_complex)
  add_cthyb_variant(complex_all ${PROJECT_NAME}_c_complex_all)
endif()
