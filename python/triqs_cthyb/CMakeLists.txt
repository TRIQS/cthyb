# Configure the version
configure_file(version.py.in version.py)

# ------- Python files ----------

# Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(file ${python_sources})
 configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py")

# ------- Python Extension ----------

find_package(Python COMPONENTS Interpreter Development NumPy)
find_package(Clair REQUIRED)

set(module_name solver_core)

Python_add_library( ${module_name} MODULE "${module_name}.${CLAIR_C2PY_MODULE_EXTENSION}")

target_include_directories(${module_name} SYSTEM PRIVATE  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(${module_name} PRIVATE ${PROJECT_NAME}_c clair::c2py Python::NumPy triqs_py)

# What is the use of this ? 
add_library(${PROJECT_NAME}::${module_name} ALIAS ${module_name})

install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST}/${module_dir})

if(ASAN)
  target_link_libraries(${module_name}  PRIVATE $<BUILD_INTERFACE:asan>)
endif()

if(UBSAN)
  target_link_libraries(${module_name}  PRIVATE $<BUILD_INTERFACE:ubsan>)
endif()


# OLD STUFF
## Build and install any python modules
#foreach(gen ${wrap_generators})
  #string(REPLACE "_desc.py" "" gen ${gen})
  #get_filename_component(module_name ${gen} NAME_WE)
  #get_filename_component(module_dir ${gen} DIRECTORY)
  #add_cpp2py_module(NAME ${module_name} DIRECTORY ${module_dir})
  #add_library(${PROJECT_NAME}::${module_name} ALIAS ${module_name})
  #target_link_libraries(${module_name} ${PROJECT_NAME}_c triqs_py)
  #install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST}/${module_dir})
#endforeach()


