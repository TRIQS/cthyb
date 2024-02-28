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

# Find python
find_package(Python COMPONENTS Interpreter Development NumPy)

# Get c2py
include(FetchContent)
FetchContent_Declare(
  c2py
  GIT_REPOSITORY  https://github.com/flatironinstitute/c2py
  GIT_TAG        unstable 
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(c2py)

# 
option(GENERATE_PYTHON_BINDINGS "Use Clair python bindings generators" OFF)
if (GENERATE_PYTHON_BINDINGS)
  find_package(Clair REQUIRED)
endif()


set(module_name solver_core)

Python_add_library( ${module_name} MODULE ${module_name}.$<IF:$<BOOL:${GENERATE_PYTHON_BINDINGS}>,cpp,wrap.cxx>)

target_include_directories(${module_name} SYSTEM PRIVATE  ${CMAKE_SOURCE_DIR}/src)

target_link_libraries(${module_name} PRIVATE ${PROJECT_NAME}_c $<$<BOOL:${GENERATE_PYTHON_BINDINGS}>:clair::c2py_plugin> c2py::c2py Python::NumPy triqs_py)

# What is the use of this ? 
add_library(${PROJECT_NAME}::${module_name} ALIAS ${module_name})

install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST}/${module_dir})

if(ASAN)
  target_link_libraries(${module_name}  PRIVATE $<BUILD_INTERFACE:asan>)
endif()

if(UBSAN)
  target_link_libraries(${module_name}  PRIVATE $<BUILD_INTERFACE:ubsan>)
endif()

#--------------------------------------------------------
#  A target to copy the wrap.cxx files back to the src
#--------------------------------------------------------

add_custom_target(copy_wrap_cxx_to_src)
add_custom_command(
  TARGET copy_wrap_cxx_to_src
  COMMAND cp ${CMAKE_CURRENT_BINARY_DIR}/*.wrap.cxx ${CMAKE_CURRENT_SOURCE_DIR}
  COMMAND echo "** Copying all wrap.cxx to src"
)

