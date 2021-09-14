# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests anderson spinless kanamori kanamori_offdiag legendre rbt impurity_trace_atomic_gf impurity_trace_bug_try_insert impurity_trace_op_insert)
if(MeasureG2)
  list(APPEND all_tests G2)
endif()
#file(GLOB_RECURSE all_tests RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} ${PROJECT_NAME}::${PROJECT_NAME}_c ${PROJECT_NAME}_warnings gtest_main)
  set_property(TARGET ${test_name} PROPERTY RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  add_test(NAME ${test_name} COMMAND ${test_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${test_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()
  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${test_name}
      COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++20
      --template=gcc
      --verbose
      --force
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${test}
    )
  endif()
endforeach()

# add_test with various -DXX -DYYY
macro(add_test_defs test)
 set(t "${test}${ARGV1}")
 add_executable(${t} ${test})
 target_link_libraries(${t} PRIVATE ${PROJECT_NAME}::${PROJECT_NAME}_c gtest_main)
 target_compile_definitions(${t} PRIVATE "${ARGV2}") 
 add_test(${t} ${t})
endmacro()

add_test_defs(anderson _qn "QN")
add_test_defs(anderson _block "BLOCK")
add_test_defs(anderson _block_qn "BLOCK;QN")
add_test_defs(spinless _qn "QN")
add_test_defs(kanamori _qn "QN")
add_test_defs(kanamori_offdiag _qn "QN")
