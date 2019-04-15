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

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} triqs_cthyb_c gtest)
  add_test(${test_name} ${test_name})
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
      --std=c++14
      --template=gcc
      --verbose
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${test}.cpp
    )
  endif()
endforeach()

# add_test with various -DXX -DYYY
macro(add_test_defs test)
 set(t "${test}${ARGV1}")
 add_executable(${t} ${test})
 target_link_libraries(${t} PRIVATE triqs_cthyb_c gtest)
 target_compile_definitions(${t} PRIVATE "${ARGV2}") 
 add_test(${t} ${t})
endmacro()

add_test_defs(anderson _qn "QN")
add_test_defs(anderson _block "BLOCK")
add_test_defs(anderson _block_qn "BLOCK;QN")
add_test_defs(spinless _qn "QN")
add_test_defs(kanamori _qn "QN")
add_test_defs(kanamori_offdiag _qn "QN")
