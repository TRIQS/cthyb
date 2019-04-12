# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests setup_Delta_tau_and_h_loc single_site_bethe atomic_observables kanamori_py slater measure_static histograms move_global h5_read_write O_tau_ins)
if(Local_hamiltonian_is_complex)
  list(APPEND all_tests atomic_gf_complex atomdiag_ed complex_bug81)
  if(Hybridisation_is_complex)
    list(APPEND all_tests complex_Gtau_ED setup_Delta_tau_and_h_loc_complex)
  endif()
endif()

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME ${test_name} COMMAND ${TRIQS_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST ${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endforeach()
