# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# List all tests
set(all_tests setup_Delta_tau_and_h_loc single_site_bethe atomic_observables kanamori slater measure_static histograms move_global h5_read_write O_tau_ins)

if(Local_hamiltonian_is_complex)
 list(APPEND all_tests atomic_gf_complex atomdiag_ed complex_bug81)
 if(Hybridisation_is_complex)
  list(APPEND all_tests complex_Gtau_ED setup_Delta_tau_and_h_loc_complex)
 endif()
endif()

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# Set the PythonPath : put the build dir first (in case there is an installed version). 
set_property(TEST ${all_tests} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )


