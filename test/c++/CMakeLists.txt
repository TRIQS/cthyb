# Copy h5 files to binary dir
file(GLOB all_h5_files *.h5)
file(COPY ${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# add_test with possibly various -DXX -DYYY
macro(add_test_defs name)
 set(t "${name}${ARGV1}")
 add_executable(${t} ${name}.cpp)
 target_link_libraries(${t} PRIVATE cthyb_c gtest)
 foreach(opt ${ARGV2})
  target_compile_definitions(${t} PRIVATE ${opt}) 
 endforeach()
 add_test(NAME ${t} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${t})
endmacro()

#
if(NFFT_FOUND)
 # Disable nfft test, scheduled for removal
 #add_test_defs(test_nfft)
 #target_link_libraries(test_nfft PRIVATE nfft)
endif(NFFT_FOUND)

add_test_defs(anderson)
add_test_defs(anderson _block "BLOCK")
add_test_defs(anderson _qn "QN")
add_test_defs(anderson _block_qn "BLOCK;QN")

add_test_defs(spinless)
add_test_defs(spinless _qn "QN")

add_test_defs(kanamori)
add_test_defs(kanamori _qn "QN")

add_test_defs(kanamori_offdiag)
add_test_defs(kanamori_offdiag _qn "QN")

add_test_defs(legendre)
add_test_defs(G2)

add_test_defs(rbt)

add_test_defs(impurity_trace)

# Not ported, should be checked by atom_diag
#add_test_defs(h_diag_test)
#add_test_defs(atomic_gf)
