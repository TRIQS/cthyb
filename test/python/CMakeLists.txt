find_package(TriqsTest)

# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

triqs_add_python_test(single_site_bethe)
triqs_add_python_test(atomic_observables)
triqs_add_python_test(kanamori)
triqs_add_python_test(slater)
triqs_add_python_test(measure_static)
triqs_add_python_test(histograms)
