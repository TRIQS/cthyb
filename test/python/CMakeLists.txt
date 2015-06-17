find_package(TriqsTest)

FILE(COPY DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

triqs_add_test_hdf(kanamori " -d 1e-5" )
triqs_add_test_hdf(single_site_bethe " -d 1e-5" )
triqs_add_test_hdf(slater " -d 1e-5" )

triqs_add_all_subdirectories_with_cmakelist()
