find_package(TriqsTest)

triqs_add_test_txt(operator_python_test)
triqs_add_test_hdf(solver_python_test " -p 1e-6" )

triqs_add_all_subdirectories_with_cmakelist()