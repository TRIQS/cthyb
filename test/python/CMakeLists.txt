find_package(TriqsTest)

triqs_add_test_hdf(solvercore_python_test " -d 1e-7" )
triqs_add_test_hdf(solver_python_test " -d 1e-7" )

triqs_add_all_subdirectories_with_cmakelist()
