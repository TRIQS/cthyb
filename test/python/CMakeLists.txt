find_package(TriqsTest)

triqs_add_test_hdf(solvercore " -d 1e-7" )
triqs_add_test_hdf(solver_params " -d 1e-7" )
triqs_add_test_hdf(solver_kwargs " -d 1e-7" )

triqs_add_all_subdirectories_with_cmakelist()
