find_package(TriqsTest)

triqs_add_test_hdf(kanamori_core.py " -d 1e-7" )
triqs_add_test_hdf(kanamori_kwargs.py " -d 1e-7" )
triqs_add_test_hdf(kanamori_params.py " -d 1e-7" )
triqs_add_test_hdf(single_site_bethe.py " -d 1e-7" )
triqs_add_test_hdf(slater.py " -d 1e-7" )

triqs_add_all_subdirectories_with_cmakelist()
