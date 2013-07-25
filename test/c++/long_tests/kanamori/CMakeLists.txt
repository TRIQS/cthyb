add_executable(kanamori kanamori.cpp)
target_link_libraries(kanamori triqs triqs_utility triqs_solver_ctqmc_krylov)

add_custom_target(kanamori_py_scripts ALL)
add_custom_command(TARGET kanamori_py_scripts PRE_BUILD
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/kanamori.hyb.py
		${CMAKE_CURRENT_SOURCE_DIR}/kanamori.plot.py
	COMMAND ${CMAKE_COMMAND} -DSRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}" -P
	${CMAKE_CURRENT_SOURCE_DIR}/../prepare_python_scripts.cmake)

add_custom_target(kanamori_job_files ALL)
add_custom_command(TARGET kanamori_job_files PRE_BUILD
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/hyb.job.in
		${CMAKE_CURRENT_SOURCE_DIR}/krylov.job.in
	COMMAND ${CMAKE_COMMAND} -DSRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}" -P
	${CMAKE_CURRENT_SOURCE_DIR}/create_job_files.cmake)

# Prepare an ED calculation (make use of pomerol ED library)
find_package(PkgConfig)
pkg_search_module(POMEROL pomerol)
if(${POMEROL_FOUND})
	add_executable(kanamori_ed kanamori_ed.cpp)
	add_definitions(${POMEROL_CFLAGS})
	target_link_libraries(kanamori_ed triqs triqs_utility ${POMEROL_LDFLAGS})
endif(${POMEROL_FOUND})

