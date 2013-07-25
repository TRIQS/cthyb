add_executable(anderson anderson.cpp)
target_link_libraries(anderson triqs triqs_utility triqs_solver_ctqmc_krylov)

add_custom_target(anderson_py_scripts ALL)
add_custom_command(TARGET anderson_py_scripts PRE_BUILD
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/anderson.hyb.py
		${CMAKE_CURRENT_SOURCE_DIR}/anderson.plot.py
	COMMAND ${CMAKE_COMMAND} -DSRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}" -P
	${CMAKE_CURRENT_SOURCE_DIR}/../prepare_python_scripts.cmake)

add_custom_target(anderson_job_files ALL)
add_custom_command(TARGET anderson_job_files PRE_BUILD
	DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/hyb.job.in
		${CMAKE_CURRENT_SOURCE_DIR}/krylov.job.in
	COMMAND ${CMAKE_COMMAND} -DSRC_DIR="${CMAKE_CURRENT_SOURCE_DIR}" -P
	${CMAKE_CURRENT_SOURCE_DIR}/create_job_files.cmake)

# Prepare an ED calculation (make use of pomerol ED library)
find_package(PkgConfig)
pkg_search_module(POMEROL pomerol)
if(${POMEROL_FOUND})
	add_executable(anderson_ed anderson_ed.cpp)
	add_definitions(${POMEROL_CFLAGS})
	target_link_libraries(anderson_ed triqs triqs_utility ${POMEROL_LDFLAGS})
endif(${POMEROL_FOUND})

