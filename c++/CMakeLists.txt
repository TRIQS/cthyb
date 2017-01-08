# The solver
set(LIBRARY_SOURCES
    solver_core.cpp
    atom_diag.cpp
    atom_diag_functions.cpp
    atom_diag_worker.cpp
    impurity_trace.cpp
    ./measures/density_matrix.cpp
    ./moves/shift.cpp
    ./moves/remove.cpp
    ./moves/insert.cpp
    ./moves/global.cpp
    ./moves/double_remove.cpp
    ./moves/double_insert.cpp
    )
if(MEASURE_G2)
 set(LIBRARY_SOURCES ${LIBRARY_SOURCES} measure_g2_inu.cpp measure_g2_legendre.cpp)
endif(MEASURE_G2)

add_library(cthyb_c ${LIBRARY_SOURCES})

target_link_libraries(cthyb_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(cthyb_c)

configure_file(config.hpp.in config.hpp)
install(TARGETS cthyb_c DESTINATION lib)

set(HEADERS_TO_INSTALL
    array_suppl.hpp # FIXME
    atom_diag.hpp
    atom_diag_functions.hpp
    solve_parameters.hpp
    solver_core.hpp
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.hpp DESTINATION include/cthyb)
install(FILES ${HEADERS_TO_INSTALL} DESTINATION include/cthyb)
