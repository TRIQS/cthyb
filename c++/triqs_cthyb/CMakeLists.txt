set(LIBRARY_SOURCES
    parameters.cpp
    container_set.cpp
    solver_core.cpp
    impurity_trace.cpp
    moves/insert.cpp
    moves/remove.cpp
    moves/shift.cpp
    moves/global.cpp
    moves/double_insert.cpp
    moves/double_remove.cpp
    measures/density_matrix.cpp
    measures/G_tau.cpp
    measures/O_tau_ins.cpp
    measures/G_l.cpp
    )

#FIXME : for cmake > 3.1, use target_sources below
if(MeasureG2)
  list(APPEND LIBRARY_SOURCES
       measures/G2_tau.cpp
       measures/G2_iw.cpp
       measures/G2_iw_acc.cpp
       measures/G2_iw_nfft.cpp
       measures/G2_iwll.cpp
       )
endif()

add_library(cthyb_c ${LIBRARY_SOURCES})

# Link against the triqs library
target_link_libraries(cthyb_c PUBLIC triqs)

target_compile_options(cthyb_c PUBLIC -std=c++17 -fPIC)
target_include_directories(cthyb_c PUBLIC $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/c++>)
target_compile_definitions(cthyb_c
			  PUBLIC
			    $<$<BOOL:${Hybridisation_is_complex}>:HYBRIDISATION_IS_COMPLEX>
			    $<$<BOOL:${Local_hamiltonian_is_complex}>:LOCAL_HAMILTONIAN_IS_COMPLEX>
			  PRIVATE
			    CTHYB_GIT_HASH=${CTHYB_GIT_HASH}
			    TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
			    $<$<CONFIG:Debug>:DEBUG_CTHYB>
			    $<$<BOOL:${MeasureG2}>:CTHYB_G2_NFFT>
			    $<$<CONFIG:Debug>:TRIQS_DEBUG>
			    $<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			  )

# FIXME : should go ?
option(EXT_DEBUG "Enable extended debugging output [developers only]" OFF)
if(EXT_DEBUG)
 target_compile_definitions(cthyb_c PRIVATE EXT_DEBUG)
endif()

# FIXME : To be simplied
option(SAVE_CONFIGS "Save visited configurations to configs.h5 [developers only]" OFF)
if(SAVE_CONFIGS)
 set(NUM_CONFIGS_TO_SAVE 50000 CACHE STRING "Number of visited configurations to save [developers only]")
 target_compile_options(cthyb_c PRIVATE -DSAVE_CONFIGS -DNUM_CONFIGS_TO_SAVE=${NUM_CONFIGS_TO_SAVE})
endif()

# 2-particle GF measurement requires NFFT
if(MeasureG2)
 target_link_libraries(cthyb_c PRIVATE nfft)
 #target_sources(cthyb_c PRIVATE measures/G2_iw.cpp measures/G2_iwll.cpp)
endif()

# Install
install(TARGETS cthyb_c DESTINATION lib)

#triqs_set_rpath_for_target(cthyb_c)

# FIXME : WHY DO WE INSTALL HEADERS HERE ??
set(HEADERS_TO_INSTALL
 #array_suppl.hpp # FIXME
    parameters.hpp
    solver_core.hpp
    config.hpp
)
install(FILES ${HEADERS_TO_INSTALL} DESTINATION include/triqs/cthyb)
