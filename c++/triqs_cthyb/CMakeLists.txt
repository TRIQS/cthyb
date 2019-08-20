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
target_link_libraries(cthyb_c PUBLIC triqs)

# Pubic options : any option affecting the headers should go here
target_compile_options(cthyb_c PUBLIC 
                     $<$<BOOL:${Hybridisation_is_complex}>:-DHYBRIDISATION_IS_COMPLEX>
                     $<$<BOOL:${Local_hamiltonian_is_complex}>:-DLOCAL_HAMILTONIAN_IS_COMPLEX>
		     )

# Private options : any option here shoud affect ONLY the cpp
target_compile_options(cthyb_c PRIVATE 
                    -DCTHYB_GIT_HASH=${CTHYB_GIT_HASH} -DTRIQS_GIT_HASH=${TRIQS_GIT_HASH}
		     #$<$<CONFIG:Debug>:-DDEBUG_CTHYB>
                     $<$<BOOL:${MeasureG2}>:-DCTHYB_G2_NFFT>
		     ) 

# FIXME : should go ?
option(EXT_DEBUG "Enable extended debugging output [developers only]" OFF)

if(EXT_DEBUG)
 target_compile_options(cthyb_c PRIVATE -DEXT_DEBUG)
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
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")
