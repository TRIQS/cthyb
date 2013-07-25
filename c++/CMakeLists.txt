triqs_add_rpath_to_triqs_install()
#triqs_add_rpath_for_cxx_code() 
add_library(cthyb_krylov ctqmc.cpp )
target_link_libraries (cthyb_krylov ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL})
install(TARGETS cthyb_krylov DESTINATION lib)

