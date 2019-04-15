configure_file(triqs_cthyb-config.cmake.in triqs_cthyb-config.cmake @ONLY)
configure_file(triqs_cthyb-config-version.cmake.in triqs_cthyb-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/triqs_cthyb-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/triqs_cthyb-config-version.cmake
  DESTINATION lib/cmake/triqs_cthyb
)

install(EXPORT triqs_cthyb-targets DESTINATION lib/cmake/triqs_cthyb)
