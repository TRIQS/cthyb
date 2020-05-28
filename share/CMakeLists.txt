add_subdirectory(cmake)

if(NOT CMAKE_INSTALL_PREFIX STREQUAL TRIQS_ROOT AND NOT IS_SUBPROJECT)

  if(PythonSupport)
    set(EXPORT_PYTHON_PATH "export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${CPP2PY_PYTHON_LIB_DEST_ROOT}:$PYTHONPATH")
    set(MODFILE_PYTHON_PATH "prepend-path    PYTHONPATH          $root/${CPP2PY_PYTHON_LIB_DEST_ROOT}")
  endif()

  configure_file(triqs_cthyb.modulefile.in triqs_cthyb.modulefile @ONLY)
  configure_file(triqs_cthybvars.sh.in triqs_cthybvars.sh @ONLY)

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/triqs_cthyb.modulefile
      ${CMAKE_CURRENT_BINARY_DIR}/triqs_cthybvars.sh
    DESTINATION share
  )
 
  message(STATUS "***************************************************************")
  message(STATUS "* Custom install Location. Use:                                ")
  message(STATUS "*                                                              ")
  message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/share/triqs_cthybvars.sh    ")
  message(STATUS "*                                                              ")
  message(STATUS "* to set up the environment variables                          ")
  message(STATUS "***************************************************************")

endif()
