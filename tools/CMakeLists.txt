add_executable(analysis analyze_static_ops.cpp)
target_link_libraries(analysis cthyb_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(analysis)

add_executable(simple test.cpp)
target_link_libraries(simple cthyb_c ${TRIQS_LIBRARY_ALL})
include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
triqs_set_rpath_for_target(simple)


#install(TARGETS analysis DESTINATION lib)
