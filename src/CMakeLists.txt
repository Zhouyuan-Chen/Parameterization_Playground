list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(libigl)
igl_include(core)

target_include_directories(Parameterization PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Parameterization PUBLIC geometry-central Eigen3::Eigen igl::core)
add_subdirectory(BFF)