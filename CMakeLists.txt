list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

cmake_minimum_required(VERSION 3.24)

project(Parameterization VERSION 0.1)

add_library(Parameterization)

include(eigen)
include(libigl)
igl_include(glfw)
igl_include(restricted triangle)
include_directories(${CMAKE_SOURCE_DIR}/external/happly)
# igl_include(imgui)
add_subdirectory(src)
# add_subdirectory(examples)
# add_subdirectory(tests)
add_subdirectory(external/polyscope)
add_subdirectory(external/geometry-central)

add_executable(main main.cpp)

target_link_libraries(main PUBLIC
    Parameterization
    Eigen3::Eigen
    igl::glfw
    geometry-central
    polyscope
    )

enable_testing()
