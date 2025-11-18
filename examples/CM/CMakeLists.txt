set(EXAMPLE_SOURCE
    CM.cpp
)

add_executable(CM ${EXAMPLE_SOURCE})

target_link_libraries(CM PUBLIC
    Parameterization
    Eigen3::Eigen
    igl::glfw
    geometry-central
    polyscope
    )