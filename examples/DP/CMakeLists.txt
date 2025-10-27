set(EXAMPLE_SOURCE
    DP.cpp
)

add_executable(DP ${EXAMPLE_SOURCE})

target_link_libraries(DP PUBLIC
    Parameterization
    Eigen3::Eigen
    igl::glfw
    geometry-central
    polyscope
    )