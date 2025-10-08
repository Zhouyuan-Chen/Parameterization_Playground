set(EXAMPLE_SOURCE
    BFF.cpp
)

add_executable(BFF ${EXAMPLE_SOURCE})

target_link_libraries(BFF PUBLIC
    Parameterization
    Eigen3::Eigen
    igl::glfw
    geometry-central
    polyscope
    )