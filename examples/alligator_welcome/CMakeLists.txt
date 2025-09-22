set(EXAMPLE_SOURCE
    alligator_welcome.cpp
)

add_executable(alligator_welcome ${EXAMPLE_SOURCE})

target_link_libraries(alligator_welcome PUBLIC
    TinyAlligator
    
    )