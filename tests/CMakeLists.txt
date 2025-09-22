list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/tests")
include(catch2)
FetchContent_GetProperties(catch2)
include(Catch)

# Sources
set(TEST_SOURCES
    tiny_alligator.cpp
)
add_executable(tests ${TEST_SOURCES})

target_link_libraries(tests PUBLIC
    Catch2::Catch2WithMain
    MedialAxis
    )

catch_discover_tests(tests)