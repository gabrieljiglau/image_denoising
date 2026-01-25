find_package(Catch2 REQUIRED)

add_executable(unit_tests
    test_algorithm.cpp
    test_noise_addition.cpp)

target_link_libraries(unit_tests PRIVATE Catch2::Catch2WithMain algorithm_library core)

include(Catch)
catch_discover_tests(unit_tests)