find_package(Catch2 REQUIRED)

add_executable(tests
    test_algorithm.cpp
    test_noise_addition.cpp)

target_link_libraries(tests PRIVATE Catch2::Catch2WithMain algorithm core)

include(CTest)
include(Catch)
catch_discover_tests(tests)