cmake_minimum_required(VERSION 4.0 FATAL_ERROR)
project(image_denoising)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# the source code itself 
add_library(core
            src/app/utils.cpp)

add_library(algorithm_library
            src/app/kernel_methods.cpp
            src/app/powerIteration.cpp
            src/app/toolkit.cpp)

# libraries used
add_library(lodepng src/libs/lodepng.cpp)
target_include_directories(lodepng PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs)

add_library(cli INTERFACE)
target_include_directories(cli INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/libs)

find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)

# linking
target_link_libraries(core PUBLIC lodepng PUBLIC fmt::fmt)
target_link_libraries(algorithm_library PUBLIC core PUBLIC Eigen3::Eigen PUBLIC lodepng)

# entrypoints
add_executable(algorithm src/app/entry_points/main.cpp)
target_link_libraries(algorithm PRIVATE algorithm_library cli)

add_executable(add_noise src/app/entry_points/add_noise.cpp)
target_link_libraries(add_noise PRIVATE cli core)

# tests
include(CTest)
enable_testing()
add_subdirectory(tests)