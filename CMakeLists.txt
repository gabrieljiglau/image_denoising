cmake_minimum_required(VERSION 4.0 FATAL_ERROR)
project(image_denoising)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(core
            src/utils.cpp)

add_library(algorithm
            src/kernel_methods.cpp
            src/powerIteration.cpp
            src/toolkit.cpp)

add_library(lodepng src/libs/lodepng.cpp)
target_include_directories(lodepng PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs)

add_library(cli INTERFACE)
target_include_directories(cli INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/libs)

find_package(fmt REQUIRED)
find_package(Eigen3 REQUIRED)

target_link_libraries(core PUBLIC lodepng PUBLIC fmt::fmt)
target_link_libraries(algorithm PUBLIC core PUBLIC Eigen3::Eigen PUBLIC lodepng)

add_executable(image_denoising src/apps/main.cpp)
target_link_libraries(image_denoising PRIVATE algorithm cli)

add_executable(add_noise src/apps/add_noise.cpp)
target_link_libraries(add_noise PRIVATE cli core)

set_property(TARGET image_denoising PROPERTY CXX_STANDARD 17)