cmake_minimum_required(VERSION 4.0 FATAL_ERROR)
project(image_denoising)

set(CMAKE_CXX_STANDARD 17)

add_library(core
            src/utils.cpp)

add_library(algorithm
            src/kernel_methods.cpp
            src/powerIteration.cpp
            src/toolkit.cpp)


target_include_directories(core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs)
target_include_directories(algorithm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libs)

find_package(Eigen3 REQUIRED)
target_link_libraries(algorithm PUBLIC Eigen3::Eigen)

add_executable(image_denoising src/apps/main.cpp)
add_executable(add_noise src/apps/add_noise.cpp)

set_property(TARGET image_denoising PROPERTY CXX_STANDARD 17)