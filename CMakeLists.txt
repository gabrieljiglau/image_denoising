cmake_minimum_required(VERSION 4.0 FATAL_ERROR)
project(image_denoising)

set(CMAKE_CXX_STANDARD 17)

add_library(image_denoising
            src/utils.cpp
            src/kernel_methods.cpp
            src/powerIteration.cpp
            src/toolkit.cpp)


find_package(Eigen3 REQUIRED)
target_include_directories(image_denoising PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(image_denoising apps/main.cpp)
add_executable(add_noise apps/add_noise.cpp)

set_property(TARGET image_denoising PROPERTY CXX_STANDARD 17)