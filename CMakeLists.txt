cmake_minimum_required(VERSION 4.0 FATAL_ERROR)
project(image_denoising)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_PREFIX_PATH "$ENV{HOME}/libs/libtorch")

find_package(Torch REQUIRED)

add_subdirectory(src/custom_network)
add_subdirectory(src/libtorch)

add_executable(fashion_MNIST src/apps/vanilla_main.cpp)
target_link_libraries(fashion_MNIST custom_network)

add_executable(torch_only src/apps/torch_main.cpp)
target_link_libraries(torch_only libtorch)

set_property(TARGET image_denoising PROPERTY CXX_STANDARD 17)