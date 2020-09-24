cmake_minimum_required(VERSION 3.5)

project(opticpp LANGUAGES CXX CUDA)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CUDA_STANDARD 20)
set(CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_BUILD_TYPE Release)

find_package(CUDA 10.2)
message(STATUS "Found CUDA ${CUDA_VERSION_STRING} at ${CUDA_TOOLKIT_ROOT_DIR}")


set(GCC_COVERAGE_COMPILE_FLAGS "-ltbb -std=c++20")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(opticpp
                main.cpp
                mainwindow.cpp
                mainwindow.h
                optics/optics.h
                optics/optics_gpu.cuh   
                mainwindow.ui
                polarization.cpp
                polarization.h
                polarization.ui
                )

target_link_libraries(opticpp PRIVATE Qt5::Widgets)
target_compile_options(opticpp PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-gencode arch=compute_30,code=sm_30>)
target_include_directories(opticpp  PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(opticpp  ${THRUST_LIBRARY})
