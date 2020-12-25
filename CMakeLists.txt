cmake_minimum_required(VERSION 3.0)

project(opticpp LANGUAGES CXX)

    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(Qt5 COMPONENTS Widgets REQUIRED)
find_package(Qt5 COMPONENTS PrintSupport REQUIRED)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)

set(GCC_COVERAGE_COMPILE_FLAGS "-ltbb -std=c++20")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS}")

add_executable(opticpp
                main.cpp
                mainwindow.cpp
                mainwindow.h
                optics/optics.h
                mainwindow.ui
                polarization.cpp
                polarization.h
                polarization.ui
                force.cpp
                force.h
                force.ui
                qcustomplot/qcustomplot.h
                qcustomplot/qcustomplot.cpp
                )

target_link_libraries(opticpp PRIVATE Qt5::Widgets)
target_link_libraries(opticpp PRIVATE Qt5::PrintSupport)
