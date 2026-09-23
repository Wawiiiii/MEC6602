cmake_minimum_required(VERSION 3.15)
project(mec6602_hw1 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Eigen is header-only. Try to find an installed package first; otherwise
# fall back to EIGEN_DIR, which you can point at your local install, e.g.:
#   cmake -B build -DEIGEN_DIR="C:/Program Files/eigen-5.0.0"
find_package(Eigen3 QUIET NO_MODULE)

if(NOT Eigen3_FOUND)
    set(EIGEN_DIR "C:/Program Files/eigen-5.0.0" CACHE PATH "Path to Eigen headers")
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN_DIR}")
endif()

add_library(schemes src/schemes.cpp)
target_include_directories(schemes PUBLIC include)
target_link_libraries(schemes PUBLIC Eigen3::Eigen)


add_executable(main_cfl_study main_cfl_study.cpp)
target_link_libraries(main_cfl_study PRIVATE schemes)

add_executable(main_convergence main_convergence.cpp)
target_link_libraries(main_convergence PRIVATE schemes)


add_executable(main_euler1d_MacCormack main_euler1d_MacCormack.cpp)
target_link_libraries(main_euler1d_MacCormack PRIVATE schemes)