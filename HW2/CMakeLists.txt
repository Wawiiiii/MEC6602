cmake_minimum_required(VERSION 3.16)
project(mec6602_hw2 CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(MSVC)
    add_compile_options(/W4 /EHsc)
else()
    add_compile_options(-Wall -Wextra -Wpedantic)
    # No -ffast-math: it lets the compiler assume NaN never occurs, which would
    # silently disable the NaN checks on unfilled ghost cells.
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
endif()

# Solver library: every src/*.cpp file gets added here as we go.
add_library(euler2d
    src/mesh.cpp
    src/output.cpp
    src/flow.cpp
    src/roe.cpp
    src/residual.cpp
    src/bc.cpp
    src/solver.cpp
)
target_include_directories(euler2d PUBLIC include)

# Step 1 check: read a plot3D mesh and print it.
add_executable(test_mesh tests/test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE euler2d)

# Step 2 check: metrics, constant-flux residual, write mesh with areas (Euler2D slides 21-22).
add_executable(test_metrics tests/test_metrics.cpp)
target_link_libraries(test_metrics PRIVATE euler2d)

# Step 3 check: Roe flux, MUSCL + Hemker-Koren, residual (Roe_scheme_accel slides 5-16).
add_executable(test_roe tests/test_roe.cpp)
target_link_libraries(test_roe PRIVATE euler2d)

# Step 4 check: boundary conditions with ghost cells (basic_algo slides 7, 13-16, 20).
add_executable(test_bc tests/test_bc.cpp)
target_link_libraries(test_bc PRIVATE euler2d)

# The solver program: euler <mesh.x> <mach> <alpha_deg> [options].
add_executable(euler apps/euler_main.cpp)
target_link_libraries(euler PRIVATE euler2d)

# Step 5 check: explicit Euler / multistage RK with local time stepping (basic_algo slides 17-18).
add_executable(test_explicit tests/test_explicit.cpp)
target_link_libraries(test_explicit PRIVATE euler2d)
