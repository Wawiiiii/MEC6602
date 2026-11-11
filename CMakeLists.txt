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
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
endif()

# Solver library: every src/*.cpp file gets added here as we go.
add_library(euler2d
    src/mesh.cpp
)
target_include_directories(euler2d PUBLIC include)

# Step 1 check: read a plot3D mesh and print it.
add_executable(test_mesh tests/test_mesh.cpp)
target_link_libraries(test_mesh PRIVATE euler2d)
