cmake_minimum_required(VERSION 3.10)
project(QuickSortMPIOpenMP VERSION 1.0)

# Check if we are on an Apple Silicon processor
if(APPLE AND CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
    # Set the C and C++ compilers to GCC installed by Homebrew
    set(CMAKE_C_COMPILER "/opt/homebrew/bin/gcc-13")
    set(CMAKE_CXX_COMPILER "/opt/homebrew/bin/g++-13")
endif()

# Require C11 standard
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED True)

# Find the MPI package
find_package(MPI REQUIRED)

# Find the OpenMP package
find_package(OpenMP REQUIRED)

# Add the executable
add_executable(QuickSort quicksort.c)

# Include MPI
target_include_directories(QuickSort PRIVATE ${MPI_INCLUDE_PATH})
target_link_libraries(QuickSort ${MPI_C_LIBRARIES})

# Set the compiler flags and link libraries
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native ${OpenMP_C_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
target_link_libraries(QuickSort OpenMP::OpenMP_C)
