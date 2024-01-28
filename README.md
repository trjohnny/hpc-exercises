# QuickSort Parallel Implementation

## Overview
This QuickSort program is an implementation of the QuickSort algorithm with parallel processing capabilities. It is designed for educational purposes as part of the exercises for the Lectures on "Foundations of High Performance Computing".

## Compilation Instructions

### For MacBook with Silicon Processor (Without MPI)
1. **Install OpenMP (if not already installed)**: Ensure that OpenMP is installed on your system. OpenMP can be installed using package managers like Homebrew. If you don't have Homebrew, install it first from [https://brew.sh/](https://brew.sh/).

2. **Install Compiler with OpenMP Support**: Install a compiler that supports OpenMP, such as GCC. You can install GCC using Homebrew:
   ```
   brew install gcc
   ```

3. **Generate Makefile and Compile**:
    - Navigate to the directory containing the `QuickSort` program.
    - Create a build directory and navigate into it:
      ```
      mkdir build && cd build
      ```
    - Generate the Makefile using CMake (set the appropriate compilers path):
      ```
      cmake -DCMAKE_C_COMPILER=/opt/homebrew/bin/gcc-13 -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/g++-13 ..
      ```
    - Compile the program:
      ```
      make
      ```

4. **Run the Program**: After successful compilation, run the program using the format:
   ```
   ./QuickSort <length_of_array> <number_of_threads>
   ```

### For Systems with MPI Support (e.g., ORFEO)
1. **Install MPI (if not already installed)**: Ensure that MPI (such as OpenMPI or MPICH) is installed on your system. You can install it using a package manager or by downloading it from the respective websites.

2. **Allocate the desired number of nodes, tasks per node and cpus per task**:
   e.g. :
   ```
   salloc -p THIN --nodes 1 --tasks-per-node 24 --cpus-per-task 1 --time=00:30:00
   ```

3. **Load the MPI module**:
   ```
   module load openMPI/4.1.5/gnu/12.2.1
   ```

4. **Compile the Program with mpicc**:
    - Navigate to the directory containing the `QuickSort` program.
    - Use `mpicc` to compile the program:
      ```
      mpicc -fopenmp -o QuickSort.x quicksort.c
      ```

5. **Run the Program**: After successful compilation, run the program using `mpiexec` or `mpirun`, depending on your MPI installation, and specify the number of processes along with the program arguments:

   ```
   mpirun QuickSort.x <length_of_array> <number_of_threads>
   ```

### Notes
- Ensure that the correct GCC version is used when compiling on a MacBook with a Silicon Processor, as the default `clang` compiler might not support OpenMP.
- For MPI compilation, the `mpicc` command is typically included with your MPI installation. Ensure that `mpicc` is in your system's PATH.
- The number of threads specified when running the program should not exceed the number of available CPU cores for optimal performance.
## Running the Program
To run the program, use the following command format:
```
./QuickSort <length_of_array_to_be_sorted> <number_of_threads_to_spawn>
```
For example:
```
./QuickSort 10000 12
```

This command will sort an array of 10,000 elements using 12 threads.

## Design Choices

### 1. Use of Tasks in the Quicksort Implementation
Tasks are used in the QuickSort implementation instead of sections for the following reasons:

- **Dynamic Nature of Tasks**: Tasks in OpenMP are more dynamic and flexible compared to sections. They allow for the creation of a new unit of work (task) that can be executed by any thread in the team. This is particularly beneficial for recursive algorithms like QuickSort, where the number of tasks to be created is not known in advance and can vary depending on the data.
- **Load Balancing**: Tasks can be executed by any available thread, leading to better load balancing. In contrast, sections are static and each section is executed by a specific thread, which might lead to inefficient utilization of resources if the workload is not evenly distributed.

### 2. Merging Algorithm and Its Complexity
The merging algorithm used in the program is designed to merge sorted chunks of an array after they have been sorted individually. The complexity of the algorithm is as follows:

- **Algorithm Description**: The merge function sequentially takes two adjacent sorted chunks and merges them into a single sorted chunk. This process is repeated in a bottom-up manner until all chunks are merged into a single sorted array.
- **Complexity**: The time complexity of the merging process is O(n log k), where n is the total number of elements in the array and k is the number of chunks. Each merge operation is linear in the size of the chunks being merged, and there are log k levels of merging since the number of chunks is halved in each level.

## Important Notes
- Ensure that the MPI library is installed on your system if you intend to compile and run the program with MPI support.
- When running on a MacBook with a Silicon Processor, MPI may not be supported or optimal; thus, it is recommended to compile and run without MPI.
- The program dynamically allocates memory for the data to be sorted; ensure your system has sufficient memory for large data sizes.

