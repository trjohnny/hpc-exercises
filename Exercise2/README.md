# QuickSort Parallel Implementation

## Compilation Instructions

### For MacBook with Silicon Processor (Without MPI)
1. **Install OpenMP (if not already installed)**: Ensure that OpenMP is installed on your system. OpenMP can be installed using package managers like Homebrew. If you don't have Homebrew, install it first from [https://brew.sh/](https://brew.sh/).

2. **Install Compiler with OpenMP Support**: Install a compiler that supports OpenMP, such as GCC. You can install GCC using Homebrew:
   ```
   brew install gcc
   ```

3. **Generate Makefile and Compile**:
    - Navigate to the ``Exercise2`` directory
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
    - Navigate to the ``Exercise2`` directory
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