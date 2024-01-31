#!/bin/bash

# Change to the directory where your OSU Micro-Benchmarks are located
# Update this path according to your actual OSU Micro-Benchmark directory path
OMB_PATH="osu-micro-benchmarks-7.3"

# Navigate to the OSU Micro-Benchmark folder
cd "$OMB_PATH"

# Check if the benchmarks are already configured (e.g., Makefile exists)
if [ ! -f Makefile ]; then
    echo "Configuring the OSU Micro-Benchmarks..."
    ./configure CC=/opt/programs/openMPI/4.1.5/bin/mpicc CXX=/opt/programs/openMPI/4.1.5/bin/mpicc && echo "Configuration completed."
fi

# Build the benchmarks
echo "Building the OSU Micro-Benchmarks..."
make && make install && echo "Build completed."

# Navigate to the specific benchmark directory
cd c/mpi/collective/blocking

# Run the benchmarks
echo "Running osu_bcast benchmarks..."

echo "Baseline:"
mpirun osu_bcast -x 100 -i 1000 -f

echo "Scatter AllGather Ring:"
mpirun --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_bcast_algorithm 9 ./osu_bcast -x 100 -i 1000 -f

echo "Pipeline:"
mpirun --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_bcast_algorithm 3 ./osu_bcast -x 100 -i 1000 -f

echo "osu_bcast benchmark runs completed."


echo "Running osu_gather benchmarks..."

echo "Baseline:"
mpirun osu_gather -x 100 -i 1000 -f

echo "Binomial:"
mpirun --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_gather_algorithm 2 osu_gather -x 100 -i 1000 -f

echo "Basin Linear:"
mpirun --mca coll_tuned_use_dynamic_rules true --mca coll_tuned_gather_algorithm 1 osu_gather -x 100 -i 1000 -f

echo "osu_gather benchmark runs completed."
