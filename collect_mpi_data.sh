#!/bin/bash

output_file="quick_sort_MPI_results.txt"

array_sizes=(1000 10000 100000 1000000 10000000)

core_counts=(1 2 4 8 12 16 20 24)

# Number of repetitions for each run
N=20

# Clear the output file
> "$output_file"

# Run QuickSort with different array sizes and thread counts
for cores in "${core_counts[@]}"; do
    
    for size in "${array_sizes[@]}"; do
        echo "Array Size: $size, Cores: $cores"
        
        for (( i=1; i<=N; i++ )); do
            echo "Repetition $i out of $N"
            mpirun -np $cores QuickSort.x $size 1 >> "$output_file"
        done
    done

done

echo "Data collection complete. Results stored in $output_file"

