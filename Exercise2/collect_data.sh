#!/bin/bash

output_file="quick_sort_results.txt"

array_sizes=(1000 10000 100000 1000000 10000000)

thread_counts=(1 2 4 8 12 16 20 24)

# Number of repetitions for each run
N=20

# Clear the output file
> "$output_file"

# Run QuickSort with different array sizes and thread counts
for size in "${array_sizes[@]}"; do
    for threads in "${thread_counts[@]}"; do
        echo "Array Size: $size, Threads: $threads"
        for (( i=1; i<=N; i++ )); do
            echo "Repetition $i out of $N"
            mpirun QuickSort.x $size $threads >> "$output_file"
        done
    done
done

echo "Data collection complete. Results stored in $output_file"

