
/* ────────────────────────────────────────────────────────────────────────── *
 │                                                                            │
 │ This file is part of the exercises for the Lectures on                     │
 │   "Foundations of High Performance Computing"                              │
 │ given at                                                                   │
 │   Master in HPC and                                                        │
 │   Master in Data Science and Scientific Computing                          │
 │ @ SISSA, ICTP and University of Trieste                                    │
 │                                                                            │
 │ contact: luca.tornatore@inaf.it                                            │
 │                                                                            │
 │     This is free software; you can redistribute it and/or modify           │
 │     it under the terms of the GNU General Public License as published by   │
 │     the Free Software Foundation; either version 3 of the License, or      │
 │     (at your option) any later version.                                    │
 │     This code is distributed in the hope that it will be useful,           │
 │     but WITHOUT ANY WARRANTY; without even the implied warranty of         │
 │     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          │
 │     GNU General Public License for more details.                           │
 │                                                                            │
 │     You should have received a copy of the GNU General Public License      │
 │     along with this program.  If not, see <http://www.gnu.org/licenses/>   │
 │                                                                            │
 * ────────────────────────────────────────────────────────────────────────── */


#if defined(__STDC__)
#  if (__STDC_VERSION__ >= 199901L)
#     define _XOPEN_SOURCE 700
#  endif
#endif
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <mpi.h>



// ================================================================
//  MACROS and DATATYPES
// ================================================================


#if defined(_OPENMP)

// measure the wall-clock time
#define CPU_TIME (clock_gettime( CLOCK_REALTIME, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)

// measure the cpu thread time
#define CPU_TIME_th (clock_gettime( CLOCK_THREAD_CPUTIME_ID, &myts ), (double)myts.tv_sec +     \
                     (double)myts.tv_nsec * 1e-9)

#else

// measure ther cpu process time
#define CPU_TIME (clock_gettime( CLOCK_PROCESS_CPUTIME_ID, &ts ), (double)ts.tv_sec + \
                  (double)ts.tv_nsec * 1e-9)
#endif


#if defined(DEBUG)
#define VERBOSE
#endif

#if defined(VERBOSE)
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

//
// The data to be sorted consists in a structure, just to mimic that you
// may have to sort not basic types, whch may have some effects on the memory efficiency.
// The structure, defined below as data_t is just an array of DATA_SIZE double,
// where DATA_SIZE is defined here below.
// However, for the sake of simplicity, we assume that we sort with respect one of
// the double in each structure, the HOTst one, where HOT is defined below (my
// choice was to set HOT to 0, so the code will use the first double to sort the
// data).
//
// Note that you can change the data structure and the ordering as you like.

#if !defined(DATA_SIZE)
#define DATA_SIZE 4
#endif
#define HOT       0

// let's define the default amount of data
//
#if (!defined(DEBUG) || defined(_OPENMP))
#define N_dflt    10000
#else
#define N_dflt    10000
#endif

#define TASK_SIZE 500

// let's define the data_t type
//
typedef struct
{
  double data[DATA_SIZE];
} data_t;


// let's defined convenient macros for max and min between
// two data_t objects
//
#define MAX( a, b ) ( (a)->data[HOT] >(b)->data[HOT]? (a) : (b) )
#define MIN( a, b ) ( (a)->data[HOT] <(b)->data[HOT]? (a) : (b) )


// ================================================================
//  PROTOTYPES
// ================================================================

// let'ìs define the compare funciton that will be used by the
// sorting routine
//
typedef int (compare_t)(const void*, const void*);

// let's define the verifying function type, used to test the
// results
//
typedef int (verify_t)(data_t *, int, int, int);


// declare the functions
//
extern inline compare_t compare;        // the compare function
extern inline compare_t compare_ge;     // the compare for "greater or equal"
verify_t  verify_partitioning;          // verification functions
verify_t  verify_sorting;
verify_t  show_array;
void merge_chunks(data_t* data, int data_size, int num_chunks);

// declare the partitioning function
//
extern inline int partitioning( data_t *, int, int, compare_t );

// declare the sorting function
//
void quicksort( data_t *, int, int, compare_t);


// ================================================================
//  CODE
// ================================================================


int main ( int argc, char **argv )
{

    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int N = N_dflt;
    if (argc > 1) N = atoi(argv[1]);

    int n_threads = 1;
    if (argc > 2) n_threads = atoi(argv[2]);

    omp_set_num_threads(n_threads);

    data_t *data = NULL;
    data_t *local_data = NULL;

    if (rank == 0) {
        // Root process generates the full data set
        data = (data_t *) malloc(N * sizeof(data_t));
        #pragma omp parallel
        {
            int me = omp_get_thread_num();
            short int seed = time(NULL) % ((1 << sizeof(short int)) - 1);
            short int seeds[3] = {seed - me, seed + me, seed + me * 2};

            #pragma omp for
            for (int i = 0; i < N; i++)
                data[i].data[HOT] = erand48(seeds);
        }

    }


    // ---------------------------------------------
    //  process
    //
    struct timespec ts;
    double tstart = CPU_TIME;

    // Calculate the number of elements each process will handle
    int *sendcounts = malloc(size * sizeof(int));
    int *displs = malloc(size * sizeof(int));
    int sum = 0;
    for (int i = 0; i < size; ++i) {
        sendcounts[i] = N / size;
        if (i < N % size) ++sendcounts[i];
        sendcounts[i] *= DATA_SIZE;
        displs[i] = sum;
        sum += sendcounts[i];
    }

    // Allocate memory for local data on each process
    int local_N = sendcounts[rank] / DATA_SIZE;
    local_data = (data_t*)malloc(local_N * sizeof (data_t));


    // Distribute the data to all processes
    MPI_Scatterv(data, sendcounts, displs, MPI_DOUBLE,
                 local_data, sendcounts[rank], MPI_DOUBLE,
                 0, MPI_COMM_WORLD);

    #pragma omp parallel
    {
        #pragma omp single
        quicksort(local_data, 0, local_N, compare_ge);
        #pragma omp taskwait
    }

    // Gather the sorted data back to the root process
    MPI_Gatherv(local_data, sendcounts[rank], MPI_DOUBLE,
                data, sendcounts, displs, MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    double tend = CPU_TIME;

    if (rank == 0) {
        // Root process merges the sorted chunks
        //show_array(data, 0, N, 0);
        merge_chunks(data, N, size);
        //show_array(data, 0, N, 0);

        // Verify and output results
        if (verify_sorting(data, 0, N, 0))
            printf("%d\t%d\t%g sec\n", N, n_threads, tend - tstart);
        else
            printf("Array is not sorted correctly\n");

        free(data);
    }

    free(sendcounts);
    free(displs);
    free(local_data);
    MPI_Finalize();
    return 0;
}

void merge_chunks(data_t* data, int data_size, int num_chunks) {
    int regular_chunk_size = data_size / num_chunks;
    int remainder = data_size % num_chunks;
    int num_larger_chunks = remainder;
    int larger_chunk_size = regular_chunk_size + (remainder > 0 ? 1 : 0);

    data_t* temp = (data_t*)malloc(data_size * sizeof(data_t));
    int* chunk_start_indices = (int*)malloc(num_chunks * sizeof(int));
    int* chunk_sizes = (int*)malloc(num_chunks * sizeof(int));

    // Calculate start indices and sizes for each chunk
    int current_start = 0;
    for (int i = 0; i < num_chunks; i++) {
        chunk_start_indices[i] = current_start;
        chunk_sizes[i] = (i < num_larger_chunks) ? larger_chunk_size : regular_chunk_size;
        current_start += chunk_sizes[i];
    }

    int i, j, k;
    for (int step = 1; step < num_chunks; step *= 2) {
        for (int chunk = 0; chunk < num_chunks; chunk += 2 * step) {
            if (chunk + step >= num_chunks) {
                break; // No more pairs of chunks to merge
            }

            int chunk1_start = chunk_start_indices[chunk];
            int chunk1_end = chunk1_start + chunk_sizes[chunk];
            int chunk2_start = chunk_start_indices[chunk + step];
            int chunk2_end = (chunk + 2 * step < num_chunks) ?
                             chunk_start_indices[chunk + 2 * step] : data_size;

            i = chunk1_start;
            j = chunk2_start;
            k = chunk1_start;

            while (i < chunk1_end && j < chunk2_end) {
                if (compare(&data[i], &data[j]) <= 0) {
                    temp[k++] = data[i++];
                } else {
                    temp[k++] = data[j++];
                }
            }

            // Copy remaining elements from the first chunk
            while (i < chunk1_end) {
                temp[k++] = data[i++];
            }

            // Copy remaining elements from the second chunk
            while (j < chunk2_end) {
                temp[k++] = data[j++];
            }

            // Update start indices for merged chunks
            chunk_start_indices[chunk + step] = chunk1_start;
            chunk_sizes[chunk] = k - chunk1_start;
        }

        // Copy the merged data back to the original array
        memcpy(data, temp, data_size * sizeof(data_t));
    }

    free(temp);
    free(chunk_start_indices);
    free(chunk_sizes);
}



#define SWAP(A,B,SIZE) do {int sz = (SIZE); char *a = (A); char *b = (B); \
    do { char _temp = *a;*a++ = *b;*b++ = _temp;} while (--sz);} while (0)

inline int partitioning( data_t *data, int start, int end, compare_t cmp_ge )
{
    // Pick up the median of [start], [mid], and [end-1] as pivot
    int mid = start + (end - start) / 2;
    data_t *median;
    if (cmp_ge(&data[start], &data[mid])) {
        median = MAX(&data[start], MIN(&data[mid], &data[end-1]));
    } else {
        median = MAX(&data[mid], MIN(&data[start], &data[end-1]));
    }

    // Swap median with the last element
    SWAP((void*)median, (void*)&data[end-1], sizeof(data_t));
    void *pivot = (void*)&data[--end];

    // partition around the pivot element
    int pointbreak = end-1;
    for ( int i = start; i <= pointbreak; i++ )
        if( cmp_ge( (void*)&data[i], pivot ) )
        {
	    while( (pointbreak > i) && cmp_ge( (void*)&data[pointbreak], pivot ) ) pointbreak--;
	    if (pointbreak > i )
	        SWAP( (void*)&data[i], (void*)&data[pointbreak--], sizeof(data_t) );
        }
    pointbreak += !cmp_ge( (void*)&data[pointbreak], pivot ) ;
    SWAP( (void*)&data[pointbreak], pivot, sizeof(data_t) );

    return pointbreak;
}

void quicksort(data_t *data, int start, int end, compare_t cmp) {

    int size = end - start;

    if (size <= 1) return;
    int pivot = partitioning(data, start, end, cmp);

    #pragma omp task shared(data) if (size > TASK_SIZE)
    quicksort(data, start, pivot, cmp);// Sort the left part up to pivot

    #pragma omp task shared(data) if (size > TASK_SIZE)
    quicksort(data, pivot + 1, end, cmp);// Sort the right part after pivot

}

int verify_sorting( data_t *data, int start, int end, int not_used )
{
    int i = start;
    while( (++i < end) && (data[i].data[HOT] >= data[i-1].data[HOT]) );
    return ( i == end );
}

int verify_partitioning( data_t *data, int start, int end, int mid )
{
    int failure = 0;
    int fail = 0;

    for( int i = start; i < mid; i++ )
        if ( compare( (void*)&data[i], (void*)&data[mid] ) >= 0 )
            fail++;

    failure += fail;
    if ( fail )
    {
        printf("failure in first half\n");
        fail = 0;
    }

    for( int i = mid+1; i < end; i++ )
        if ( compare( (void*)&data[i], (void*)&data[mid] ) < 0 )
            fail++;

    failure += fail;
    if ( fail )
        printf("failure in second half\n");

    return failure;
}

int show_array( data_t *data, int start, int end, int not_used )
{
    for ( int i = start; i < end; i++ )
        printf( "%f ", data[i].data[HOT] );
    printf("\n");
    return 0;
}

inline int compare( const void *A, const void *B )
{
    data_t *a = (data_t*)A;
    data_t *b = (data_t*)B;

    double diff = a->data[HOT] - b->data[HOT];
    return ( (diff > 0) - (diff < 0) );
}

inline int compare_ge( const void *A, const void *B )
{
    data_t *a = (data_t*)A;
    data_t *b = (data_t*)B;

    return (a->data[HOT] >= b->data[HOT]);
}
