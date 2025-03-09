/**
 * Copyright 2025 Jim Haslett
 *
 * This work part of a university assignment.  If you are taking the course
 * this work was assigned for, do the right thing, and solve the assignment
 * yourself!
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "mpi.h"

#define OUTFILENAME "result.txt"
#define REPORTFILENAME "report.txt"


/**
 * This is the CompareSplit function from the textbook p.250, it is horribly unreadable code
 * 
 * This function takes in two arrays, compares the combined elements of both, and keeps either
 * the smaller half or the larger half of the combined set of the two arrays
 * 
 * Parameters:
 *      nlocal      - The number of elements in the local array
 *      elmnts      - The array of elements local to this process
 *      relmnts     - The array fo elements received from the remote process
 *      wspace      - Array space used for temporary storage
 *      keepsmall   - Used as a boolean to tell the function if it should keep the smaller or larger elements
 */
void CompareSplit(int nlocal, int *elmnts, int *relmnts, int *wspace, int keepsmall)
{
    int i, j, k;
    
    for (i=0; i<nlocal; i++)
        wspace[i] = elmnts[i]; /* Copy the elmnts array into the wspace array */
    
    if (keepsmall) { /* Keep the nlocal smaller elements */
        for (i=j=k=0; k<nlocal; k++) {
            if (j == nlocal || (i < nlocal && wspace[i] < relmnts[j]))
                elmnts[k] = wspace[i++];
            else
                elmnts[k] = relmnts[j++];
        }
    }
    else { /* Keep the nlocal larger elements */
        for (i=k=nlocal-1, j=nlocal-1; k>=0; k--) {
            // if (j == 0 || (i >= 0 && wspace[i] >= relmnts[j]))  /* This is a typo in the textbook! should be j<0 */ 
            if (j < 0 || (i >= 0 && wspace[i] >= relmnts[j]))
                elmnts[k] = wspace[i--];
            else
                elmnts[k] = relmnts[j--];
        }
    }
}

/**
 * This function generates the initial array of values 1 to nlocal * number_of_processes.  It then does a
 * random shuffle using the Fischer-Yates algorigthm.  This is the initial random array.  This array is then
 * printed to the screen, and scattered to the processes.
 * 
 * Parameters:
 *      rank                - MPI rank of the current process
 *      nlocal              - number of elements in the local array
 *      local_array         - pointer to the local array storing the n/p chunk assigned to the current process
 *      number_of_processes - number of MPI processes running
 */
void set_up_initial_array(int rank, int nlocal, int *local_array, int number_of_processes){
    int *initial_array;  /* pointer to starting array, only allocated on process 0 */
    if (rank == 0){
        /* This process should generate and transmit the initial array */

        int init_array_size = nlocal * number_of_processes;  /* size of initial array */
        initial_array = (int *)malloc(init_array_size * sizeof(int)); /* allocate memory for initial array */
        srand(time(NULL));  /* initial random seed from clock */

        /* initialize starting array */
        for (int i = 0; i < init_array_size; ++i){
            initial_array[i] = i+1;
        }

        /* random shuffle starting array (Fisher-Yates) */
        int j, t;
        for (int i = init_array_size-1; i >= 0; --i){
            j = rand() % (i+1);
            t = initial_array[i];
            initial_array[i] = initial_array[j];
            initial_array[j] = t;
        }
        
        /* Print unsorted starting array */
        printf("Starting Array:\n");
        for (int i = 0; i < init_array_size; ++i){
            printf("%d ", initial_array[i]);
        }
        printf("\n");


        /* Distribute start_array to each process */
        MPI_Scatter(initial_array, nlocal, MPI_INT, local_array, nlocal, MPI_INT, 0, MPI_COMM_WORLD);

        free(initial_array); /* free initial array memory */
    } else {
        /* Distribute start_array to each process */
        MPI_Scatter(initial_array, nlocal, MPI_INT, local_array, nlocal, MPI_INT, 0, MPI_COMM_WORLD);        
    }
}

/**
 * This function gathers the results from all processes back into process 0, which prints the result.
 * 
 * Parameters:
 *      rank                - MPI rank of the current process
 *      nlocal              - number of elements in the local array
 *      local_array         - pointer to the local array storing the n/p chunk assigned to the current process
 *      number_of_processes - number of MPI processes running
 */
void collect_results(int rank, int nlocal, int *local_array, int number_of_processes){
    int *result;  /* result array, only allocated on process 0 */
    if (rank == 0){
        /* This process should gather and report the results */
        int init_array_size = nlocal * number_of_processes;  /* size of initial array */
        result = (int *)malloc(init_array_size * sizeof(int)); /* allocate memory for initial array */

        /* Gather result from each process */
        MPI_Gather(local_array, nlocal, MPI_INT, result, nlocal, MPI_INT, 0, MPI_COMM_WORLD);
        
        /* Print Results */
        printf("Sorted Array:\n");
        for (int i = 0; i < init_array_size; ++i){
            printf("%d ", result[i]);
        }
        printf("\n");

        free(result);  /* free allocated memory */
    } else {
        /* Gather result from each process */
        MPI_Gather(local_array, nlocal, MPI_INT, result, nlocal, MPI_INT, 0, MPI_COMM_WORLD);
    }
}

/**
 * comparer function used for build in quicksort.  Simple integer comparsion.
 */
int comparer(const void *e1, const void *e2)
{
    return (*((int *)e1) - *((int *)e2));
}

/**
 * This is the actual bitonic sort function.  Prior to this running, the initial array has been generated and
 * distributed to each process.
 * 
 * Parameters:
 *      rank                - MPI rank of the current process
 *      nlocal              - number of elements in the local array
 *      local_array         - pointer to the local array storing the n/p chunk assigned to the current process
 *      number_of_processes - number of MPI processes running
 */
void parallel_bitonic_sort(int rank, int nlocal, int *local_array, int number_of_processes){
    MPI_Status status;  /* MPI_Recv status variable */
    int *recieved_array;  /* Received array for each process */
    int *working_space;  /* Working space for sort operation */

    /** Allocate memory for working arrays */
    recieved_array = (int *)malloc(nlocal * sizeof(int));
    working_space = (int *)malloc(nlocal * sizeof(int));
    /** */

    /* Step 1:  Sort local elements */
    qsort(local_array, nlocal, sizeof(int), comparer);

    /* Step 2: CompareSplit with neighbors */
    int d = (int)log2(number_of_processes); /* dimension of they hypercube */
    int parnter_rank;  /* used to store rank of parnter this process will communicate with at each step */
    int partner_mask;  /* bitmask used at each step in the calculation of partner_rank */
    int increasing_mask;  /* bitmask used to determine if this step is an increasing or decreasing comparator operation */
    int take_lower_value;  /* boolean calculated based on the parnter and the increasing_mask to tell if this process   
                                should keep the lower or higher values in the comparesplit*/
    for (int i = 0; i < d; ++i){
        partner_mask = (int)pow(2, i); /* Initial comm parter partner_mask for this round */
        increasing_mask = partner_mask << 1;  /* Calculate mask for increasing or decreasing comparator operation */
        for (int j = i; j >= 0; --j){
            parnter_rank = rank ^ partner_mask;  /* Calculate communication partner this round */
            MPI_Sendrecv(local_array, nlocal, MPI_INT, parnter_rank, 1, recieved_array, nlocal, MPI_INT, parnter_rank, 1, MPI_COMM_WORLD, &status);
            
            /**
             * This is the magic sauce that determines for each process at each step if it should take the lower
             * or higher values in the CompareSpit function.d
             * 
             * The way this works:
             *      1. increasing_mask is calculated for each i value.  It is one bit more significant than the partner
             *          mask when i = j.
             *      2. take_lower_value is true if the current process's rank is lower than the partner's rank and
             *          the bit of the current rank value corresponding to the increasing_rank bit value is NOT set.  
             *          OR, the opposite, if the current process rank is higher than the partner's rank the bit of 
             *          current rank value corresponding ot the increasing_rank bit value IS set.
             * 
             * if take_lower_value is true, this process keeps the lower half of the CompareSplit otherwise, it keeps
             * the higher value half.
             */
            take_lower_value = ((rank & increasing_mask) > 0) != (rank < status.MPI_SOURCE);
            
            CompareSplit(nlocal, local_array, recieved_array, working_space, take_lower_value);
            partner_mask = partner_mask >> 1;  /* rightshift partner_mask to get next comm partner */
        }
    }

    /** Free allocated memory */
    free(recieved_array);
    free(working_space);
    /** */
}


/**
 * Main function.
 */
int main(int argc,char* argv[]) {

    /** Set up MPI */
    int rank;  /* MPI rank */
    MPI_Init(&argc, &argv);  /* Initialize MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  /* Fetch rank */
    int number_of_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &number_of_processes);  /* Fetch number of processes */
    /** */
    
    /** Set up local storage, define local array size */
    int nlocal = number_of_processes;  /* use number fo processes as the number of elements per process, this ensures n^2 array */
    int *local_array;  /* Local array for each process */
    /** */

    /** Allocate memory for arrays */
    local_array = (int *)malloc(nlocal * sizeof(int));
    /** */


    /* Initialize starting array, and distribute */
    set_up_initial_array(rank, nlocal, local_array, number_of_processes);  
    
    
    /* Perform Parallel Bitonic Sort */
    parallel_bitonic_sort(rank, nlocal, local_array, number_of_processes);  


    /* Initialize starting array, and distribute */
    collect_results(rank, nlocal, local_array, number_of_processes);  


    /** Free allocated memory */
    free(local_array);
    /** */

    /** Shut down MPI and exit */
    MPI_Finalize();
    exit(0);
    /** */
}