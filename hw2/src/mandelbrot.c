/* Jerry Chee
 * MPCS 51087 hw 2
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>
#include "complex.h"
#include "decomp.h"
#include "mpi.h"

/* Global var for criterion for mandelbrot divergenc 
 *  * and maximum number of iterations*/
int r_max = 2;
int r_max_sq = 4;
int iter_max = 1000; 
int x_min = -2; 
int y_min = -2; 

/* I am looking in a fixed W = (-2,2)x(-2,2) in R^2 window */
int main(int argc, char** argv)
{
	/* MPI initialize */
	int nprocs; 
	int mype; 
	int stat;

	MPI_Init(&argc, &argv); /* return number of procs */
	stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	assert(stat == MPI_SUCCESS);

	/* return interger proc id */
	stat = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	assert(stat == MPI_SUCCESS);

	/*
	printf("mype:%d   nprocs:%d\n",mype,nprocs);
	fflush(stdout);
	*/

	/* ======================================================================= */

	/* read in number of times to partion window along an axis */ 
	int n = atoi(argv[1]);	
	
	/* read in dynamic chunk size */
	int k = atoi(argv[2]);

	/* timing argument */
	int timing = atoi(argv[3]);

	time_master(static_decomp, "static", n, k, nprocs, mype, timing);

	//time_master(dynamic_decomp, "dynamic",n, k, nprocs, mype, timing);

	/* ======================================================================= */
	MPI_Finalize();
}
