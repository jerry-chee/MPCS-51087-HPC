/* Jerr Chee
 * MPCS 51087 hw2
 * dynamic_decomp.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include "complex.h"
#include "decomp.h"
#include "mpi.h"
#include <limits.h>


/* initialize reArr */
void initArr(reArr* arr, int initial_size)
{
	arr->array = (int*) malloc(initial_size * sizeof(int));
	assert(arr->array != NULL);
	arr->used = 0;
	arr->size = initial_size;
}

/* insert element, resize of reach size limits */
void insertArr(reArr* arr, int x)
{
	//printf("max_int: %d\n", INT_MAX);
	
	if (arr->used == arr->size)
	{
		arr->size *= 2;
		arr->array = (int*) realloc(arr->array, arr->size * sizeof(int));
	}
	//int i = arr->used;
	//if (arr->array != NULL) printf("not NULL\n");
	arr->array[arr->used++] = x;
	//arr->used = (i+1);
}

/* free the reArr */
void freeArr(reArr* arr)
{
	free(arr->array);
	arr->array = NULL;
	arr->used = arr->size = 0;
}

/* prints reArr out to file, every n elements per line */
void printArr(FILE* fp, reArr* arr, int n)
{
	for (int i=0; i<(arr->used); i++)
	{
		if (i % n == 0) fprintf(fp, "\n");
		fprintf(fp, "%d ", arr->array[i]);
	}
}

/* a job, which the master assigns to a slave at a time */
void dynamic_job(int mype, reArr* data, double step, 
		int x_start, int x_end,
		int y_start, int y_end)
{
	double x_coord;
	double y_coord;

	for (int y = y_start; y < y_end; y++)
	{
		for (int x = x_start; x < x_end; x++)
		{
			x_coord = x_min + x*step;
			y_coord = y_min + y*step;

			mandelbrot_save(data, x_coord, y_coord);
			//count = mandelbrot_count(x_coord, y_coord);
			//fprintf(fp, "%lf,%lf,%d\n", x_coord, y_coord, count);
		}
	}	
}

void dynamic_decomp(int n, int k, int nprocs, int mype, reArr* data, int timing)
{
	/* MPI variables */
	MPI_Status status;
	int stat;
	int msg = 1;
	int requestor;
	int exit = -10; // can't be confused with any y values
													// outside the window

	assert(n % k == 0); // divisiblity requirment
	double step = 4 / ((double) n);
	int x_start = 0;
	int x_end = n;

	/* master-slave algorithm */	
	if (mype == 0) // master
	{
		int y_start = 0;
		int y_end = n;

		/* send out all the jobs */
		for (int y = y_start; y < y_end; y += k)
		{
			/* Listen for request */
			stat = MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, 0, 
					MPI_COMM_WORLD, &status);
			assert(stat == MPI_SUCCESS);

			/* find out who the message came from */
			requestor = status.MPI_SOURCE;
			//printf("0proc received job request from %dproc\n", requestor);

			/* send them a job */
			stat = MPI_Send(&y, 1, MPI_INT, requestor, 0, MPI_COMM_WORLD);
			assert(stat == MPI_SUCCESS);
			//printf("0proc sent job y:%ld to %dproc\n", y, requestor);
		}

		/* tell all workers to exit function */
		for (int i = 0; i < nprocs-1; i++)
		{
			/* Listen for request */
			stat = MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, 0, 
					MPI_COMM_WORLD, &status);
			assert(stat == MPI_SUCCESS);

			/* find out who the message came from */
			requestor = status.MPI_SOURCE;
			//printf("0proc received job request from %dproc\n", requestor);

			/* tell them to exit the function*/
			stat = MPI_Send (&exit, 1, MPI_INT, requestor, 0, MPI_COMM_WORLD);
			assert(stat == MPI_SUCCESS);
			//printf("0proc sent exit status to %dproc\n", requestor);
		}
	} else // slaves
	{
		int y_job_start;
		int y_job_end;

		while (1)
		{
			/* Ask master for a job */
			stat = MPI_Send(&msg, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			assert(stat == MPI_SUCCESS);
			//printf("%dproc sent job requesto to 0proc\n", mype);

			/* Wait for response */
			stat = MPI_Recv(&y_job_start, 1, MPI_INT, 0, 0, 
					MPI_COMM_WORLD, &status);
			assert(stat == MPI_SUCCESS);
			//printf("%dproc received msg from 0proc\n", mype);

			/* check message */
			if (y_job_start == exit)
			{
				//printf("%dproc received exit status from 0proc\n", mype);
				break;
			}

			/* do a portion of the mandelbrot calculation */
			y_job_end = y_job_start + k;
			dynamic_job(mype, data, step, x_start, x_end, 
					y_job_start, y_job_end); 			
		
			/* print out to file if timing == 0 */
			if (timing == 0)
			{
				char* dyn = "../data/mandelbrot/dyanmic_decomposition";
				char buf[512];
				sprintf(buf, "%s%dproc_%dk_%dy_start_runtime.csv", 
						dyn, mype, k, y_job_start);
				FILE* fp = fopen(buf, "w");
				assert(fp != NULL);
				printArr(fp, data, n);
				fclose(fp);
			}

			/* free reArr, reinit */
			freeArr(data);
			initArr(data, (n*k));
		}
	}

}

