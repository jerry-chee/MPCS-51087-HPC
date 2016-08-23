/* Jerry Chee
 * MPCS 51087 hw 2
 * static_decomp.c
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "complex.h"
#include "decomp.h"
#include "mpi.h"

/* returns iterations needed to break mandelbrot condition or iter_max */
int mandelbrot_count(double x_coord, double y_coord)
{
	/* keep track of iterations */                
	int count = 0;

	/* initialize complex structs */  
	struct Complex c = init(x_coord, y_coord);  
	struct Complex z = init(0,0);
	double tmp;
	//print_complex(c, "c");
	//print_complex(z, "z0"); 

	/* perform iterations, check if z exceeds condition */
	for (int i = 1; i <= iter_max; i++)
	{
		/* perform next iteration */
		tmp = z.x * z.x - z.y * z.y + c.x;
		z.y = 2 * z.x * z.y + c.y;
		z.x = tmp;
		
		//z_tmp = add(sq(z), c);
		//z = z_tmp; 
		//print_complex(z, "z");

		if (dist_sq(z) > r_max_sq)
		{   
			break;
		}   

		/* iterate count */   
		count++;
	}

	return count;
}

/* prints mandelbrot_count to file pointer */
void mandelbrot_save(reArr* data, double x_coord, double y_coord)
{
	//printf("calling mandelbrot_count\n");
	int count = mandelbrot_count(x_coord, y_coord);

	//printf("calling insertArr\n");
	insertArr(data, count);		
	//printf("insertArr success\n");
}

/* runs mandelbrot decomp functions, aggregates times, prints to file */
void time_master(void (*f)(int, int, int, int, reArr*, int), char* type, 
		int n, int k, int nprocs, int mype, int timing)
{
	printf("%dproc running decomp fn\n", mype);

	/* initialize struct */
	reArr data;
	initArr(&data, (n*k));
	printf("%dproc initialized reArr\n", mype);

	/* 0proc is the timer */
	clock_t start, end;
	double cpu_time_used;

	/* run function on each proc */
	printf("%dproc running function\n",mype);
	if ((timing == 1) && (mype == 0)) start = clock();
	f(n, k, nprocs, mype, &data, timing);	

	printf("%dproc waiting at MPI_Barrier()\n", mype);
	/* wait for all procs to finish */
	MPI_Barrier(MPI_COMM_WORLD);

	/* print out times */
	if ((timing == 1) && (mype == 0)) // master 
	{
		/* timing end */
		end = clock();
		cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;	

		printf("0proc printing out to file, %lf\n", cpu_time_used);
		/* file out for time */
		char* path_time = "../data/mandelbrot/times/";
		char path_buffer[512];
		sprintf(path_buffer, "%s%dn_%dk_time.csv", path_time, n, k);
		FILE* fp = fopen(path_buffer, "a");
		assert(fp != NULL);
		fprintf(fp, "%s,%d,%d,%d,%lf\n", type, n, nprocs, k, cpu_time_used);
		fclose(fp);
		printf("0proc finished\n");
	}

	/* wait for 0proc to finish printing times */
	MPI_Barrier(MPI_COMM_WORLD);

	/* print results out to one file, static decomp only */
	if ( (timing == 0) && (strcmp(type, "static") == 0) )
	{
		/* open file pointer */
		char* path_data = "../data/mandelbrot/";
		char buffer_data[512];
		sprintf(buffer_data, "%s%dprocs_data.csv", path_data, nprocs);
		FILE* fp = NULL;

		/* 0proc opens file pointer first to clear*/
		if (mype == 0)
		{
			fp = fopen(buffer_data, "w");
			assert(fp != NULL);
		}

		/* wait for 0proc to finish */
		MPI_Barrier(MPI_COMM_WORLD);

		if (mype != 0)
		{
			fp = fopen(buffer_data, "a");
			assert(fp != NULL);
		}

		/* wait for rest of procs to finish */
		MPI_Barrier(MPI_COMM_WORLD);

		for (int i=0; i<nprocs; i++)
		{
			if (mype == i) 
			{
				printArr(fp, &data, n);
				printf("%dproc printed to file, used:%d\n", mype, data.used);
			}
			/* wait for this proc to finish */
			MPI_Barrier(MPI_COMM_WORLD);
		}

		/* all procs close  file pointer */
		fclose(fp);
	}

	/* free reArr */
	freeArr(&data);
}

/* compute mandelbrot set with static data decomposition */
void static_decomp(int n, int k, int nprocs, int mype, reArr* data, int timing)
{

	/* divide the window W into nprocs regions 
	 * all procs get xrange (-2,2), but yrange (i,j).
	 * assigning ranges by indices.
	 */
	double step = 4 / ((double) n);
	/* NOTE: n >= must be guaranteed in order for 
	 * code to work
	 * decomposition by taking the floor of n/procs 
	 * for all block sizes except the last, for which
	 * the n % nprocs is added.*/
	int block = 0; 
	int block_floor = n / nprocs;
	int block_mod = n % nprocs;

	if ((mype == nprocs-1) && (n % nprocs != 0))
	{
		block = block_mod + block_floor;
	} else if ((mype == nprocs-1) && (n % nprocs == 0))
	{
		block = block_floor;
	}	else 
	{
		block = block_floor; /* -1 */
	}
	int y_start = mype * block_floor;
	int y_end = y_start + block;
	int x_start = 0;
	int x_end = n;

	/* testing index assignment */
	printf("mype: %d \nstep: %lf \t block: %d \ny_start: %d \t y_end: %d \nx_start: %d \t x_end: %d \n\n", mype, step, block, y_start, y_end, x_start, x_end); 
	fflush(stdout); 

	double x_coord;
	double y_coord;

	for (int y = y_start; y < y_end; y++)
	{
		for (int x = x_start; x < x_end; x++)
		{
			x_coord = x_min + x*step;
			y_coord = y_min + y*step;		

			/* here we calculate the number of iterations it takes to
			 * escape mandelbrot set */
			//printf("%dproc calling mandelbrot_save\n", mype);
			mandelbrot_save(data, x_coord, y_coord);
		}
	}
}

