#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include "mpi.h"

int main(int argc, char** argv)
{
	/* max number iterations */
	int max = 10000000;
	int nreps = 1000;
	//int max = 1000000;
	//int step = 100000;

	/* MPI Init variables */
	int nprocs;
	int mype;
	int stat;
	MPI_Status status;

	MPI_Init(&argc, &argv);

	/* return number of procs */
	stat = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	assert(stat == MPI_SUCCESS);

	/* return interger proc id */
	stat = MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	assert(stat == MPI_SUCCESS);

	/* wait for everyone to get here */
	MPI_Barrier(MPI_COMM_WORLD);

	/* ======================================================================= */
	if (mype == 0) // master
	{
		printf("master starting\n");

		/* file pointer */
		char* path = "../data/latency_bandwidth/lt_bd_times.csv";
		FILE* fp = fopen(path, "w");
		assert(fp != NULL);
		fprintf(fp, "bytes,time,bytes_per_sec\n");
		long int bytes = 0;

		/* timing variables */
		double start, end, time;

		printf("master mallocing\n");
		/* malloc array which will be sent and received in */
		double* m_send = (double*) malloc(max * sizeof(double));
		double* m_recv = (double*) malloc(max * sizeof(double));

		printf("master entering computation loop\n");
		/* timing loop */
		double t_sum = 0;
		double t_avg = 0;
		long unsigned int i;
		for (i = 1; i < max; i*=2)
		{
			for (int j=0; j < nreps; j++)
			{
				/* start clock */
				start = MPI_Wtime();

				/* send message */
				stat = MPI_Send(m_send, i, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
				assert(stat == MPI_SUCCESS);

				/* wait for message to get sent back */
				stat = MPI_Recv(m_recv, i, MPI_DOUBLE, 1, 0, 
						MPI_COMM_WORLD, &status);
				assert(stat == MPI_SUCCESS);

				/* end clock */
				end = MPI_Wtime();
				time = end - start;

				/* add to t_sum */
				t_sum += time;
			}
			
			/* take the average */
			t_avg = t_sum / nreps;		

			bytes = i*sizeof(double);
			/* print to file */
			fprintf(fp, "%ld,%lf,%lf\n", bytes, t_avg, (bytes/t_avg));
			fflush(stdout);

			/* zero out sum and avg */
			t_sum = 0;
			t_avg = 0;
		}

		printf("master exited computation loop\n");
		/* close file pointer */
		fclose(fp);	

		/* free arrays */
		free(m_send);
		free(m_recv);

	}else if (mype == 1) // slave
	{
		printf("slave running\n");

		/* malloc array recved in and sent from */
		double* s_arr = (double*) malloc(max * sizeof(double));

		printf("slave entering computation loop\n");
		long unsigned int i;
		for (i = 1; i < max; i*=2)
		{
			for (int j=0; j < nreps; j++)
			{
				/* wait for message from 0 */
				stat = MPI_Recv(s_arr, i, MPI_DOUBLE, 0, 0, 
						MPI_COMM_WORLD, &status);
				assert(stat == MPI_SUCCESS);

				/* send back message */
				stat = MPI_Send(s_arr, i, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
				assert(stat == MPI_SUCCESS);
			}
		}

		printf("slave exited computation loop\n");
		/* free array */
		free(s_arr);

	} else // do nothing
	{
	}

	MPI_Finalize();
}
