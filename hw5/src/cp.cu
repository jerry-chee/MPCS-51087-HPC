#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#define NDIM 5

__global__ void init_rand(curandState *state)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(1, index, 0, &state[index]);
}

/* I believe curand_uniform gen between 0 and 1 */
/* this code isn't code right, yet */
__global__ void gen_rand(curandState *state, double *randArray)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	
	int i;
	for (i = index; i < 5*n; i += nblocks*nthreads)
		randArray[i] = curand_uniform(&state[index]);
}

/* implement local array to save */
__global__ void eval_function(double *randArray, double *int_f, int *n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int start = 5 * index;

	double f = exp(-randArray[start]*randArray[start] 
			-randArray[start+1]*randArray[start+1] 
			-randArray[start+2]*randArray[start+2]  
			-randArray[start+3]*randArray[start+3] 
			-randArray[start+4]*randArray[start+4]);
	int_f[index] = f / (*n);
}

__global__ void reduction(double *int_f, double *int_f_master, int *n)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	if (index == 0)
	{
		*int_f_master = 0.0;
		int i;
		for (i = 0; i < *n; i++)
			*int_f_master += int_f[i];
	}
}

int main(int argc, char **argv)
{
	int n, *d_n;
	double t0, t1;
	double vol=1.0;

	int nBlocks  = 2;
	int nThreads = 64;

	n = atoi(argv[1]);

	cudaMalloc(&d_n, sizeof(int));
	cudaMemcpy(d_n, &n, sizeof(int), cudaMemcpyHostToDevice);

	curandState *d_state;
	cudaMalloc(&d_state, nThreads * nBlocks);

	double *randArray;
	cudaMalloc(&randArray, 5 * n * sizeof(double));

	double *int_f;
	cudaMalloc(&int_f, n * sizeof(double));

	double *int_f_master;
	cudaMalloc(&int_f_master, sizeof(double));

	init_rand    <<<nThreads, nBlocks>>>( d_state );
	gen_rand     <<<nThreads, nBlocks>>>( d_state, randArray);
	eval_function<<<nThreads, nBlocks>>>( randArray, int_f, d_n);
	reduction    <<<nThreads, nBlocks>>>( int_f, int_f_master, d_n);

	double integral;
	cudaMemcpy(&integral, int_f_master, sizeof(double), 
		cudaMemcpyDeviceToHost);

	//int i;
	//for (i = 0; i < NDIM; i++)
	//	vol *= (b-a);

	printf("%lf\n", integral);

	cudaFree(d_state);
	cudaFree(randArray);
	cudaFree(int_f);
	cudaFree(int_f_master);	
}
