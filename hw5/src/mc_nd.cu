#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>

#define NDIM 5
#define THREADS_PER_BLOCK 32

float uran(float a, float b, unsigned int *seed);
__host__ __device__ float f(float x[]);

__global__ void mc_integration(float *rand_device, float *int_f, int *n_device)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	int start = NDIM * index;
	__shared__ float tmp[THREADS_PER_BLOCK];

	tmp[threadIdx.x] = f(&rand_device[start]) / *n_device;

	__syncthreads();

	if (0 == threadIdx.x)
	{
		float int_f_block = 0.0;
		for (int i = 0; i < blockDim.x; i++)
			int_f_block += tmp[i];
			
		atomicAdd(int_f, int_f_block);
	}
}

int main(int argc, char **argv)
{
	// Host variables
	int n;
	float a, b;
	//float vol=1.0;
	cudaEvent_t start, stop;
	float milliseconds = 0;
	float integral;
	unsigned int seed = 1;
	float *rand_host;

	a = 0,0; b = 1.0;
	n = atoi(argv[1]); // number of points to evaluate at

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	rand_host = (float *) malloc(NDIM * n * sizeof(float));	

	// Device variables
	float *rand_device;
	int *n_device;
	float *int_f;
	float zero = 0.0;

	cudaMalloc(&rand_device, NDIM * n * sizeof(float));
	
	for (int i = 0; i < NDIM*n; i++)
	{
		rand_host[i] = uran(a, b, &seed);
		//printf("%f\n", rand_host[i]);
	}
	//putchar('\n');

	cudaMemcpy(rand_device, rand_host, NDIM * n * sizeof(float), 
			cudaMemcpyHostToDevice);

	cudaMalloc(&n_device, sizeof(int));
	cudaMemcpy(n_device, &n, sizeof(int), cudaMemcpyHostToDevice);

	cudaMalloc(&int_f, sizeof(double));
	cudaMemcpy(int_f, &zero, sizeof(float), cudaMemcpyHostToDevice);

	// run on GPU
	assert(n % THREADS_PER_BLOCK == 0);
	cudaEventRecord(start);
	mc_integration<<< n/THREADS_PER_BLOCK, THREADS_PER_BLOCK >>>( 
			rand_device, int_f, n_device );
	cudaEventRecord(stop);

	//printf("n/NblOCKS: %d\n", n/NBLOCKS);
	//for (int i = 0; i <n; i++)
	//	printf("integral: %f\n", (f(&rand_host[NDIM*i])/n)); 


	// bring back integral
	cudaMemcpy(&integral, int_f, sizeof(float), cudaMemcpyDeviceToHost);

	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&milliseconds, start, stop);
	
	printf("mc cuda integral: %f, n:%d, time:%f\n", integral, n, milliseconds/1000);

	// free
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	free(rand_host);
	cudaFree(rand_device);
	cudaFree(n_device);
	cudaFree(int_f);	

}

float uran(float a, float b, unsigned int *seed){
	  return rand_r(seed) / (RAND_MAX + 1.0) * (b - a) + a;
}

__host__ __device__ float f(float x[]){ 
	  return exp(-x[0]*x[0] -x[1]*x[1] -x[2]*x[2]  -x[3]*x[3] -x[4]*x[4]); 
}

