#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

__global__ void random_test(curandState *state, float *randArray)
{
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	curand_init(1, index, 0, &state[index]);

	randArray[index] = curand_uniform(&state[index]);
}

int main(int argc, char* argv[])
{
	int nblocks;
	int nthreads;

	nblocks = atoi(argv[1]);
	nthreads = atoi(argv[2]);

	curandState *d_state;
	cudaMalloc(&d_state, nthreads * nblocks);

	float *randArray;
	cudaMalloc(&randArray, nblocks*nthreads*sizeof(float));

	random_test<<<nthreads, nblocks>>>( d_state, randArray);

	float *rand_loc;
	rand_loc = (float *) malloc(nblocks * nthreads * sizeof(float));

	cudaMemcpy(rand_loc, randArray, nthreads*nblocks*sizeof(float),
			cudaMemcpyDeviceToHost);

	int i;
	for (i = 0; i < nblocks*nthreads; i++)
		printf("%f\n", rand_loc[i]);

	cudaFree(d_state);
	cudaFree(randArray);
	free(rand_loc);

}
