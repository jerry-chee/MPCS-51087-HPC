#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <cuda.h>

#define NDIM 5
#define THREADS_PER_BLOCK 243 // make sure a power of 3

__host__ __device__ float f( float x[] )
{
	return exp(-x[0]*x[0] -x[1]*x[1] -x[2]*x[2]  -x[3]*x[3] -x[4]*x[4]);
}

__global__ void simpson_integration(int *mesh_device, float *x_device,
		float *w_device, float *delta_device, float *h_device, int *n_device, 
		float *integral_device)
{
	float a = 0.0;
	float b = 1.0;
	long long int size = (*n_device)*(*n_device)*
		(*n_device)*(*n_device)*(*n_device)
		* 3*3*3*3*3;
	int index = blockIdx.x * blockDim.x + threadIdx.x;
	__shared__ float tmp[THREADS_PER_BLOCK];

	assert(blockDim.x * gridDim.x == size);
	for (long long int i = index; i < 5*size; i+=size)
	{
		x_device[i] = a + mesh_device[2*i] * (*delta_device) + 
			mesh_device[2*i+1] * (*h_device);	
	}

	__syncthreads();

	int start = 10*index;
	tmp[threadIdx.x] = w_device[mesh_device[start+1]] * 
		w_device[mesh_device[start+3]] * w_device[mesh_device[start+5]] * 
		w_device[mesh_device[start+7]] * w_device[mesh_device[start+9]] *
		f(&x_device[5*index]);
	
	__syncthreads();

	if (0 == threadIdx.x)
	{
		float int_block = 0.0;
		assert(blockDim.x == THREADS_PER_BLOCK);
		for (int i = 0; i < blockDim.x; i++)
			int_block += tmp[i];

		atomicAdd(integral_device, int_block);
	}
}

int main(int argc, char **argv)
{
	// Host variables
	float integral_host, delta, h;
	int n;
	float a, b;
	int i,ii,j,jj,k,kk,l,ll,m,mm;
	cudaEvent_t start, stop;
	float milliseconds;
	int *mesh_host;
	float w_host[3];

	a = 0.0; b = 1.0;
	n = atoi(argv[1]);

	delta = (b-a)/n;
	h = delta / 2.0;
	w_host[0] = h/3; w_host[1] = 4.*w_host[0]; w_host[2] = w_host[0];

	milliseconds = 0.0;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);


	long long int size = n*n*n*n*n * 3*3*3*3*3; 
	mesh_host = (int *) malloc( 10 * size * sizeof(int));

	long long int index_host = 0;
	for (i = 0; i < n; i++){
		for(ii = 0; ii < 3; ii++){
			for (j = 0; j < n; j++){
				for(jj = 0; jj < 3; jj++){
					for (k = 0; k < n; k++){
						for(kk = 0; kk < 3; kk++){
							for (l = 0; l < n; l++){
								for(ll = 0; ll < 3; ll++){
									for (m = 0; m < n; m++){
										for(mm = 0; mm < 3; mm++){
											mesh_host[index_host] = i; index_host++;

											mesh_host[index_host] = ii; index_host++;

											mesh_host[index_host] = j; index_host++;

											mesh_host[index_host] = jj; index_host++;

											mesh_host[index_host] = k; index_host++;

											mesh_host[index_host] = kk; index_host++;

											mesh_host[index_host] = l; index_host++;

											mesh_host[index_host] = ll; index_host++;

											mesh_host[index_host] = m; index_host++;

											mesh_host[index_host] = mm; index_host++;
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	// Device varibles
	int *mesh_device;
	int *n_device;
	float *integral_device, *x_device, *w_device, 
				*delta_device, *h_device, zero;

	cudaMalloc(&mesh_device, 10 * size * sizeof(int));
	cudaMemcpy(mesh_device, mesh_host, 10 * size * sizeof(int), 
			cudaMemcpyHostToDevice);

	cudaMalloc(&x_device, 5 * size * sizeof(float));

	cudaMalloc(&w_device, 3 * sizeof(float));
	cudaMemcpy(w_device, w_host, 3 * sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&delta_device, sizeof(float));
	cudaMemcpy(delta_device, &delta, sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&h_device, sizeof(float));
	cudaMemcpy(h_device, &h, sizeof(float), cudaMemcpyHostToDevice);

	cudaMalloc(&n_device, sizeof(int));
	cudaMemcpy(n_device, &n, sizeof(int), cudaMemcpyHostToDevice);

	zero = 0.0;
	cudaMalloc(&integral_device, sizeof(float));
	cudaMemcpy(integral_device, &zero, sizeof(float), cudaMemcpyHostToDevice);

	// Run on device
	assert(size % THREADS_PER_BLOCK == 0);
	cudaEventRecord(start);
	simpson_integration<<< size/THREADS_PER_BLOCK, THREADS_PER_BLOCK >>> (
			mesh_device, x_device, w_device, 
			delta_device, h_device, n_device, integral_device);
	cudaEventRecord(stop);

	// bring back
	cudaMemcpy(&integral_host, integral_device, sizeof(float), 
			cudaMemcpyDeviceToHost);

	// Finish timing
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&milliseconds, start, stop);

	printf("simpson cuda integral: %f, n:%d, time:%f\n", 
			integral_host, n, milliseconds/1000);

	// Free
	free(mesh_host);
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	cudaFree(mesh_device);
	cudaFree(x_device);
	cudaFree(w_device);
	cudaFree(delta_device);
	cudaFree(h_device);
	cudaFree(n_device);
	cudaFree(integral_device);
}
