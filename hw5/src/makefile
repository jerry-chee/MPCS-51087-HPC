FLAGS = -lm -g

mc_openmp: mc_nd.c
	gcc $(FLAGS) -Wall mc_nd.c -o mc_openmp -fopenmp

mc_cuda: mc_nd.cu
	nvcc $(FLAGS) mc_nd.cu -o mc_cuda 

simpson_openmp: simpson_nd.c
	gcc $(FLAGS) -Wall simpson_nd.c -o simpson_openmp -fopenmp

simpson_cuda: simpson_nd.cu
	nvcc $(FLAGS) -G simpson_nd.cu -o simpson_cuda
