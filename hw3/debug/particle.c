#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "mpi.h"

//#pragma pack(push,1)
typedef struct Particle_{
	int alive;
	double mu; // cosine of polar angle
	double phi; // azimuthal angle
	double u; // direction
	double v;
	double w;
	double x; // position
	double y;
	double z;
	int surface_crossed;
} Particle;
//#pragma pack(pop)

void print_particle(Particle *p) 
{
	printf("p->alive:%d \t p->mu:%lf \t p->phi:%lf \n p->u:%lf \t p->v:%lf \t p->w:%lf \n p->x:%lf \t p->y:%lf \t p->z:%lf \n p->surface_crossed:%d\n", 
			p->alive, p->mu, p->phi, p->u, p->v, p->w, p->x, p->y, p->z,
			p->surface_crossed);
}
int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

	int mype, nprocs;
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);

	// Type Commit
	int status;
	int blocks[10] = {1,1,1,1,1,1,1,1,1,1};
//#pragma pack(push,1)
	MPI_Datatype types[10] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
//#pragma pack(pop)
	MPI_Aint displacements[10];

	MPI_Get_address(&types[0], &displacements[0]);     
	MPI_Get_address(&types[1], &displacements[1]);     
	MPI_Get_address(&types[2], &displacements[2]);
	MPI_Get_address(&types[3], &displacements[3]);
	MPI_Get_address(&types[4], &displacements[4]);
	MPI_Get_address(&types[5], &displacements[5]);
	MPI_Get_address(&types[6], &displacements[6]);
	MPI_Get_address(&types[7], &displacements[7]);
	MPI_Get_address(&types[8], &displacements[8]);
	MPI_Get_address(&types[9], &displacements[9]);

	for (int i = 1; i < 10; i++)
		displacements[i] = displacements[i] - displacements[0];
	displacements[0] = 0;

	/*
	MPI_Aint displacements[10];
	displacements[0] = 0;
	for (int i=1; i<10; i++)
	{
		displacements[i] = mem_locations[i] - mem_locations[i-1];                        
		//printf("displacements[i]:%d\n", displacements[i]);
	}   
	*/

	MPI_Datatype MPI_Particle;
	status = MPI_Type_create_struct(10, blocks, displacements,
			types, &MPI_Particle);
	assert(status == MPI_SUCCESS);
	MPI_Type_commit(&MPI_Particle);
	
	Particle p;
	int n;

	if (mype == 0)
	{
		p.alive           = 1;
		p.mu              = 1;
		p.phi             = 1;
		p.u               = 1;
		p.v               = 1;
		p.w               = 1;
		p.x               = 1;
		p.y               = 1;
		p.z               = 1;
		p.surface_crossed = 1;

		MPI_Send(&p, 1, MPI_Particle, 1, 0, MPI_COMM_WORLD);
		//n = 1;
		//MPI_Send(&n, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
	}
	if (mype == 1)
	{
		MPI_Status status;
		//MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
		//printf("%d\n", n);
		MPI_Recv(&p, 1, MPI_Particle, 0, 0, MPI_COMM_WORLD, &status);
		print_particle(&p);
	}

	MPI_Finalize();
}
