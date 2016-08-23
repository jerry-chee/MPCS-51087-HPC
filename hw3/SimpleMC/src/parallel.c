#include "header.h"

// MPI_Particle struct
void declare_MPI_Particle(void)
{	
	int status;
	int blocks[10] = {1,1,1,1,1,1,1,1,1,1};
	MPI_Datatype types[10] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
	MPI_Aint displacements[10];
	MPI_Get_address(&types[1], &displacements[1]);     
	MPI_Get_address(&types[2], &displacements[2]);
	MPI_Get_address(&types[3], &displacements[3]);
	MPI_Get_address(&types[4], &displacements[4]);
	MPI_Get_address(&types[5], &displacements[5]);
	MPI_Get_address(&types[6], &displacements[6]);
	MPI_Get_address(&types[7], &displacements[7]);
	MPI_Get_address(&types[8], &displacements[8]);
	MPI_Get_address(&types[9], &displacements[9]);

	MPI_Datatype MPI_Particle;
	status = MPI_Type_struct(10, blocks, displacements, types, &MPI_Particle);
	assert(status == MPI_SUCCESS);
	MPI_Type_commit(&MPI_Particle);
}

void add_buffer(Buffer *buff, int mype_dest)
{
	assert(buff != NULL);
	Buffer *last = buff->last;
	last->next = init_buffer(mype_dest);
	buff->last = last->next;
}	

Buffer *search_buffer(Buffer *buff, int mype_dest)
{
	assert(buff != NULL);
	Buffer *curr, *next;
	curr = buff;
	next = curr->next;
	while(next != NULL)
	{
		if (curr->mype_dest == mype_dest)
			return curr;
		curr = next;
		next = next->next;
	}
	if (curr->mype_dest == mype_dest)
		return curr;
	else
		return NULL;
}

void particle_deep_copy_coord(double x_loc_dest, double y_loc_dest,
		 double z_loc_dest, Particle *p_to, Particle *p_from)
{
	p_to->alive						= p_from->alive;
	p_to->mu							= p_from->mu;
	p_to->phi							= p_from->phi;
	p_to->u								= p_from->u;
	p_to->v								= p_from->v;
	p_to->w								= p_from->w;
	p_to->x								= x_loc_dest;
	p_to->y								= y_loc_dest;
	p_to->z								= z_loc_dest;
	p_to->surface_crossed	= p_from->surface_crossed;
}

void particle_deep_copy(Particle *p_to, Particle *p_from)
{
	assert(p_to != NULL);
	assert(p_from != NULL);
	p_to->alive						= p_from->alive;
	p_to->mu							= p_from->mu;
	p_to->phi							= p_from->phi;
	p_to->u								= p_from->u;
	p_to->v								= p_from->v;
	p_to->w								= p_from->w;
	p_to->x								= p_from->x;
	p_to->y								= p_from->y;
	p_to->z								= p_from->z;
	p_to->surface_crossed	= p_from->surface_crossed;
}

void add_particle_to_buffer(Buffer *buff, double x_loc_dest, double y_loc_dest,
		double z_loc_dest, Particle *p, int mype_dest)
{
	assert(buff != NULL);
	Buffer *dest;
	dest = search_buffer(buff, mype_dest);

	if (dest == NULL)
	{
		add_buffer(buff, mype_dest);
		dest = buff->last;
	}

	long int index = dest->bank->n;	
	if (index+1 > dest->bank->sz)
		dest->bank->resize(dest->bank, 2*(index + 1));

	particle_deep_copy_coord(x_loc_dest, y_loc_dest, z_loc_dest, 
			&(dest->bank->p[index]), p);
	dest->bank->n++;
}

void print_buffer(Buffer *buff)
{
	assert(buff != NULL);
	Buffer *curr, *next;
	curr = buff;
	next = curr->next;
	while (next != NULL)
	{
		printf("mype_dest:%d, bank.n:%ld, bank,sz:%ld\n", 
				curr->mype_dest, curr->bank->n, curr->bank->sz);
		curr = next;
		next = next->next;
	}
	printf("mype_dest:%d, bank.n:%ld, bank,sz:%ld\n", 
			curr->mype_dest, curr->bank->n, curr->bank->sz);
}                      

void print_single_buffer(Buffer *buff)
{
	assert(buff != NULL);
	printf("mype_dest:%d, bank.n:%ld, bank,sz:%ld\n", 
			buff->mype_dest, buff->bank->n, buff->bank->sz);
}

void add_particles_to_bank(Bank *bank, Particle *ps, long int count)
{
	assert(bank != NULL);
	assert(ps		!= NULL);
	
	long int start = bank->n;
	long int end   = bank->n + count;
	//printf("add_particles_to_bank_count:%ld, start:%ld, end:%ld, n:%ld, sz:%ld\n",
	//		count, start, end, bank->n, bank->sz);
	if (end		> bank->sz)
		bank->resize(bank, 2*end);	
	for (long int i_p = start; i_p < end; i_p++)
	{
		//printf("i_p:%ld\n", i_p);
		particle_deep_copy(&(bank->p[i_p]), &(ps[i_p - start]));
	}
	bank->n = end;
}

// writes over source_bank for workers
void source_scatter(int mype, int n_workers, double d_proc, 
		Parameters *parameters, Bank *source_bank, int buff_lim)
{
	//Type Commit
	int status;
	int blocks[10] = {1,1,1,1,1,1,1,1,1,1};
	MPI_Datatype types[10] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
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
	{
		displacements[i] = displacements[i] - displacements[0];
		//printf("displacements[i]:%d\n", displacements[i]);
	}
	displacements[0] = 0;
	MPI_Datatype MPI_Particle;
	status = MPI_Type_create_struct(10, blocks, displacements, 
			types, &MPI_Particle);
	assert(status == MPI_SUCCESS);
	MPI_Type_commit(&MPI_Particle);

	int count;
	int tag = 0;
	int master = n_workers;

	if (mype == master) // master
	{
		int mype_dest;
		int x_lvl, y_lvl, z_lvl;
		double x_loc, y_loc, z_loc;
		double delta = 0.00000001;

		Particle *p_tmp;
		Buffer *master_buff = init_buffer(-1);
		Buffer *search_buff;
		Particle *p_bundle;

		// loop over source_bank particles
		for (int i_p = 0; i_p < source_bank->n; i_p++)
		{
			if (0 > source_bank->p[i_p].x)
				source_bank->p[i_p].x = 0;
			if (0 > source_bank->p[i_p].y)
				source_bank->p[i_p].y = 0;
			if (0 > source_bank->p[i_p].z)
				source_bank->p[i_p].z = 0;

			if (source_bank->p[i_p].x >= parameters->gx)
			 source_bank->p[i_p].x = parameters->gx - delta;	
			if (source_bank->p[i_p].y >= parameters->gy)
			 source_bank->p[i_p].y = parameters->gy - delta;	
			if (source_bank->p[i_p].z >= parameters->gz)
			 source_bank->p[i_p].z = parameters->gz - delta;	
		
			//if (p_tmp->y >= parameters->gy)
			//	printf("%lf\n", p_tmp->y-delta);

			p_tmp = &(source_bank->p[i_p]);
			assert(p_tmp != NULL);
			
			assert(0 <= p_tmp->x && p_tmp->x < parameters->gx);
			assert(0 <= p_tmp->y);
			assert(p_tmp->y < parameters->gy);
			assert(0 <= p_tmp->z && p_tmp->z < parameters->gz);
		
			x_lvl = lvl_from_global(p_tmp->x, d_proc);
			y_lvl = lvl_from_global(p_tmp->y, d_proc);
			z_lvl = lvl_from_global(p_tmp->z, d_proc);

			mype_dest = mype_from_lvls(x_lvl, y_lvl, z_lvl, n_workers);
			if (mype_dest < 0 || mype_dest >= n_workers)
				printf("%d,%d,%d,mype_dest:%d\n", x_lvl, y_lvl, z_lvl, mype_dest);
			assert(0 <= mype_dest && mype_dest < n_workers);

			x_loc = coord_loc(p_tmp->x, d_proc, x_lvl);
			y_loc = coord_loc(p_tmp->y, d_proc, y_lvl);
			z_loc = coord_loc(p_tmp->z, d_proc, z_lvl);
			
			assert(0 <= x_loc && x_loc < d_proc); 
			assert(0 <= y_loc && y_loc < d_proc); 
			assert(0 <= z_loc && z_loc < d_proc); 
			
			//printf("p->x:%lf, x_loc:%lf\n", p_tmp->x, x_loc);
			//printf("p->y:%lf, y_loc:%lf\n", p_tmp->y, y_loc);
			//printf("p->z:%lf, z_loc:%lf\n", p_tmp->z, z_loc);
			//printf("mype_dest:%d\n", mype_dest);

			add_particle_to_buffer(master_buff, x_loc, y_loc, z_loc,
					p_tmp, mype_dest);
			search_buff = search_buffer(master_buff, mype_dest);
			assert(search_buff != NULL);

			// if buffer exceeds limit
			if (search_buff->bank->n >= buff_lim)
			{
				p_bundle = search_buff->bank->p;
				count = (int) search_buff->bank->n; //buff_lim an int so wont be too big
				//printf("master count:%d\n", count);
				assert(p_bundle != NULL);
				
				MPI_Send(&count, 1, MPI_INT, mype_dest, 
						tag, MPI_COMM_WORLD);
				MPI_Send(p_bundle, count, MPI_Particle, mype_dest,
						tag, MPI_COMM_WORLD);
				
				search_buff->bank->n = 0;
			}
		}
	
		// send off the rest of buffer
		Buffer *curr, *next;
		curr = master_buff->next; //cause head is mype_dest=-1

		if (curr == NULL)
			next = NULL;
		else
			next = curr->next;

		while (next != NULL)
		{
			mype_dest = curr->mype_dest;
			assert(0 <= mype_dest && mype_dest < n_workers);
			p_bundle = curr->bank->p;
			count = (int) curr->bank->n; //buff_lim an int so wont be too big
			//printf("master count:%d\n", count);
			assert(p_bundle != NULL);

			MPI_Send(&count, 1, MPI_INT, mype_dest, 
					tag, MPI_COMM_WORLD);
			MPI_Send(p_bundle, count, MPI_Particle, mype_dest,
					tag, MPI_COMM_WORLD);

			curr = next;
			next = next->next;
		}

		// last buff
		assert(curr != NULL);

		mype_dest = curr->mype_dest;
		assert(0 <= mype_dest && mype_dest < n_workers);
		p_bundle = curr->bank->p;
		count = (int) curr->bank->n; //buff_lim an int so wont be too big
		//printf("last in buff master count:%d\n", count);
		assert(p_bundle != NULL);

		MPI_Send(&count, 1, MPI_INT, mype_dest, 
				tag, MPI_COMM_WORLD);
		MPI_Send(p_bundle, count, MPI_Particle, mype_dest,
				tag, MPI_COMM_WORLD);

		// free buff
		free_buffer(master_buff);


		// exit sending
		int exit = -1;
		for (int dest = 0; dest < n_workers; dest++)
		{
			MPI_Send(&exit, 1, MPI_INT, dest,
					tag, MPI_COMM_WORLD);	
		}
	}
	else // workers
	{
		// empty source bank
		source_bank->n = 0;

		MPI_Status status;
		int stat;

		int p_sz = 1;
		Particle *p_buff = malloc(sizeof(Particle) * p_sz);

		while(1)
		{
			stat = MPI_Recv(&count, 1, MPI_INT, master,
					tag, MPI_COMM_WORLD, &status);
			assert(stat == MPI_SUCCESS);
			//printf("%dproc count:%d\n", mype, count);

			if (count == -1)
				break;

			if (count > p_sz)
			{
				p_sz = 2*count;
				p_buff = realloc(p_buff, sizeof(Particle)*p_sz);
				assert(p_buff != NULL);
			}

			stat = MPI_Recv(p_buff, count, MPI_Particle, master,
					tag, MPI_COMM_WORLD, &status);
			assert(stat == MPI_SUCCESS);

			//print_particle(p_buff);

			//printf("%dproc, add_particles_to_bank count:%d\n", mype, count);
			add_particles_to_bank(source_bank, p_buff, count);
		}

		for (int i = 0; i < source_bank->n; i++)
		{
			assert(0 <= source_bank->p[i].x && source_bank->p[i].x <= d_proc);
			assert(0 <= source_bank->p[i].y && source_bank->p[i].y <= d_proc);
			assert(0 <= source_bank->p[i].z && source_bank->p[i].z <= d_proc);
		}
		free(p_buff);
	}

	MPI_Type_free(&MPI_Particle);
}	

// writes over fission_bank in master
void fission_gather(int mype, int n_workers, double d_proc,
		Parameters *parameters, Parameters *parameters_loc, Bank *fission_bank, int buff_lim)
{
	//Type Commit
	int status;
	int blocks[10] = {1,1,1,1,1,1,1,1,1,1};
	MPI_Datatype types[10] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
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

	MPI_Datatype MPI_Particle;
	status = MPI_Type_create_struct(10, blocks, displacements, 
			types, &MPI_Particle);
	assert(status == MPI_SUCCESS);
	MPI_Type_commit(&MPI_Particle);

	// coordinate transform to global
	if (mype != n_workers)
	{
		int x_lvl, y_lvl, z_lvl;
		double x_global, y_global, z_global;
		int n_workers_axis = (int) cbrt(n_workers);
		Particle *p;

		double delta = 0.0000001;
		for (long int i_p = 0; i_p < fission_bank->n; i_p++)
		{
			if (0 > fission_bank->p[i_p].x)
				fission_bank->p[i_p].x = 0;
			if (0 > fission_bank->p[i_p].y)
				fission_bank->p[i_p].y = 0;
			if (0 > fission_bank->p[i_p].z)
				fission_bank->p[i_p].z = 0;

			if (fission_bank->p[i_p].x >= parameters_loc->gx)
			 fission_bank->p[i_p].x = parameters_loc->gx - delta;	
			if (fission_bank->p[i_p].y >= parameters_loc->gy)
			 fission_bank->p[i_p].y = parameters_loc->gy - delta;	
			if (fission_bank->p[i_p].z >= parameters_loc->gz)
			 fission_bank->p[i_p].z = parameters_loc->gz - delta;	

			x_lvl = x_lvl_from_mype(mype, n_workers);
			y_lvl = y_lvl_from_mype(mype, n_workers);
			z_lvl = z_lvl_from_mype(mype, n_workers);

			assert(0 <= x_lvl && x_lvl < n_workers_axis);
			assert(0 <= y_lvl && y_lvl < n_workers_axis);
			assert(0 <= z_lvl && z_lvl < n_workers_axis);

			p = &(fission_bank->p[i_p]);
			x_global = coord_global(p->x, d_proc, x_lvl);
			y_global = coord_global(p->y, d_proc, y_lvl);
			z_global = coord_global(p->z, d_proc, z_lvl);

			assert(0 <= p->x && p->x <= d_proc);
			assert(0 <= p->y && p->y <= d_proc);
			assert(0 <= p->z && p->z <= d_proc);

			// reasign global coord
			p->x = x_global;
			p->y = y_global;
			p->z = z_global;

			assert(0 <= p->x && p->x <= parameters->gx);
			assert(0 <= p->y && p->y <= parameters->gy);
			assert(0 <= p->z && p->z <= parameters->gz);
		}
	}

	// communicate to master how much being sent over
	// do it in rounds
	int sendcount = 1;
	int recvcount = 1;
	int root = n_workers;
	int round_counts[n_workers+1];
	long int sum_tot = 0;
	int round_sum = 0;
	int round_disp[n_workers+1];

	// I cannot send all of them at once I think
	int round_size = 0;
	long int curr_index = 0;
	int running_sum = 0;
	int i;


	Particle *p_sendbuf;
	int p_sendcount;
	Particle *p_recvbuf;
	int *p_recvcounts;

	// if curr_index < loc_count for a proc then send nothing
	while (1)
	{
		if (fission_bank->n < curr_index+1)
			round_size = 0;
		else if (buff_lim > (fission_bank->n - curr_index+1))
			round_size = fission_bank->n - curr_index+1;
		else
			round_size = buff_lim;

		if (mype == n_workers)
			round_size = 0;
		assert(round_size >= 0);
		//printf("mype:%d, round_size:%d\n", mype, round_size);

		MPI_Allgather(&round_size, sendcount, MPI_INT,
				round_counts, recvcount, MPI_INT,
				MPI_COMM_WORLD);

		// check that not all round_counts are zero except master
		round_sum = 0;
		for (i = 0; i < n_workers; i++)
		{
			round_sum += round_counts[i];
			//printf("mype:%d, i:%d, round_counts[i]:%d\n",
			//		mype, i, round_counts[i]);
		}
		if (round_sum == 0)
			break;

		if (mype == n_workers) // master
		{
			sum_tot += round_sum;
			//printf("round_sum:%d\n", round_sum);
			round_disp[0] = 0;
			running_sum = 0;
			for (i = 1; i < n_workers; i++)
			{
				running_sum += round_counts[i-1];
				round_disp[i] = running_sum;
			}
			round_disp[n_workers] = round_disp[n_workers-1];

			// resize fission_bank if necessary
			if (sum_tot > fission_bank->sz)
				fission_bank->resize(fission_bank, 2*sum_tot);
		}

		// send over particles
		p_sendbuf = &(fission_bank->p[curr_index]);
		p_sendcount = round_size;
		p_recvbuf = &(fission_bank->p[sum_tot-round_sum]);
		p_recvcounts = round_counts;

		//printf("round_sum:%d, sum_tot:%ld\n", round_sum, sum_tot);
		MPI_Gatherv(p_sendbuf, p_sendcount, MPI_Particle,
				p_recvbuf, p_recvcounts, round_disp,
				MPI_Particle, root, MPI_COMM_WORLD); 

		if (mype == n_workers)
			fission_bank->n = sum_tot;

		if (mype != n_workers) 
			curr_index += buff_lim;
	}

	/*
		 if (mype == n_workers)
		 {
		 printf("fission_bank->n:%ld\n", fission_bank->n);
		 printf("fission_bank->sz:%ld\n", fission_bank->sz);
		 }
		 */
	MPI_Type_free(&MPI_Particle);
}

// returns whether or not to do another round of buffer exchanges
int comm_ready(int mype, int n_workers, Buffer *buff)
{
	int buff_content;
	if (mype != n_workers) // workers
		if (buff->next != NULL)
			buff_content = 1;

	if (mype == n_workers)
		buff_content = 0;

	int sendcount = 1;
	int recvcount = 1;
	int buff_tally[n_workers+1];

	MPI_Allgather(&buff_content, sendcount, MPI_INT,
			buff_tally, recvcount, MPI_INT,
			MPI_COMM_WORLD);

	// tally
	int buff_sum = 0;
	for (int i = 0; i < n_workers; i++)
		buff_sum += buff_tally[i];

	return buff_sum;
}

// overwrites incoming particles into source_bank. Does not free buff
void comm_routine(int mype, int n_workers, Buffer *buff,
		Bank *source_bank)
{
	//Type Commit
	int stat;
	int blocks[10] = {1,1,1,1,1,1,1,1,1,1};
	MPI_Datatype types[10] = {MPI_INT, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, 
		MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
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

	MPI_Datatype MPI_Particle;
	stat = MPI_Type_create_struct(10, blocks, displacements, 
			types, &MPI_Particle);
	assert(stat == MPI_SUCCESS);
	MPI_Type_commit(&MPI_Particle);

	if (mype != n_workers)
	{
		// overwrite source_bank
		source_bank->n = 0;
	}
	Buffer *curr, *next;
	Particle *p_buff;

	MPI_Status status;
	int sendcount;
	int dest;
	int tag = 0;
	int exit = -1;

	// Each proc shares once at a time
	for (int i = 0; i < n_workers; i++)
	{
		if (i == mype)
		{
			// if buff non-empty
			curr = buff->next;

			if (curr != NULL)
			{
				next = curr->next;
				assert(curr != NULL);

				while (next != NULL)
				{
					p_buff = curr->bank->p;
					sendcount = (int) curr->bank->n; // shouldn't be a type problem
					dest = curr->mype_dest;

					MPI_Send(&sendcount, 1, MPI_INT, dest, 
							tag, MPI_COMM_WORLD);
					MPI_Send(p_buff, sendcount, MPI_Particle, dest, 
							tag, MPI_COMM_WORLD);

					curr = next;
					next = next->next;
				}

				assert(curr != NULL);
				p_buff = curr->bank->p;
				sendcount = (int) curr->bank->n; // shouldn't be a type problem
				dest = curr->mype_dest;

				MPI_Send(&sendcount, 1, MPI_INT, dest, 
						tag, MPI_COMM_WORLD);
				MPI_Send(p_buff, sendcount, MPI_Particle,
						dest, tag, MPI_COMM_WORLD);
			}

			// free buffer
			free_buffer(buff);

			// exit sending
			for (dest = 0; dest < n_workers; dest++)
			{
				if (dest != i)
				{
					MPI_Send(&exit, 1, MPI_INT, dest,
							tag, MPI_COMM_WORLD);
				}
			}
		}	
		else if ((i != mype) && (n_workers != mype))
		{
			while (1)
			{
				stat = MPI_Recv(&sendcount, 1, MPI_INT, i,
						tag, MPI_COMM_WORLD, &status);
				assert(stat == MPI_SUCCESS);

				if (sendcount == -1)
					break;

				if (source_bank->n + sendcount > source_bank->sz)
					source_bank->resize(source_bank, 2*(source_bank->n + sendcount));

				p_buff = &(source_bank->p[source_bank->n]);
				stat = MPI_Recv(p_buff, sendcount, MPI_Particle, i,
						tag, MPI_COMM_WORLD, &status);
				assert(stat == MPI_SUCCESS);

				source_bank->n += sendcount;
			}
		}
		else
		{
		}

		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Type_free(&MPI_Particle);
}
