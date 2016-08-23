#include "debug_buffer.c"

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

void add_particle_to_buff(Buffer *buff, double x_loc_dest, double y_loc_dest,
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

	int index = dest->bank->n;	
	if (index+1 > dest->bank->sz)
		dest->bank->resize(dest->bank);

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

void add_particles_to_bank(Bank *bank, Particle *ps, int count)
{
	assert(bank != NULL);
	assert(ps		!= NULL);
	
	int start = bank->n;
	int end   = bank->n + count;
	if (end		> bank->sz)
		bank->resize(bank);	

	for (int i_p = start; i_p < end; i_p++)
	{
		particle_deep_copy(&(bank->p[i_p]), &(ps[i_p - start]));
	}
}

void source_scatter(int mype, int n_workers, double d_proc, 
		Bank *source_bank, int buff_lim)
{
	// TEMP	
	if (mype == 0)
	{
		int mype_dest;
		int x_lvl_dest, y_lvl_dest, z_lvl_dest; 
		double x_loc_dest, y_loc_dest, z_loc_dest;
		
		Buffer *master_buff = init_buffer(-1);
		Buffer *search_buff;
		int count;
		int tag = 0;

		// loop over source_bank
		for (int i_p = 0; i_p < source_bank->n; i_p++)
		{
			// find loc coords and mype_dest
			p = &(source_bank->p[i_p]);
			
			x_lvl_dest = lvl_from_global(p->x, d_proc);
			y_lvl_dest = lvl_from_global(p->y, d_proc);
			z_lvl_dest = lvl_from_global(p->z, d_proc);

			mype_dest = mype_from_lvls(x_lvl_dest, y_lvl_dest, 
					z_lvl_dest, n_workers);

			x_loc_dest = coord_loc(p->x, d_proc, x_lvl_dest);
			y_loc_dest = coord_loc(p->y, d_proc, y_lvl_dest);
			z_loc_dest = coord_loc(p->z, d_proc, z_lvl_dest);

			// deep copy particle to buff
			add_particle_to_buff(master_buff, x_loc_dest, y_loc_dest, 
					z_loc_dest, p, mype_dest);
			
			search_buff = search_buffer(master_buff, mype_dest);
			assert(search_buff != NULL);
			// if search_buff too big send it off
			if (search_buff->bank->n > buff_lim)
			{
				p_bundle = search_buff->bank->p;
				count = search_buff->bank->n;				
				MPI_Send(&count, 1, MPI_LONG_INT, 
						mype_dest, tag, MPI_COMM_WORLD);
				MPI_Send(p_bundle, count, MPI_Particle, 
						mype_dest, tag, MPI_COMM_WORLD);
				
				// empty the *p in search_buff
				free(search_buff->bank->p);
				search_buff->bank->n = 0;
			}	
		}

		// send off remaing buffer
		Buffer *curr, *next;
		curr = master_buff->next;
		if (curr == NULL)
			next = NULL;
		else 
			next = curr->next;
		while (next != NULL)
		{
			p_bundle = curr->bank->p;
			count = curr->bank->n;				
			MPI_Send(&count, 1, MPI_LONG_INT, 
					mype_dest, tag, MPI_COMM_WORLD);
			MPI_Send(p_bundle, count, MPI_Particle, 
					mype_dest, tag, MPI_COMM_WORLD);
		
			curr = next;
			next = next->next;			
		}
		assert(curr != NULL);
		p_bundle = curr->bank->p;
		count = curr->bank->n;				
		MPI_Send(&count, 1, MPI_LONG_INT, 
				mype_dest, tag, MPI_COMM_WORLD);
		MPI_Send(p_bundle, count, MPI_Particle, 
				mype_dest, tag, MPI_COMM_WORLD);

		// tell all workers done sending src_bank
		for (int dest = 1; dest <= n_workers; dest++)
		{
			int exit = -1;
			MPI_Send(&exit, 1, MPI_Particle, 
					dest, tag, MPI_COMM_WORLD);		 	
		}

		free_buffer(master_buff);
	} 
	else // workers
	{
		long int p_count = 1;
		Particle *p_tmp = (Particle *) malloc(p_count * sizeof(Particle));
		int count;
		int tag = 0;
		int src = 0;
		MPI_Status status;
		int stat;

		while (1)
		{
			stat = MPI_Recv(&count, 1, MPI_LONG_INT, src, tag, 
					MPI_COMM_WORLD, &status);
			assert(stat == MPI_SUCCESS);

			// check exist status
			if (count == -1)
				break;

			// resize p_tmp if necessary
			if(count > p_count)
			{
				free(p_tmp);
				p_count = count;
				p_tmp = (Particle *) malloc(p_count * sizeof(Particle));
			}
			
			stat = MPI_Recv(&p_tmp, count, MPI_Particle, src, tag,
					MPI_COMM_WORLD, &status);
			assert(stat == MPI_SUCCESS);

			add_particles_to_bank(source_bank, p_tmp, count);
		}

		free(p_tmp);
	}
}

