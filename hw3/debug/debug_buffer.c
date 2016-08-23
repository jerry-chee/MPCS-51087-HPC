#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

typedef struct Bank_{
	unsigned long n; // number of particles
	unsigned long sz; // size of bank
	int *p; // particle array
	void (*resize)(struct Bank_ *b);
} Bank;

typedef struct Buffer_{         
	int mype_dest;
	Bank *bank;
	struct Buffer_ *next;
	struct Buffer_ *last;
} Buffer;

void resize_particles(Bank *b);
Bank *init_bank(unsigned long n_particles)
{
	Bank *b = malloc(sizeof(Bank));
	b->p = malloc(n_particles*sizeof(int));
	b->sz = n_particles;
	b->n = 0;
	b->resize = resize_particles;

	return b;
}

Buffer *init_buffer(int mype_dest)
{                                                                                    
	Buffer *buff = malloc(sizeof(Buffer));
	buff->mype_dest = mype_dest;
	buff->bank = init_bank(1);
	buff->next = NULL;
	buff->last = buff;

	return buff;
}

void resize_particles(Bank *b)
{
	b->p = realloc(b->p, sizeof(int)*2*b->sz);                                    
	b->sz = 2*b->sz;

	return;
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

void particle_deep_copy(int *p_to, int *p_from)
{
	*p_to = *p_from;
	/*
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
		 */
}

void add_particle(Buffer *buff, int *p, int mype_dest)
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

	particle_deep_copy(&(dest->bank->p[index]), p);
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
		printf("mype_dest:%d, bank.n:%d, bank,sz:%d\n", 
				curr->mype_dest, curr->bank->n, curr->bank->sz);
		curr = next;
		next = next->next;
	}
	printf("mype_dest:%d, bank.n:%d, bank,sz:%d\n", 
			curr->mype_dest, curr->bank->n, curr->bank->sz);
}

