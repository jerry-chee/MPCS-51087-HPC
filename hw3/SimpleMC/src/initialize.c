#include "header.h"

Parameters *init_parameters(void)
{
  Parameters *p = malloc(sizeof(Parameters));

  p->n_particles = 1000000;
  p->n_batches = 10;
  p->n_generations = 1;
  p->n_active = 10;
  p->bc = REFLECT;
  p->n_nuclides = 1;
  p->tally = TRUE;
  p->n_bins = 16;
  p->seed = 1;
  p->nu = 2.5;
  p->xs_f = 0.012;
  p->xs_a = 0.03;
  p->xs_s = 0.27;
  p->gx = 400;
  p->gy = 400;
  p->gz = 400;
  p->write_tally = FALSE;
  p->write_keff = FALSE;
  p->tally_file = NULL;
  p->keff_file = NULL;

  return p;
}

Parameters *init_parameters_loc(Parameters *parameters, int n_workers)
{
	Parameters *p_loc = malloc(sizeof(Parameters));

	int n_workers_axis = (int) cbrt(n_workers);
	assert(parameters->n_bins % n_workers_axis == 0);

	p_loc->n_particles = parameters->n_particles; 
	p_loc->n_batches = parameters->n_batches;
	p_loc->n_generations = parameters->n_generations;
	p_loc->n_active = parameters->n_active;
	p_loc->bc = parameters->bc;
	p_loc->tally = parameters->tally;
	p_loc->n_bins = parameters->n_bins / n_workers_axis; // different
	p_loc->seed = parameters->seed;
	p_loc->nu = parameters->nu;
	p_loc->xs_f = parameters->xs_f;
	p_loc->xs_a = parameters->xs_a;
	p_loc->xs_s = parameters->xs_s;
	p_loc->gx = parameters->gx / n_workers_axis; // different
	p_loc->gy = parameters->gy / n_workers_axis; // different
	p_loc->gz = parameters->gz / n_workers_axis; // different
	p_loc->write_tally = parameters->write_tally;
	p_loc->write_keff = parameters->write_keff;
	p_loc->tally_file = parameters->tally_file;
	p_loc->keff_file = parameters->keff_file;

	return p_loc;
}

Geometry *init_geometry(Parameters *parameters)
{
  Geometry *g = malloc(sizeof(Geometry));

  g->x = parameters->gx;
  g->y = parameters->gy;
  g->z = parameters->gz;
  g->bc = parameters->bc;

  return g;
}

Tally *init_tally(Parameters *parameters)
{
  Tally *t = malloc(sizeof(Tally));

  t->tallies_on = FALSE;
  t->n = parameters->n_bins;
  t->dx = parameters->gx/t->n;
  t->dy = parameters->gy/t->n;
  t->dz = parameters->gz/t->n;
  t->flux = calloc(t->n*t->n*t->n, sizeof(double));

  return t;
}

Material *init_material(Parameters *parameters)
{
  int i;
  Nuclide sum = {0, 0, 0, 0, 0};

  // Hardwire the material macroscopic cross sections for now to produce a keff
  // close to 1 (fission, absorption, scattering, total, atomic density)
  Nuclide macro = {parameters->xs_f, parameters->xs_a, parameters->xs_s,
     parameters->xs_f + parameters->xs_a + parameters->xs_s, 1.0};

  Material *m = malloc(sizeof(Material));
  m->n_nuclides = parameters->n_nuclides;
  m->nuclides = malloc(m->n_nuclides*sizeof(Nuclide));

  // Generate some arbitrary microscopic cross section values and atomic
  // densities for each nuclide in the material such that the total macroscopic
  // cross sections evaluate to what is hardwired above
  for(i=0; i<m->n_nuclides; i++){
    if(i<m->n_nuclides-1){
      m->nuclides[i].atom_density = rn()*macro.atom_density;
      macro.atom_density -= m->nuclides[i].atom_density;
    }
    else{
      m->nuclides[i].atom_density = macro.atom_density;
    }
    m->nuclides[i].xs_a = rn();
    sum.xs_a += m->nuclides[i].xs_a * m->nuclides[i].atom_density;
    m->nuclides[i].xs_f = rn();
    sum.xs_f += m->nuclides[i].xs_f * m->nuclides[i].atom_density;
    m->nuclides[i].xs_s = rn();
    sum.xs_s += m->nuclides[i].xs_s * m->nuclides[i].atom_density;
  }
  for(i=0; i<m->n_nuclides; i++){
    m->nuclides[i].xs_a /= sum.xs_a/macro.xs_a;
    m->nuclides[i].xs_f /= sum.xs_f/macro.xs_f;
    m->nuclides[i].xs_s /= sum.xs_s/macro.xs_s;
    m->nuclides[i].xs_t = m->nuclides[i].xs_a + m->nuclides[i].xs_s;
  }

  m->xs_f = parameters->xs_f;
  m->xs_a = parameters->xs_a;
  m->xs_s = parameters->xs_s;
  m->xs_t = parameters->xs_a + parameters->xs_s;

  return m;
}

Bank *init_source_bank(Parameters *parameters, Geometry *geometry)
{
  unsigned long i_p; // index over particles
  Bank *source_bank;

  // Initialize source bank
  source_bank = init_bank(parameters->n_particles);

  // Sample source particles
  for(i_p=0; i_p<parameters->n_particles; i_p++){
    sample_source_particle(geometry, &(source_bank->p[i_p]));
    source_bank->n++;
  }

  return source_bank;
}

Bank *init_fission_bank(Parameters *parameters)
{
  Bank *fission_bank;
  fission_bank = init_bank(2*parameters->n_particles);

  return fission_bank;
}

Bank *init_bank(unsigned long n_particles)
{
  Bank *b = malloc(sizeof(Bank));
  b->p = malloc(n_particles*sizeof(Particle));
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

void sample_source_particle(Geometry *geometry, Particle *p)
{
  p->alive = TRUE;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = rn()*geometry->x;
  p->y = rn()*geometry->y;
  p->z = rn()*geometry->z;

  return;
}

void sample_fission_particle(Particle *p, Particle *p_old)
{
  p->alive = TRUE;
  p->mu = rn()*2 - 1;
  p->phi = rn()*2*PI;
  p->u = p->mu;
  p->v = sqrt(1 - p->mu*p->mu)*cos(p->phi);
  p->w = sqrt(1 - p->mu*p->mu)*sin(p->phi);
  p->x = p_old->x;
  p->y = p_old->y;
  p->z = p_old->z;

  return;
}

void resize_particles(Bank *b, long int newsize)
{
	assert(newsize > b->sz);

	b->p = realloc(b->p, sizeof(Particle)*newsize);
	b->sz = newsize;

  //b->p = realloc(b->p, sizeof(Particle)*2*b->sz);
  //b->sz = 2*b->sz;

  return;
}

void print_particle(Particle *p)
{
	printf("p->alive:%d \t p->mu:%lf \t p->phi:%lf \n p->u:%lf \t p->v:%lf \t p->w:%lf \n p->x:%lf \t p->y:%lf \t p->z:%lf \n p->surface_crossed:%d\n", 
			p->alive, p->mu, p->phi, p->u, p->v, p->w, p->x, p->y, p->z,
			p->surface_crossed);
}

void print_bank(Bank *b, int p_true)
{
	printf("bank->n:%ld, bank->sz:%ld\n", b->n, b->sz);

	if (p_true)
		for (int i = 0; i < b->n; i++)
			print_particle(&(b->p[i]));
}

void free_bank(Bank *b)
{
	free(b->p);
	b->p = NULL;
	free(b);
	b = NULL;

	return;
}

void free_buffer(Buffer *buff)
{
	if (buff != NULL)
	{
		Buffer *curr, *next;
		curr = buff;
		next = curr->next;
		while(next != NULL)
		{
			free_bank(curr->bank);
			free(curr);
			
			curr = next;
			next = next->next;
		}
		free_bank(curr->bank);
		free(curr);
	}
}

void free_material(Material *material)
{
	free(material->nuclides);
	material->nuclides = NULL;
	free(material);
	material = NULL;

	return;
}

void free_tally(Tally *tally)
{
	free(tally->flux);
	tally->flux = NULL;
	free(tally);
	tally = NULL;

	return;
}
