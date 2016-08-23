#include "header.h"

// Main logic to move particle
void transport(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, Particle *p)
{
  double d_b;
  double d_c;
  double d;

  while(p->alive){

    // Find distance to boundary
    d_b = distance_to_boundary(geometry, p);

    // Find distance to collision
    d_c = distance_to_collision(material);

    // Take smaller of two distances
    d = d_b < d_c ? d_b : d_c;

    // Advance particle
    p->x = p->x + d*p->u;
    p->y = p->y + d*p->v;
    p->z = p->z + d*p->w;

    // Case where particle crosses boundary
    if(d_b < d_c){
      cross_surface(geometry, p);
    }
    // Case where particle has collision
    else{
      collision(material, fission_bank, parameters->nu, p);

      // Score tallies
      if(tally->tallies_on == TRUE){
        score_tally(parameters, material, tally, p);
      }
    }
  }
  return;
}

void transport_loc(int mype, int n_workers, double d_proc, Parameters *parameters_loc, Geometry *geometry, Geometry *geometry_local, Material *material, Bank *source_bank, Bank *fission_bank, Buffer *buff, Tally *tally, Particle *p)
{
  double d_b;
  double d_c;
	//double d;
	double delta = 0.0000001;

	if (mype != n_workers) // workers
	{
		while(p->alive){
			// Find distance to boundary
			d_b = distance_to_boundary(geometry_local, p);

			// Find distance to collision
			d_c = distance_to_collision(material);

			// Take smaller of two distances
			//d = d_b < d_c ? d_b : d_c;

			// Advance particle by d
			p->x = p->x + d_c*p->u;
			p->y = p->y + d_c*p->v;
			p->z = p->z + d_c*p->w;

			// if particles out of bounds, place onto boundary
			if (d_b > d_c)
			{
				if (p->x > geometry_local->x)
					p->x = geometry_local->x - delta;
				if (p->y > geometry_local->y)
					p->y = geometry_local->y - delta;
				if (p->z > geometry_local->z)
					p->z = geometry_local->z - delta;

				if (0 > p->x)
					p->x = 0;
				if (0 > p->y)
					p->y = 0;
				if (0 > p->z)
					p->z = 0;
			}

			// Case where particle crosses boundary
			if(d_b < d_c){
				//printf("%dproc cross_surface\n", mype);
				cross_surface_loc(mype, n_workers, d_proc, buff, 
						geometry, geometry_local, p);
			}
			// Case where particle has collision
			else{

				assert(0 <= p->x && p->x < geometry_local->x);
				assert(0 <= p->y && p->y < geometry_local->y);
				assert(0 <= p->z && p->z < geometry_local->z);
				
				collision(material, fission_bank, parameters_loc->nu, p);

				// Score tallies
				if(tally->tallies_on == TRUE){
					score_tally(parameters_loc, material, tally, p);
				}
			}
		}
	}
	return;
}

// Returns the distance to the nearest boundary for a particle traveling in a
// certain direction
double distance_to_boundary(Geometry *geometry, Particle *p)
{
	int i;
	double dist;
	double d = D_INF;
	int    surfaces[6] = {X0, X1, Y0, Y1, Z0, Z1};
	double p_angles[6] = {p->u, p->u, p->v, p->v, p->w, p->w};
	double p_coords[6] = {p->x, p->x, p->y, p->y, p->z, p->z};
	double s_coords[6] = {0, geometry->x, 0, geometry->y, 0, geometry->z};

	for(i=0; i<6; i++){
		if(p_angles[i] == 0){
			dist = D_INF;
		}
		else{
			dist = (s_coords[i] - p_coords[i])/p_angles[i];
			if(dist <= 0){
				dist = D_INF;
			}
		}
		if(dist < d){
			d = dist;
			p->surface_crossed = surfaces[i];
		}
	}

	return d;
}

// Returns the distance to the next collision for a particle
double distance_to_collision(Material *material)
{
	double d;

	if(material->xs_t == 0){
		d = D_INF;
	}
	else{
		d = -log(rn())/material->xs_t;
	}

	return d;
}

// Handles a particle crossing a surface in the geometry
void cross_surface(Geometry *geometry, Particle *p)
{
	// Handle vacuum boundary conditions (particle leaks out)
	if(geometry->bc == VACUUM){
		p->alive = FALSE;
	}

	// Handle reflective boundary conditions
	else if(geometry->bc == REFLECT){
		if(p->surface_crossed == X0){
			p->u = -p->u;
			p->x = 0.0;
		}
		else if(p->surface_crossed == X1){
			p->u = -p->u;
			p->x = geometry->x;
		}
		else if(p->surface_crossed == Y0){
			p->v = -p->v;
			p->y = 0.0;
		}
		else if(p->surface_crossed == Y1){
			p->v = -p->v;
			p->y = geometry->y;
		}
		else if(p->surface_crossed == Z0){
			p->w = -p->w;
			p->z = 0.0;
		}
		else if(p->surface_crossed == Z1){
			p->w = -p->w;
			p->z = geometry->z;
		}
	}

	// Handle periodic boundary conditions
	else if(geometry->bc == PERIODIC){
		if(p->surface_crossed == X0){
			p->x = geometry->x;
		}
		else if(p->surface_crossed == X1){
			p->x = 0;
		}
		else if(p->surface_crossed == Y0){
			p->y = geometry->y;
		}
		else if(p->surface_crossed == Y1){
			p->y = 0;
		}
		else if(p->surface_crossed == Z0){
			p->z = geometry->z;
		}
		else if(p->surface_crossed == Z1){
			p->z = 0;
		}
	}

	return;
}

void cross_surface_loc(int mype, int n_workers, double d_proc,
		Buffer *buff, Geometry *geometry, 
		Geometry *geometry_local, Particle *p)
{
	// test if particle exited external surface
	if (exterior_surface(mype, n_workers, p->surface_crossed))
	{
		//printf("%dproc exterior external surface\n", mype);
		cross_surface(geometry_local, p);
	}
	else
	{
		// convert loc coord to global
		int mype_dest;
		int x_lvl_src, y_lvl_src, z_lvl_src;
		int x_lvl_dest, y_lvl_dest, z_lvl_dest;	
		double x_global, y_global, z_global;
		double x_loc_dest, y_loc_dest, z_loc_dest;

		x_lvl_src = x_lvl_from_mype(mype, n_workers);
		y_lvl_src = y_lvl_from_mype(mype, n_workers);
		z_lvl_src = z_lvl_from_mype(mype, n_workers);

		x_global = coord_global(p->x, d_proc, x_lvl_src);
		y_global = coord_global(p->y, d_proc, y_lvl_src);
		z_global = coord_global(p->z, d_proc, z_lvl_src);

		// check that coord_global are not out of bounds
		if (out_of_bounds(x_global, y_global, z_global, geometry))
		{
			//printf("%dproc out of bounds\n", mype);
			cross_surface(geometry_local, p);
		}
		else 
		{
			// convert global to loc of dest
			x_lvl_dest = lvl_from_global(x_global, d_proc);
			y_lvl_dest = lvl_from_global(y_global, d_proc);
			z_lvl_dest = lvl_from_global(z_global, d_proc);

			mype_dest = mype_from_lvls(x_lvl_dest, y_lvl_dest, 
					z_lvl_dest, n_workers);

			x_loc_dest = coord_loc(x_global, d_proc, x_lvl_dest);
			y_loc_dest = coord_loc(y_global, d_proc, y_lvl_dest);
			z_loc_dest = coord_loc(z_global, d_proc, z_lvl_dest);

			assert(0 <= x_loc_dest && x_loc_dest <= geometry_local->x);
			assert(0 <= y_loc_dest && y_loc_dest <= geometry_local->y);
			assert(0 <= z_loc_dest && z_loc_dest <= geometry_local->z);

			add_particle_to_buffer(buff, x_loc_dest, y_loc_dest, 
					z_loc_dest, p, mype_dest);

			// set particle off. only on this proc
			p->alive = FALSE;
		}
	}
}

void collision(Material *material, Bank *fission_bank, double nu, Particle *p)
{
	int n;
	int i = 0;
	double prob = 0.0;
	double cutoff;
	Nuclide nuc = {0, 0, 0, 0, 0};

	// Cutoff for sampling nuclide
	cutoff = rn()*material->xs_t;

	// Sample which nuclide particle has collision with
	while(prob < cutoff){
		nuc = material->nuclides[i];
		prob += nuc.atom_density*nuc.xs_t;
		i++;
	}

	// Cutoff for sampling reaction
	cutoff = rn()*nuc.xs_t;

	// Sample fission
	if(nuc.xs_f > cutoff){

		// Sample number of fission neutrons produced
		if(rn() > nu - (int)nu){
			n = nu;
		}
		else{
			n = nu + 1;
		}

		// Sample n new particles from the source distribution but at the current
		// particle's location
		if(fission_bank->n+n >= fission_bank->sz){
			fission_bank->resize(fission_bank, 2*(fission_bank->n+n));
		}
		for(i=0; i<n; i++){
			sample_fission_particle(&(fission_bank->p[fission_bank->n]), p);
			fission_bank->n++;
		}
		p->alive = FALSE;
	}

	// Sample absorption (disappearance)
	else if(nuc.xs_a > cutoff){
		p->alive = FALSE;
	}

	// Sample scattering
	else{
		p->mu = rn()*2 - 1;
		p->phi = rn()*2*PI;
		p->u = p->mu;
		p->v = sqrt(1 - p->mu*p->mu) * cos(p->phi);
		p->w = sqrt(1 - p->mu*p->mu) * sin(p->phi);
	}

	return;
}

// lvls from mype
int x_lvl_from_mype(int mype, int n_workers)
{
	return mype % (int) cbrt(n_workers);  
}

int y_lvl_from_mype(int mype, int n_workers)
{
	int z = z_lvl_from_mype(mype, n_workers);  
	int x = x_lvl_from_mype(mype, n_workers);
	int n_workers_axis = (int) cbrt(n_workers);
	int inv_y_lvl = (mype - x - z*n_workers_axis*n_workers_axis) / n_workers_axis;

	return n_workers_axis - 1 - inv_y_lvl;
}

int z_lvl_from_mype(int mype, int n_workers)
{
	int n_workers_axis = (int) cbrt(n_workers);     
	int n2 = n_workers_axis * n_workers_axis;

	return (int) floor(mype / n2);
}

// lvls from global coord
int lvl_from_global(double coord_global, double d_proc)
{
	return (int) floor(coord_global / d_proc);
}

// mype from lvls
int mype_from_lvls(int x_lvl, int y_lvl, int z_lvl, int n_workers)
{
	int n_workers_axis = (int) cbrt(n_workers);
	int n2 = n_workers_axis * n_workers_axis;
	int base = z_lvl * n2;
	int y_inv = (n_workers_axis-1) - y_lvl;
	return base + x_lvl + n_workers_axis*y_inv; 
}

// coordinate transforms
double coord_global(double coord_loc, double d_proc, int coord_lvl)
{
	return coord_lvl * d_proc + coord_loc;
}

double coord_loc(double coord_global, double d_proc, int coord_lvl)
{
	return coord_global - coord_lvl * d_proc;
}

int out_of_bounds(double x_global, double y_global, double z_global,
		Geometry *geometry_global)
{
	if (x_global < 0 || x_global > geometry_global->x)
		return TRUE;
	else if (y_global < 0 || y_global > geometry_global->y)
		return TRUE;
	else if (z_global < 0 || z_global > geometry_global->z)
		return TRUE;
	else 
		return FALSE;
}
