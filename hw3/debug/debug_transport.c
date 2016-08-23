#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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
int mype_from_lvl(int x_lvl, int y_lvl, int z_lvl, int n_workers)
{
	int n_workers_axis = (int) cbrt(n_workers);
	int n2 = n_workers_axis * n_workers_axis;
	int base = z_lvl * n2;
	int y_inv = (n_workers_axis-1) - y_lvl;
	printf("n_workers_axis:%d, base:%d, y_inv:%d\n", 
			n_workers_axis, base, y_inv);

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

/*
	 double coord_tranform()
	 {
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

	 x_lvl_dest = lvl_from_global(x_global, d_proc);
	 y_lvl_dest = lvl_from_global(y_global, d_proc);
	 z_lvl_dest = lvl_from_global(z_global, d_proc);

	 mype_dest = mype_from_lvls(x_lvl_dest, y_lvl_dest,
	 z_lvl_dest, n_workers);

	 x_loc_dest = coord_loc(x_global, d_proc, x_lvl_dest);
	 y_loc_dest = coord_loc(y_global, d_proc, y_lvl_dest);
	 z_loc_dest = coord_loc(z_global, d_proc, z_lvl_dest);

	 add_particle(buff, x_loc_dest, y_loc_dest, z_loc_dest, p, mype_dest);
	 }
	 */
