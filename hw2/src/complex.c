/* Jerry Chee
 * MPCS 51807 hw2
 * Complex add, sq, etc functions
 */

#include <stdlib.h>
#include <stdio.h>
#include "complex.h"


struct Complex init(double x_coord, double y_coord)
{
	struct Complex z_out;
	z_out.x = x_coord;
	z_out.y = y_coord;

	return z_out;
} 

/* length of vector from origin */
double dist_sq(struct Complex z)
{
	return z.x*z.x + z.y*z.y;
}

/* square of a complex number */
struct Complex sq(struct Complex z)
{
	double a = z.x;
	double b = z.y;

	struct Complex z_out;
	z_out.x = a*a - b*b;
	z_out.y = 2 * a * b;

	return z_out;
}

/* sum of 2 complex numbers */
struct Complex add(struct Complex z1, struct Complex z2)
{
	struct Complex z_out;
	z_out.x = z1.x + z2.x;
	z_out.y = z1.y + z2.y;

	return z_out;
}

/* prints out a complex number */
void print_complex(struct Complex z, char* name)
{
	printf("%s(%lf, %lf)\n", name, z.x, z.y);
}
