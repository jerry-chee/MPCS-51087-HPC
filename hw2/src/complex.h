/* Jerry Chee
 * MPCS 51807 hw2
 * complex.h
 */

#ifndef COMPLEX_H
#define COMPLEX_H

#include <stdio.h>
#include <stdlib.h>

/* struct for each complex value represented as a pair of doubles */
struct Complex {
	double x; // 'real'
	double y; // 'complex'
};


/* creates instance of struct */
struct Complex init(double x_coord, double y_coord);

/* length of vector from origin */
double dist_sq(struct Complex z);

/* square of a complex number */
struct Complex sq(struct Complex z);

/* sum of 2 complex numbers */ 
struct Complex add (struct Complex z1, struct Complex z2);

/* prints out a complex number */ 
void print_complex(struct Complex z, char* name);

#endif // COMPLEX_H
