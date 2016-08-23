/*
an alternative 2d simpson's rule that calculates weights explicitly.
this may be easier to extend to n-dimenions.
 */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<omp.h>

#define NDIM 5

double f( double x[] )
{
	return exp(-x[0]*x[0] -x[1]*x[1] -x[2]*x[2]  -x[3]*x[3] -x[4]*x[4]);
}

int main(int argc, char **argv)
{
  // inputs
  int i,ii,j,jj,k,kk,l,ll,m,mm;
	int n;
	double start, stop;
  double a, b;
  double delta, h, integral = 0.0;
  double x[NDIM];
  double w[3];
  assert (argc == 2);
  n = atoi(argv[1]); /* "resolution", ie number of points at which to sample f */

  a = 0.0; b = 1.0; /* upper and lower bounds of the integral */

  delta =(b-a)/n;   /* "grid spacing" -- fixed interval between function sample points */

  h = delta / 2.0;  /* h is used for convenience to find half distance between adjacent samples */
  integral = 0.0;   /* the accumulated integral estimate */

  /* three point weights that define Simpson's rule */
  w[0] = h/3.; w[1] = 4.*w[0]; w[2] = w[0];

	start = omp_get_wtime();
	omp_set_nested(1);
#pragma omp parallel for private(x,ii,j,jj,k,kk,l,ll,m,mm) schedule(static) reduction(+:integral) 
	for (i = 0; i < n; i++){
		for (ii=0;ii<3;ii++){
			x[0] = a + i * delta + ii * h;
			for( j = 0; j < n; j++ ){
				for (jj=0;jj<3;jj++){
					x[1] = a + j * delta + jj * h;
					for (k = 0; k < n; k++){
						for(kk=0;kk<3;kk++){
							x[2] = a + k * delta + kk * h;
							for (l = 0; l < n; l++){
								for(ll=0;ll<3;ll++){
									x[3] = a + l * delta + ll * h;
									for (m = 0; m < n; m++){
										for(mm=0;mm<3;mm++){
											x[4] = a + m * delta + mm * h;
											integral += w[ii]*w[jj]*w[kk]*w[ll]*w[mm]*f(x);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	stop = omp_get_wtime();
	printf("simpson openmp integral:%lf, n:%d, time:%lf\n", 
			integral, n, stop-start);

	return 0;
}
