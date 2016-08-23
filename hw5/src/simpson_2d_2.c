/*
an alternative 2d simpson's rule that calculates weights explicitly.
this may be easier to extend to n-dimenions.
 */


#include<stdio.h>
#include<math.h>
#include<assert.h>

double f( double x, double y )
{
  return sin(100*x)/x*sin(y);
}

int main(int argc, char **argv)
{
  // inputs
  int i,ii,j,jj,n;
  double a, b;
  double delta, h, integral = 0.0;
  double x,y;
  double w[3];
  assert (argc == 2);
  n = atoi(argv[1]); /* "resolution", ie number of points at which to sample f */

  a = .5; b = 1.5; /* upper and lower bounds of the integral */

  delta =(b-a)/n;   /* "grid spacing" -- fixed interval between function sample points */

  h = delta / 2.0;  /* h is used for convenience to find half distance between adjacent samples */
  integral = 0.0;   /* the accumulated integral estimate */

  /* three point weights that define Simpson's rule */
  w[0] = h/3.; w[1] = 4.*w[0]; w[2] = w[0];

	for (j = 0; j < n; j++){
		for (jj=0;jj<3;++jj){
			y = a + j * delta + jj * h;
			for( i = 0; i < n; i++ ){
				for (ii=0;ii<3;++ii){
					x = a + i * delta + ii * h;
					integral += w[ii]*w[jj]*f(x,y);
				}
			}
		}
	}

	printf("integral is: %lf\n", integral);

	return 0;
}
