#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>

#define NDIM 5

double uran(double a, double b, unsigned int *seed);
double f(double x[]);

int main(int argc, char **argv){
  int i,j,n;
  double a, b;
	double t0, t1;
  double vol=1.0;
	double int_f = 0.0;
  double x[NDIM];

  a = 0.0; b = 1.0;
  n = atoi(argv[1]);

	t0 = omp_get_wtime();
#pragma omp parallel for private(j, x) schedule(static) reduction(+:int_f)
//#pragma ordered reduction(+:int_f) 
 	for (i = 1; i <= n; ++i){

		//if (omp_get_thread_num() == 0)
		//	printf("i:%d\n", i);

    /* get one uniformly distributed random number per coordinate */
		unsigned int seed = 1;
		for (j = 0; j < NDIM; ++j){
			//*seed = omp_get_thread_num();
			x[j] = uran(a,b, &seed);    /* NDIM random numbers at which to evaluate function */
		}

		int_f +=  f(x)/n;      /* estimate of average value of function */
	}

	for (i=0;i<NDIM;++i)
		vol *= (b-a);         /* compute n-dimensional volume */

	//int_f *= vol;           /* scale by n-dimensional volume */
	t1 = omp_get_wtime();

	printf("mc openmp integral:%lf, n:%d time:%lf\n", int_f, n, t1-t0);
	return 0;
}

double uran(double a, double b, unsigned int *seed){
	return rand_r(seed) / (RAND_MAX + 1.0) * (b - a) + a;
}

double f(double x[]){
	return exp(-x[0]*x[0] -x[1]*x[1] -x[2]*x[2]  -x[3]*x[3] -x[4]*x[4]); 
}

