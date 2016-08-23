/* Jerr Chee
 * hw1
 * Problem 2
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* high q function */
void high_q(long int n)
{
	long int data = 100;
	for (int i=0; i<n; i++)
	{
		data = (data * data) / data;
	}
	//printf("data: %ld\n", data);
}

/* timing function */
double timer(void (*f)(long int), long int n)
{
	clock_t start, end;
	double cpu_time_used;

	start = clock();
	f(n);
	end = clock();
	cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;

	return cpu_time_used;
}

int main(int argc, char** argv)
{
	int exp;
	long int coeff;
	long n;
	double t;

	/* open file pointer */
	FILE* fp = fopen("./p2.csv", "w");

	/* print header */
	fprintf(fp, "n,time\n"); 
	
	for (exp=6; exp<=8; exp++)
	{
		for (coeff=0; coeff<=8; coeff+=2)
		{
			/* set the the exponenet of n */
			n = (long int) pow(10, exp); 

			/* if coeff != 0, then multiply coeff by n */
			if (coeff != 0)
			{
				n = coeff * n;
			}			
			
			/* call the timing function */
			t = timer(high_q, n);

			/* print to csv */
			fprintf(fp, "%ld,%lf\n", n, t);
		}
	}

	fclose(fp);
}
