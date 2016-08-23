/* Jerry Chee
 * hw1
 * Problem 1
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mkl.h>
#include <math.h>

/* global variable matric size */
int size = 256; // default

/* matrix allocation and free functions */
double** matrix_allocate(int m, int n)
{
	double* data = (double*) malloc(n*m * sizeof(double));
	double** M = (double**) malloc(m * sizeof(double));

	for (int i=0; i<m; i++)
		M[i] = &data[i*n];

	return M;
}

void matrix_free(double** M)
{
	free(M[0]);
	free(M);
}

/* matrix fill with a constant */
void matrix_fill(double** M, double value)
{
	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			M[i][j] = value;
		}
	}
}

/* matrix fill with continuous integers */
void matrix_cont_fill(double** M, int start)
{
	int val = start;
	
	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			M[i][j] = val;
			val++;
		}
	}
}

/* matrix compaire. assumes they are of same size */
int matrix_compare(double** A, double** B)
{
	int same = 1; // true
	
	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			if (A[i][j] != B[i][j])
			{
				same = 0; // false
				printf("(%d,%d) don't match\n",i,j);
			}			
		}
	}
	return same;
}

/* matrix print function */
void matrix_print(double** M)
{
	for (int i=0; i<size; i++)
	{
		for (int j=0; j<size; j++)
		{
			printf("%lf ", M[i][j]);
		}
		printf("\n");
	}
}

/* naive matrix multiplication of C = C+A*B*/
void mmx_naive(double** A, double** B, double** C, int n, int bs)
{
	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
		{
			for (int k=0; k<n; k++)
			{
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}		
}

/* blocking matrix multiplication */
void mmx_blocking(double** A, double** B, double** C, int n, int bs)
{
	/* number of blocks */
	int num_blocks = n / bs;		

	/* For moving over blocks */
	int i,j,k = 0;

	/* For moving over elements when multiplying blocks */
	int ib,jb,kb = 0;

	/* block offset in matrix */
	int i_offset = 0;
	int j_offset = 0;
	int k_offset = 0;

	for (i=0; i<num_blocks; i++)
	{
		/* find the i (col) offset */
		i_offset = i * bs;
		//printf("i_offset: %d\n", i_offset);

		for (j=0; j<num_blocks; j++)
		{
			/* find j (row) offset */
			j_offset = j * bs;
			//printf("\tj_offset: %d\n", j_offset);

			for (k=0; k<num_blocks; k++)
			{
				k_offset = k * bs; 
				/* now we multiply each matrix block */
				for (ib=i_offset; ib < (i_offset+bs); ib++)
				{
					for (jb=j_offset; jb < (j_offset+bs); jb++)
					{
						for (kb=k_offset; kb < (k_offset+bs); kb++)
						{	
							C[ib][jb] += A[ib][kb] * B[kb][jb];
						}
					}
				}
			}
		}
	}
}

/* matrix multiplication of submatrices within larger */
void mmx_subm(double** A, int i_a, int j_a,		
		double** B, int i_b, int j_b,
		double** C, int i_c, int j_c,
		int n)
{
	for (int i=i_a; i < (i_a+n); i++)
	{
		for (int j=j_b; j < (j_b+n); j++)
		{
			for (int k=0; k<n; k++)
			{
				C[i][j] += A[i][k+j_a] * B[k+i_b][j];
			}
		}
	}		
}

/* recursive helper */
void mmx_rhelper(double** A, int i_a, int j_a, 
		double** B, int i_b, int j_b,
		double** C, int i_c, int j_c,
		int n, int bs)
{
	if (n == bs)
	{
		mmx_subm(A, i_a, j_a,
			B, i_b, j_b,
			C, i_c, j_c,
			n);	
		/* print statements for debugging
		printf("multiplying\n");
		printf("n: %d \n\ti_a: %d, j_a %d\n", n, i_a, j_a);
		printf("\ti_b: %d, j_b: %d \n\ti_c:%d, j_c: %d\n", i_b, j_b, i_c, j_c);
		printf("n: %d\n", n);
		*/
	} else 
	{
		/*
		printf("n: %d \n\ti_a: %d, j_a %d\n", n, i_a, j_a);
		printf("\ti_b: %d, j_b: %d \n\ti_c:%d, j_c: %d\n", i_b, j_b, i_c, j_c);
		*/
		mmx_rhelper(A, i_a, j_a, 
				B, i_b, j_b, 
				C, i_c, j_c, 
				n/2, bs);
		mmx_rhelper(A, i_a, (j_a + n/2), 
				B, (i_b +n/2), j_b, 
				C, i_c, j_c, 
				n/2, bs);
		
		mmx_rhelper(A, i_a, j_a, 
				B, i_b, (j_b + n/2), 
				C, i_c, (j_c + n/2), 
				n/2,bs);
		mmx_rhelper(A, i_a, (j_a + n/2), 
				B, (i_b + n/2), (j_b + n/2), 
				C, i_c, (j_c + n/2), 
				n/2, bs);
		
		mmx_rhelper(A, (i_a + n/2), j_a, 
				B, i_b, j_b, 
				C, (i_c + n/2), j_c, 
				n/2, bs);
		mmx_rhelper(A, (i_a + n/2), (j_a + n/2), 
				B, (i_b + n/2), j_b, 
				C, (i_c + n/2), j_c, 
				n/2, bs);

		mmx_rhelper(A, (i_a + n/2), j_a, 
				B, i_b, (j_b + n/2), 
				C, (i_c + n/2), (j_c + n/2), 
				n/2, bs);
		mmx_rhelper(A, (i_a + n/2), (j_a + n/2), 
				B, (i_b + n/2), (j_b + n/2), 
				C, (i_c + n/2), (j_c + n/2), 
				n/2, bs);
	}
}

/* recursive matrix multiplication */
void mmx_recursive(double** A, double** B, double** C, int n, int bs)
{
	mmx_rhelper(A, 0, 0, B, 0, 0, C, 0, 0, n, bs);
}

/* BLAS matrix multiplication */
void mmx_blas(double** A, double** B, double** C, int n, int bs)
{
	/* make cblas function call */
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		n, n, n, 1, *A, n, *B, n, 1, *C, n);
}

/* timing function */
double mmx_time(void (*f)(double** , double**, double**, int, int),
		double** A, double** B, double** C, int n, int bs)
{
	clock_t start, end;
	double cpu_time_used;

	start = clock();
	f(A, B, C, n, bs); 
	end = clock();
	cpu_time_used = ((double) (end-start)) / CLOCKS_PER_SEC;

	return cpu_time_used;
}

/* testing */
void testing()
{
	/* test run */
	/*
	if (argc == 2) 
	{
		size = atoi(argv[1]);
		printf("Input matrix size: %d\n", size);  
	} else
	{
		printf("Default matrix size: %d\n", size);
	}		
	*/

	/* allocate matrices */
	double** A_n = matrix_allocate(size,size);
	double** B_n = matrix_allocate(size,size);
	double** C_n = matrix_allocate(size,size);

	double** A_b = matrix_allocate(size,size);
	double** B_b = matrix_allocate(size,size);
	double** C_b = matrix_allocate(size,size);

	double** A_r = matrix_allocate(size,size);
	double** B_r = matrix_allocate(size,size);
	double** C_r = matrix_allocate(size,size);

	/* fill with values */	
	matrix_cont_fill(A_n,1);
	matrix_cont_fill(B_n,2);

	matrix_cont_fill(A_b,1);
	matrix_cont_fill(B_b,2);

	matrix_cont_fill(A_r,1);
	matrix_cont_fill(B_r,2);

	/* matrix multiply functions */
	//mmx_naive(A_n, B_n, C_n, size, 0);
	//matrix_print(C_n);

	//putchar('\n');

	double t = mmx_time(mmx_blocking, A_b, B_b, C_b, size, size/4);
	printf("t: %lf\n", t);
	//matrix_print(C_b);

	//putchar('\n');

	//mmx_recursive(A_r, B_r, C_r, size, size/4);
	//matrix_print(C_r);

	int equality = matrix_compare(C_n, C_r);
	printf("equality=%d\n", equality);

	double t_b = mmx_time(mmx_blas, A_r, B_r, C_r, size, 0);
	printf("t_b: %lf\n", t_b);

	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	//	size, size, size, 1, *A_r, size, *B_r, size, 1, *C_r, size);

	/* free */
	matrix_free(A_n);
	matrix_free(B_n);
	matrix_free(C_n);

	matrix_free(A_b);
	matrix_free(B_b);
	matrix_free(C_b);
}

/* 1st plot, broad range of n */
void run_1()
{
	/* matrix variable */
	long int n; // size
	int bsB = 32; // block size 
	int bsR = 32; // min recursive problem size

	FILE* fp = fopen("plot1.csv", "w");

	/* matrix timing variables */
	double tN, tB, tR, tBLAS;

	/* print header */
	fprintf(fp, "n,naive,blocking,recursive,BLAS,\n");

	for (int i=5; i<=11; i++)
	{
		/* set size of n */
		n = pow(2,i);	

		/* allocate matrices */
		double** A = matrix_allocate(n,n);
		double** B = matrix_allocate(n,n);
		double** C = matrix_allocate(n,n);

		/* matrix  multiply and time */
		tN = mmx_time(mmx_naive, A, B, C, n, 0);
		tB = mmx_time(mmx_blocking, A, B, C, n, bsB);
		tR = mmx_time(mmx_recursive, A, B, C, n, bsR);		
		tBLAS = mmx_time(mmx_blas, A, B, C, n, 0);

		/* print this line in the csv */
		fprintf(fp, "%ld,%lf,%lf,%lf,%lf,\n", n, tN, tB, tR, tBLAS);

		/* matrix free */
		matrix_free(A);
		matrix_free(B);
		matrix_free(C);
	}
}

/* 2nd plot, blocking method fixed n and variable block size */
void run_2()
{
	int exp = 11;
	int n = pow(2, exp);
	int bs; 
	double t;

	FILE* fp = fopen("plot2.csv", "w");

	/* print header */
	fprintf(fp, "bs,blocking,\n");

	for (int i=1; i<=exp; i++)
	{
		/* set block size */
		bs = pow(2,i);

		/* allocate matrices */
		double** A = matrix_allocate(n,n);
		double** B = matrix_allocate(n,n);
		double** C = matrix_allocate(n,n);

		/* matrix multiply and time */
		t = mmx_time(mmx_blocking, A, B, C, n, bs);

		/* print line in csv */
		fprintf(fp, "%d,%lf,\n", bs, t);

		/* matrix free */
		matrix_free(A);
		matrix_free(B);
		matrix_free(C);
	}
}

/* 3rd plot, recursive method fixed n and variable min problem size */
void run_3()
{
	int exp = 11;
	int n = pow(2, exp);
	int bs; 
	double t;

	FILE* fp = fopen("plot3.csv", "w");

	/* print header */
	fprintf(fp, "bs,recursive,\n");

	for (int i=1; i<=exp; i++)
	{
		/* set block size */
		bs = pow(2,i);

		/* allocate matrices */
		double** A = matrix_allocate(n,n);
		double** B = matrix_allocate(n,n);
		double** C = matrix_allocate(n,n);

		/* matrix multiply and time */
		t = mmx_time(mmx_recursive, A, B, C, n, bs);

		/* print line in csv */
		fprintf(fp, "%d,%lf,\n", bs, t);

		/* matrix free */
		matrix_free(A);
		matrix_free(B);
		matrix_free(C);
	}
}


int main(int argc, char** argv)
{
	//run_1();
	//run_2();
	//run_3();
	/* test run */
	
	if (argc == 2) 
	{
		size = atoi(argv[1]);
		printf("Input matrix size: %d\n", size);  
	} else
	{
		printf("Default matrix size: %d\n", size);
	}		
	

	/* allocate matrices */
	double** A_n = matrix_allocate(size,size);
	double** B_n = matrix_allocate(size,size);
	double** C_n = matrix_allocate(size,size);

	double** A_b = matrix_allocate(size,size);
	double** B_b = matrix_allocate(size,size);
	double** C_b = matrix_allocate(size,size);

	double** A_r = matrix_allocate(size,size);
	double** B_r = matrix_allocate(size,size);
	double** C_r = matrix_allocate(size,size);

	/* fill with values */	
	matrix_cont_fill(A_n,1);
	matrix_cont_fill(B_n,2);

	matrix_cont_fill(A_b,1);
	matrix_cont_fill(B_b,2);

	matrix_cont_fill(A_r,1);
	matrix_cont_fill(B_r,2);

	/* matrix multiply functions */
	mmx_naive(A_n, B_n, C_n, size, 0);
	matrix_print(C_n);

	//putchar('\n');

	mmx_blocking(A_b, B_b, C_b, size, size/4);
	//double t = mmx_time(mmx_blocking, A_b, B_b, C_b, size, size/4);
	//printf("t: %lf\n", t);
	matrix_print(C_b);

	//putchar('\n');

	mmx_recursive(A_r, B_r, C_r, size, size/4);
	matrix_print(C_r);

	//int equality = matrix_compare(C_n, C_r);
	//printf("equality=%d\n", equality);

	//double t_b = mmx_time(mmx_blas, A_r, B_r, C_r, size, 0);
	//printf("t_b: %lf\n", t_b);

	//cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
	//	size, size, size, 1, *A_r, size, *B_r, size, 1, *C_r, size);

	/* free */
	matrix_free(A_n);
	matrix_free(B_n);
	matrix_free(C_n);

	matrix_free(A_b);
	matrix_free(B_b);
	matrix_free(C_b);
}
