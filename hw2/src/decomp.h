/* Jerry Chee
 * MPCS 51807 hw2
 * static_decomp.h
 */

#ifndef STATIC_DECOMP
#define STATIC_DECOMP

/* Global var for criterion for mandelbrot divergence 
 * and maximum number of iterations*/
extern int r_max;
extern int r_max_sq;
extern int iter_max; 
extern int x_min; 
extern int y_min; 

/* struct for resizable array (reArr) in C */
typedef struct {  
	int* array;
	int used;
	int size;
} reArr;

/* init function for resizable array (reArr) */
void initArr(reArr* arr, int initial_size);

/* insert element, resize if reach size limit */
void insertArr(reArr* arr, int x);

/* free reArr */
void freeArr(reArr* arr);
/* returns iterations needed to break mandelbrot condition or iter_max */

/* prints reArr out to file, every n elements per line */
void printArr(FILE* fp, reArr* arr, int n);

int mandelbrot_count(double x_coord, double y_coord); 

/* prints mandelbrot_count to file pointer */
/* timing = 1 to time and not print out to file */
void mandelbrot_save(reArr* data, double x_coord, double y_coord);

/* timing function aggregates all times across procs */
void time_master(void (*f)(int, int, int, int, reArr*, int), char* type, 
		int n, int k, int nprocs, int mype, int timing);

/* compute mandelbrot set with static data decomposition */
/* n: number of ways to divide window, */
/* k: number of rows in window in a dynamic job */
/* returns length of int* where data stored */
void static_decomp(int n, int k, int nprocs, int mype, reArr* data, int timing);

/* function which each worker rungs */
void dynamic_job(int mype, reArr* data, double step,
		int x_start, int x_end,
		int y_start, int y_end);

/* compute mandelbrot set with dynamic decomposition */
void dynamic_decomp(int n, int k, int nprocs, int mype, reArr* data, int timing);

#endif
