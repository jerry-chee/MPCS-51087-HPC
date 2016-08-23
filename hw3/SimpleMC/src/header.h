#ifndef HEADER
#define HEADER

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include<sys/time.h>
#include<math.h>
#include<float.h>
#include<unistd.h>
#include<string.h>
#include<strings.h>
#include<assert.h>
#include"mpi.h"

#define TRUE 1
#define FALSE 0

// Constants
#define PI 3.1415926535898
#define D_INF DBL_MAX

// Geometry boundary conditions
#define VACUUM 0
#define REFLECT 1
#define PERIODIC 2

// Reaction types
#define TOTAL 0
#define ABSORPTION 1
#define SCATTER 2
#define FISSION 3

// Surfaces
#define X0 0
#define X1 1
#define Y0 2
#define Y1 3
#define Z0 4
#define Z1 5

// RNG streams
#define N_STREAMS 3
#define STREAM_INIT 0
#define STREAM_TRACK 1
#define STREAM_OTHER 2

typedef struct Parameters_{
  unsigned long long seed; // RNG seed
  unsigned long n_particles; // number of particles
  int n_batches; // number of batches
  int n_generations; // number of generations per batch
  int n_active; // number of active batches
  int bc; // boundary conditions
  int n_nuclides; // number of nuclides in material
  int tally; // whether to tally
  int n_bins; // number of bins in each dimension of mesh
  double nu; // average number of fission neutrons produced
  double xs_a; // absorption macro xs
  double xs_s; // scattering macro xs
  double xs_f; // fission macro xs
  double gx; // geometry size in x
  double gy; // geometry size in y
  double gz; // geometry size in z
  int write_tally; // whether to output tallies
  int write_keff; // whether to output keff
  char *tally_file; // path to write tallies to
  char *keff_file; // path to write keff to
} Parameters;

typedef struct Particle_{
  int alive;
  double mu; // cosine of polar angle
  double phi; // azimuthal angle
  double u; // direction
  double v;
  double w;
  double x; // position
  double y;
  double z;
  int surface_crossed;
} Particle;

typedef struct Geometry_{
	int bc;
	double x;
	double y;
	double z;
} Geometry;

typedef struct Nuclide_{
	double xs_f; // fission micro xs
	double xs_a; // absorption micro xs
	double xs_s; // scattering micro xs
	double xs_t; // total micro xs
	double atom_density; // atomic density of nuclide in material
} Nuclide;

typedef struct Material_{
	double xs_f; // fission macro xs
	double xs_a; // absorption macro xs
	double xs_s; // scattering macro xs
	double xs_t; // total macro xs
	int n_nuclides;
	Nuclide *nuclides;
} Material;

typedef struct Tally_{
	int tallies_on; // whether tallying is currently turned on
	int n; // mumber of grid boxes in each dimension 
	double dx; // grid spacing
	double dy;
	double dz;
	double *flux;
} Tally;

typedef struct Bank_{
	unsigned long n; // number of particles
	unsigned long sz; // size of bank
	Particle *p; // particle array
	void (*resize)(struct Bank_ *b, long int newsize);
} Bank;

typedef struct Buffer_{
	int mype_dest;
	Bank *bank;
	struct Buffer_ *next;
	struct Buffer_ *last;
} Buffer;

// io.c function prototypes
void parse_parameters(Parameters *parameters);
void parse_parameters_loc(Parameters *parameters_loc);
void read_CLI(int argc, char *argv[], Parameters *parameters);
void print_error(char *message);
void print_parameters(Parameters *parameters);
void border_print(void);
void fancy_int(long a);
void center_print(const char *s, int width);
void print_status(int i_a, int i_b, double keff_batch, double keff_mean, double keff_std);
void init_output(Parameters *parameters);
void write_tally(Tally *t, char *filename);
void write_tally_parallel(int mype, int n_workers, 
		Parameters *parameters, Tally *taly_loc,
		char *filename);
void write_tally_loc(int mype, Tally *tally_loc,
		int i_loc, int j_loc,
		int mype_dest, char *filename);
void write_keff(double *keff, int n, char *filename);

// utils.c funtion prototypes
double timer(void);
int interior(int mype, int n_workers);
int exterior_surface(int mype, int n_workers, int surface);

// prng.c function prototypes
double rn(void);
int rni(int min, int max);
void set_stream(int rn_stream);
void set_initial_seed(unsigned long long rn_seed0);
void rn_skip(long long n);

// initialize.c function prototypes
Parameters *init_parameters(void);
Parameters *init_parameters_loc(Parameters *parameters, int n_workers);
Geometry *init_geometry(Parameters *parameters);
Tally *init_tally(Parameters *parameters);
Material *init_material(Parameters *parameters);
Bank *init_fission_bank(Parameters *parameters);
Bank *init_source_bank(Parameters *parameters, Geometry *geometry);
Bank *init_bank(unsigned long n_particles);
Buffer *init_buffer(int mype_dest);
void sample_source_particle(Geometry *geometry, Particle *p);
void sample_fission_particle(Particle *p, Particle *p_old);
void resize_particles(Bank *b, long int newsize);
void print_particle(Particle *p);
void print_bank(Bank *b, int p_true);
void free_bank(Bank *b);
void free_buffer(Buffer* buff);
void free_material(Material *material);
void free_tally(Tally *tally);

// transport.c function prototypes
void transport(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, Particle *p);
void transport_loc(int mype, int n_workers, double d_proc, Parameters *parameters_loc, Geometry *geometry, Geometry *geometry_local, Material *material, Bank *source_bank, Bank *fission_bank, Buffer *buff, Tally *tally, Particle *p);
double distance_to_boundary(Geometry *geometry, Particle *p);
double distance_to_collision(Material *material);
void cross_surface(Geometry *geometry, Particle *p);
void cross_surface_loc(int mype, int n_workers, double d_proc, Buffer *buff, Geometry *geometry, Geometry *geometry_local, Particle *p);
void collision(Material *material, Bank *fission_bank, double nu, Particle *p);
int x_lvl_from_mype(int mype, int n_workers);
int y_lvl_from_mype(int mype, int n_workers);
int z_lvl_from_mype(int mype, int n_workers);
int lvl_from_global(double coord_global, double d_proc);
int mype_from_lvls(int x_lvl, int y_lvl, int z_lvl, int n_workers);
double coord_global(double coord_loc, double d_proc, int coord_lvl);
double coord_loc(double coord_global, double d_proc, int coord_lvl);
int out_of_bounds(double x_global, double y_global, double z_global,
		Geometry *geometry_global);

// eigenvalue.c function prototypes
void run_eigenvalue(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, double *keff);
void run_eigenvalue_parallel(int mype, int n_workers, double d_proc, int buff_lim, Parameters *parameters, Parameters *parameters_loc, Geometry *geometry, Geometry *geometry_loc, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally_loc, double *keff);

void synchronize_bank(Bank *source_bank, Bank *fission_bank);
void calculate_keff(double *keff, double *mean, double *std, int n);

// tally.c function prototypes
void score_tally(Parameters *parameters, Material *material, Tally *tally, Particle *p);
void reset_tally(Tally *tally);

// parallel.c function prototypes
void declare_MPI_Particle(void);
void add_buffer(Buffer *buff, int mype_dest);
Buffer *search_buffer(Buffer *buff, int mype_dest);
void particle_deep_copy_coord(double x_loc_dest, double y_loc_dest,
		double z_loc_dest, Particle *p_to, Particle *p_from);
void particle_deep_copy(Particle *p_to, Particle *p_from);
void add_particle_to_buffer(Buffer *buff, double x_loc_dest, double y_loc_dest, 
		double z_loc_dest, Particle *p, int mype_dest);
void add_particles_to_bank(Bank *bank, Particle *ps, long int count);
void print_buffer(Buffer *buff);
void print_single_buffer(Buffer *buff);
void source_scatter(int mype, int n_workers, double d_proc, 
		Parameters *parameters, Bank *source_bank, int buff_lim);
void fission_gather(int mype, int n_workers, double d_proc, 
		Parameters *parameters, Parameters *parameters_loc, 
		Bank *fission_bank, int buff_lim);
int comm_ready(int mype, int n_workers, Buffer *buff);
void comm_routine(int mype, int n_workers, Buffer *buff, 
		Bank *source_bank);
#endif
