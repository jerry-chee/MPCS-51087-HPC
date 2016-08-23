#include "header.h"

int main(int argc, char *argv[])
{
	// MPI Stuff
	MPI_Init(&argc, &argv);
	int mype;
	int nprocs;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	int master = nprocs - 1;;
	int n_workers = nprocs - 1;
	int n_workers_axis = (int) cbrt(n_workers);

	double epsilon = 0.0000000001;
	// nprocs 3rd deg power 
	assert(abs(floor(cbrt(n_workers)) - cbrt(n_workers)) < epsilon);; 
	assert(n_workers > 1);

// ============================================================================

	Parameters *parameters; // user defined parameters
  Parameters *parameters_loc; // local parameters
	Geometry *geometry; // homogenous cube geometry
	Geometry *geometry_loc; // local geometry
  Material *material; // problem material
  Bank *source_bank; // array for particle source sites
  Bank *fission_bank; // array for particle fission sites
  Tally *tally_loc; // scalar flux tally
	Tally *tally;
  double *keff; // effective multiplication factor
  double t1, t2; // timers

  // Get inputs: set parameters to default values, parse parameter file,
  // override with any command line inputs, and print parameters
  parameters = init_parameters();
	parse_parameters(parameters);
	parameters_loc = init_parameters_loc(parameters, n_workers);
	read_CLI(argc, argv, parameters);
  if (mype == master) 
		print_parameters(parameters);

  // Set initial RNG seed
  set_initial_seed(parameters->seed);
  set_stream(STREAM_INIT);

  // Create files for writing results to
  init_output(parameters);

  // Set up geometry
  geometry = init_geometry(parameters);
	geometry_loc = init_geometry(parameters_loc);

  // Set up material
  material = init_material(parameters);

  // Set up tallies
	//if (mype != master)
	tally_loc = init_tally(parameters_loc);
	// for serial version
	if (mype == master)
		tally = init_tally(parameters);

	// Create source bank and initial source distribution
	if (mype == master)
		source_bank = init_source_bank(parameters, geometry);
	else
		source_bank = init_bank((long int) (parameters->n_particles / n_workers));

	// Create fission bank
	if (mype == master)
		fission_bank = init_fission_bank(parameters);
	else 
		fission_bank = init_bank((long int) 2*parameters->n_particles / n_workers);

	// Set up array for k effective
	if (mype == master)
		keff = calloc(parameters->n_active, sizeof(double));

	// RUN TESTS
	if (mype == master)
	{
		center_print("SIMULATION", 79);
		border_print();
		printf("%-15s %-15s %-15s\n", "BATCH", "KEFF", "MEAN KEFF");
	}

	// Start time
	if (mype == master)
		t1 = timer();

	double d_proc = geometry->x / n_workers_axis;
	int buff_lim = 1000;
	
	run_eigenvalue_parallel(mype, n_workers, d_proc, buff_lim,
			parameters, parameters_loc, geometry, geometry_loc, material, 
			source_bank, fission_bank, tally_loc , keff);

	if (mype == master)
	{
		//run_eigenvalue(parameters, geometry, material, 
		//		source_bank, fission_bank, tally, keff);
	}

	// Stop time
	if (mype == master)
	{
		printf("master end timing\n");
		t2 = timer();
		printf("Simulation time: %f secs\n", t2-t1);
	}

	// Free memory
	if (mype == master)
		free(keff);
	if (mype == master)
		free_tally(tally);
	free_tally(tally_loc);
	free_bank(fission_bank);
	free_bank(source_bank);
	free_material(material);
	free(geometry);
	free(geometry_loc);
	free(parameters);
	free(parameters_loc);

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}

