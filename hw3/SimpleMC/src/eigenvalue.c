#include "header.h"

// given properly sized things
void run_eigenvalue(Parameters *parameters, Geometry *geometry, Material *material, Bank *source_bank, Bank *fission_bank, Tally *tally, double *keff)
{
  int i_b; // index over batches
  int i_a = -1; // index over active batches
  int i_g; // index over generations
  unsigned long i_p; // index over particles
  double keff_gen = 1; // keff of generation
  double keff_batch; // keff of batch
  double keff_mean; // keff mean over active batches
  double keff_std; // keff standard deviation over active batches

  // Loop over batches
  for(i_b=0; i_b<parameters->n_batches; i_b++){

    keff_batch = 0;

    // Turn on tallying and increment index in active batches
    if(i_b >= parameters->n_batches - parameters->n_active){
      i_a++;
      if(parameters->tally == TRUE){
        tally->tallies_on = TRUE;
      }
    }

    // Loop over generations
    for(i_g=0; i_g<parameters->n_generations; i_g++){

      // Set RNG stream for tracking
      set_stream(STREAM_TRACK);

      // Loop over particles
      for(i_p=0; i_p<parameters->n_particles; i_p++){

	// Set seed for particle i_p by skipping ahead in the random number
	// sequence stride*(total particles simulated) numbers from the initial
	// seed. This allows for reproducibility of the particle history.
        rn_skip((i_b*parameters->n_generations + i_g)*parameters->n_particles + i_p);

        // Transport the next particle
        transport(parameters, geometry, material, source_bank, fission_bank, tally, &(source_bank->p[i_p]));
      }

      // Switch RNG stream off tracking
      set_stream(STREAM_OTHER);
      rn_skip(i_b*parameters->n_generations + i_g);

      // Calculate generation k_effective and accumulate batch k_effective
      keff_gen = (double) fission_bank->n / source_bank->n;
      keff_batch += keff_gen;

      // Sample new source particles from the particles that were added to the
      // fission bank during this generation
      synchronize_bank(source_bank, fission_bank);
    }

    // Calculate k effective
    keff_batch /= parameters->n_generations;
    if(i_a >= 0){
      keff[i_a] = keff_batch;
    }
    calculate_keff(keff, &keff_mean, &keff_std, i_a+1);

    // Tallies for this realization
    if(tally->tallies_on == TRUE){
      if(parameters->write_tally == TRUE){
        write_tally(tally, parameters->tally_file);
      }
      reset_tally(tally);
    }

    // Status text
    print_status(i_a, i_b, keff_batch, keff_mean, keff_std);
  }

  // Write out keff
  if(parameters->write_keff == TRUE){
    write_keff(keff, parameters->n_active, parameters->keff_file);
  }

  return;
}

void run_eigenvalue_parallel(int mype, int n_workers, 
		double d_proc, int buff_lim, 
		Parameters *parameters, Parameters *parameters_loc, 
		Geometry *geometry, Geometry *geometry_local, 
		Material *material, Bank *source_bank, Bank *fission_bank, 
		Tally *tally_loc, double *keff)
{
	Buffer *buff;
  int i_b; // index over batches
  int i_a = -1; // index over active batches
  int i_g; // index over generations
  unsigned long i_p; // index over particles
  double keff_gen = 1; // keff of generation
  double keff_batch; // keff of batch
  double keff_mean; // keff mean over active batches
  double keff_std; // keff standard deviation over active batches

  // Loop over batches
  for(i_b=0; i_b<parameters->n_batches; i_b++){

    keff_batch = 0;

    // Turn on tally_locing and increment index in active batches
		//if (mype != n_workers) {
		if(i_b >= parameters_loc->n_batches - parameters_loc->n_active){
			i_a++;
			if(parameters_loc->tally == TRUE){
				tally_loc->tallies_on = TRUE;
			}
		}
		//}

		// Loop over generations
		for(i_g=0; i_g<parameters_loc->n_generations; i_g++){
			// printf("%dproc, i_g:%d\n", mype, i_g);

			// Set RNG stream for tracking
			set_stream(STREAM_TRACK);

			// Distribute particles
			//printf("%dproc, i_g:%d, source_scatter\n", mype, i_g);
			MPI_Barrier(MPI_COMM_WORLD);
			source_scatter(mype, n_workers, d_proc, 
					parameters, source_bank, buff_lim);
			MPI_Barrier(MPI_COMM_WORLD);

			while (1)
			{
				//printf("%dproc, i_g:%d, init_buffer\n", mype, i_g);
				buff = init_buffer(-1);

				// only transport particles if have them
				// Loop over particles
				for(i_p=0; i_p<source_bank->n; i_p++){
					if (mype != n_workers)
					{
						assert(0 <= source_bank->p[i_p].x && 
								source_bank->p[i_p].x <= geometry_local->x);
						assert(0 <= source_bank->p[i_p].y && 
								source_bank->p[i_p].y <= geometry_local->y);
						assert(0 <= source_bank->p[i_p].z && 
								source_bank->p[i_p].z <= geometry_local->z);
					}

					// Set seed for particle i_p by skipping ahead in the random number
					// sequence stride*(total particles simulated) numbers from the initial
					// seed. This allows for reproducibility of the particle history.
					rn_skip((i_b*parameters_loc->n_generations + i_g)*
							parameters_loc->n_particles + i_p);

					// Transport the next particle
					//printf("%dproc, i_g:%d, i_p:%ld, transport_loc\n", mype, i_g, i_p);
					transport_loc(mype, n_workers, d_proc, parameters_loc, 
							geometry, geometry_local, material, 
							source_bank, fission_bank, buff, 
							tally_loc, &(source_bank->p[i_p]));
				}
				MPI_Barrier(MPI_COMM_WORLD);					

				// Particle swap among procs
				if (comm_ready(mype, n_workers, buff))	
				{
					//printf("%dproc, i_g:%d, comm_routine\n", mype, i_g);
					comm_routine(mype, n_workers, buff, source_bank);
				}
				else
				{
					//printf("%dproc, i_g:%d, free_buffer and break\n", mype, i_g);
					free_buffer(buff);
					break;
				}
			}

			// Gather particles
			//printf("%dproc, i_g:%d, fission_gather\n", mype, i_g);
			MPI_Barrier(MPI_COMM_WORLD);
			fission_gather(mype, n_workers, d_proc, 
					parameters, parameters_loc, fission_bank, buff_lim);
			MPI_Barrier(MPI_COMM_WORLD);

			// Switch RNG stream off tracking
			set_stream(STREAM_OTHER);
			rn_skip(i_b*parameters_loc->n_generations + i_g);

			if (mype == n_workers)
			{
				// Calculate generation k_effective and accumulate batch k_effective
				//printf("master, i_g:%d calculating gen k_eff & etc\n", i_g);
				keff_gen = (double) fission_bank->n / source_bank->n;
				keff_batch += keff_gen;

				// Sample new source particles from the particles that were added to the
				// fission bank during this generation
				//printf("master, i_g;%d synchronize_bank\n", i_g);
				synchronize_bank(source_bank, fission_bank);
			}
		}

		// Calculate k effective
		if (mype == n_workers)
		{
			//printf("master, calc k_eff\n");
			keff_batch /= parameters->n_generations;
			if(i_a >= 0){
				keff[i_a] = keff_batch;
			}
			calculate_keff(keff, &keff_mean, &keff_std, i_a+1);
		}

		// Tallies for this realization
		//printf("%dproc print tallies\n", mype);
		if(tally_loc->tallies_on == TRUE){
			if(parameters_loc->write_tally == TRUE){
				write_tally_parallel(mype, n_workers, 
						parameters, tally_loc, parameters->tally_file);
			}
			reset_tally(tally_loc);
		}

		// Status text
		if (mype == n_workers)
		{
			//printf("master print status\n"); 
			print_status(i_a, i_b, keff_batch, keff_mean, keff_std);
		}
	}

	// Write out keff
	if (mype == n_workers)
	{
		if(parameters->write_keff == TRUE){
			write_keff(keff, parameters->n_active, parameters->keff_file);
		}
	}
	return;
}

void synchronize_bank(Bank *source_bank, Bank *fission_bank)
{
	unsigned long i, j;
	unsigned long n_s = source_bank->n;
	unsigned long n_f = fission_bank->n;

	// If the fission bank is larger than the source bank, randomly select
	// n_particles sites from the fission bank to create the new source bank
	if(n_f >= n_s){

		// Copy first n_particles sites from fission bank to source bank
		memcpy(source_bank->p, fission_bank->p, n_s*sizeof(Particle));

		// Replace elements with decreasing probability, such that after final
		// iteration each particle in fission bank will have equal probability of
		// being selected for source bank
		for(i=n_s; i<n_f; i++){
			j = rni(0, i+1);
			if(j<n_s){
				memcpy(&(source_bank->p[j]), &(fission_bank->p[i]), sizeof(Particle));
			}
		}
	}

	// If the fission bank is smaller than the source bank, use all fission bank
	// sites for the source bank and randomly sample remaining particles from
	// fission bank
	else{

		// First randomly sample particles from fission bank
		for(i=0; i<(n_s-n_f); i++){
			j = rni(0, n_f);
			memcpy(&(source_bank->p[i]), &(fission_bank->p[j]), sizeof(Particle));
		}

		// Fill remaining source bank sites with fission bank
		memcpy(&(source_bank->p[n_s-n_f]), fission_bank->p, n_f*sizeof(Particle));
	}

	fission_bank->n = 0;

	return;
}

void calculate_keff(double *keff, double *mean, double *std, int n)
{
	int i;

	*mean = 0;
	*std = 0;

	// Calculate mean
	for(i=0; i<n; i++){
		*mean += keff[i];
	}
	*mean /= n;

	// Calculate standard deviation
	for(i=0; i<n; i++){
		*std += pow(keff[i] - *mean, 2);
	}
	*std = sqrt(*std/(n-1));

	return;
}
