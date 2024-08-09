#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "parameters.h"
#include "lib_sim.h"

#define RND gsl_rng_uniform_pos(rng) // generate number from Unif(0,1)

// Define propensity updating rule for SSD simulation
void propensity_update(double** propensity, double* propensity_sums, int cell_idx, int* wildtype_populations, int* ssd_populations, double degradation_rate, double diffusion_rate, double nucleus_control_factor, double target_population, double density) {
	/* Updates the propensity of reactions in a cell
    
    Inputs
    ------
	propensity, 2d array of propensity vectors for each cell
	propensity_sums, 1d array of propensity sums for each cell
	cell_idx, index of cell to update propensity array of
	wildtype_populations, 1d array of wildtype population size of each cell
	ssd_populations, 1d array of SSD population size of each cell
    degradation_rate, degradation rate
	diffusion_rate, diffusion rate
	nucleus_control_factor, nucleus conrol factor
    target_population, target population
	density, density */

	int cell_wildtype_population = wildtype_populations[cell_idx];
	int cell_ssd_population = ssd_populations[cell_idx];
	double replication_rate = fmax(0.0, degradation_rate + nucleus_control_factor * (target_population - cell_wildtype_population - density*cell_ssd_population));
	
	propensity[cell_idx][0] = degradation_rate * cell_wildtype_population;
	propensity[cell_idx][1] = replication_rate * cell_wildtype_population;
	propensity[cell_idx][2] = 2 * diffusion_rate * cell_wildtype_population;
	propensity[cell_idx][3] = degradation_rate * cell_ssd_population;
	propensity[cell_idx][4] = replication_rate * cell_ssd_population;
	propensity[cell_idx][5] = 2 * diffusion_rate * cell_ssd_population;

	propensity_sums[cell_idx] = 0.00;
	for (int i=0; i<6; ++i) {propensity_sums[cell_idx] += propensity[cell_idx][i];}
	return;
}

int main(int argc, char *argv[]) {
    int seed;
    const double log_site_std_mutation_rate = LOG_SITE_STD_MUTATION_RATE;
    const double degradation_rate = DEGRADATION_RATE;
	const double diffusion_rate = DIFFUSION_RATE;
    const double nucleus_control_factor = NUCLEUS_CONTROL_FACTOR;
    const double target_population = TARGET_POP;
	const double density = DENSITY;

	if (argc==2) {
		seed = atoi(argv[1]);
	} else {
		printf("argc = %d\n", argc);
		printf("Usage: %s  seed\n", argv[0]);
		return 0;
	}

	/* set up GSL RNG */
	gsl_rng* rng = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(rng, seed);
	/* end of GSL setup */

	// Generate replication rate and mutation rate for this run
    const long double site_std_mutation_rate = pow(10, - log_site_std_mutation_rate);

    // Change to directory which stores simulation results. If directory doesn't exist, create one. 
    if (chdir(DIR_NAME)) {
		mkdir(DIR_NAME, 0700); // Linux
		//mkdir(DIR_NAME); // Windows
		chdir(DIR_NAME);
    }

	// Set up file to save parameter values
	FILE *fp_parameters = fopen("parameters.txt", "w");
	fprintf(fp_parameters, "parameter,value\n");
	fprintf(fp_parameters, "cells,%d\n", CELLS);
	fprintf(fp_parameters, "log_site_std_mutation_rate,%e\n", LOG_SITE_STD_MUTATION_RATE);
	fprintf(fp_parameters, "degradation_rate,%e\n", DEGRADATION_RATE);
	fprintf(fp_parameters, "diffusion_rate,%e\n", DIFFUSION_RATE);
	fprintf(fp_parameters, "nucleus_control_factor,%e\n", NUCLEUS_CONTROL_FACTOR);
	fprintf(fp_parameters, "target_pop,%e\n", TARGET_POP);
	fprintf(fp_parameters, "density,%e\n", DENSITY);
	fprintf(fp_parameters, "len_genome,%d\n", LEN_GENOME);

	fprintf(fp_parameters, "sim_length,%e\n", SIM_LENGTH);
	fprintf(fp_parameters, "introduce_after,%e\n", INTRODUCE_AFTER);
	fprintf(fp_parameters, "recording_space,%e\n", RECORDING_SPACE);

	fprintf(fp_parameters, "max_n_events,%d\n", MAX_N_EVENTS);
	fprintf(fp_parameters, "max_mutants,%d\n", MAX_MUTANTS);
	fclose(fp_parameters);

	// Set up files to write population data in
	char ssd_population_filename[50];
	sprintf(ssd_population_filename, "ssd_sim_populations%d.txt", seed);
	FILE *fp_ssd_population = fopen(ssd_population_filename, "w");
	fprintf(fp_ssd_population, "t,cell,w,m\n");
	fclose(fp_ssd_population);

	// Set up files to write site frequency spectrum data in
	char ssd_sfs_filename[50];
	sprintf(ssd_sfs_filename, "ssd_sim_site_frequency_spectrum%d.txt", seed);
	FILE *fp_ssd_sfs = fopen(ssd_sfs_filename, "w");
	fprintf(fp_ssd_sfs, "t,cell,sfs\n");
	fclose(fp_ssd_sfs);

	// Declare variables to store standard mutational information of wildtype individuals
	/* wildtype_state contains mutational information about wildtype individuals
	wildtype_state[k][i][0] is the number of standard mutations which individual i possesses in cell k
	wildtype_state[k][i][1:] are the identities of the standard mutations which individual i possesses in cell k */
	int*** wildtype_state = malloc(CELLS * sizeof(int**));
	/* wildtype_populations[k] is the number of wildtype individuals in cell k */
	int* wildtype_populations = malloc(CELLS * sizeof(int));

	/* mutant_counts[0][0] is the latest mutation identity across the whole system
	mutant_counts[1:][0] = 0 are dummy entries
	mutant_counts[k][i] is the number of individuals in cell k with mutation identity i
	mutant_counts has number of columns mutant_counts[0][0] + 1 */
	int** mutant_counts = malloc(CELLS * sizeof(int*));

	/* wildtype_propensity[k][0] is the propensity of degradation of a wildtype individual in cell k
	wildtype_propensity[k][1] is the propensity of replication of a wildtype individual in cell k
	wildtype_propensity[k][2] is the propensity of diffusion of a  wildtype in cell k */
	double** wildtype_propensity = malloc(CELLS * sizeof(double*));
	/* wildtype_propensity_sums[k] is the sum of the wildtype propensity vector of cell k */
	double* wildtype_propensity_sums = malloc(CELLS * sizeof(double));

	for (int k=0; k<CELLS; ++k) {
		// Initialise wildtype population steady state
		wildtype_state[k] = malloc((int) target_population * sizeof(int*));
		for (int i=0; i<(int) target_population; ++i) {wildtype_state[k][i] = calloc(1, sizeof(int));}
		wildtype_populations[k] = (int) target_population;
	
		// Initialise mutant counts
		mutant_counts[k] = calloc(1, sizeof(int));

		// Compute initial wildtype propensity
		wildtype_propensity[k] = malloc(3 * sizeof(double));
		wildtype_propensity_update(wildtype_propensity, wildtype_propensity_sums, k, wildtype_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
	}

	double wildtype_propensity_sum_across_cells;
	double current_time = 0.0;
	double recording_time = 0.0;
	int n_event = 0;

	// Record initial data
	write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, recording_time, ssd_population_filename, ssd_sfs_filename);
	recording_time += RECORDING_SPACE;

	// Gillespie algorithm until time threshold reached
	while (current_time<INTRODUCE_AFTER && n_event < MAX_N_EVENTS) {
		
		// Realise event according to propensity
		wildtype_propensity_sum_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {wildtype_propensity_sum_across_cells += wildtype_propensity_sums[k];}
		current_time += -log(RND) / wildtype_propensity_sum_across_cells;
		wildtype_gillespie_event(rng, wildtype_propensity, wildtype_propensity_sums, wildtype_state, wildtype_populations, mutant_counts, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
		n_event++;

		// Record data at recording time
		if (current_time>=recording_time) {
			compact_relabel_wildtype_mutations(mutant_counts, wildtype_state, wildtype_populations);
			write_data_to_file_pre_introduce(wildtype_populations, mutant_counts, recording_time, ssd_population_filename, ssd_sfs_filename);

			recording_time += RECORDING_SPACE;
		}

		// Debug for when mutant counts > population
		int* args = malloc(2 * sizeof(int));
		int mutant_counts_max = argmax_mutant_counts(args, mutant_counts);
		int total_wildtype_populations;
		for (int k=0; k<CELLS; ++k) {total_wildtype_populations += wildtype_populations[k];}
		if (mutant_counts_max>total_wildtype_populations) {
			printf("n_event = %d\nmutant_counts_max = %d\nargs = (%d, %d)\ntotal_wildtype_populations = %d", n_event, mutant_counts_max, args[0], args[1], total_wildtype_populations);
			exit(99);
		}
		free(args);
	}

	int introduce_cell_idx = CELLS / 2;
	while (wildtype_populations[introduce_cell_idx]==0) {

		// Realise event according to propensity
		wildtype_propensity_sum_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {wildtype_propensity_sum_across_cells += wildtype_propensity_sums[k];}
		current_time += -log(RND) / wildtype_propensity_sum_across_cells;
		wildtype_gillespie_event(rng, wildtype_propensity, wildtype_propensity_sums, wildtype_state, wildtype_populations, mutant_counts, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population);
		n_event++;
	}

	// Free memory
	for (int k=0; k<CELLS; ++k) {free(wildtype_propensity[k]);}
	free(wildtype_propensity);
	free(wildtype_propensity_sums);

	/* ssd_state contains mutational information about RA individuals
	ssd_state[k][i][0] is the number of standard mutations which individual i possesses in cell k
	ssd_state[k][i][1:] are the identities of the standard mutations which individual i possesses in cell k */
	int*** ssd_state = malloc(CELLS * sizeof(int**));
	/* ssd_populations[k] is the population size of RA individuals in cell k */
	int* ssd_populations = malloc(CELLS * sizeof(int));

	// Introduce 1 SSD individual at central cell
	introduce_ra_or_ssd(rng, wildtype_state, ssd_state, wildtype_populations, ssd_populations, introduce_cell_idx);			

	// Initialise propensity
	/* propensity[k][0] is the propensity of degradation of a wildtype individual in cell k
	propensity[k][1] is the propensity of replication of a wildtype individual in cell k
	propensity[k][2] is the propensity of diffusion of a wildtype individual in cell k
	propensity[k][3] is the propensity of degradation of an RA/SSD individual in cell k
	propensity[k][4] is the propensity of replication of an RA/SSD individual in cell k
	propensity[k][5] is the propensity of diffusion of an RA/SSD individual in cell k*/
	double** propensity = malloc(CELLS * sizeof(double*));
	/* propensity_sums[k] is the sum of the propensity vector of cell k */
	double* propensity_sums = malloc(CELLS * sizeof(double));
	for (int k=0; k<CELLS; ++k) {
		propensity[k] = malloc(6 * sizeof(double));
		propensity_update(propensity, propensity_sums, k, wildtype_populations, ssd_populations, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, density);
	}

	double propensity_sum_across_cells;
	
	// Simulate
	while (current_time<SIM_LENGTH && n_event<MAX_N_EVENTS) {
		// if (n_event%100==0) {printf("n_event = %d\n", n_event);}
		propensity_sum_across_cells = 0;
		for (int k=0; k<CELLS; ++k) {propensity_sum_across_cells += propensity_sums[k];}
		current_time += -log(RND) / propensity_sum_across_cells;
		gillespie_event(rng, propensity, propensity_sums, wildtype_state, ssd_state, mutant_counts, wildtype_populations, ssd_populations, site_std_mutation_rate, degradation_rate, diffusion_rate, nucleus_control_factor, target_population, density);
		n_event++;

		// In case of extinction of SSD individuals, record data before beginning new simulation
		int total_ssd_population = 0;
		for (int k=0; k<CELLS; ++k) {total_ssd_population += ssd_populations[k];}
		if (total_ssd_population==0) {
			compact_relabel_mutations(mutant_counts, wildtype_state, ssd_state, wildtype_populations, ssd_populations);
			write_data_to_file(wildtype_populations, ssd_populations, mutant_counts, recording_time, ssd_population_filename, ssd_sfs_filename);
			break;
		}

		// Record data at recording time
		if (current_time>=recording_time){
			// printf("t = %.1f\n", current_time);
			compact_relabel_mutations(mutant_counts, wildtype_state, ssd_state, wildtype_populations, ssd_populations);
			write_data_to_file(wildtype_populations, ssd_populations, mutant_counts, recording_time, ssd_population_filename, ssd_sfs_filename);

			recording_time += RECORDING_SPACE;
		}
	}

	// Free memory
	for (int k=0; k<CELLS; ++k){
		for (int i=0; i<wildtype_populations[k]; ++i){free(wildtype_state[k][i]);}
		free(wildtype_state[k]);
		for (int i=0; i<ssd_populations[k]; ++i){free(ssd_state[k][i]);}
		free(ssd_state[k]);

		free(mutant_counts[k]);

		free(propensity[k]);
	}
	free(wildtype_state);
	free(wildtype_populations);
	free(ssd_state);
	free(ssd_populations);

	free(mutant_counts);

	free(propensity);
	free(propensity_sums);

    gsl_rng_free(rng);
	return 0;
}