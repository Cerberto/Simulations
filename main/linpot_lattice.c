/****
 *
 *	Variational Monte Carlo (VMC) for the single particle in a linear potential
 *
 ***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "random.h"
#include "simulations.h"
#include "cluster.h"

#include "extras.h"

#define delta	1.49	/* hypercube side */
#define NTH		200		/* number of sweeps needed for thermalization */
#define NSW		1000000	/* number of "effective" sweeps */
#define NBIN	100		/* number of sweeps skipped in order to avoid correlation */
#define TMAX	50		/* maximum time delay for autocorrelation */

double varpar;	/* variational parameter of the trial wave funct */
double T;		/* hopping amplitude */
double V;		/* strength of the local potential

/* (not normalized) probability distribution */
double probability (double *x);


int main (int argc, char *argv[]) {	
	int sw, i;
	double site;

	int NDAT = NSW/NBIN;	/* size of the data sample */
	printf("Number of data points : %d", NDAT)

	/* jackknife structure for the energy estimator */
	cluster energy;
		cluster_init(&energy, NDAT);

	/* array needed to compute the autocorrelation */
	double *autocorr;
		autocorr = malloc(NSW*sizeof(double));
	
	/* Name of the output files */
	char *therm_name, *integral_name, *autocorr_name;
			= malloc(10*sizeof(char));
		therm_name		= malloc(100*sizeof(char));
		integral_name	= malloc(100*sizeof(char));
		autocorr_name	= malloc(100*sizeof(char));
	sprintf(therm_name, "lp_output/thermalization.dat");
	sprintf(integral_name, "lp_output/expectationvalues.dat");
	sprintf(autocorr_name, "lp_output/autocorrelation.dat");
	
	/* Output streams */
	FILE *distrib_file, *therm_file, *integral_file, *autocorr_file;
		therm_file		= fopen(therm_name, "w");
		integral_file	= fopen(integral_name, "w");
		autocorr_file	= fopen(autocorr_name, "w");

	printf("\nData will be saved in the following files:\n");
	printf("\t\"%s\"\t-> estimations of energy;\n", integral_name);
	printf("\t\"%s\"\t-> values of x at each sweep of the metropolis;\n", therm_name);
	printf("\t\"%s\"\t-> autocorrelation of data.\n\n", autocorr_name);
			
	/* Initialization of the randomizer */
	srand(time(NULL));
	rlxd_init(1,rand());
		
	/* first point is picked randomly ("hot" initialization) */
	hot_init(&site,1);

	/* Process is left free for a certain number NTH of sweeps:
	 * no data are used for estimating the integral
	 **/
	for(sw=0; sw<NTH; sw++) {
		metropolis(probability, site, 1, delta);
		fprintf(therm_file, "%d\t%.10e\n", sw, site);
	}
	
	/* From now on data are collected and used for the estimation */
	i=0;
	for(sw=0; sw<NSW; sw++) {
		metropolis(probability, &site, 1, delta);
		
		/* All data corresponding to each sweep of the metropolis algorithm
		* (apart from those needed for thermalization) are stored in the
		* autocorr array (in order to compute the autocorrelation), just
		* for a particular value of the parameter varpar
		**/
		autocorr[sw] = site;
		fprintf(distrib_file, "%.10e\n", site);
		fprintf(therm_file, "%d\t%.10e\n", sw+NTH, site);
		
		/* In order to get independent samples we take into account just
		 * the configuration every NBIN sweeps; data are stored in an array
		 * contained in a jackknife structure.
		 **/
		if((sw+1)%NBIN==0) {
			energy.Vec[i] = localenergy(site);
			i++;
		}
	}
			
	/* Expectation value and variance of the integral are computed
	 * with the jackknife re-sampling method
	 **/
	clusterJK(&energy);
	fprintf(integral_file, "%lf\t%.10e\t%.10e\n", varpar, energy.Mean, energy.Sigma);
	
	for(i=0; i<TMAX+1; i++)
			fprintf(autocorr_file, "%d\t%e\n", i, autocorrelation(autocorr, i, NDAT));
	
	fclose(therm_file);
	fclose(distrib_file);
	fclose(autocorr_file);
	fclose(integral_file);
	free(autocorr);
	
	exit(EXIT_SUCCESS);
}


double localenergy (double *x) {
	return -2.0*T*exp(varpar) + V*x
}

double trialWF (double *x) {
	return exp(-varpar*x);
}

double probability (double *x) {
	return trialWF(x)*trialWF(x);
}
