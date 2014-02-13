/****
 *
 *	Variational Monte Carlo (VMC) for the single particle in a linear potential
 *
 ***/

#define MAIN_PROGRAM


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "random.h"
#include "simulations.h"
#include "cluster.h"

#include "extras.h"


int NTH;		/* number of sweeps needed for thermalization */
int NSW;		/* number of "effective" sweeps */
int NBIN;		/* number of sweeps skipped in order to avoid correlation */
int TMAX;		/* maximum time delay for autocorrelation */

double delta;	/* hypercube side */
double varpar;	/* variational parameter of the trial wave funct */
double T;		/* hopping amplitude */
double V;		/* strength of the local potential */
double L;			/* number of lattice sites */

double probability (double *x);
double trialWF (double x);
double localenergy (double x);

void pluto () {
	printf("\nPluto!");
	fflush(stdout);
}

int main (int argc, char *argv[]) {	

	int NDAT;	/* size of the data sample */
	int sw, i;
	double en_sum, site;
	double *site_p;
		site_p = &site;

	/* jackknife structure for the energy estimator */
	cluster energy;

	/* array needed to compute the autocorrelation */
	double *autocorr;
	
	/* Name of the output files */
	char *therm_name, *integral_name, *autocorr_name;
		therm_name		= malloc(100*sizeof(char));
		integral_name	= malloc(100*sizeof(char));
		autocorr_name	= malloc(100*sizeof(char));
	sprintf(therm_name, "lp_output/thermalization.dat");
	sprintf(integral_name, "lp_output/expectationvalues.dat");
	sprintf(autocorr_name, "lp_output/autocorrelation.dat");
	
	/* Output streams */
	FILE *therm_file, *integral_file, *autocorr_file;
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

	scanf("%d",&NTH);
	scanf("%d",&NSW);
	scanf("%d",&NBIN);
	scanf("%d",&TMAX);
		
	scanf("%lf",&delta);
	scanf("%lf",&L);
	scanf("%lf",&site);
	scanf("%lf",&T);
	scanf("%lf",&V);
	scanf("%lf",&varpar);

	NDAT = NSW/NBIN;
	printf("Number of data points : %d", NDAT);
	JKinit(&energy, NDAT);
	autocorr = malloc(NSW*sizeof(double));

do {

	/* Process is left free for a certain number NTH of sweeps:
	 * no data are used for estimating the integral
	 **/
	for(sw=0; sw<NTH; sw++) {
		L1metropolis(probability, site_p, delta);
		site = *site_p;
//		fprintf(therm_file, "%d\t%.10e\n", sw, site);
	}
	
	/* From now on data are collected and used for the estimation */
	i=0;
	en_sum = 0;
	for(sw=0; sw<NSW; sw++) {
		L1metropolis(probability, site_p, delta);
		site = *site_p;

		/* store the value of 'site' to compute the autocorrelation */
		autocorr[sw] = site;
		// fprintf(therm_file, "%d\t%.10e\n", sw+NTH, site);
		en_sum += localenergy(site);

		/* Applying binning method */
		if((sw+1)%NBIN==0) {
			energy.Vec[i] = en_sum;
			en_sum = 0;
			i++;
		}
	}
			
	/* Expectation value and variance of the integral are computed
	 * with the jackknife re-sampling method
	 **/
	JKcluster(&energy);
	fprintf(integral_file, "%lf\t%.10e\t%.10e\n", varpar, energy.Mean, energy.Sigma);
	
	for(i=0; i<TMAX+1; i++)
		fprintf(autocorr_file, "%d\t%e\n", i, autocorrelation(autocorr, i, NDAT));

	varpar += 0.05;
} while (varpar<1.0);


	fclose(therm_file);
	fclose(autocorr_file);
	fclose(integral_file);
	free(autocorr);
	
	exit(EXIT_SUCCESS);
}


double localenergy (double x) {
	return -2.0*T*exp(varpar) + V*x;
}

double trialWF (double x) {
	return exp(-varpar*x);
}

double probability (double *x) {
	if (*x < 0.5 || *x > L+.5) return 0;
	else return trialWF(*x)*trialWF(*x);
}
