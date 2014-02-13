
/*******************************************************************************
 *
 *	Variational Monte Carlo (VMC) for the single particle in a linear potential
 *	on a 1D lattice.
 *	Optimization performed with the steepest descent.
 *
 *******************************************************************************/

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
double vpar;	/* variational parameter of the trial wave funct */
double T;		/* hopping amplitude */
double V;		/* strength of the local potential */
double L;			/* number of lattice sites */

double probability (double *x);
double trialWF (double x);
double localenergy (double x);


int main (int argc, char *argv[]) {	

	int NDAT;	/* size of the data sample */
	int sw, i;
	double en_sum, enld_sum, ld_sum, site, vpar_init, vpar_old, Dvpar, accuracy;
	double *site_p;
		site_p = &site;

	/* jackknife structure for the energy estimator */
	cluster energy, enldcorr, ld;

	/* array needed to compute the autocorrelation */
	double *autocorr;
	
	/* Name of the output files */
	char *therm_name, *integral_name, *autocorr_name, *distrib_name;
		therm_name		= malloc(100*sizeof(char));
		integral_name	= malloc(100*sizeof(char));
		autocorr_name	= malloc(100*sizeof(char));
		distrib_name	= malloc(100*sizeof(char));
	sprintf(therm_name, "lp_output/thermalization.dat");
	sprintf(integral_name, "lp_output/expectationvalues.dat");
	sprintf(autocorr_name, "lp_output/autocorrelation.dat");
	sprintf(distrib_name, "lp_output/distribution.dat");
	
	/* Output streams */
	FILE *therm_file, *integral_file, *autocorr_file, *distrib_file;
		therm_file		= fopen(therm_name, "w");
		integral_file	= fopen(integral_name, "w");
		autocorr_file	= fopen(autocorr_name, "w");
		distrib_file	= fopen(distrib_name, "w");
			
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
	scanf("%lf",&vpar_init);
	scanf("%lf",&accuracy);
	scanf("%lf",&Dvpar);

	NDAT = NSW/NBIN;
	printf("Number of lattice sites : %d \n", (int)L);
	printf("Number of data points   : %d \n", NDAT);
	JKinit(&energy, NDAT);
	JKinit(&enldcorr, NDAT);
	JKinit(&ld, NDAT);
	autocorr = malloc(NSW*sizeof(double));

	vpar = vpar_init;
	do {

		/* Process is left free for a certain number NTH of sweeps:
	 	* no data are used for estimating the integral
	 	**/
		for(sw=0; sw<NTH; sw++) {
			L1metropolis(probability, site_p, delta);
			site = *site_p;
			fprintf(therm_file, "%d\t%3.1lf\n", sw, site);
		}
		
		/* From now on data are collected and used for the estimation */
		i=0;
		en_sum = ld_sum = enld_sum = 0;
		for(sw=0; sw<NSW; sw++) {
			L1metropolis(probability, site_p, delta);
			site = *site_p;

			/* Print site to see the distribution */
			if (vpar == vpar_init)
				fprintf(distrib_file, "%3.1lf\n", site);

			/* store the value of 'site' to compute the autocorrelation */
			autocorr[sw] = site;
			en_sum += localenergy(site)/NBIN;
			ld_sum += -site/NBIN;
			enld_sum += -site*localenergy(site)/NBIN;

			/* Applying binning method */
			if((sw+1)%NBIN==0) {
				energy.Vec[i] = en_sum;
				enldcorr.Vec[i] = enld_sum;
				ld.Vec[i] = ld_sum;
				en_sum = enld_sum = ld_sum = 0;
				i++;
			}
		}
				
		/* Expectation value and variance computed
	 	* with the jackknife re-sampling method
	 	**/
		JKcluster(&energy);
		JKcluster(&enldcorr);
		JKcluster(&ld);
		fprintf(integral_file, "%lf\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\t%.10e\n", \
			vpar, energy.Mean, sqrt(energy.Var), \
			enldcorr.Mean, sqrt(enldcorr.Var), \
			ld.Mean, sqrt(ld.Var) );

		vpar_old = vpar;
		
		/* New variational parameter calculated via Steepest Descent */
		vpar = vpar - Dvpar*2.0*(enldcorr.Mean - energy.Mean * ld.Mean - 2*T);
		
	} while (fabs(vpar - vpar_old) > accuracy);

	for(i=0; i<TMAX+1; i++)
		fprintf(autocorr_file, "%d\t%e\n", i, autocorrelation(autocorr, i, NDAT));
	printf("Optimal variational parameter : %lf +- %lf \n\n", vpar, accuracy);

	fclose(therm_file);
	fclose(autocorr_file);
	fclose(integral_file);
	free(autocorr);
	
	exit(EXIT_SUCCESS);
}


double localenergy (double x) {
	return -2.0*T*exp(vpar) + V*x;
}

double trialWF (double x) {
	return exp(-vpar*x);
}

double probability (double *x) {
	if (*x < 0.5 || *x > L+.5) return 0;
	else return trialWF(*x)*trialWF(*x);
}
