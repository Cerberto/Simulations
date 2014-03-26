
/*******************************************************************************
 *
 *	Variational Monte Carlo (VMC) for the single particle in a linear potential
 *	on a 1D lattice.
 *	Optimization performed with the steepest descent.
 *
 ******************************************************************************/

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

double delta;	/* hypercube side */
double vpar;	/* variational parameter of the trial wave funct */
double T;		/* hopping amplitude */
double V;		/* strength of the local potential */
double L;		/* number of lattice sites */

double probability (double *x);
double trialWF (double x);
double localenergy (double x);


int main (int argc, char *argv[]) {

	int NDAT;	/* size of the data sample */
	int sw, i;
	double en_sum, enld_sum, ld_sum, site;
	double en_mean, enld_mean, ld_mean;
	double *en_vec, *enld_vec, *ld_vec;
	double *site_p;
		site_p = &site;
	
	/* Output file name */
	char *integral_name;
		integral_name	= malloc(100*sizeof(char));
	sprintf(integral_name, "lp_output/expectationvalues_bf.dat");
	
	/* Output stream */
	FILE *integral_file;
		integral_file	= fopen(integral_name, "a");
			
	/* Initialization of the randomizer */
	srand(time(NULL));
	rlxd_init(1,rand());

	scanf("%d",&NTH);
	scanf("%d",&NSW);
	scanf("%d",&NBIN);
		
	//scanf("%lf",&delta);	/* ATTENTION!! */
	scanf("%lf",&L);
	scanf("%lf",&site);
	scanf("%lf",&T);
	scanf("%lf",&V);
	scanf("%lf",&delta);

	vpar  = atof(argv[1]);

	NDAT = NSW/NBIN;
	en_vec   = malloc(NDAT*sizeof(double));
	ld_vec   = malloc(NDAT*sizeof(double));
	enld_vec = malloc(NDAT*sizeof(double));

	/* Process is left free for a certain number NTH of sweeps:
 	* no data are used for estimating the integral
 	**/
	for(sw=0; sw<NTH; sw++) {
		L1metropolis(probability, site_p, delta);
		site = *site_p;
	}
	
	/* From now on data are collected and used for the estimation */
	i=0;
	en_sum = ld_sum = enld_sum = 0;
	for(sw=0; sw<NSW; sw++) {
		L1metropolis(probability, site_p, delta);
		site = *site_p;
		en_sum   += localenergy(site)/NBIN;
		ld_sum   += -site/NBIN;
		enld_sum += -site*localenergy(site)/NBIN;
	
		/* Applying binning method */
		if((sw+1)%NBIN==0) {
			en_vec[i]   = en_sum;
			ld_vec[i]   = ld_sum;
			enld_vec[i] = enld_sum;
			en_sum = enld_sum = ld_sum = 0;
			i++;
		}
	}

	en_mean = ld_mean = enld_mean = 0;
	for (i=0; i<NDAT; i++) {
		en_mean   += en_vec[i]/NDAT;
		ld_mean	  += ld_vec[i]/NDAT;
		enld_mean += enld_vec[i]/NDAT;
	}

	/* print variational energy and its derivative on file */
	fprintf(integral_file, "%lf\t%.10e\t%.10e\n", \
		vpar, en_mean, 2.0*(enld_mean - en_mean*ld_mean));

	fclose(integral_file);
	
	exit(EXIT_SUCCESS);
}


double localenergy (double x) {
	return -T*(trialWF(x+1.0)+trialWF(x-1.0))/trialWF(x) + V*x;
}

double trialWF (double x) {
	if (x < 0.5 || x > L+.5) return 0.;
	else return exp(-vpar*x);
}

double probability (double *x) {
	return trialWF(*x)*trialWF(*x);
}
