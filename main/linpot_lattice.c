
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
int TMAX;		/* maximum time delay for autocorrelation */
int NCONV;		/* number of iterations needed to convergence */
int NAV;		/* number of iterations to average on after convergence */

double delta;	/* hypercube side */
double vpar;	/* variational parameter of the trial wave funct */
double T;		/* hopping amplitude */
double V;		/* strength of the local potential */
double L;		/* number of lattice sites */

double probability (double *x);
double trialWF (double x);
double localenergy (double x);
//double ender (double *x);


int main (int argc, char *argv[]) {

	int NDAT;	/* size of the data sample */
	int sw, i, counter;
	double en_sum, enld_sum, ld_sum, en_wsum, en_varsum, \
		site, vpar_init, Dvpar, acptrate;
	double *site_p;
		site_p = &site;

	/* jackknife structure for the energy estimator */
	cluster vpar_av, enld, en, ld, ender;

	/* array needed to compute the autocorrelation */
	double *autocorr;
	
	/* Name of the output files */
	char *therm_name, *integral_name, *autocorr_name, *distrib_name, *vpar_name;
		therm_name		= malloc(100*sizeof(char));
		integral_name	= malloc(100*sizeof(char));
		autocorr_name	= malloc(100*sizeof(char));
		distrib_name	= malloc(100*sizeof(char));
		vpar_name		= malloc(100*sizeof(char));
	sprintf(therm_name, "lp_output/thermalization.dat");
	sprintf(integral_name, "lp_output/expectationvalues.dat");
	sprintf(autocorr_name, "lp_output/autocorrelation.dat");
	sprintf(distrib_name, "lp_output/distribution.dat");
	sprintf(vpar_name, "lp_output/variationalparameter.dat");
	
	/* Output streams */
	FILE *therm_file, *integral_file, *autocorr_file, *distrib_file, *vpar_file;
		therm_file		= fopen(therm_name, "w");
		integral_file	= fopen(integral_name, "w");
		autocorr_file	= fopen(autocorr_name, "w");
		distrib_file	= fopen(distrib_name, "w");
		vpar_file		= fopen(vpar_name, "w");
			
	/* Initialization of the randomizer */
	srand(time(NULL));
	rlxd_init(1,rand());

	scanf("%d",&NTH);
	scanf("%d",&NSW);
	scanf("%d",&NBIN);
	scanf("%d",&TMAX);
	scanf("%d",&NCONV);
	scanf("%d",&NAV);
		
	scanf("%lf",&delta);
	scanf("%lf",&L);
	scanf("%lf",&site);
	scanf("%lf",&T);
	scanf("%lf",&V);
	scanf("%lf",&vpar_init);
	scanf("%lf",&Dvpar);

	NDAT = NSW/NBIN;
	printf("Number of lattice sites : %d \n", (int)L);
	printf("Number of data points   : %d \n", NDAT);
	autocorr = malloc(NSW*sizeof(double));
	JKinit(&enld, NDAT);
	JKinit(&en, NDAT);
	JKinit(&ld, NDAT);
	JKinit(&vpar_av, NAV);

	vpar = vpar_init;
	en_wsum = en_varsum = 0;
	for(counter=-NCONV; counter<NAV; counter++) {

		/* Process is left free for a certain number NTH of sweeps:
	 	* no data are used for estimating the integral
	 	**/
	 	acptrate = 0;
		for (sw=0; sw<NTH; sw++) {
			acptrate += L1metropolis (probability, site_p, delta)/NTH;
			site = *site_p;
			fprintf(therm_file, "%d\t%3.1lf\n", sw, site);
		}
		// printf("Acceptance rate in thermalization : %lf\n", acptrate);
		
		/* From now on data are collected and used for the estimation */
		i=0;
		en_sum = ld_sum = enld_sum = 0;
		acptrate = 0;
		for(sw=0; sw<NSW; sw++) {
			acptrate += L1metropolis (probability, site_p, delta)/NSW;
			site = *site_p;

			/* Print site to see the distribution */
			if (counter == 0) {
				fprintf(distrib_file, "%3.1lf\n", site);
				/* store the value of 'site' to compute the autocorrelation */
				autocorr[sw] = site;
			}

			en_sum   += localenergy(site)/NBIN;
			ld_sum   += -site/NBIN;
			enld_sum += -site*localenergy(site)/NBIN;

			/* Applying binning method */
			if((sw+1)%NBIN==0) {
				enld.Vec[i] = enld_sum;
				en.Vec[i] = en_sum;
				ld.Vec[i] = ld_sum;
				en_sum = enld_sum = ld_sum = 0;
				i++;
			}
		}
		// printf("Acceptance rate when thermalized  : %lf\n", acptrate);

		/* Average and variance calculated with the jackknife re-sampling */
		JKcluster(&enld);	/* covariance energy-logarithmic derivative */
		JKcluster(&en);		/* energy */
		JKcluster(&ld);		/* logarithmic derivative */

		/* Average and variance of the variational energy derivative */
		ender = energy_derivative (&enld, &en, &ld);
		printf("vpar = %lf,\tenergy mean / derivative = %lf / %.5e\n", \
			vpar, en.Mean, ender.Mean);

		/* Print mean and variance of energy and its derivative on file */
		if (counter < 0)
			fprintf(integral_file, "%lf\t%.10e\t%.10e\t%.10e\t%.10e\n", \
				vpar, en.Mean, sqrt(en.Var), ender.Mean, sqrt(ender.Var));

		/* New variational parameter calculated via Steepest Descent */
		fprintf(vpar_file, "%d\t%.10e\n", counter, vpar);

		/* Average over several iterations
		 * when (assuming to have) converged to the minimum */
		if (counter >= 0) {
			vpar_av.Vec[counter] = vpar;
			en_wsum		+= en.Mean/en.Var;
			en_varsum 	+= 1.0/en.Var;
		}
		vpar = vpar - Dvpar*ender.Mean;
		// vpar = vpar - Dvpar*2.0*correlation(en.Vec, ld.Vec, NDAT);
	}

	/* Calculate the mean variational parameter and its variance */
	JKcluster(&vpar_av);
	printf("Optimal variational parameter      : %lf +- %lf \n", \
		vpar_av.Mean, sqrt(vpar_av.Var));

	/* Calculate the minimum variational energy as a weighted average */
	printf("Average minimum variational energy : %lf +- %.10e \n\n", \
		en_wsum/en_varsum, 1.0/en_varsum);

	/* Print the autocorrelation of the process 'site' for the initial vpar */
	for(i=0; i<TMAX+1; i++)
		fprintf(autocorr_file, "%d\t%e\n", i, autocorrelation(autocorr,i,NDAT));

	fclose(therm_file);
	fclose(autocorr_file);
	fclose(integral_file);
	free(autocorr);
	JKdelete(&enld);
	JKdelete(&en);
	JKdelete(&ld);
	JKdelete(&ender);
	
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
