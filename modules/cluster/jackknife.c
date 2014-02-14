
/*******************************************************************************
 * 
 * 		File "jackknife.c"
 * 
 * File containing routines useful to handle cluster jackknife to compute
 * averages and variances of primary and secondary random variables
 * 
 * Routines are:
 * _ JKinit -> initialization of a cluster jackknife;
 * 
 * _ JKcluster -> assignment of mean and variance (of the mean) of the sample
 *		stored in the array contained in the cluster.
 *
 * _ JKfunction -> returns a cluster jackknife with samples, mean and variance
 *		(of the mean) of a secondary random variable defined as a function
 *		(passed as argument) of several primary random variables.
 *		For the sake of generality this routine requires a POINTER to an ARRAY
 *		of jackknife clusters as argument.
 * 
 ******************************************************************************/

#define JACKKNIFE_C

#include <stdio.h>
#include <stdlib.h>
#include "cluster.h"


/* Assignment of mean value and variance in a cluster structure */
void JKcluster (cluster *C)
{
	int i;
	int dim = C->Dim;
	
	C->Mean = 0;
	for(i=0; i<dim; i++)
		C->Mean += (C->Vec[i])/((double)dim);
	
	for(i=0; i<dim; i++)
		C->Vec[i] = C->Mean + (C->Mean - C->Vec[i])/((double)(dim - 1));
	
	C->Var = 0;
	for(i=0; i<dim; i++)
	{
		C->Var += (C->Vec[i] - C->Mean)*(C->Vec[i] - C->Mean);
		C->Var *= ((double)(dim - 1))/((double)dim);
	}
}


/* Cluster jackknive initialization */
void JKinit (cluster *C, int dim)
{
	C->Dim	= dim;
	C->Mean	= 0;
	C->Var= 0;
	C->Vec	= malloc(dim*sizeof(double));
}


/* Compute the JK of a secondary random variable */
cluster JKfunction (double (*f)(double *), int narg, cluster *X)
{
	int i, j;
	double temp = 0;
	double *args;

	cluster result;
	int dim = X->Dim;
	JKinit(&result,dim);

	args = malloc(narg*sizeof(double));
	for (j=0; j<narg; j++)
		args[j] = X[j].Mean;
	result.Mean = f(args);
	for(i=0; i<dim; i++) {
		for (j=0; j<narg; j++)
			args[j] = X[j].Vec[i];
		result.Vec[i] = f(args);
		temp += (result.Vec[i] - result.Mean)*(result.Vec[i] - result.Mean);
		temp *= ((double)(dim - 1)/(double)dim);
	}
	result.Var = temp;
	free(args);
	
	return result;
}
