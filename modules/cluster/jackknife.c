
/*******************************************************************************
 * 
 * 		File "jackknife.c"
 * 
 * File contenente le definizioni delle routines necessarie per la manipolazione
 * di cluster jackknife.
 * 
 * Le routines sono:
 * _ cluster_init -> inizializzazione di un cluster jackknife;
 * 
 * _ clusterJK -> assegnamento di media e varianza della media nel cluster
 * 		jackknife passato ad argomento.
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
cluster JKfunction (double (*f)(double), cluster *X)
{
	int i;
	double temp = 0;
	int dim = X->Dim;
	cluster result;
	JKinit(&result,dim);
	result.Mean = f(X->Mean);
	for(i=0; i<dim; i++)
	{
		result.Vec[i] = f(X->Vec[i]);
		temp += (result.Vec[i] - result.Mean)*(result.Vec[i] - result.Mean);
		temp *= ((double)(dim - 1)/(double)dim);
	}
	result.Var = temp;
	
	return result;
}
