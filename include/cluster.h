
#ifndef CLUSTER_H
#define CLUSTER_H


typedef struct
{
	double *Vec;	/* Array of data */
	int Dim;		/* Array length */
	double Mean;	/* Mean value on the array */
	double Var;		/* Variance of the mean */
} cluster;


#ifndef JACKKNIFE_C
extern void JKcluster(cluster *C);
extern void JKinit(cluster *C, int dim);
extern cluster JKfunction (double (*f)(double *), int narg, cluster *X);
extern void JKdelete (cluster *C);
#endif

#endif

