
#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include "cluster.h"

#ifndef HARMOSC_C
extern double HOmetropolis (double *state, int state_dim);
extern double HOeAction (double *state, int state_dim);
extern void HOautocorrelation (double *x, int t, double *v, int x_dim, int steps);
extern cluster DeltaE (cluster *A, cluster *B, cluster *C);
extern cluster MatrixElementX (cluster *DE, cluster *Corr, int t, int N);
extern cluster sqrt_jk (cluster *X);
#endif

#ifndef METROPOLIS_C
extern void cold_init (double *v, int dim);
extern void hot_init (double *v, int dim);
extern double metropolis (double (*P)(double *), double *state, int state_dim, double delta);
#endif

#ifndef ONEDLATTICE_C
extern double L1metropolis (double (*P)(double *), double *state, double delta);
extern cluster energy_derivative (cluster *enld, cluster *en, cluster *ld);
#endif

#ifndef STATISTICS_C
extern double autocorrelation (double *x, int t, int dim);
extern double correlation (double *x, double *y, int dim);
#endif

#endif
