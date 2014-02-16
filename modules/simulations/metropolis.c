

#define METROPOLIS_C

#include <stdlib.h>
#include <math.h>
#include "simulations.h"
#include "random.h"

/* Cold initialization */
void cold_init(double *v, int dim) {
	int i;
	for(i=0; i<dim; i++)
		v[i] = 0;
}


/* Hot initialization */
void hot_init(double *v, int dim) {
	int i;
	double *temp;
	temp = malloc(dim*sizeof(double));
	ranlxd(temp,dim);
	for(i=0; i<dim; i++)
		v[i] = 10.0*(2*temp[i] - 1);
	
	free(temp);
}


/* Routine that executes a sweep of the Metropolis algorithm */
double metropolis (double (*P)(double *), double *state, int state_dim, double delta) {
	int i;
	double swap, x_new, acceptance;
	double chosen = 0;
	double u[2];

	for(i=0; i<state_dim; i++) {
		ranlxd(u,2);
		x_new = state[i] + delta*(u[0] - 0.5);
		swap = state[i];
		state[i] = x_new;
			acceptance = P(state);
		state[i] = swap;
			acceptance /= P(state);
		
		if(acceptance >= u[1]) {
			state[i] = x_new;
			chosen += 1.0/state_dim;
		}
	}
	return chosen;
}


/* Same as 'metropolis' but configuration space on a 1D lattice */
double L1metropolis (double (*P)(double *), double *state, double delta) {
	double swap, x_new, acceptance;
	double chosen = 0;
	double u[2];

	ranlxd(u,2);
	x_new = *state + delta*(2.0*u[0] - 1.0);
	x_new = rint(x_new);
	swap = *state;
	*state = x_new;
		acceptance = P(state);
	*state = swap;
		acceptance /= P(state);
	
    if(acceptance >= u[1]) {
    	*state = x_new;
    	chosen++;
	}
    return chosen;
}
