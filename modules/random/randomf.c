/*
 *      Routine tu be used in Fortran to generate random numbers with the
 *      Martin Luescher algorithm
 **/

#include "random.h"

void frand_ (double x[], int n) {
    ranlxd(x,n);
}
