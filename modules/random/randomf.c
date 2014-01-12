/*
 *      Routine tu be used in Fortran to generate random numbers with the
 *      Martin Luscher algorithm
 **/

#include "random.h"

void randf_ (double x[], int n) {
    ranlxd(x,n);
}
