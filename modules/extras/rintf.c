/*******************************************************************************
 * 
 *      File 'rintf.c'
 * 
 *  rint wrapper for Fortran
 * 
 ******************************************************************************/

#include <tgmath.h>

void rintf_ (double * val, double * rval) {
  *rval = rint(*val);
}
