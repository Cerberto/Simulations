#include <tgmath.h>

void myrint_ (float * val, float * rval) {
  *rval = rint(*val);
}


double frint_ (double x) {
    return rint(x);
}
