#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "random.h"

void rlxdinit_(int *lux,int *seed){
    int lux1,seed1;
    lux1=*lux;seed1=*seed;
    rlxd_init(lux1,seed1);
}

void ranlxdf_(double vec[],int *lvec){
    int lvec1; 
    lvec1=*lvec;
    ranlxd(vec,lvec1);
}

void rlxdgetf_(int *state){
    rlxd_get(state);
}

void rlxdresetf_(int *state1){
    rlxd_reset(state1);
}

void rlxd_sizef_(int *n){
    int n1;
    n1=rlxd_size();
    *n=n1;
}
