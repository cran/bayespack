#include<R.h>
#include<Rmath.h>

void F77_SUB(rndstart)(void) { GetRNGstate(); }
void F77_SUB(rndend)(void) { PutRNGstate(); }

double F77_SUB(uniran)(void){
  return unif_rand();
}

double F77_SUB(rnruni)(void){
  return unif_rand();
}

double F77_SUB(norran)(void){
  return norm_rand();
}

double F77_SUB(rnrnor)(void){
  return norm_rand();
}

double F77_SUB(gamran)(double *a){
  return rgamma(*a, 1);
}

double F77_SUB(betran)(double *a, double *b){
  return rbeta(*a, *b);
}

