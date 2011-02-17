#include <R.h>
#include <Rdefines.h>

/*

This is heavily inspired by the Fortran interface
in the ucminf package (by Stig Bousgaard Mortensen).
Author: Bjoern Bornkamp

*/

SEXP logPost;
SEXP mnFns;
SEXP env;

extern void F77_NAME(banint)(int *m, double *relreq, int *maxvls,
			     int *rs,  int *mn, char *problm, 
			     double *mu, double *C, int *pru,
			     int *numtrn, int *method, double *nrmcon, 
			     double *means, double *errors, 
			     double *covrnc, int *inform);

double F77_SUB(usrlgpc)(double *par){
  SEXP val,m;
  int i,m2;
  double out;
  PROTECT(m = findVarInFrame(env, install(".m")));
  m2 = asInteger(m);
  PROTECT(val = findVarInFrame(env, install(".par")));
  for (i = 0; i < m2; i++) 
    REAL(val)[i] = par[i];
  out = asReal(eval(logPost, env));
  if (!R_FINITE(out))
      error("non-finite function value from user function");
  UNPROTECT(2);
  return out;
}

void F77_SUB(usrmnsc)(double *par, double *mns) {
  SEXP val,out,m;
  int i,m2;
  PROTECT(m = findVarInFrame(env, install(".m")));
  m2 = asInteger(m);
  PROTECT(out = allocVector(REALSXP, m2));
  PROTECT(val = findVarInFrame(env, install(".par")));
  for (i = 0; i < m2; i++) 
    REAL(val)[i] = par[i];
  out = eval(mnFns, env);
  for (i = 0; i < m2; i++) 
    mns[i] = REAL(out)[i];
  UNPROTECT(3);
}

SEXP BAYPOST(SEXP RlogPost, SEXP RmnFns, SEXP rho){
  /* fixed parameters */
  char problm[7] = "generic";
  int rs = 0;
  /* R objects */
  SEXP m,mn,mu,means,errors,covrnc,covmod;
  SEXP relreq,maxvls,numtrn,method,pru;
  SEXP nrmcon,inform;

  logPost = RlogPost;
  mnFns = RmnFns;
  env = rho;
  PROTECT(m      = findVarInFrame(rho, install(".m")));
  PROTECT(mn     = findVarInFrame(rho, install(".mn")));
  PROTECT(mu     = findVarInFrame(rho, install(".mu")));
  PROTECT(means  = findVarInFrame(rho, install(".means")));
  PROTECT(errors = findVarInFrame(rho, install(".errors")));
  PROTECT(covrnc = findVarInFrame(rho, install(".covrnc")));
  PROTECT(covmod = findVarInFrame(rho, install(".covmod")));
  PROTECT(relreq = findVarInFrame(rho, install(".relreq")));
  PROTECT(maxvls = findVarInFrame(rho, install(".maxvls")));
  PROTECT(numtrn = findVarInFrame(rho, install(".numtrn")));
  PROTECT(method = findVarInFrame(rho, install(".method")));
  PROTECT(nrmcon = findVarInFrame(rho, install(".nrmcon")));
  PROTECT(pru = findVarInFrame(rho, install(".pru")));
  PROTECT(inform = findVarInFrame(rho, install(".inform")));

  F77_CALL(banint)(INTEGER(m), REAL(relreq), INTEGER(maxvls), &rs, 
  		   INTEGER(mn), problm, REAL(mu), REAL(covmod), INTEGER(pru), 
  		   INTEGER(numtrn), INTEGER(method), REAL(nrmcon), 
		   REAL(means), REAL(errors), REAL(covrnc), 
		   INTEGER(inform));

  UNPROTECT(14);
  return R_NilValue;
}
