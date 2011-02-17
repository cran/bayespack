static double ((*cfun)());
static void ((*cmean)());
double usrfnc_( x ) double *x; { return (*cfun)(x); }
void meansb_( x, fvl ) double *x, *fvl; { cmean( x, fvl ); }
int bzlpnt_( cfunp, cmeanp ) void *cfunp,*cmeanp; { cfun = cfunp; cmean = cmeanp; }
