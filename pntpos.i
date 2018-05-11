/* pntpos.i */

%module pntpos

%{

	/* Put header files here or function declarations like below */

	

	extern double varerr(const prcopt_t *opt, double el, int sys);

	extern int rescode(int iter, const obsd_t *obs, int n, const double *rs,

                   const double *dts, const double *vare, const int *svh,

                   const nav_t *nav, const double *x, const prcopt_t *opt,

                   double *v, double *H, double *var, double *azel, int *vsat,

                   double *resp, int *ns);

	extern int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,

                  const double *vare, const int *svh, const nav_t *nav,

                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,

                  double *resp, char *msg);

%}



extern double varerr(const prcopt_t *opt, double el, int sys);

	extern int rescode(int iter, const obsd_t *obs, int n, const double *rs,

                   const double *dts, const double *vare, const int *svh,

                   const nav_t *nav, const double *x, const prcopt_t *opt,

                   double *v, double *H, double *var, double *azel, int *vsat,

                   double *resp, int *ns);

	extern int estpos(const obsd_t *obs, int n, const double *rs, const double *dts,

                  const double *vare, const int *svh, const nav_t *nav,

                  const prcopt_t *opt, sol_t *sol, double *azel, int *vsat,

                  double *resp, char *msg);
