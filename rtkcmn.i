/* rtkcmn.i */

%module rtkcmn



%{

static void fatalerr(const char *format, ...);

extern int satno(int sys, int prn);

extern int satsys(int sat, int *prn);

extern double *mat(int n, int m);

extern int *imat(int n, int m);

extern double dot(const double *a, const double *b, int n);

extern double norm(const double *a, int n);

extern void matcpy(double *A, const double *B, int n, int m);



extern void matmul(const char *tr, int n, int k, int m, double alpha,

                   const double *A, const double *B, double beta, double *C);

static int ludcmp(double *A, int n, int *indx, double *d);



static void lubksb(const double *A, int n, const int *indx, double *b);

extern int matinv(double *A, int n);

extern int solve(const char *tr, const double *A, const double *Y, int n,

                 int m, double *X);

extern int lsq(const double *A, const double *y, int n, int m, double *x,

               double *Q);

extern void time2epoch(gtime_t t, double *ep);

extern gtime_t timeadd(gtime_t t, double sec);



extern void time2str(gtime_t t, char *s, int n);



extern char *time_str(gtime_t t, int n);

extern void ecef2pos(const double *r, double *pos);

%}



static void fatalerr(const char *format, ...);

extern int satno(int sys, int prn);

extern int satsys(int sat, int *prn);

extern double *mat(int n, int m);

extern int *imat(int n, int m);

extern double dot(const double *a, const double *b, int n);

extern double norm(const double *a, int n);

extern void matcpy(double *A, const double *B, int n, int m);



extern void matmul(const char *tr, int n, int k, int m, double alpha,

                   const double *A, const double *B, double beta, double *C);

static int ludcmp(double *A, int n, int *indx, double *d);



static void lubksb(const double *A, int n, const int *indx, double *b);

extern int matinv(double *A, int n);

extern int solve(const char *tr, const double *A, const double *Y, int n,

                 int m, double *X);

extern int lsq(const double *A, const double *y, int n, int m, double *x,

               double *Q);

extern void time2epoch(gtime_t t, double *ep);

extern gtime_t timeadd(gtime_t t, double sec);



extern void time2str(gtime_t t, char *s, int n);



extern char *time_str(gtime_t t, int n);

extern void ecef2pos(const double *r, double *pos);
