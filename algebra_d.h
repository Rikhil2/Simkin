#include <stdio.h>
#include <stdlib.h>
#include <mem.h>
#include <math.h>
#include <string.h>

#define EPS 1.0e-10

/* Definitions for v_err error handler (in ALGEBRA.C) */
#define MATH_ALLOC_ERR		-1000
#define FILE_ERR			-2000
#define CHOL				-3000
#define SVD_ERR				-4000

#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))


/******************** from ALGEBRA_D.C **********************/
#if defined(__cplusplus)
    extern "C" {
#endif
    int *ic(int size);
    void id(int *v);
    int **imc(int xsize, int ysize);
    void imd(int **ma, int ysize);
    unsigned char **bmc(int xsize, int ysize);
    void bmd(unsigned char **ma, int ysize);
    double *vc(int size);
    void vd(double *v);
    double **mc(int xsize, int ysize);
    void md(double **ma, int ysize);
    void copyvect(double *dest, double *src, int size);
    void matcopy(double **src, double **dest, int xsize, int ysize);
    double scal_prod(double *v1,  double  *v2, int size);
    double scal_prod_w(double *v1,  double  *v2, double *ww, int size);
    double scal_prodcol(double **ma, int m,  int n, int size);
    double scal_prodcolve(double **ma, int m, double *v, int size, int *status);
    double scal_prodcolv(double **ma, int m,  double *v, int size);
    double *vect_add(double *v1, double  *v2, double *v_res, int size);
    double *const_mul(double  *v1, double p, double *v_res, int size);
    double *lin_com(double *v1, double *v2, double c1, double c2,
                                                    double *v_res, int size);
    double *mat_vect(double  **mat, double  *vect, double *v_res, int size);
    double *matvectstat(double **mat, double *vect, int *status,
                                                    double *v_res, int size);
    double **matrmul(double **res, double **a, double **b, int size);
    void zeromat(double **a, int size);
    int matprint(double **mat, int size, const char *name);
    int vectprint(double *vec_x, double *vec_y, int size, const char *name);
    int vectprintInt(double *vec_x, int *vec_y, int size, const char *name);
    int vectprintShort(double *vec_x, short *vec_y, int size, const char *name);

    int findmax(double *v, int num);
    int findmin(double *v, int num);
    int matprint(double **mat, int size, const char *name);
    double **sqmatrix (double **a, double **res, int m, int n);
    void transmatvect(double **a, double *b, double *c, int m, int n);
    void vectmat(double **a, double *b, double *c, int m, int n);
    void fit_straight_line(double *x, double *y, int n, double *a,double *b);
    void fit_weighted_straight_line(double *x, double *y, int n, double *a, double *b);
    void draw_straight_line(double *p1, double *p2, double *a, double *b);
    double average_ampl(double *v, int n, double dc);
    void run_average(int pass, int order, double *v, int num);
    void div_fract(int *numer, int *denom);
    int whole_int_divisor(int number, int divisor);
    int whole_dbl_divisor(double number, double divisor);

/*************************** from CHOLESK_D.C ***********************/
    void svdcmp(double **a, int m, int n, double *w, double **v, double *rv1);
    void svbksb(double **u, double *w, double **v, int m, int n, double *b,
                                                        double *x, double *vec);
    int checkeigen(double *w, double minimum, int size);
    double **invsvd(double **q, double **zeta, double *eig, int parnum);
    int choldc(double **a, int n, double *p);
    int cholsl(double **a, int n, double *p, double *b, double *x);
    double **invchol(double **a, int n, double *p, double **aa);

     #if TWOFREQCORR
     void lpf(double *X, int numpoints, double f);
     #endif

#if defined(__cplusplus)
    }
#endif
