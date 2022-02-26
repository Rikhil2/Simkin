// SV
// File ALGEBRA_D.C - basic matrix algebra. Definitions and prototypes are in ALGEBRA_D.H

// Last revision:	Apr.17, 2007 SV


#include "algebra_d.h"

extern void err_hand(int mode);

// err_hand is a generic error handler, which must be present in every
// project to handle exceptions resulted from the errors in this file.
// The example of function of this type is given below:

// Error handler
void err_hand(int mode)
{
switch(mode)
        {
        // Functions in ALGEBRA.C - this file
        case MATH_ALLOC_ERR:
                puts("ERROR: memory allocation error\n");
                exit(0);
        case FILE_ERR:
                puts("ERROR: cannot open file\n");
                exit(0);
	     // Functions in ALGEBRA.C - this file
        case CHOL:
                puts("ERROR: Cholesky - input matrix must be positive definite\n");
                exit(0);
        }
}


/* Integer vector constructor */
int *ic(int size)
{
int *v;

if (!(v = (int *)malloc(size*sizeof(int)))) err_hand(MATH_ALLOC_ERR);
return(v);
}


/* Integer destructor */
void   id(int *v)
{
free(v);
}


// Integer matrix constructor */
int **imc(int xsize, int ysize)
{
int **ma;
int i;

if (!(ma = (int **)malloc(ysize*sizeof(int*)))) return(NULL);
for(i=0; i<ysize; i++)
        if (!(ma[i] = (int *)malloc(xsize*sizeof(int)))) return(NULL);
return(ma);
}


// Integer matrix destructor
void imd(int **ma, int ysize)
{
int i;

for(i=0; i<ysize; i++) free(ma[i]);
free(ma);
}


// Byte (unsigned char) matrix constructor */
unsigned char **bmc(int xsize, int ysize)
{
unsigned char **ma;
int i;

if (!(ma = (unsigned char **)malloc(ysize*sizeof(unsigned char*)))) return(NULL);
for(i=0; i<ysize; i++)
        if (!(ma[i] = (unsigned char *)malloc(xsize*sizeof(unsigned char)))) return(NULL);
return(ma);
}


// Byte (unsigned char) matrix destructor
void bmd(unsigned char **ma, int ysize)
{
int i;

for(i=0; i<ysize; i++) free(ma[i]);
free(ma);
}


/* Vector constructor */
double *vc(int size)
{
double *v;

if (!(v = (double *)malloc(size*sizeof(double)))) err_hand(MATH_ALLOC_ERR);
return(v);
}


/* Vector destructor */
void   vd(double *v)
{
free(v);
}


/* Matrix constructor. ysize - number of rows, xsize - number of columns */
double **mc(int xsize, int ysize)
{
double **ma;
int i;

if (!(ma = (double **)malloc(ysize*sizeof(double*)))) err_hand(MATH_ALLOC_ERR);
for(i=0; i<ysize; i++)
        if (!(ma[i] = (double *)malloc(xsize*sizeof(double)))) err_hand(MATH_ALLOC_ERR);
return(ma);
}


/* Matrix destructor. ysize - number of elements in the column (number of rows */
void md(double **ma, int ysize)
{
int i;

for(i=0; i<ysize; i++) free(ma[i]);
free(ma);
}


/* Copying contents of one vector of size 'size' to another */
void copyvect(double *dest, double *src, int size)
{
memcpy(dest, src, (unsigned)(size*sizeof(double)));
}


// Running-average filter of order 'order' on vector 'v' with 'num' number of
// elements performs 'pass' passes and returns result in the same location
// order must be 3, 5, 7, etc.
void run_average(int pass, int order, double *v, int num)
{
int i, j, k, n, p;
double *u;
double x;

u = vc(num);
n = order/2;
for (p=0; p<pass; p++)		// pass counter
	{
	for (i=0; i<num; i++)
		{
	   x = 0;
   	if (i<n) j=0;		// starting elements
	   else if (i>num-1-n) j=num-order;		// end elements
   	else j=i-n;		// normal elements
	   for (k=0; k<order; k++, j++) x += v[j];
   	x /= order;
	   u[i] = x;
   	}
	copyvect(v, u, num);
   }
}


/* Copying contents of one matrix (xsize x ysize) to another */
void matcopy(double **src, double **dest, int xsize, int ysize)
{
int i;

for (i=0; i<ysize; i++)
	memcpy(dest[i], src[i], (unsigned)(xsize*sizeof(double)));
}


/* Scalar product of two vectors: v1 and v2. p - result */
double   scal_prod(double *v1,  double  *v2, int size)
{
int i;
double p = 0.0;

for(i=0; i<size; i++)
        p += v1[i] * v2[i];

return(p);
}


/* WEIGHTED scalar product of two vectors: v1 and v2. ww containes weights
   p - result */
double   scal_prod_w(double *v1,  double  *v2, double *ww, int size)
{
int i;
double p = 0.0;

for(i=0; i<size; i++)
        p += (v1[i] * v2[i]) / ww[i];

return(p);
}


/* Scalar product of two vectors-columns (m and n) of matrix ma */
double   scal_prodcol(double **ma, int m, int n, int size)
{
int i;
double p = 0.0;

for(i=0; i<size; i++)
        p += ma[i][m] * ma[i][n];

return(p);
}


/* Scalar product of vector-column (m) of matrix 'ma' and vector 'v'.
   int vector 'status' describes active rows of matrix: non-active
   elements marked with 0  */
double   scal_prodcolve(double **ma, int m,  double *v, int size, int *status)
{
int i, l;
double p = 0.0;

for(i=0, l=0; i<size; i++, l++)
        {
        while (!status[l]) l++;
        p += ma[l][m] * v[i];
        };
return(p);
}


/* Scalar product of vector-column (m) of matrix 'ma' and vector 'v' */
double   scal_prodcolv(double **ma, int m,  double *v, int size)
{
int i;
double p = 0.0;

for(i=0; i<size; i++)
        p += ma[i][m] * v[i];

return(p);
}


/* Addition of v2 and v2. v_re - result */
double   *vect_add(double *v1, double  *v2, double *v_res, int size)
{
int i;

for(i=0; i<size; i++)
        v_res[i] = v1[i] + v2[i];

return(v_res);
}


/* Multiplication of vector v1 by constant value p. v_res - result */
double   *const_mul(double  *v1, double p, double *v_res, int size)
{
int i;

for(i=0; i<size; i++)
        v_res[i] = v1[i] * p;

return(v_res);
}


/* Linear combination: v_res=v1*c1 + v2*c2 */
double   *lin_com(double  *v1, double  *v2, double c1, double c2,
                  double *v_res, int size)
{
int i;

for(i=0; i<size; i++)
        v_res[i] = (v1[i]*c1 + v2[i]*c2);

return(v_res);
}


/* Effect of linear operator MAT (matrix size*size) on the
   vector vect. v_res = MAT * vect */
double   *mat_vect(double  **mat, double  *vect, double *v_res, int size)
{
int i;

for (i=0; i<size; i++)
        v_res[i] = scal_prod(mat[i], vect, size);

return(v_res);
}


/* Effect of the linear operator 'mat' (matrix size*size) on
   vector 'vect'. Integer vector 'status' shows which components of 'mat' and
   'vect' should be ommited. v_res = MAT * vect. Elements of v_res
   for which status[i]=0 remain unchanged. */
double *matvectstat(double **mat, double *vect, int *status, double *v_res, int size)
{
int i, j;
double x;

for (i=0; i<size; i++)
        {
        if (status[i])
        			for (j=0, x=0.0; j<size; j++)
               	if (status[j]) x += mat[i][j]*vect[j];
        v_res[i] =x;
        };
return(v_res);
}



/* Performs left-to-right multiplication of square matrixes. RES = A*B */
double **matrmul(double **res, double **a, double **b, int size)
{
double ge, p;
int i, j, k;

for (i=0; i<size; i++)
        for (j=0; j<size; j++)
                {
                ge = 0.0;
                for (k=0; k<size; k++)
                        {
                        p = a[i][k]*b[k][j];
                        ge = ge + p;
                        }
                res[i][j] = ge;
                };
return(res);
}



/* Computation of square matrix: A*A(trans).
   a(m*n) - input matrix, res(m*m) - resulting square matrix,
   m - number of rows, n - number of columns. Both matrices allocated elsewhere */
double **sqmatrix (double **a, double **res, int m, int n)
{
double ge;
int i, j, k;

for (i=0; i<m; i++)
	for (j=i; j<m; j++)
		{
		ge = 0.0;
		for (k=0; k<n; k++)
			ge += a[i][k]*a[j][k];
		res[i][j] = res[j][i] = ge;
		};
return(res);
}



/* Matrix-vector operation for General Linear Least Squares.
   b(row-vector) * A(trans).  a(m*n) - input matrix (not transposed), b(n)-input vector.
   c(m)-output vector.  m-number of rows, n-number of columns.
   Vectors and matrices allocated elsewhere */
void transmatvect(double **a, double *b, double *c, int m, int n)
{
int i, j;

for (i=0; i<m; i++)
	{
	c[i] = 0.0;
	for (j=0; j<n; j++)
		c[i] +=  a[i][j] * b[j] ;
	}
}



/* Another matrix-vector operation for General Linear Least Squares.
   b(m) * A(m*n) = c(n). Used to compute predicted vector c(n) using design
   matrix A(m*n) and parameter vector b(m).
   a(m*n) - input matrix (not transposed), b(m)-input vector.
   c(n)-output vector.  m-number of rows, n-number of columns.
*/
void vectmat(double **a, double *b, double *c, int m, int n)
{
int i, j;

for (i=0; i<n; i++)
	{
	c[i] = 0.0;
	for (j=0; j<m; j++)
		c[i] +=  b[j] * a[j][i];
	}
}

/* Replaces matrix elements less then EPS with 0.0 */
void zeromat(double **a, int size)
{
int i, j;

for (i=0; i<size; i++)
        for (j=0; j<size; j++)
                if (fabs(a[i][j])<EPS) a[i][j] = 0.0;
}



/* Prints matrix into the ASCII file. Debugging function. Prints lines, until
   1950 character limit is exceeded, and then ignores the rest. If matrix has
   a very small values, use 'zeromat' is highly recommended before
   this function */
int matprint(double **mat, int size, const char *name)
{
FILE *output;
int i, j;
char *li, *st;
char line[2000];

if (!(output = fopen(name, "wt"))) err_hand(FILE_ERR);
for (i=0; i<size; i++)
        {
        li = line;
        for (j=0; j<size; j++)
                {
                if ((li - line) > 1950) break;
                gcvt((double)mat[i][j], 10, li);
                st = li;
                while (*st != '.') st++; st+=6; *st='\0';
                li += strlen(li);
                strcpy(li, "   "); li += 3;
                };
        *li = '\0';
        fprintf(output, "%s\n", line);
        };
fclose(output);
return(1);
}

/* Prints vectors vec_x and vec_y into ASCII file as a two column. */
int vectprint(double *vec_x, double *vec_y, int size, const char *name)
{
FILE *output;
int i;

if (!(output = fopen(name, "wt"))) err_hand(FILE_ERR);
if (vec_x == NULL)
    for (i=0; i<size; i++) fprintf(output, "%d\t%e\n", i, vec_y[i]);
else
    for (i=0; i<size; i++) fprintf(output, "%e\t%e\n", vec_x[i], vec_y[i]);
fclose(output);

return(1);
}

//---------------- DB -----------------------------------------
// Same as the above, but for 32-bit integer arrays (with x-valus double, or NULL)
int vectprintInt(double *vec_x, int *vec_y, int size, const char *name)
{
FILE *output;
int i;

if (!(output = fopen(name, "wt"))) err_hand(FILE_ERR);
if (vec_x == NULL)
	for (i=0; i<size; i++) fprintf(output, "%d\t%d\n", i, vec_y[i]);
else
    for (i=0; i<size; i++) fprintf(output, "%e\t%d\n", vec_x[i], vec_y[i]);
fclose(output);

return(1);
}
//---------------------------------------------------------------------

//---------------- DB -----------------------------------------
// Same as the above, but for short integer arrays (with x-valus double, or NULL)
int vectprintShort(double *vec_x, short *vec_y, int size, const char *name)
{
FILE *output;
int i;

if (!(output = fopen(name, "wt"))) err_hand(FILE_ERR);
if (vec_x == NULL)
	for (i=0; i<size; i++) fprintf(output, "%d\t%hd\n", i, vec_y[i]);
else
	for (i=0; i<size; i++) fprintf(output, "%g\t%hd\n", vec_x[i], vec_y[i]);
fclose(output);

return(1);
}
//---------------------------------------------------------------------

/* Returns the number of the maximal element in array v */
int findmax(double *v, int num)
{
int i, max=0;
double m;

m = v[max];
for(i=1; i<num; i++)
        if (v[i] > m)
                {
                m = v[i];
                max = i;
                }
return(max);
}

/* Returns number of the minimal element of vector v */
int findmin(double *v, int num)
{
int i, min=0;
double m;

m = v[min];
for(i=1; i<num; i++)
        if (v[i] < m)
                {
                m = v[i];
                min = i;
                }
return(min);
}

// draws straight line: y = a*x + b, through n points (x, y), supplied in
// the vectors x and y. Places coefficients into a and b
void fit_straight_line(double *x, double *y, int n, double *a, double *b)
{
int i;
double S, Sx, Sy, Sxx, Sxy, delta;

S = (double)n;
Sx = Sy = Sxx = Sxy = 0.0;
for (i=0; i<n; i++)
	{
   Sx += x[i];
   Sy += y[i];
   Sxx += (x[i]*x[i]);
   Sxy += (x[i]*y[i]);
   }
delta = S*Sxx - Sx*Sx;

*b = (Sxx*Sy - Sx*Sxy)/delta;
*a = (S*Sxy - Sx*Sy)/delta;
}


// draws y2 weighted straight line: y = a*x + b, through n points (x, y), supplied in
// the vectors x and y. Places coefficients into a and b
void fit_weighted_straight_line(double *x, double *y, int n, double *a, double *b)
{
int i;
double s, sx, sy, sxx, sxy, delta;

*a=*b=0.0;

// fit straight line
s=sx=sxx=sy=sxy=delta=0.0;
for (i=0; i<n; i++)
        {
        s += (y[i]*y[i]);
        sx += (x[i]*y[i]*y[i]);
        sxx += (x[i]*x[i]*y[i]*y[i]);
        sy += (y[i]*y[i]*y[i]);
        sxy += (x[i]*y[i]*y[i]*y[i]);
        };

delta = s*sxx - sx*sx;
if (delta==0.0) return;

*b = (sxx*sy - sx*sxy)/delta;
*a = (s*sxy - sx*sy)/delta;
}


// draws straight line y = a*x + b through two points: p1(x1, y1) and p2(x2, y2)
// Places coefficients into a and b.
void draw_straight_line(double *p1, double *p2, double *a, double *b)
{
double delta_x;

if (p1[0]==p2[0])
        {
        *a = *b = 0.0;
        return;
        }
delta_x = p1[0]-p2[0];
*a = (p1[1]-p2[1])/delta_x;
*b = (p1[0]*p2[1]-p2[0]*p1[1])/delta_x;
}


// finds common dividers in a fraction numer/denom
void div_fract(int *numer, int *denom)
{
int i, smaller;

smaller = *denom;
if ((*numer) < (*denom)) smaller = *numer;
for (i=2; i<=smaller; )
	{
	if (whole_int_divisor(*denom, i) && whole_int_divisor(*numer, i))
   	{
      (*numer) = (*numer)/i;
      (*denom) = (*denom)/i;
      smaller /= i;
      }
   else i++;
   }
}


// returns 1 if 'divisor' is the whole integer divisor of 'number'
// both values must be positive integers
int whole_int_divisor(int number, int divisor)
{
double x;
int intdiv;

if (divisor<=0 || number<=0) return(0);
if (divisor > number) return(0);
if (divisor==1 || divisor==number) return(1);

x = (double)number/(double)divisor;
intdiv = number/divisor;
if ((x-(double)intdiv)==0.0) return(1);
else return(0);
}


// returns 1 if 'divisor' is the whole integer divisor of 'number'
// both values must be positive doubles
int whole_dbl_divisor(double number, double divisor)
{
double x, intdiv;

if (divisor<=0.0 || number<=0.0) return(0);
if (divisor > number) return(0);
if (divisor==1.0 || divisor==number) return(1);

x = number/divisor;
intdiv = ceil(x);
if ((x-intdiv)==0.0) return(1);
else return(0);
}


// finds an average amplitude of a signal
double average_ampl(double *v, int n, double dc)
{
int i;
double a;

a = 0.0;
for (i=0; i<n; i++) a+=fabs(v[i]-dc);
a /= n;

return(a);
}

