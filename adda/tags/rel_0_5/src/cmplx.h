#ifndef __cmplx_h
#define __cmplx_h

#define MIN(A,B) (((A) > (B)) ? (B) : (A))
#define MAX(A,B) (((A) < (B)) ? (B) : (A))

#define FAST inline  /* should work faster with inline */


typedef struct {
        double r;              /* real part      */
        double i;              /* imaginary part */
}       doublecomplex;                /* double-precision Complex */

/*some basic complex operations*/

double cAbs2(doublecomplex a);
void cAdd(doublecomplex a,doublecomplex b,doublecomplex *c);
void cSquare(doublecomplex a,doublecomplex *b);
void cMultReal(double a,doublecomplex b,doublecomplex *c);
void cMult_i(doublecomplex *c);
void cMult(doublecomplex a,doublecomplex b,doublecomplex *c);
double cMultConRe(doublecomplex a,doublecomplex b);
double cMultConIm(doublecomplex a,doublecomplex b);
void cLinComb(doublecomplex a,doublecomplex b, double  c1, double c2, doublecomplex *c);
void cInv(doublecomplex a,doublecomplex *b);
void cDiv(doublecomplex a,doublecomplex b,doublecomplex *c);
void cExp(double arg, doublecomplex *c);



/*some basic vector operations*/

void MultScal(double a,double *b,double *c);
void vMult(double *a,double *b,double *c);
void cScalMultRVec(double *ComplexVector, doublecomplex ComplexScalar, doublecomplex *Result);
double cvNorm2(doublecomplex *a);
double DotProd(double *a, double *b);
void crDotProd(doublecomplex *ComplexVector, double *RealVector, doublecomplex *Result);  
void LinComb(double *Vector1,double *Vector2, double Constant1, double Constant2, double *Result);
double TrSym(double *a);
double QuadForm(double *matr, double *vec);
void MatrVec(double matr[3][3], double *vec, double *res);
void permutate(double *vec, int *ord);
void permutate_i(int *vec, int *ord);

/* all vectors should be of length 3 */

/* auxillary functions */
double deg2rad(double deg);
double rad2deg(double rad);

#endif /*__cmplx_h*/