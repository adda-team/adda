/* FILE: types.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of various structures
 */
#ifndef __types_h
#define __types_h

typedef struct    /* integration parameters */
{
  double INT_EPS; /* convergence criterium */
  int JMAX;       /* maximal number of refinements */
  int K;          /* number of points to use for extra-polation */
  double min;     /* minimum */
  double max;     /* maximum */
  int Grid_size;  /* number of gridpoints */
} Parms_1D;

typedef struct    /* integration parameters */
{                 /* !!! All angles are in degrees */
  double min;     /* minimum */
  double max;     /* maximum */
  int N;          /* number of points */
  double *val;    /* values of points*/
} integr_parms;

#endif /*__types_h*/
