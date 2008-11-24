/* FILE: types.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of various structures
 */
#ifndef __types_h
#define __types_h

typedef struct    /* integration parameters */
{
  double eps;     /* convergence criterium */
  int Jmax;       /* maximal number of refinements */
  int K;          /* number of points to use for extra-polation */
  double min;     /* minimum */
  double max;     /* maximum */
  int Grid_size;  /* number of gridpoints */
  int equival;    /* whether max and min points are equivalent */
} Parms_1D;

typedef struct    /* values of angles */
{                 /* !!! All angles are in degrees */
  double min;     /* minimum; for convenience (not really needed) */
  double max;     /* maximum; for convenience (not really needed) */
  int N;          /* number of points */
  double *val;    /* values of points*/
} angle_set;

typedef struct     /* integration parameters */
{                  /* !!! All angles are in degrees */
  int type;        /* if pairs are used or grid */
  int N;           /* total number of pairs (grid points) */
  angle_set theta; /* values of theta */
  angle_set phi;   /* values of phi */
} scat_grid_angles;

#endif /*__types_h*/
