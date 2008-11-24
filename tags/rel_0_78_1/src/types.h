/* FILE: types.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of various structures
 *
 * Copyright (C) 2006-2007 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __types_h
#define __types_h

/* complex numbers */
typedef double doublecomplex[2];    /* complies with FFTW3 definition */
#define RE 0
#define IM 1

typedef struct      /* integration parameters */
{
  double eps;       /* convergence criterium */
  int Jmin;         /* minimal number of refinements */
  int Jmax;         /* maximal number of refinements */
  double min;       /* minimum */
  double max;       /* maximum */
  size_t Grid_size; /* number of gridpoints */
  int equival;      /* whether max and min points are equivalent */
  int periodic;     /* whether integrated function is periodic */
} Parms_1D;

typedef struct      /* values of angles */
{                   /* !!! All angles are in degrees */
  double min;       /* minimum; for convenience (not really needed) */
  double max;       /* maximum; for convenience (not really needed) */
  size_t N;         /* number of points */
  double *val;      /* values of points*/
} angle_set;

typedef struct       /* integration parameters */
{                    /* !!! All angles are in degrees */
  int type;          /* if pairs are used or grid */
  size_t N;          /* total number of pairs (grid points) */
  angle_set theta;   /* values of theta */
  angle_set phi;     /* values of phi */
} scat_grid_angles;

#endif /*__types_h*/
