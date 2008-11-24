/*
 *	types.h
 *      contains the data-types used by the Romberg-routine
 *      for internal use and communication of data
 */

typedef struct    /* integration parameters */
{
  REAL INT_EPS;   /* convergence criterium */
  int JMAX;       /* maximal number of refinements */
  int K;          /* number of points to use for extra-polation */
  REAL min;       /* minimum */
  REAL max;       /* maximum */
  int Grid_size;  /* number of gridpoints */
} Parms_1D;

typedef struct    /* multiple dimensional vector (crude polymorphism) */
{
  REAL *x;
} Rvector;
