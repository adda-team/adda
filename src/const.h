/* FILE: const.h
 * AUTH: Maxim Yurkin
 * DESCR: All the constants used by ADDA code
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __const_h
#define __const_h

/* version number (string) */
#define ADDA_VERSION "0.73"

/* basic constants */
#define FALSE        0
#define TRUE         1    /* better not use for comparison, since true may be != TRUE
                                use only for assignment */
#define UNDEF        -1   /* should be used only for variables,
                             which are naturally non-negative */

/* specify to inline some functions */
/* if problems with compiler change to "static" */
#define INLINE static __inline

/* simple functions */
#define MIN(A,B) (((A) > (B)) ? (B) : (A))
#define MAX(A,B) (((A) < (B)) ? (B) : (A))

/* complex numbers */
typedef double doublecomplex[2];    /* complies with FFTW3 definition */
#define RE 0
#define IM 1

/* math constants */
#define PI           3.141592653589793
#define TWO_PI       6.283185307179586
#define ONE_THIRD    0.333333333333333

/* sizes of some arrays */
#define MAX_NMAT      10   /* maximum different materials (<128) */
#define MAX_N_SH_PARMS 25  /* maximum shape parameters */

/* shape types */
#define SH_SPHERE     0        /* sphere */
#define SH_BOX        1        /* box (may be rectangular) */
#define SH_PRISMA     3        /* prisma (triangular) */
#define SH_LINE       4        /* line with width of one dipole */
#define SH_COATED     5        /* coated sphere */
#define SH_SPHEREBOX  6        /* sphere in a box */
#define SH_RBC	      8        /* Red Blood Cell */
#define SH_ELLIPSOID  11       /* general ellipsoid */
#define SH_SDISK_ROT  14       /* disc cut of a sphere */
#define SH_CYLINDER   18       /* cylinder */
#define SH_READ       20       /* read from file */

/* which way to calculate coupleconstant */
#define POL_CM       0     /* Clausius Mossotti */
#define POL_RR       1     /* Radiative Reaction correction */
#define POL_LDR      2     /* Lattice Dispersion Relation */
#define POL_CLDR     3     /* Corrected Lattice Dispersion Relation */
#define POL_SO	     4     /* Second Order formulation */

/* how to calculate scattering quantities */
#define SQ_DRAINE    0     /* classical, as Draine */
#define SQ_SO        1     /* Second Order formulation */

/* how to calculate interaction term */
#define G_POINT_DIP  0     /* as point dipoles */
#define G_SO         1     /* Second Order formulation */

/* ldr constants */
#define LDR_B1       1.8915316
#define LDR_B2      -0.1648469
#define LDR_B3       1.7700004

/* 2nd_order constants */
#define SO_B1        1.5867182
#define SO_B2        0.13488017
#define SO_B3        0.11895826

/* iterative methods; see iterative.c for info */
#define IT_CGNR      0
#define IT_BICGSTAB  1
#define IT_BICG_CS   2
#define IT_QMR_CS    3

/* type of E field calculation */
#define CE_NORMAL    0    /* normal */
#define CE_PARPER    1    /* use symmetry to calculate both incident polarizations
                             from one calculation of internal fields */

/* path and size of tables */
#define TAB_PATH     "tables"
#define TAB_SIZE     142
#define TAB_RMAX     10

/* beam types */
#define B_PLANE      0
#define B_LMINUS     1
#define B_DAVIS1     2
#define B_DAVIS3     3
#define B_DAVIS5     4
#define B_BARTON1    5
#define B_BARTON3    6
#define B_BARTON5    7
#define B_WEIRD      8
#define B_BUGGY      9

/* types of scattering grid */
#define SG_GRID      0      /* grid of angles */
#define SG_PAIRS     1      /* set of independent pairs */
/* types of angles set */
#define SG_RANGE     0      /* range with uniformly spaced points */
#define SG_VALUES    1      /* any set of values */

/* types of phi_integr (should be different one-bit numbers) */
#define PHI_UNITY    1     /* just integrate */
#define PHI_COS2     2     /* integrate with cos(2*phi) */
#define PHI_SIN2     4     /* integrate with sin(2*phi) */
#define PHI_COS4     8     /* integrate with cos(4*phi) */
#define PHI_SIN4     16    /* integrate with sin(4*phi) */

/* types of checkpoint (to save) */
#define CHP_NONE     0  /* do not save checkpoint */
#define CHP_NORMAL   1  /* save checkpoint if not finished in time and exit */
#define CHP_REGULAR  2  /* save checkpoints in regular time intervals
                             (until finished or halted) */
#define CHP_ALWAYS   3  /* save checkpoint either if finished or time elapsed
                             and calculate all the scattering quantities */

/* return values for functions */
#define CHP_EXIT    -2        /* exit after saving checkpoint */

/**************   Global Defines and Data structures (all for LogError) *****************/

#define POSIT __FILE__, __LINE__    /* position of the error in source code */

#define ALL 0       /* each processor may report this error */
#define ONE 1       /* only root processor reports an error */

/* error codes; only 2 of them are actually used */
#define EC_MASK  0xF0000000
#define EC_FATAL 0xE0000000
#define EC_CRIT  0xD0000000
#define EC_ERROR 0xC0000000    /* error */
#define EC_WARN  0xB0000000    /* warning */
#define EC_DEBUG 0xA0000000
#define EC_INFO  0x90000000
#define EC_MESS  0x80000000

#endif /*__const_h*/
