/* FILE: vars.c
 * AUTH: Maxim Yurkin
 * DESCR: All the global variables are defined here
 *        Global means: used in three or more source files.
 *        Variables that are used in only two source files are calles 'semi-global'
 *           and not listed here. They are defined in one file and referenced with
 *           'extern' in another one.
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <time.h>
#include "cmplx.h"
#include "const.h"
#include "types.h"

/* basic variables */
int boxX,boxY,boxZ;  /* sizes of box enclosing the particle */
double gridspace;    /* inter-dipole distance */
double kd;           /* k*d=2*PI/dpl */
double ka_eq;        /* volume-equivalent size parameter */
double inv_G;        /* inverse of equivalent cross section */
double WaveNum;      /* wavenumber of incident light */
double *DipoleCoord; /* vector to hold the coordinates of the dipoles */
short int *position; /* position of the dipoles */
size_t memory;       /* total memory usage in bytes */
int IntRelation;     /* type of formula for interaction term */

/* symmetries */
int symX,symY,symZ;   /* symmetries of reflection relative to the planes
                         perpendicular to x, y, and z axes. Only Y is actually used */
int symR;             /* symmetry of 90 deg. rotation about z axes */

/* flags (TRUE or FALSE) */
int prognose;     /* make a prognose about needed ram */
int yzplane;      /* Calculate the field in the yz-plane */
int all_dir;      /* Calculate the field for all directions on a theta-phi grid (internal
                     parameter - initialized by other options: calculation of Csca and asym) */
int scat_grid;    /* calculate field on a grid of scattering angles */
int phi_integr;   /* integrate over the phi angle */
int reduced_FFT;  /* reduced number of storage for FFT, when matrix is symmetric */
int orient_avg;   /* whether to use orientation averaging*/
int load_chpoint; /* whether to load checkpoint */

/* 3D vectors */
double prop[3];                /* incident direction (in particle reference frame) */
double incPolX[3],incPolY[3];  /* incident polariztions (in particle RF) */

/* file info */
char directory[MAX_DIRNAME];   /* directory to save data in */
FILE *logfile;                 /* file where all the information about the run is saved */

/* refractive index */
int Nmat;                           /* number of different refractive indexes */
doublecomplex ref_index[MAX_NMAT];  /* a set of refractive indexes */
doublecomplex cc_sqrt[MAX_NMAT][3]; /* sqrt of couple constants */
unsigned char *material;            /* material: index for cc */

/* iterative solver */
int IterMethod;           /* iterative method to use */
int maxiter;              /* maximum number of iterations */
doublecomplex *pvec;      /* polarization of dipoles */
doublecomplex *Einc;      /* incident field on dipoles */

/* scattering at different angles */
int nTheta;	                   /* number of angles in scattering profile */
double alph_deg, bet_deg, gam_deg; /* Euler angles of particle orientation in degrees */
angle_set alpha_int;               /* sets of angles */
scat_grid_angles angles;       /* angle sets for scat_grid */
doublecomplex *EgridX,*EgridY; /* E calculated on a grid for many different directions
                                 (holds Eper and Epar) for two incident polarizations */
double *Egrid_buffer;          /* buffer to accumulate Egrid */

/* checkpoint */
int chp_type;              /* type of checkpoint (to save) */
time_t chp_time;           /* time of checkpoint (in sec) */
char chp_dir[MAX_DIRNAME]; /* directory name to save/load checkpoint */

/* auxillary grids and their partition over processors */
int gridX,gridY,gridZ,gridYZ;   /* sizes of the 'matrix' X */

int smallY,smallZ;              /* the size of the reduced matrix X */
int local_Nsmall;               /* number of  points of expanded grid per one processor */

int nprocs;                     /* total number of processes */
int ringid;                     /* id of current process */
int local_z0,local_z1,local_Nz; /* starting, ending z for current processor and
                                      number of z layers (based on the division of smallZ) */
int local_z1_coer;              /* ending z, coerced to be not greater than boxZ */
int local_x0,local_x1,local_Nx; /* starting, ending x for current processor and
                                      number of x layers (based on the division of smallX) */
int local_d0,local_d1,local_Ndip; /* starting, ending realdipoles and local quantity */
int local_nvoid_Ndip,nvoid_Ndip;  /* number of local and total non-void dipoles */
size_t nlocalRows;   /* number of local rows of decomposition (only real dipoles) */

/* timing */
time_t wt_start,             /* starting wall time */
       last_chp_wt;          /* wall time of the last checkpoint */
clock_t Timing_OneIterComm,  /* communication time during one iteration */
        Timing_EField,       /* time for calculating scattered fields */
        Timing_FileIO,       /* time for input and output */
        Timing_Integration;  /* time for all integrations (with precomputed values) */

