/* FILE: vars.h
 * AUTH: Maxim Yurkin
 * DESCR: All the global variables are declared here.
 *        Global means: used in three or more source files.
 *        Variables that are used in only two source files are calles 'semi-global'
 *           and not listed here. They are defined in one file and referenced with
 *           'extern' in another one.
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __vars_h
#define __vars_h

#include <stdio.h>    /* for FILE and size_t */
#include <time.h>     /* for time_t */
#include "const.h"    /* for MAX_NMAT */
#include "types.h"    /* for doublecomplex, angle_set, scat_grid_angles */
#include "timing.h"   /* for TIME_TYPE */

/* basic variables */
extern int boxX,boxY,boxZ;
extern double gridspace,kd,ka_eq,inv_G,WaveNum;
extern double *DipoleCoord;
extern unsigned short *position;
extern size_t memory;
extern int IntRelation;
extern int beamtype;

/* symmetries */
extern int symX,symY,symZ,symR;

/* flags (TRUE or FALSE) */
extern int prognose,yzplane,all_dir,scat_grid,phi_integr,sh_granul,
           reduced_FFT,orient_avg,load_chpoint,NoSymmetry,beam_asym,anisotropy;
/* 3D vectors */
extern double prop[3],incPolX[3],incPolY[3],beam_center[3];

/* file info */
extern char directory[];
extern FILE *logfile;

/* refractive index */
extern int Nmat,Ncomp;
extern doublecomplex ref_index[MAX_NMAT];
extern doublecomplex cc_sqrt[MAX_NMAT][3];
extern unsigned char *material;

/* iterative solver */
extern int IterMethod,maxiter;
extern doublecomplex *xvec,*pvec,*Einc;

/* scattering at different angles */
extern int nTheta;
extern double alph_deg, bet_deg, gam_deg;
extern angle_set alpha_int;
extern scat_grid_angles angles;
extern doublecomplex *EgridX,*EgridY;
extern double *Egrid_buffer;

/* checkpoint */
extern int chp_type;
extern time_t chp_time;
extern char chp_dir[];

/* auxillary grids and their partition over processors */
extern int gridX,gridY,gridZ,gridYZ;

extern int smallY,smallZ,local_Nsmall;
extern int nprocs,ringid;
extern int local_z0,local_z1,local_Nz,local_z1_coer,
           local_x0,local_x1,local_Nx,
           local_d0,local_d1,local_Ndip,
           local_nvoid_Ndip,nvoid_Ndip;
extern size_t nlocalRows;

/* timing */
extern time_t wt_start,last_chp_wt;
extern TIME_TYPE Timing_EField,Timing_FileIO,Timing_Integration,Timing_OneIterComm,tstart_main;

#endif /*__vars_h*/
