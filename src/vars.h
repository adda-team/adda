/* File: vars.h
 * $Date::                            $
 * Descr: all the global variables are declared here
 *
 *        'Global' means used in three or more source files. Variables that are used in only two source files are called
 *        'semi-global' and not listed here. They are defined in one file and referenced with 'extern' in another one.
 *
 * Copyright (C) 2006-2014 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __vars_h
#define __vars_h

// project headers
#include "const.h"   // for MAX_NMAT
#include "timing.h"  // for TIME_TYPE
#include "types.h"   // for doublecomplex, angle_set, scat_grid_angles
// system headers
#include <stdbool.h> // for bool
#include <stdio.h>   // for FILE and size_t
#include <time.h>    // for time_t

// basic variables

extern int boxX,boxY,boxZ;
extern size_t boxXY;
extern double gridspace,dipvol,kd,ka_eq,inv_G,WaveNum;
extern double * restrict DipoleCoord;

extern double memory,memPeak;
extern enum inter IntRelation;
extern enum pol PolRelation;
extern enum beam beamtype;

// symmetries
extern bool symX,symY,symZ,symR;

// flags
extern bool prognosis,yzplane,scat_plane,store_mueller,all_dir,scat_grid,phi_integr,sh_granul,reduced_FFT,orient_avg,
	load_chpoint,beam_asym,anisotropy,save_memory,ipr_required;
extern double propAlongZ;

// 3D vectors
extern double prop_0[3],prop[3],incPolX[3],incPolY[3],beam_center[3],box_origin_unif[3];

// file info
extern const char * restrict directory;
extern FILE * restrict logfile;
extern int term_width;

// refractive index
extern int Nmat,Ncomp;
extern doublecomplex ref_index[MAX_NMAT];
extern doublecomplex cc_sqrt[MAX_NMAT][3];
extern doublecomplex chi_inv[MAX_NMAT][3];
extern unsigned char * restrict material;

// iterative solver
extern enum iter IterMethod;
extern int maxiter;
extern doublecomplex *xvec,*pvec,* restrict Einc;

// scattering at different angles
extern int nTheta;
extern double alph_deg, bet_deg, gam_deg;
extern angle_set alpha_int;
extern scat_grid_angles angles;
extern doublecomplex * restrict EgridX,* restrict EgridY;

extern int nprocs,ringid;

extern size_t local_Ndip,local_nvoid_Ndip,local_nRows,local_nvoid_d0,local_nvoid_d1,nvoid_Ndip;

// timing
extern TIME_TYPE Timing_EField,Timing_FileIO,Timing_Integration,tstart_main;

// related to a nearby surface
extern bool surface,msubInf;
extern enum refl ReflRelation;
extern doublecomplex msub;
extern double inc_scale,hsub,prIncRefl[3],prIncTran[3];

#ifndef SPARSE //These variables are exclusive to the FFT mode

extern unsigned short * restrict position;

extern doublecomplex * restrict Xmatrix;

// auxiliary grids and their partition over processors
extern size_t gridX,gridY,gridZ;
extern size_t gridYZ;
extern size_t smallY,smallZ;
extern size_t local_Nsmall;
extern int local_z0,local_z1,local_z1_coer,local_Nz_unif;
extern size_t local_Nz,local_x0,local_x1,local_Nx;

#else //These variables are exclusive to the sparse mode

extern int *position;
extern int * restrict position_full;

#endif //SPARSE

#endif // __vars_h
