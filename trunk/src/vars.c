/* File: vars.c
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
#include "vars.h" // corresponding header

// basic variables
int boxX,boxY,boxZ;       // sizes of box enclosing the particle
size_t boxXY;             // boxX*boxY, used for indexing
double gridspace;         // dipole size (d)
double dipvol;            // dipole volume
double kd;                // k*d=2*PI/dpl
double ka_eq;             // volume-equivalent size parameter
double inv_G;             // inverse of equivalent cross section
double WaveNum;           // wavenumber of incident light
double * restrict DipoleCoord;      // vector to hold the coordinates of the dipoles
double memory;            // total memory usage in bytes
double memPeak;           // peak memory usage in bytes
enum inter IntRelation;   // type of formula for interaction term
enum pol PolRelation;     // type of formula for self-term (polarizability relation)
enum beam beamtype;       // type of incident beam

// symmetries (in particle reference frame)
	// symmetries of reflection relative to the planes perpendicular to x, y, and z axes
bool symX,symY,symZ;
bool symR;         // symmetry of 90-degrees rotation about z axes

// flags (true or false)
bool prognosis;     // make a prognosis about needed ram
bool yzplane;       // Calculate the field in the yz-plane
bool scat_plane;    // Calculate the scattering in plane through ez (in lab RF), prop, incPolX
bool store_mueller; // Calculate and write Mueller matrix to file
bool all_dir;       /* Calculate the field for all directions on a theta-phi grid (internal parameter - initialized by
                       other options: calculation of Csca and asym) */
bool scat_grid;     // calculate field on a grid of scattering angles
bool phi_integr;    // integrate over the phi angle
bool reduced_FFT;   // reduced number of storage for FFT, when matrix is symmetric
bool orient_avg;    // whether to use orientation averaging
bool load_chpoint;  // whether to load checkpoint
bool beam_asym;     // whether the beam center is shifted relative to the origin
bool sh_granul;     // whether to fill one domain with granules
bool anisotropy;    // whether the scattering medium is anisotropic
bool save_memory;   // whether to sacrifice some speed for memory
bool ipr_required;  /* whether inner product in MatVec will be used by iterative solver (causes additional
                       initialization, e.g., for OpenCL) */
double propAlongZ;  // equal 0 for general incidence, and +-1 for incidence along the z-axis (can be used as flag)

// 3D vectors (in particle reference frame)
double prop_0[3],prop[3];     // incident direction (in laboratory and particle reference frame)
double incPolX[3],incPolY[3]; // incident polarizations (in particle RF)
double beam_center[3];        // coordinates of the beam center
double box_origin_unif[3];    /* coordinates of the center of the first dipole in the local computational box (after
                                 uniform distribution of non-void dipoles among all processors) */

// file info
const char * restrict directory; // directory to save data in
FILE * restrict logfile;         // file where all the information about the run is saved
int term_width;                  // width of the terminal to which ADDA produces output

// refractive index
int Nmat;  // number of different domains (for each either scalar or tensor refractive index is specified
int Ncomp; // number of components of each refractive index (1 or 3)
doublecomplex ref_index[MAX_NMAT];  // a set of refractive indexes
doublecomplex cc_sqrt[MAX_NMAT][3]; // sqrt of couple constants
doublecomplex chi_inv[MAX_NMAT][3]; // normalized inverse susceptibility: = 1/(V*chi)
unsigned char * restrict material;  // material: index for cc

// iterative solver
enum iter IterMethod; // iterative method to use
int maxiter;          // maximum number of iterations
	// the following two can't be declared restrict due to SwapPointers
doublecomplex *xvec;  // total electric field on the dipoles
doublecomplex *pvec;  // polarization of dipoles, also an auxiliary vector in iterative solvers
doublecomplex * restrict Einc;    // incident field on dipoles

// scattering at different angles
int nTheta;                        // number of angles in scattering profile
double alph_deg, bet_deg, gam_deg; // Euler angles of particle orientation in degrees
angle_set alpha_int;               // sets of angles
scat_grid_angles angles;           // angle sets for scat_grid
	// E calculated on a grid for many different directions (holds Eper and Epar) for two incident polarizations
doublecomplex * restrict EgridX,* restrict EgridY;

int nprocs;                        // total number of processes
int ringid;                        // ID of current process

size_t local_Ndip;                 // number of local total dipoles
size_t local_nvoid_Ndip;           // number of local and ...
size_t nvoid_Ndip;                 // ... total non-void dipoles
size_t local_nvoid_d0,local_nvoid_d1; // starting and ending non-void dipole for current processor
/* By defining nvoid_Ndip, local_nvoid_d0, and local_nvoid_d1 as size_t we limit the possible number of dipoles in
 * 32-bit version by 4*10^9. This can be restricting for such huge runs distributed among more than 1000 processors. But
 * we assume that using such a large number of processors implies modern cluster and hence 64-bit compilation of ADDA.
 * Anyway, a direct test for Ndip larger than the limit is made and a meaningful error message is produced if needed.
 *
 * The same limitation is implied in a few other places (like number of lines in dipole file, etc.). Definitions of
 * mat_count[] and Ndip are made as size_t due to the same reasoning.
 */
size_t local_nRows;                 // number of local rows of decomposition (only real dipoles)

// timing
TIME_TYPE Timing_EField,      // time for calculating scattered fields
          Timing_FileIO,      // time for input and output
          Timing_Integration, // time for all integrations (with precomputed values)
          tstart_main;        // starting time of the program (after MPI_Init in parallel)
          

// related to a nearby surface
bool surface;           // whether nearby surface is present
enum refl ReflRelation; // method to calculate reflected Green's tensor
doublecomplex msub;     // complex refractive index of the substrate
double inc_scale;       // scale to account for irradiance of the incident beam - 1/Re(msub)
bool msubInf;           // whether msub is infinite (perfectly reflecting surface)
double hsub;            // height of particle center above surface
/* Propagation (phase) directions of secondary incident beams above (A) and below (B) the surface (unit vectors)
 * When msub is complex, one of this doesn't tell the complete story, since the corresponding wave is inhomogeneous,
 * given by the complex wavenumber ktVec
 */
double prIncRefl[3],prIncTran[3];

#ifndef SPARSE //These variables are exclusive to the FFT mode

// position of the dipoles; in the very end of make_particle() z-components are adjusted to be relative to the local_z0
unsigned short * restrict position;

/* holds input vector (on expanded grid) to matvec. Also used as buffer in certain algorithms, that do not call MatVec
 * (this should be strictly ensured !!!)
 */
doublecomplex * restrict Xmatrix;

// auxiliary grids and their partition over processors
size_t gridX,gridY,gridZ; /* sizes of the 'matrix' X, size_t - to remove type conversions we assume that 'int' is enough
                             for it, but this declaration is to avoid type casting in calculations */
size_t gridYZ;            // gridY*gridZ
size_t smallY,smallZ;     // the size of the reduced matrix X
size_t local_Nsmall;      // number of  points of expanded grid per one processor

int local_z0,local_z1;    // starting and ending z for current processor
size_t local_Nz;          // number of z layers (based on the division of smallZ)
int local_Nz_unif;        /* number of z layers (distance between max and min values), belonging to this processor,
                             after all non_void dipoles are uniformly distributed between all processors */
int local_z1_coer;        // ending z, coerced to be not greater than boxZ (and not smaller than local_z0)
	// starting, ending x for current processor and number of x layers (based on the division of smallX)
size_t local_x0,local_x1,local_Nx;

#else //These variables are exclusive to the sparse mode

int *position; // no reason to restrict this to short in sparse mode; actually it points to a part of position_full
//in sparse mode, all coordinates must be available to each process
int * restrict position_full;

#endif //SPARSE

