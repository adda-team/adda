/* File: param.c
 * $Date::                            $
 * Descr: initialization, parsing and handling of input parameters; also printout general information; contains file
 *        locking routines
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
#include "const.h" // keep this first
#include "param.h" // corresponding header
// project headers
#include "cmplx.h"
#include "comm.h"
#include "crosssec.h"
#include "debug.h"
#include "fft.h"
#include "function.h"
#include "io.h"
#include "oclcore.h"
#include "os.h"
#include "parbas.h"
#include "vars.h"
// system headers
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef CLFFT_AMD
/* One can also include clFFT.h (the only recommended public header), which can be redundant, but more portable.
 * However, version macros are not documented anyway (in the manual), so there seem to be no perfectly portable way to
 * obtain them.
 */
#	include <clFFT.version.h>
#endif

#ifdef OCL_BLAS
/* In contrast to clFFT, here including main header (clBLAS.h) is not an option, since it doesn't include the
 * following header.
 */
#	include <clBLAS.version.h>
#endif

#ifndef NO_SVNREV
#	include "svnrev.h" // for SVNREV, this file is automatically created during compilation
#endif

// definitions for file locking
#ifdef USE_LOCK
#	ifdef WINDOWS
#		define FILEHANDLE HANDLE
#	elif defined(POSIX)
#		include <unistd.h>
#		include <fcntl.h>
#		ifdef LOCK_FOR_NFS
#			include <errno.h> // for error handling of fcntl call
#		endif
#		define FILEHANDLE int
#	else
#		error "Unknown operation system. Creation of lock files is not supported."
#	endif
#	define LOCK_WAIT 1 // in seconds
#	define MAX_LOCK_WAIT_CYCLES 60
#	define ONLY_FOR_LOCK
#else
#	define FILEHANDLE int
#	define ONLY_FOR_LOCK ATT_UNUSED
#endif

// GLOBAL VARIABLES

opt_index opt; // main option index; it is also defined as extern in param.h

// SEMI-GLOBAL VARIABLES

// defined and initialized in crosssec.c
extern const char avg_string[];
// defined and initialized in GenerateB.c
extern const char beam_descr[];
// defined and initialized in make_particle.c
extern const bool volcor_used;
extern const char *sh_form_str1,*sh_form_str2;
#ifndef SPARSE
extern const int gr_N;
extern const double gr_vf_real;
#endif //SPARSE
extern const size_t mat_count[];

// used in CalculateE.c
bool store_int_field; // save full internal fields to text file
bool store_dip_pol;   // save dipole polarizations to text file
bool store_beam;      // save incident beam to file
bool store_scat_grid; // Store the scattered field for grid of angles
bool calc_Cext;       // Calculate the extinction cross-section - always do
bool calc_Cabs;       // Calculate the absorption cross-section - always do
bool calc_Csca;       // Calculate the scattering cross-section by integration
bool calc_vec;        // Calculate the unnormalized asymmetry-parameter
bool calc_asym;       // Calculate the asymmetry-parameter
bool calc_mat_force;  // Calculate the scattering force by matrix-evaluation
bool store_force;     // Write radiation pressure per dipole to file
bool store_ampl;      // Write amplitude matrix to file
int phi_int_type; // type of phi integration (each bit determines whether to calculate with different multipliers)
// used in calculator.c
bool avg_inc_pol;                 // whether to average CC over incident polarization
double polNlocRp;                 // Gaussian width for non-local polarizability
const char *alldir_parms;         // name of file with alldir parameters
const char *scat_grid_parms;      // name of file with parameters of scattering grid
// used in crosssec.c
double incPolX_0[3],incPolY_0[3]; // initial incident polarizations (in lab RF)
enum scat ScatRelation;           // type of formulae for scattering quantities
// used in GenerateB.c
int beam_Npars;
double beam_pars[MAX_N_BEAM_PARMS]; // beam parameters
opt_index opt_beam;                 // option index of beam option used
const char *beam_fnameY;            // names of files, defining the beam (for two polarizations)
const char *beam_fnameX;
// used in interaction.c
double igt_lim; // limit (threshold) for integration in IGT
double igt_eps; // relative error of integration in IGT
double nloc_Rp; // Gaussian width for non-local interaction
bool InteractionRealArgs; // whether interaction (or reflection) routines can be called with real arguments
// used in io.c
char logfname[MAX_FNAME]=""; // name of logfile
// used in iterative.c
double iter_eps;           // relative error to reach
enum init_field InitField; // how to calculate initial field for the iterative solver
const char *infi_fnameY;   // names of files, defining the initial field (for two polarizations)
const char *infi_fnameX;
bool recalc_resid;         // whether to recalculate residual at the end of iterative solver
enum chpoint chp_type;     // type of checkpoint (to save)
time_t chp_time;           // time of checkpoint (in sec)
char const *chp_dir;       // directory name to save/load checkpoint
// used in make_particle.c
enum sh shape;                   // particle shape definition
int sh_Npars;                    // number of shape parameters
double sh_pars[MAX_N_SH_PARMS];  // storage for shape parameters
double sizeX;                    // size of particle along x-axis
double dpl;                      // number of dipoles per lambda (wavelength)
double lambda;                   // incident wavelength (in vacuum)
int jagged;                      // size of big dipoles, used to construct a particle
const char *shape_fname;         // name of file, defining the shape
const char *save_geom_fname;     // geometry file name to save dipole configuration
const char *shapename;           // name of the used shape
bool volcor;                     // whether to use volume correction
bool save_geom;                  // whether to save dipole configuration in .geom file
opt_index opt_sh;                // option index of shape option used
double gr_vf;                    // granules volume fraction
double gr_d;                     // granules diameter
int gr_mat;                      // domain number to granulate
double a_eq;                     // volume-equivalent radius of the particle
enum shform sg_format;           // format for saving geometry files
bool store_grans;                // whether to save granule positions to file

// LOCAL VARIABLES

#define GFORM_RI_DIRNAME "%.4g" // format for refractive index in directory name

static const char *run_name;    // first part of the dir name ('run' or 'test')
static const char *avg_parms;   // name of file with orientation averaging parameters
static const char *exename;     // name of executable (adda, adda.exe, adda_mpi,...)
static int Nmat_given;          // number of refractive indices given in the command line
static enum sym sym_type;       // how to treat particle symmetries
/* The following '..._used' flags are, in principle, redundant, since the structure 'options' contains the same flags.
 * However, the latter can't be easily addressed by the option name (a search over the whole options is required).
 */
static bool prop_used;          // whether '-prop ...' was used in the command line
static bool orient_used;        // whether '-orient ...' was used in the command line
static bool yz_used;            // whether '-yz ...' was used in the command line
static bool scat_plane_used;    // whether '-scat_plane ...' was used in the command line
static bool int_surf_used;      // whether '-int_surf ...' was used in the command line

/* TO ADD NEW COMMAND LINE OPTION
 * If you need new variables or flags to implement effect of the new command line option, define them here. If a
 * variable is used only in this source file, define it as local (static). If it is used in one another file, define
 * them as semi-global, i.e. define them here under corresponding 'used in ...' line and put 'extern' declaration in
 * corresponding source file under 'defined and initialized in param.c' line. If it is used in this and two or more
 * files, define it as global, by putting definition in vars.c and 'extern' declaration in vars.h.
 */

// structure definitions
struct subopt_struct {
	const char *name;  // name of option
	const char *usage; // how to use (argument list)
	const char *help;  // help string
	const int narg;    /* possible number of arguments; UNDEF -> should not be checked; may contain also some special
	                      negative codes, like FNAME_ARG. Currently it is assumed that all arguments are float (if no
	                      special argument is specified), but it can be changed if will become a limitation. */
	const int type;    // type of suboption
};
struct opt_struct {
	const char *name;                // name of option
	void (*func)(int Narg,char **argv); // pointer to a function, parsing this parameter
	bool used;                       // flag to indicate, if the option was already used
	const char *usage;               // how to use (argument list)
	const char *help;                // help string
	const int narg;                  // possible number of arguments; UNDEF -> should not be checked
	const struct subopt_struct *sub; // suboptions
};
// const string for usage of ADDA
static const char exeusage[]="[-<opt1> [<args1>] [-<opt2> <args2>]...]]";

/* initializations of suboptions (comments to elements of subopt_struct are above); Contrary to 'options', suboptions
 * are defined as null-terminated array, because they may be referenced not directly by their names but also as
 * options[i].sub. In the latter case macro LENGTH can't be used to estimate the length of the array. So getting this
 * length is at least nontrivial (e.g. can be done with some intermediate variables). So using NULL-termination seems
 * to be the easiest.
 */
static const struct subopt_struct beam_opt[]={
	{"barton5","<width> [<x> <y> <z>]","5th order approximation of the Gaussian beam (by Barton). The beam width is "
		"obligatory and x, y, z coordinates of the center of the beam (in laboratory reference frame) are optional "
		"(zero, by default). All arguments are in um. This is recommended option for simulation of the Gaussian beam.",
		UNDEF,B_BARTON5},
	{"davis3","<width> [<x> <y> <z>]","3rd order approximation of the Gaussian beam (by Davis). The beam width is "
		"obligatory and x, y, z coordinates of the center of the beam (in laboratory reference frame) are optional "
		"(zero, by default). All arguments are in um.",UNDEF,B_DAVIS3},
	{"dipole","<x> <y> <z>","Field of a unit point dipole placed at x, y, z coordinates (in laboratory reference "
		"frame). All arguments are in um. Orientation of the dipole is determined by -prop command line option."
		"Implies '-scat_matr none'. If '-surf' is used, dipole position should be above the surface.",3,B_DIPOLE},
	{"lminus","<width> [<x> <y> <z>]","Simplest approximation of the Gaussian beam. The beam width is obligatory and "
		"x, y, z coordinates of the center of the beam (in laboratory reference frame) are optional (zero, by"
		" default). All arguments are in um.",UNDEF,B_LMINUS},
	{"plane","","Infinite plane wave",0,B_PLANE},
	{"read","<filenameY> [<filenameX>]","Defined by separate files, which names are given as arguments. Normally two "
		"files are required for Y- and X-polarizations respectively, but a single filename is sufficient if only "
		"Y-polarization is used (e.g. due to symmetry). Incident field should be specified in a particle reference "
		"frame in the same format as used by '-store_beam'.",FNAME_ARG_1_2,B_READ},
	/* TO ADD NEW BEAM
	 * add a row to this list in alphabetical order. It contains: beam name (used in command line), usage string, help
	 * string, possible number of float parameters, beam identifier (defined inside 'enum beam' in const.h). Usage and
	 * help string are shown either when -h option is invoked or error is found in input command line. Usage string is
	 * one line giving a list of possible arguments or argument combinations. Do not include beam name in it, use
	 * <arg_name> to denote argument, [...] for optional arguments, and {...|...|...} for multiple options of an
	 * argument. Help string should contain general description of the beam type and its arguments. Instead of number
	 * of parameters UNDEF can be used (if beam can accept variable number of parameters, then check it explicitly in
	 * function PARSE_FUNC(beam) below) or one of FNAME_ARG family (see explanation in const.h) if beam accepts a
	 * filenames as arguments. Number of parameters should not be greater than MAX_N_BEAM_PARMS (defined in const.h).
	 */
	{NULL,NULL,NULL,0,0}
};
static const struct subopt_struct shape_opt[]={
#ifndef SPARSE
	{"axisymmetric","<filename>","Axisymmetric homogeneous shape, defined by its contour in ro-z plane of the "
		"cylindrical coordinate system. Its symmetry axis coincides with the z-axis, and the contour is read from "
		"file.",FNAME_ARG,SH_AXISYMMETRIC},
	{"bicoated","<R_cc/d> <d_in/d>","Two identical concentric coated spheres with outer diameter d (first domain), "
		"inner diameter d_in, and center-to-center distance R_cc (along the z-axis). It describes both separate and "
		"sintered coated spheres. In the latter case sintering is considered symmetrically for cores and shells.",
		2,SH_BICOATED},
	{"biellipsoid","<y1/x1> <z1/x1> <x2/x1> <y2/x2> <z2/x2>","Two general ellipsoids in default orientations with "
		"centers on the z-axis, touching each other. Their semi-axes are x1,y1,z1 (lower one, first domain) and "
		"x2,y2,z2 (upper one, second domain).",5,SH_BIELLIPSOID},
	{"bisphere","<R_cc/d> ","Two identical spheres with diameter d and center-to-center distance R_cc (along the "
		"z-axis). It describe both separate and sintered spheres.",1,SH_BISPHERE},
	{"box","[<y/x> <z/x>]","Homogeneous cube (if no arguments are given) or a rectangular parallelepiped with edges "
		"x,y,z.",UNDEF,SH_BOX},
	{"capsule","<h/d>","Homogeneous capsule (cylinder with half-spherical end caps) with cylinder height h and "
		"diameter d (its axis of symmetry coincides with the z-axis).",1,SH_CAPSULE},
	{"chebyshev","<eps> <n>","Axisymmetric Chebyshev particle of amplitude eps and order n, r=r_0[1+eps*cos(n*theta)]. "
		"eps is a real number, such that |eps|<=1, while n is a natural number",2,SH_CHEBYSHEV},
	{"coated","<d_in/d> [<x/d> <y/d> <z/d>]","Sphere with a spherical inclusion; outer sphere has a diameter d (first "
		"domain). The included sphere has a diameter d_in (optional position of the center: x,y,z).",UNDEF,SH_COATED},
	{"cylinder","<h/d>","Homogeneous cylinder with height (length) h and diameter d (its axis of symmetry coincides "
		"with the z-axis).",1,SH_CYLINDER},
	{"egg","<eps> <nu>","Axisymmetric egg shape given by a^2=r^2+nu*r*z-(1-eps)z^2, where 'a' is scaling factor. "
		"Parameters must satisfy 0<eps<=1, 0<=nu<eps.",2,SH_EGG},
	{"ellipsoid","<y/x> <z/x>","Homogeneous general ellipsoid with semi-axes x,y,z",2,SH_ELLIPSOID},
	{"line","","Line along the x-axis with the width of one dipole",0,SH_LINE},
	{"plate", "<h/d>","Homogeneous plate (cylinder with rounded side) with cylinder height h and full diameter d (i.e. "
		"diameter of the constituent cylinder is d-h). Its axis of symmetry coincides with the z-axis.",1,SH_PLATE},
	{"prism","<n> <h/Dx>","Homogeneous right prism with height (length along the z-axis) h based on a regular polygon "
		"with n sides of size 'a'. The polygon is oriented so that positive x-axis is a middle perpendicular for one "
		"of its sides. Dx is size of the polygon along the x-axis, equal to 2Ri and Rc+Ri for even and odd n "
		"respectively. Rc=a/[2sin(pi/n)] and Ri=Rc*cos(pi/n) are radii of circumscribed and inscribed circles "
		"respectively.",2,SH_PRISM},
	{"rbc","<h/d> <b/d> <c/d>","Red Blood Cell, an axisymmetric (over z-axis) biconcave homogeneous particle, which is "
		"characterized by diameter d, maximum and minimum width h, b, and diameter at the position of the maximum "
		"width c. The surface is described by ro^4+2S*ro^2*z^2+z^4+P*ro^2+Q*z^2+R=0, ro^2=x^2+y^2, P,Q,R,S are "
		"determined by the described parameters.",3,SH_RBC},
#endif // !SPARSE
	{"read","<filename>","Read a particle geometry from file <filename>",FNAME_ARG,SH_READ},
#ifndef SPARSE
	{"sphere","","Homogeneous sphere",0,SH_SPHERE},
	{"spherebox","<d_sph/Dx>","Sphere (diameter d_sph) in a cube (size Dx, first domain)",1,SH_SPHEREBOX},
#endif // !SPARSE
	/* TO ADD NEW SHAPE
	 * add a row to this list in alphabetical order. It contains: shape name (used in command line), usage string, help
	 * string, possible number of float parameters, shape identifier (defined inside 'enum sh' in const.h). Usage and
	 * help string are shown either when -h option is invoked or error is found in input command line. Usage string is
	 * one line giving a list of possible arguments or argument combinations. Do not include shape name in it, use
	 * <arg_name> to denote argument, [...] for optional arguments, and {...|...|...} for multiple options of an
	 * argument. Help string should contain general description of the shape and its arguments. Instead of number of
	 * parameters UNDEF can be used (if shape can accept variable number of parameters, then check it explicitly in
	 * function PARSE_FUNC(shape) below) or FNAME_ARG (if the shape accepts a single string argument with file name).
	 * Number of parameters should not be greater than MAX_N_SH_PARMS (defined in const.h). It is recommended to use
	 * dimensionless shape parameters, e.g. aspect ratios.
	 */
	{NULL,NULL,NULL,0,0}
};
/* TO ADD NEW COMMAND LINE OPTION
 * If a new option requires separate description of suboptions, add a static structure here (similar to already existing
 * ones). It should contain a number of rows corresponding to different suboptions (see definition of subopt_struct
 * above for details) and terminate with a NULL row.
 */

// EXTERNAL FUNCTIONS

// GenerateB.c
void InitBeam(void);

//======================================================================================================================
/* Prototypes of parsing functions; definitions are given below. defines are for conciseness. Since we use one common
 * prototype style for all parsing functions and many of the latter do not actually use the passed parameters, these
 * parameters have 'unused' attribute to eliminate spurious warnings.
 */
#define PARSE_NAME(a) parse_##a
#define PARSE_FUNC(a) static void PARSE_NAME(a)(int Narg ATT_UNUSED,char **argv ATT_UNUSED)
#define PAR(a) #a,PARSE_NAME(a),false
PARSE_FUNC(alldir_inp);
PARSE_FUNC(anisotr);
PARSE_FUNC(asym);
PARSE_FUNC(beam);
PARSE_FUNC(chp_dir);
PARSE_FUNC(chp_load);
PARSE_FUNC(chp_type);
PARSE_FUNC(chpoint);
PARSE_FUNC(Cpr);
PARSE_FUNC(Csca);
PARSE_FUNC(dir);
PARSE_FUNC(dpl);
PARSE_FUNC(eps);
PARSE_FUNC(eq_rad);
#ifdef OPENCL
PARSE_FUNC(gpu);
#endif
#ifndef SPARSE
PARSE_FUNC(granul);
#endif
PARSE_FUNC(grid);
PARSE_FUNC(h) ATT_NORETURN;
PARSE_FUNC(init_field);
PARSE_FUNC(int);
PARSE_FUNC(int_surf);
PARSE_FUNC(iter);
PARSE_FUNC(jagged);
PARSE_FUNC(lambda);
PARSE_FUNC(m);
PARSE_FUNC(maxiter);
PARSE_FUNC(no_reduced_fft);
PARSE_FUNC(no_vol_cor);
PARSE_FUNC(ntheta);
PARSE_FUNC(opt);
PARSE_FUNC(orient);
PARSE_FUNC(phi_integr);
PARSE_FUNC(pol);
PARSE_FUNC(prognosis);
PARSE_FUNC(prop);
PARSE_FUNC(recalc_resid);
#ifndef SPARSE
PARSE_FUNC(save_geom);
#endif
PARSE_FUNC(scat);
PARSE_FUNC(scat_grid_inp);
PARSE_FUNC(scat_matr);
PARSE_FUNC(scat_plane);
#ifndef SPARSE
PARSE_FUNC(sg_format);
#endif
PARSE_FUNC(shape);
PARSE_FUNC(size);
PARSE_FUNC(store_beam);
PARSE_FUNC(store_dip_pol);
PARSE_FUNC(store_force);
#ifndef SPARSE
PARSE_FUNC(store_grans);
#endif
PARSE_FUNC(store_int_field);
PARSE_FUNC(store_scat_grid);
PARSE_FUNC(surf);
PARSE_FUNC(sym);
PARSE_FUNC(test);
PARSE_FUNC(V) ATT_NORETURN;
PARSE_FUNC(vec);
PARSE_FUNC(yz);
/* TO ADD NEW COMMAND LINE OPTION
 * add a function prototype to this list. Add a line 'PARSE_FUNC(option_name);' in alphabetical order. It will be
 * expanded automatically using defines specified above.
 */

static struct opt_struct options[]={
	{PAR(alldir_inp),"<filename>","Specifies a file with parameters of the grid of scattering angles for calculating "
		"integral scattering quantities.\n"
		"Default: "FD_ALLDIR_PARMS,1,NULL},
	{PAR(anisotr),"","Specifies that refractive index is anisotropic (its tensor is limited to be diagonal in particle "
		"reference frame). '-m' then accepts 6 arguments per each domain. Can not be used with CLDR polarizability and "
		"all SO formulations.",0,NULL},
	{PAR(asym),"","Calculate the asymmetry vector. Implies '-Csca' and '-vec'",0,NULL},
	{PAR(beam),"<type> [<args>]","Sets the incident beam, either predefined or 'read' from file. All parameters of "
		"predefined beam types (if present) are floats.\n"
		"Default: plane",UNDEF,beam_opt},
	{PAR(chp_dir),"<dirname>","Sets directory for the checkpoint (both for saving and loading).\n"
		"Default: "FD_CHP_DIR,1,NULL},
	{PAR(chp_load),"","Restart a simulation from a checkpoint",0,NULL},
	{PAR(chp_type),"{normal|regular|always}",
		"Sets type of the checkpoint. All types, except 'always', require '-chpoint'.\n"
		"Default: normal",1,NULL},
	{PAR(chpoint),"<time>","Specifies the time for checkpoints in format '#d#h#m#s'. All fields are optional, numbers "
		"are integers, 's' can be omitted, the format is not case sensitive.\n"
		"Examples: 12h30M, 1D10s, 3600",1,NULL},
	{PAR(Cpr),"","Calculate the total radiation force, expressed as cross section.",0,NULL},
	{PAR(Csca),"","Calculate scattering cross section (by integrating the scattered field)",0,NULL},
	{PAR(dir),"<dirname>","Sets directory for output files.\n"
		"Default: constructed automatically",1,NULL},
	{PAR(dpl),"<arg>","Sets parameter 'dipoles per lambda', float.\n"
		"Default: 10|m|, where |m| is the maximum of all given refractive indices.",1,NULL},
	{PAR(eps),"<arg>","Specifies the stopping criterion for the iterative solver by setting the relative norm of the "
		"residual 'epsilon' to reach. <arg> is an exponent of base 10 (float), i.e. epsilon=10^(-<arg>).\n"
		"Default: 5 (epsilon=1E-5)",1,NULL},
	{PAR(eq_rad),"<arg>","Sets volume-equivalent radius of the particle in um, float. If default wavelength is used, "
		"this option specifies the volume-equivalent size parameter. Can not be used together with '-size'. Size is "
		"defined by some shapes themselves, then this option can be used to override the internal specification and "
		"scale the shape.\n"
		"Default: determined by the value of '-size' or by '-grid', '-dpl', and '-lambda'.",1,NULL},
#ifdef OPENCL
	{PAR(gpu),"<index>","Specifies index of GPU that should be used (starting from 0). Relevant only for OpenCL "
		"version of ADDA, running on a system with several GPUs.\n"
		"Default: 0",1,NULL},
#endif
#ifndef SPARSE
	{PAR(granul),"<vol_frac> <diam> [<dom_number>]","Specifies that one particle domain should be randomly filled with "
		"spherical granules with specified diameter <diam> and volume fraction <vol_frac>. Domain number to fill is "
		"given by the last optional argument. Algorithm may fail for volume fractions > 30-50%.\n"
		"Default <dom_number>: 1",UNDEF,NULL},
#endif // !SPARSE
	{PAR(grid),"<nx> [<ny> <nz>]","Sets dimensions of the computation grid. Arguments should be even integers. In most "
		"cases <ny> and <nz> can be omitted (they are automatically determined by <nx> based on the proportions of the "
		"scatterer). This command line option is not relevant when particle geometry is read from a file ('-shape "
		"read'). If '-jagged' option is used the grid dimension is effectively multiplied by the specified number.\n"
		"Default: 16 (if neither '-size' nor '-eq_rad' are specified) or defined by\n"
		"         '-size' or '-eq_rad', '-lambda', and '-dpl'.",UNDEF,NULL},
	{PAR(h),"[<opt> [<subopt>]]","Shows help. If used without arguments, ADDA shows a list of all available command "
		"line options. If first argument is specified, help on specific command line option <opt> is shown (only the "
		"name of the option should be given without preceding dash). For some options (e.g. '-beam' or '-shape') "
		"specific help on a particular suboption <subopt> may be shown.\n"
		"Example: shape coated",UNDEF,NULL},
	{PAR(init_field),"{auto|inc|read <filenameY> [<filenameX>]|wkb|zero}",
		"Sets prescription to calculate initial (starting) field for the iterative solver.\n"
		"'auto' - automatically choose from 'zero' and 'inc' based on the lower residual value.\n"
		"'inc' - derived from the incident field,\n"
		"'read' - defined by separate files, which names are given as arguments. Normally two files are required for "
		"Y- and X-polarizations respectively, but a single filename is sufficient if only Y-polarization is used (e.g. "
		"due to symmetry). Initial field should be specified in a particle reference frame in the same format as used "
		"by '-store_int_field',\n"
		"'wkb' - from Wentzel-Kramers-Brillouin approximation,\n"
#ifdef SPARSE
		"!!! 'wkb' is not operational in sparse mode\n"
#endif
		"'zero' is a zero vector,\n"
		"Default: auto",UNDEF,NULL},
	{PAR(int),"{fcd|fcd_st|igt [<lim> [<prec>]]|igt_so|nloc <Rp>|nloc_av <Rp>|poi|so}",
		"Sets prescription to calculate the interaction term.\n"
		"'fcd' - Filtered Coupled Dipoles - requires dpl to be larger than 2.\n"
		"'fcd_st' - static (long-wavelength limit) version of FCD.\n"
		"'igt' - Integration of Green's Tensor. Its parameters are: <lim> - maximum distance (in dipole sizes), for "
		"which integration is used, (default: infinity); <prec> - minus decimal logarithm of relative error of the "
		"integration, i.e. epsilon=10^(-<prec>) (default: same as argument of '-eps' command line option).\n"
#ifdef NO_FORTRAN
		"!!! 'igt' relies on Fortran sources that were disabled at compile time.\n"
#endif
		"'igt_so' - approximate evaluation of IGT using second order of kd approximation.\n"
		"'nloc' - non-local interaction of two Gaussian dipole densities (based on point value of Gh), <Rp> is the "
		"width of the latter in um (must be non-negative).\n"
		"'nloc_av' - same as 'nloc' but based on averaging over the cube volume.\n"
		"'poi' - (the simplest) interaction between point dipoles.\n"
		"'so' - under development and incompatible with '-anisotr'.\n"
#ifdef SPARSE
		"!!! All options except 'poi' incur a significant slowing down in sparse mode.\n"
#endif
		"Default: poi",UNDEF,NULL},
		/* TO ADD NEW INTERACTION FORMULATION
		 * Modify string constants after 'PAR(int)': add new argument (possibly with additional sub-arguments) to list
		 * {...} and its description to the next string. If the new interaction formulation requires unusually large
		 * computational time, add a special note for sparse mode (after '#ifdef SPARSE').
		 */
	{PAR(int_surf),"{img|som}",
		"Sets prescription to calculate the reflection term.\n"
		"'img' - approximation based on a single image dipole (fast but inaccurate).\n"
		"'som' - direct evaluation of Sommerfeld integrals.\n"
	#ifdef SPARSE
		"!!! In sparse mode 'som' is expected to be very slow.\n"
	#endif
		"Default: som (but 'img' if surface is perfectly reflecting)",1,NULL},
		/* TO ADD NEW REFLECTION FORMULATION
		 * Modify string constants after 'PAR(int_surf)': add new argument (possibly with additional sub-arguments) to
		 * list {...} and its description to the next string. If the new reflection formulation requires unusually
		 * large computational time, add a special note for sparse mode (after '#ifdef SPARSE').
		 * !!! If subarguments are added, second-to-last argument should be changed from 1 to UNDEF, and consistency
		 * test for number of arguments should be implemented in PARSE_FUNC(int_surf) below.
		 */
	{PAR(iter),"{bcgs2|bicg|bicgstab|cgnr|csym|qmr|qmr2}","Sets the iterative solver.\n"
		"Default: qmr",1,NULL},
		/* TO ADD NEW ITERATIVE SOLVER
		 * add the short name, used to define the new iterative solver in the command line, to the list "{...}" in the
		 * alphabetical order.
		 */
	{PAR(jagged),"<arg>","Sets a size of a big dipole in units of small dipoles, integer. It is used to improve the "
		"discretization of the particle without changing the shape.\n"
		"Default: 1",1,NULL},
	{PAR(lambda),"<arg>","Sets incident wavelength in um, float.\n"
		"Default: 2*pi",1,NULL},
	{PAR(m),"{<m1Re> <m1Im> [...]|<m1xxRe> <m1xxIm> <m1yyRe> <m1yyIm> <m1zzRe> <m1zzIm> [...]}","Sets refractive "
		"indices, float. Each pair of arguments specifies real and imaginary part of the refractive index of one of "
		"the domains. If '-anisotr' is specified, three refractive indices correspond to one domain (diagonal elements "
		"of refractive index tensor in particle reference frame). Maximum number of different refractive indices is "
		"defined at compilation time by the parameter MAX_NMAT in file const.h (by default, 15). None of the "
		"refractive indices can be equal to 1+0i.\n"
		"Default: 1.5 0",UNDEF,NULL},
	{PAR(maxiter),"<arg>","Sets the maximum number of iterations of the iterative solver, integer.\n"
		"Default: very large, not realistic value",1,NULL},
	{PAR(no_reduced_fft),"","Do not use symmetry of the interaction matrix to reduce the storage space for the "
		"Fourier-transformed matrix.",0,NULL},
	{PAR(no_vol_cor),"","Do not use 'dpl (volume) correction'. If this option is given, ADDA will try to match size of "
		"the dipole grid along x-axis to that of the particle, either given by '-size' or calculated analytically from "
		"'-eq_rad'. Otherwise (by default) ADDA will try to match the volumes, using either '-eq_rad' or the value "
		"calculated analytically from '-size'.",0,NULL},
	{PAR(ntheta),"<arg>","Sets the number of intervals, into which the range of scattering angles [0,180] (degrees) is "
		"equally divided, integer. This is used for scattering angles in yz-plane. If particle is not symmetric and "
		"orientation averaging is not used, the range is extended to 360 degrees (with the same length of elementary "
		"interval, i.e. number of intervals is doubled).\n"
		"Default: from 90 to 720 depending on the size of the computational grid.",1,NULL},
	{PAR(opt),"{speed|mem}",
		"Sets whether ADDA should optimize itself for maximum speed or for minimum memory usage.\n"
		"Default: speed",1,NULL},
	{PAR(orient),"{<alpha> <beta> <gamma>|avg [<filename>]}","Either sets an orientation of the particle by three "
		"Euler angles 'alpha','beta','gamma' (in degrees) or specifies that orientation averaging should be "
		"performed. <filename> sets a file with parameters for orientation averaging. Here zyz-notation (or "
		"y-convention) is used for Euler angles.\n"
		"Default orientation: 0 0 0\n"
		"Default <filename>: "FD_AVG_PARMS,UNDEF,NULL},
	{PAR(phi_integr),"<arg>","Turns on and specifies the type of Mueller matrix integration over azimuthal angle "
		"'phi'. <arg> is an integer from 1 to 31, each bit of which, from lowest to highest, indicates whether the "
		"integration should be performed with multipliers 1, cos(2*phi), sin(2*phi), cos(4*phi), and sin(4*phi) "
		"respectively.\n"
		"Examples: 1 (one integration with no multipliers),\n"
		"          6 (two integration with cos(2*phi) and sin(2*phi) multipliers).",1,NULL},
	{PAR(pol),"{cldr|cm|dgf|fcd|igt_so|lak|ldr [avgpol]|nloc <Rp>|nloc_av <Rp>|rrc|so}",
		"Sets prescription to calculate the dipole polarizability.\n"
		"'cldr' - Corrected LDR (see below), incompatible with '-anisotr'.\n"
		"'cm' - (the simplest) Clausius-Mossotti.\n"
		"'dgf' - Digitized Green's Function (second order approximation to LAK).\n"
		"'fcd' - Filtered Coupled Dipoles (requires dpl to be larger than 2).\n"
		"'igt_so' - Integration of Green's Tensor over a cube (second order approximation).\n"
		"'lak' - (by Lakhtakia) exact integration of Green's Tensor over a sphere.\n"
		"'ldr' - Lattice Dispersion Relation, optional flag 'avgpol' can be added to average polarizability over "
		"incident polarizations.\n"
		"'nloc' - non-local (Gaussian dipole density, based on lattice sums), <Rp> is the width of the latter in um "
		"(must be non-negative).\n"
		"'nloc_av' - same as 'nloc' but based on averaging of Gh over the dipole volume.\n"
		"'rrc' - Radiative Reaction Correction (added to CM).\n"
		"'so' - under development and incompatible with '-anisotr'.\n"
		"Default: ldr (without averaging).",UNDEF,NULL},
		/* TO ADD NEW POLARIZABILITY FORMULATION
		 * Modify string constants after 'PAR(pol)': add new argument (possibly with additional sub-arguments) to list
		 * {...} and its description to the next string.
		 */
	{PAR(prognosis),"","Do not actually perform simulation (not even memory allocation) but only estimate the required "
		"RAM. Implies '-test'.",0,NULL},
	{PAR(prop),"<x> <y> <z>","Sets propagation direction of incident radiation, float. Normalization (to the unity "
		"vector) is performed automatically. For point-dipole incident beam this determines its direction.\n"
		"Default: 0 0 1",3,NULL},
	{PAR(recalc_resid),"","Recalculate residual at the end of iterative solver.",0,NULL},
#ifndef SPARSE
	{PAR(save_geom),"[<filename>]","Save dipole configuration to a file <filename> (a path relative to the output "
		"directory). Can be used with '-prognosis'.\n"
		"Default: <type>.geom \n"
		"(<type> is a first argument to the '-shape' option; '_gran' is added if '-granul' option is used; file "
		"extension can differ depending on argument of '-sg_format' option).",
		UNDEF,NULL},
#endif // !SPARSE
	{PAR(scat),"{dr|fin|igt_so|so}","Sets prescription to calculate scattering quantities.\n"
		"'dr' - (by Draine) standard formulation for point dipoles\n"
		"'fin' - slightly different one, based on a radiative correction for a finite dipole.\n"
		"'igt_so' - second order in kd approximation to Integration of Green's Tensor.\n"
		"'so' - under development and incompatible with '-anisotr'.\n"
		"Default: dr",1,NULL},
	{PAR(scat_grid_inp),"<filename>","Specifies a file with parameters of the grid of scattering angles for "
		"calculating Mueller matrix (possibly integrated over 'phi').\n"
		"Default: "FD_SCAT_PARMS,1,NULL},
	{PAR(scat_matr),"{muel|ampl|both|none}","Specifies which scattering matrices (from Mueller and amplitude) should "
		"be saved to file. Amplitude matrix is never integrated (in combination with '-orient avg' or '-phi_integr').\n"
		"Default: muel",1,NULL},
	{PAR(scat_plane),"","Explicitly enables calculation of the scattering in the plane through e_z (in laboratory "
		"reference frame) and propagation direction of incident wave. For default incidence this is the xz-plane. It "
		"can also be implicitly enabled by other options.",0,NULL},
#ifndef SPARSE
	{PAR(sg_format),"{text|text_ext|ddscat6|ddscat7}","Specifies format for saving geometry files. First two are ADDA "
		"default formats for single- and multi-domain particles respectively. 'text' is automatically changed to "
		"'text_ext' for multi-domain particles. Two DDSCAT formats correspond to its shape options 'FRMFIL' (version "
		"6) and 'FROM_FILE' (version 7) and output of 'calltarget' utility.\n"
		"Default: text",1,NULL},
#endif // !SPARSE
		/* TO ADD NEW FORMAT OF SHAPE FILE
		 * Modify string constants after 'PAR(sg_format)': add new argument to list {...} and add its description to
		 * the next string.
		 */
	{PAR(shape),"<type> [<args>]","Sets shape of the particle, either predefined or 'read' from file. All parameters "
		"of predefined shapes are floats except for filenames.\n"
		"Default: sphere",UNDEF,shape_opt},
	{PAR(size),"<arg>","Sets the size of the computational grid along the x-axis in um, float. If default wavelength "
		"is used, this option specifies the 'size parameter' of the computational grid. Can not be used together with "
		"'-eq_rad'. Size is defined by some shapes themselves, then this option can be used to override the internal "
		"specification and scale the shape.\n"
		"Default: determined by the value of '-eq_rad' or by '-grid', '-dpl', and '-lambda'.",1,NULL},
	{PAR(store_beam),"","Save incident beam to a file",0,NULL},
	{PAR(store_dip_pol),"","Save dipole polarizations to a file",0,NULL},
	{PAR(store_force),"","Calculate the radiation force on each dipole. Implies '-Cpr'",0,NULL},
#ifndef SPARSE
	{PAR(store_grans),"","Save granule coordinates (placed by '-granul' option) to a file",0,NULL},
#endif
	{PAR(store_int_field),"","Save internal fields to a file",0,NULL},
	{PAR(store_scat_grid),"","Calculate Mueller matrix for a grid of scattering angles and save it to a file.",0,NULL},
	{PAR(surf),"<h> {<mre> <mim>|inf}","Specifies that scatterer is located above the plane surface, parallel to the "
		"xy-plane. <h> specifies the height of particle center above the surface (along the z-axis, in um). Particle "
		"must be entirely above the substrate. Following argument(s) specify the refractive index of the substrate "
		"(below the surface), assuming that the vacuum is above the surface. It is done either by two values (real and "
		"imaginary parts of the complex value) or as effectively infinite 'inf' which corresponds to perfectly"
		"reflective surface. The latter implies certain simplifications during calculations.",UNDEF,NULL},
	{PAR(sym),"{auto|no|enf}","Automatically determine particle symmetries ('auto'), do not take them into account "
		"('no'), or enforce them ('enf').\n"
		"Default: auto",1,NULL},
	{PAR(test),"","Begin name of the output directory with 'test' instead of 'run'",0,NULL},
	{PAR(V),"","Show ADDA version, compiler used to build this executable, build options, and copyright information",
		0,NULL},
	{PAR(vec),"","Calculate the not-normalized asymmetry vector",0,NULL},
	{PAR(yz),"","Explicitly enables calculation of the scattering in the yz-plane (in incident-wave reference frame). "
		"It can also be implicitly enabled by other options.",0,NULL}
	/* TO ADD NEW COMMAND LINE OPTION
	 * add a row to this list, initializing an option. It should contain: PAR(option_name), usage string, help string,
	 * number of arguments, pointer to suboption (if exist, NULL otherwise). Usage and help string are shown either
	 * when -h option is invoked or error is found in input command line. Usage string is one line giving a list of
	 * possible arguments or argument combinations. Do not include option name in it, use <arg_name> to denote argument,
	 * [...] for optional arguments, and {...|...|...} for multiple options of an argument. Help string should contain
	 * general description of the option and its arguments, and provide default value for the arguments (if applicable).
	 * UNDEF can be used instead of number of arguments to avoid automatic checking, e.g. when command line option can
	 * accept variable number of arguments. Then check it explicitly in function PARSE_FUNC(option_name) below. A
	 * separate suboption structure is recommended only for large options such as -beam and -shape, and it should be
	 * defined above.
	 */
};

// auxiliary functions
//======================================================================================================================

static const char *OptionName(void)
// produces full option name for error messages
{
	static char buf[MAX_LINE];

	if (opt.l2==UNDEF) return options[opt.l1].name;
	else {
		sprintf(buf,"%s %s",options[opt.l1].name,options[opt.l1].sub[opt.l2].name);
		return buf;
	}
}

//======================================================================================================================

void PrintErrorHelp(const char * restrict fmt, ... )
/* print anything to stderr (on root processor), then help on the arguments used, and stop; assumes that all processors
 * call it; has line wrapping. It is designed to be relatively safe (using snprintf and vsnprintf), but do not produce
 * any additional errors in case of buffer overflows, etc. (not to distract from the main error itself). However, for
 * uncontrolled (e.g. user-input) arguments, it is recommended to use PrintErrorHelpSafe instead.
 */
{
	va_list args;
	const char * restrict optname,* restrict use;
	char msg[MAX_MESSAGE]="ERROR: ";
	int shift,tmp;

	if (IFROOT) {
		// produce error message
		va_start(args,fmt);
		shift=strlen(msg);
		VSNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE,fmt,args);
		va_end(args);
		// add help message
		if (opt.l1==UNDEF) { // no option is found
			SNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE,"\n"
				"Usage: %s %s\n"
				"Type '%s -h' for help\n",exename,exeusage,exename);
		}
		else { // at least option is found
			if (opt.l2==UNDEF) use=options[opt.l1].usage;
			else use=options[opt.l1].sub[opt.l2].usage;
			optname=OptionName();
			SNPRINTF_SHIFT_ROBUST(shift,tmp,msg,MAX_MESSAGE,"\n"
				"Usage: -%s %s\n"
				"Type '%s -h %s' for details\n",optname,use,exename,optname);
		}
		WrapLines(msg);
		fprintf(stderr,"%s",msg);
		fflush(stderr);
	}
	// wait for root to generate an error message
	Synchronize();
	Stop(EXIT_FAILURE);
}

//======================================================================================================================

static void ATT_NORETURN ATT_PRINTF(1,2) PrintErrorHelpSafe(const char * restrict fmt, ... )
/* print anything to stderr (on root processor), then help on the arguments used, and stop; assumes that all processors
 * call it; same as PrintErrorHelp but uses no internal buffers to be safe for any input parameters, which may come
 * from a command line, at a cost of lacking line wrapping.
 */
{
	va_list args;
	const char * restrict optname,* restrict use;

	if (IFROOT) {
		// produce error message
		va_start(args,fmt);
		fprintf(stderr,"ERROR: ");
		vfprintf(stderr,fmt,args);
		fprintf(stderr,"\n");
		va_end(args);
		// add help message
		if (opt.l1==UNDEF) // no option is found
			fprintf(stderr,"Usage: %s %s\n"
			               "Type '%s -h' for help\n",exename,exeusage,exename);
		else { // at least option is found
			if (opt.l2==UNDEF) use=options[opt.l1].usage;
			else use=options[opt.l1].sub[opt.l2].usage;
			optname=OptionName();
			fprintf(stderr,"Usage: -%s %s\n"
			               "Type '%s -h %s' for details\n",optname,use,exename,optname);
		}
		fflush(stderr);
	}
	// wait for root to generate an error message
	Synchronize();
	Stop(EXIT_FAILURE);
}

//======================================================================================================================

static void NargError(const int Narg,const char *expec)
// Print error of illegal number of arguments to an option (suboption) and display correct usage information
{
	PrintErrorHelp("Illegal number of arguments (%d) to '-%s' option (%s expected)",Narg,OptionName(),expec);
}


//======================================================================================================================

static void NargErrorSub(const int Narg,const char *option,const char *expec)
/* Print error of illegal number of arguments to an option (usually with suboption) and display correct usage
 * information. Narg is decremented before printing.
 */
{
	PrintErrorHelp("Illegal number of arguments (%d) to '-%s' option (%s expected)",Narg-1,option,expec);
}

//======================================================================================================================

static void TestExtraNarg(const int Narg, const bool notAllowed,const char *subopt)
// check if Narg is larger than 1 for suboption, which doesn't allow that
{
	if (Narg>1 && notAllowed) PrintErrorHelp("Additional arguments are not allowed for '-%s %s'",OptionName(),subopt);
}

//======================================================================================================================

static void TestNarg(const int Narg,const int need)
// check if Narg given to an option (or suboption) is correct; interface to NargError
{
	if (need>=0) { // usual case
		if (Narg!=need) {
			char buf[MAX_WORD];
			snprintf(buf,MAX_WORD,"%d",need);
			NargError(Narg,buf);
		}
	} // otherwise special cases are considered, encoded by negative values
	else switch (need) {
		case UNDEF: break; // do nothing
		case FNAME_ARG: if (Narg!=1) NargError(Narg,"1"); break;
		case FNAME_ARG_2: if (Narg!=2) NargError(Narg,"2"); break;
		case FNAME_ARG_1_2: if (Narg!=1 && Narg!=2) NargError(Narg,"1 or 2"); break;
		default: // rigorous test that every possible special case is taken care of
			PrintError("Critical error in TestNarg function (unknown argument 'need'=%d). Probably this comes from "
				"value of 'narg' in one of static arrays, describing command line options.",need);
			// no break
	}
}

//======================================================================================================================

static void ATT_NORETURN NotSupported(const char * restrict type,const char * restrict given)
// print error message that "type 'given' is not supported"; type should start with a capital letter
{
	PrintErrorHelpSafe("%s '%s' is not supported",type,given);
}

//======================================================================================================================

static const char *ScanStrError(const char * restrict str,const unsigned int size)
// check if string fits in buffer of size 'size', otherwise produces error message; returns the passed str (redirects)
{
	if (strlen(str)>=size) PrintErrorHelp("Too long argument to '-%s' option (only %ud chars allowed). If you really "
		"need it you may increase MAX_DIRNAME in const.h and recompile",OptionName(),size-1);
	return str;
}

//======================================================================================================================

static void ScanDoubleError(const char * restrict str,double *res)
// scanf an option argument and checks for errors
{
	if (sscanf(str,"%lf",res)!=1)
		PrintErrorHelpSafe("Non-numeric argument (%s) is given to the option '-%s'",str,OptionName());
}

//======================================================================================================================

static void ScanIntError(const char * restrict str,int *res)
// scanf an option argument and checks for errors
{
	double tmp;

	if (sscanf(str,"%lf",&tmp)!=1)
		PrintErrorHelpSafe("Non-numeric argument (%s) is given to the option '-%s'",str,OptionName());
	// The following can be rewritten by call to ConvertToInteger(), but it is complicated due to presence of OptionName
	if (tmp != floor(tmp))
		PrintErrorHelpSafe("Argument value (%s) of the option '-%s' is not an integer",str,OptionName());
	if (tmp <INT_MIN || tmp>INT_MAX)
		PrintErrorHelpSafe("Argument value (%s) of the option '-%s' is out of integer bounds",str,OptionName());
	*res=(int)tmp;
}

//======================================================================================================================

static bool ScanFnamesError(const int Narg,const int need,char **argv,const char **fname1,const char **fname2)
/* If 'need' corresponds to one of FNAME_ARG, scan Narg<=2 filenames from argv into fname1 and fname2. All consistency
 * checks are left to the caller (in particular, whether Narg corresponds to need). argv should be shifted to contain
 * only filenames. fname2 can be NULL, but it will produce an error in combination with Narg=2)
 *
 * Returns whether filenames has been scanned.
 */
{
	bool res=false;

	if (IS_FNAME_ARG(need)) {
		*fname1=ScanStrError(argv[0],MAX_FNAME);
		if (Narg==2) {
			if (fname2==NULL) // consistency check e.g. for reading shape filename
				PrintError("Failed to store the second filename in function ScanFnamesError");
			else *fname2=ScanStrError(argv[1],MAX_FNAME);
		}
		res=true;
	}
	return res;
}

//======================================================================================================================

static bool IsOption(const char * restrict str)
/* checks if string is an option. First should be '-' and then letter (any case); it enables use of negative numbers
 * as sub-parameters
 */
{
	// conversion to int is needed to remove warnings caused by the fact that str[1] is _signed_ char
	return (str[0]=='-' && isalpha((int)(str[1])));
}
//======================================================================================================================

static int TimeField(const char c)
// analyze one time multiplier
{
	switch (c) {
		case 'd':
		case 'D': return 86400;
		case 'h':
		case 'H': return 3600;
		case 'm':
		case 'M': return 60;
		case 's':
		case 'S': return 1;
	}
	PrintErrorHelp("Illegal time format specifier (%c)",c);
}

//======================================================================================================================

static int ScanTime(const char * restrict str)
// scans time in seconds from a string "%d[d,D[%d]][h,H[%d]][m,M[%d]][s,S]
{
#define TIME_N_TYPES 4 // not so easy to change
	int tim,t[TIME_N_TYPES],n,i;
	char c[TIME_N_TYPES];

	for (i=0;i<TIME_N_TYPES;i++) c[i]=0;
	n=sscanf(str,"%d%c%d%c%d%c%d%c",t,c,t+1,c+1,t+2,c+2,t+3,c+3);
	if (n<1) PrintErrorHelpSafe("Wrong time format '%s'",str);
	tim=0;
	i=0;
	while (n>0) {
		tim+=t[i]*TimeField(c[i]);
		n-=2;
		i++;
	}
	return tim;
#undef TIME_N_TYPES
}

//======================================================================================================================

static void PrintTime(char * restrict s,const time_t *time_ptr)
{
	struct tm *t;

	t=gmtime(time_ptr);
	s[0]=0; // initialize string
	if (t->tm_yday>0) sprintf(s,"%dd ",t->tm_yday);
	if (t->tm_hour>0) sprintf(s+strlen(s),"%dh ",t->tm_hour);
	if (t->tm_min>0) sprintf(s+strlen(s),"%dm ",t->tm_min);
	if (t->tm_sec>0) sprintf(s+strlen(s),"%ds ",t->tm_sec);
}

//======================================================================================================================
// parsing functions definitions

PARSE_FUNC(alldir_inp)
{
	alldir_parms=ScanStrError(argv[1],MAX_FNAME);
}
PARSE_FUNC(anisotr)
{
	anisotropy = true;
	Ncomp=3;
}
PARSE_FUNC(asym)
{
	calc_asym = true;
	calc_vec = true;
	calc_Csca = true;
}
PARSE_FUNC(beam)
{
	int i,j,need;
	bool found;

	if (Narg==0) PrintErrorHelp("Suboption (argument) is missing for '-beam' option");
	Narg--;
	found=false;
	i=-1;
	while (beam_opt[++i].name!=NULL) if (strcmp(argv[1],beam_opt[i].name)==0) {
		// set suboption and beam type
		opt.l2=i;
		beamtype=(enum beam)beam_opt[i].type;
		beam_Npars=Narg;
		opt_beam=opt;
		need=beam_opt[i].narg;
		// check number of arguments
		switch (beamtype) {
			case B_LMINUS:
			case B_DAVIS3:
			case B_BARTON5: if (Narg!=1 && Narg!=4) NargError(Narg,"1 or 4"); break;
			default: TestNarg(Narg,need); break;
		}
		/* TO ADD NEW BEAM
		 * If the beam accepts variable number of arguments (UNDEF was used in beam definition above) add a check of
		 * number of received arguments as a new case above. Use NargError function similarly as done in existing tests.
		 */
		// either parse filename or parse all parameters as float; consistency is checked later
		if (!ScanFnamesError(Narg,need,argv+2,&beam_fnameY,&beam_fnameX))
			for (j=0;j<Narg;j++) ScanDoubleError(argv[j+2],beam_pars+j);
		// stop search
		found=true;
		break;
	}
	if(!found) NotSupported("Beam type",argv[1]);
}
PARSE_FUNC(chp_dir)
{
	chp_dir=ScanStrError(argv[1],MAX_DIRNAME);
}
PARSE_FUNC(chp_load)
{
	load_chpoint = true;
}
PARSE_FUNC(chp_type)
{
	if (strcmp(argv[1],"normal")==0) chp_type=CHP_NORMAL;
	else if (strcmp(argv[1],"regular")==0) chp_type=CHP_REGULAR;
	else if (strcmp(argv[1],"always")==0) chp_type=CHP_ALWAYS;
	else NotSupported("Checkpoint type",argv[1]);
}
PARSE_FUNC(chpoint)
{
	chp_time=ScanTime(argv[1]);
	if (chp_time<=0) {
		chp_time=UNDEF;
		if (chp_type==CHP_NONE) chp_type=CHP_ALWAYS;
	}
	else if (chp_type==CHP_NONE) chp_type=CHP_NORMAL;
}
PARSE_FUNC(Cpr)
{
	calc_mat_force = true;
}
PARSE_FUNC(Csca)
{
	calc_Csca = true;
}
PARSE_FUNC(dir)
{
	directory=ScanStrError(argv[1],MAX_DIRNAME);
}
PARSE_FUNC(dpl)
{
	ScanDoubleError(argv[1],&dpl);
	TestPositive(dpl,"dpl");
}
PARSE_FUNC(eps)
{
	double tmp;

	ScanDoubleError(argv[1],&tmp);
	TestPositive(tmp,"eps exponent");
	iter_eps=pow(10,-tmp);
}
PARSE_FUNC(eq_rad)
{
	ScanDoubleError(argv[1],&a_eq);
	TestPositive(a_eq,"dpl");
}
#ifdef OPENCL
PARSE_FUNC(gpu)
{
	ScanIntError(argv[1],&gpuInd);
	TestNonNegative_i(gpuInd,"GPU index");
}
#endif
#ifndef SPARSE
PARSE_FUNC(granul)
{
	if (Narg!=2 && Narg!=3) NargError(Narg,"2 or 3");
	ScanDoubleError(argv[1],&gr_vf);
	TestRangeII(gr_vf,"volume fraction",0,PI_OVER_SIX);
	ScanDoubleError(argv[2],&gr_d);
	TestPositive(gr_d,"diameter");
	if (Narg==3) {
		ScanIntError(argv[3],&gr_mat);
		TestPositive_i(gr_mat,"domain number");
	}
	else gr_mat=1;
	gr_mat--; // converted to usual indexing starting from 0
	sh_granul=true;
}
#endif
PARSE_FUNC(grid)
{
	if (Narg!=1 && Narg!=3) NargError(Narg,"1 or 3");
	ScanIntError(argv[1],&boxX); // boxes are further multiplied by jagged if needed
	TestRange_i(boxX,"gridX",1,BOX_MAX);
	if (Narg==3) {
		ScanIntError(argv[2],&boxY);
		TestRange_i(boxY,"gridY",1,BOX_MAX);
		ScanIntError(argv[3],&boxZ);
		TestRange_i(boxY,"gridY",1,BOX_MAX);
	}
}
PARSE_FUNC(h)
{
	int i,j;
	bool found;

	if (Narg>2) NargError(Narg,"not more than 2");
	// do all output on root processor
	if (IFROOT) {
		found=false;
		if (Narg>=1) {
			for(i=0;i<LENGTH(options);i++) if (strcmp(argv[1],options[i].name)==0) {
				if (Narg==2) {
					if (options[i].sub==NULL) printf("No specific help is available for suboptions of this option\n\n");
					else {
						j=-1;
						while (options[i].sub[++j].name!=NULL)
							if (strcmp(argv[2],options[i].sub[j].name)==0) {
								printf("  -%s %s %s\n%s\n",options[i].name,options[i].sub[j].name,
									options[i].sub[j].usage,WrapLinesCopy(options[i].sub[j].help));
								found=true;
								break;
							}
						if (!found) printf("Unknown suboption '%s'\n\n",argv[2]);
					}
				}
				if (!found) {
					printf("  -%s %s\n%s\n",options[i].name,options[i].usage,WrapLinesCopy(options[i].help));
					if (options[i].sub!=NULL) {
						printf("Available suboptions:\n");
						j=-1;
						while (options[i].sub[++j].name!=NULL)
							printf("  %s %s\n",options[i].sub[j].name,options[i].sub[j].usage);
						printf("Type '%s -h %s <subopt>' for details\n",exename,options[i].name);
#ifdef SPARSE
						if (strcmp(options[i].name,"shape")==0)
							printf("!!! Most of the shape options are disabled in sparse mode\n");
#endif
					}
				}
				found=true;
				break;
			}
			if (!found) printf("Unknown option '%s'\n\n",argv[1]);
		}
		if (!found) {
			printf("Usage: '%s %s'\n"
			       "Available options:\n",exename,exeusage);
			for (i=0;i<LENGTH(options);i++) printf("  -%s %s\n",options[i].name,options[i].usage);
			printf("Type '%s -h <opt>' for details\n",exename);
#ifdef SPARSE
			printf("!!! A number of options are disabled in sparse mode\n");
#endif
		}
	}
	// exit
	Stop(EXIT_SUCCESS);
}
PARSE_FUNC(init_field)
{
	bool noExtraArgs=true;

	if (Narg<1 || Narg>3) NargError(Narg,"from 1 to 3");
	if (strcmp(argv[1],"auto")==0) InitField=IF_AUTO;
	else if (strcmp(argv[1],"inc")==0) InitField=IF_INC;
	else if (strcmp(argv[1],"read")==0) {
		if (Narg!=2 && Narg!=3) NargErrorSub(Narg,"init_field read","1 or 2");
		ScanFnamesError(Narg-1,FNAME_ARG_1_2,argv+2,&infi_fnameY,&infi_fnameX);
		InitField=IF_READ;
		noExtraArgs=false;
	}
	else if (strcmp(argv[1],"wkb")==0) {
#ifdef SPARSE
		PrintErrorHelp("Initial field 'wkb' is not supported in sparse mode");
#else
		InitField=IF_WKB;
#endif
	}
	else if (strcmp(argv[1],"zero")==0) InitField=IF_ZERO;
	else NotSupported("Initial field prescription",argv[1]);
	TestExtraNarg(Narg,noExtraArgs,argv[1]);
}
PARSE_FUNC(int)
{
	double tmp;
	bool noExtraArgs=true;
	
	if (Narg<1 || Narg>3) NargError(Narg,"from 1 to 3");
	if (strcmp(argv[1],"fcd")==0) IntRelation=G_FCD;
	else if (strcmp(argv[1],"fcd_st")==0) IntRelation=G_FCD_ST;
	else if (strcmp(argv[1],"igt")==0) {
#ifdef NO_FORTRAN
		PrintErrorHelp("To use IGT compiling of Fortran sources must be enabled (comment line "
			           "'CFLAGS += -DNO_FORTRAN' in Makefile and recompile");
#endif
		IntRelation=G_IGT;
		if (Narg<1 || Narg>3) NargErrorSub(Narg,"int igt","from 0 to 2");
		if (Narg>=2) {
			ScanDoubleError(argv[2],&igt_lim);
			TestNonNegative(igt_lim,"distance limit for IGT");
			if (Narg==3) {
				ScanDoubleError(argv[3],&tmp);
				TestPositive(tmp,"IGT precision");
				igt_eps=pow(10,-tmp);
			}
		}
		noExtraArgs=false;
	}
	else if (strcmp(argv[1],"igt_so")==0) IntRelation=G_IGT_SO;
	else if (strcmp(argv[1],"nloc")==0) {
		IntRelation=G_NLOC;
		if (Narg!=2) NargErrorSub(Narg,"int nloc","1");
		ScanDoubleError(argv[2],&nloc_Rp);
		TestNonNegative(nloc_Rp,"Gaussian width");
		noExtraArgs=false;
	}
	else if (strcmp(argv[1],"nloc_av")==0) {
		IntRelation=G_NLOC_AV;
		if (Narg!=2) NargErrorSub(Narg,"int nloc_av","1");
		ScanDoubleError(argv[2],&nloc_Rp);
		TestNonNegative(nloc_Rp,"Gaussian width");
		noExtraArgs=false;
	}
	else if (strcmp(argv[1],"poi")==0) IntRelation=G_POINT_DIP;
	else if (strcmp(argv[1],"so")==0) IntRelation=G_SO;
	/* TO ADD NEW INTERACTION FORMULATION
	 * add the line to else-if sequence above in the alphabetical order, analogous to the ones already present. The
	 * variable parts of the line are its name used in command line and its descriptor, defined in const.h. If
	 * subarguments are used, test their quantity explicitly, process (scan) subarguments, and set noExtraArgs to false.
	 * See "igt" for example. You may also need to change the test for Narg in the beginning of this function.
	 */
	else NotSupported("Interaction term prescription",argv[1]);
	TestExtraNarg(Narg,noExtraArgs,argv[1]);
}
PARSE_FUNC(int_surf)
{
	if (strcmp(argv[1],"img")==0) ReflRelation=GR_IMG;
	else if (strcmp(argv[1],"som")==0) ReflRelation=GR_SOM;
	else NotSupported("Interaction term prescription",argv[1]);
	int_surf_used=true;
	/* TO ADD NEW REFLECTION FORMULATION
	 * add the line to else-if sequence above in the alphabetical order, analogous to the ones already present. The
	 * variable parts of the line are its name used in command line and its descriptor, defined in const.h.
	 * !!! If subarguments need to be used a test for them should be implemented throughout the function, similar to
	 * PARSE_FUNC(int) above
	 */
}
PARSE_FUNC(iter)
{
	if (strcmp(argv[1],"bcgs2")==0) IterMethod=IT_BCGS2;
	else if (strcmp(argv[1],"bicg")==0) IterMethod=IT_BICG_CS;
	else if (strcmp(argv[1],"bicgstab")==0) IterMethod=IT_BICGSTAB;
	else if (strcmp(argv[1],"cgnr")==0) IterMethod=IT_CGNR;
	else if (strcmp(argv[1],"csym")==0) IterMethod=IT_CSYM;
	else if (strcmp(argv[1],"qmr")==0) IterMethod=IT_QMR_CS;
	else if (strcmp(argv[1],"qmr2")==0) IterMethod=IT_QMR_CS_2;
	/* TO ADD NEW ITERATIVE SOLVER
	 * add the line to else-if sequence above in the alphabetical order, analogous to the ones already present. The
	 * variable parts of the line are its name used in command line and its descriptor, defined in const.h
	 */
	else NotSupported("Iterative method",argv[1]);
}
PARSE_FUNC(jagged)
{
	ScanIntError(argv[1],&jagged);
	TestRange_i(jagged,"jagged",1,BOX_MAX);
}
PARSE_FUNC(lambda)
{
	ScanDoubleError(argv[1],&lambda);
	TestPositive(lambda,"wavelength");
}
PARSE_FUNC(m)
{
	int i;
	double mre,mim;

	if (!IS_EVEN(Narg) || Narg==0) NargError(Narg,"even");
	Nmat=Nmat_given=Narg/2;
	if (Nmat>MAX_NMAT) PrintErrorHelp("Too many materials (%d), maximum %d are supported. You may increase parameter "
		"MAX_NMAT in const.h and recompile.",Nmat,MAX_NMAT);
	for (i=0;i<Nmat;i++) {
		ScanDoubleError(argv[2*i+1],&mre);
		ScanDoubleError(argv[2*i+2],&mim);
		ref_index[i] = mre + I*mim;
		if (ref_index[i]==1) PrintErrorHelp("Given refractive index #%d is that of vacuum, which is not supported. "
			"Consider using, for instance, 1.0001 instead.",i+1);
	}
}
PARSE_FUNC(maxiter)
{
	ScanIntError(argv[1],&maxiter);
	TestPositive_i(maxiter,"maximum number of iterations");
}
PARSE_FUNC(no_reduced_fft)
{
	reduced_FFT=false;
}
PARSE_FUNC(no_vol_cor)
{
	volcor=false;
}
PARSE_FUNC(ntheta)
{
	ScanIntError(argv[1],&nTheta);
	TestPositive_i(nTheta,"number of theta intervals");
	nTheta++;
}
PARSE_FUNC(opt)
{
	if (strcmp(argv[1],"speed")==0) save_memory=false;
	else if (strcmp(argv[1],"mem")==0) save_memory=true;
	else NotSupported("Optimization method",argv[1]);
}
PARSE_FUNC(orient)
{
	if (Narg==0) NargError(Narg,"at least 1");
	if (strcmp(argv[1],"avg")==0) {
		if (Narg>2) NargErrorSub(Narg,"orient avg","0 or 1");
		orient_avg=true;
		if (Narg==2) avg_parms=ScanStrError(argv[2],MAX_FNAME);
	}
	else {
		if (Narg!=3) NargError(Narg,"3 or 'avg ...'");
		ScanDoubleError(argv[1],&alph_deg);
		ScanDoubleError(argv[2],&bet_deg);
		ScanDoubleError(argv[3],&gam_deg);
	}
	/* this will catch even trivial cases like '-orient 0 0 0', but we still want to point the user to existent
	 * incompatibilities. This may help to prevent misunderstanding on the user's side.
	 */
	orient_used=true;
}
PARSE_FUNC(phi_integr)
{
	phi_integr = true;
	ScanIntError(argv[1],&phi_int_type);
	TestRange_i(phi_int_type,"type of integration over phi",1,31);
}
PARSE_FUNC(pol)
{
	bool noExtraArgs=true;

	if (Narg!=1 && Narg!=2) NargError(Narg,"1 or 2");
	if (strcmp(argv[1],"cldr")==0) PolRelation=POL_CLDR;
	else if (strcmp(argv[1],"cm")==0) PolRelation=POL_CM;
	else if (strcmp(argv[1],"dgf")==0) PolRelation=POL_DGF;
	else if (strcmp(argv[1],"fcd")==0) PolRelation=POL_FCD;
	else if (strcmp(argv[1],"igt_so")==0) PolRelation=POL_IGT_SO;
	else if (strcmp(argv[1],"lak")==0) PolRelation=POL_LAK;
	else if (strcmp(argv[1],"ldr")==0) {
		PolRelation=POL_LDR;
		if (Narg!=1 && Narg!=2) NargErrorSub(Narg,"pol ldr","from 0 to 1");
		if (Narg==2) {
			if (strcmp(argv[2],"avgpol")==0) avg_inc_pol=true;
			else PrintErrorHelpSafe("Unknown argument '%s' to '-pol ldr' option",argv[2]);
		}
		noExtraArgs=false;
	}
	else if (strcmp(argv[1],"nloc")==0) {
		PolRelation=POL_NLOC;
		if (Narg!=2) NargErrorSub(Narg,"pol nloc","1");
		ScanDoubleError(argv[2],&polNlocRp);
		TestNonNegative(polNlocRp,"Gaussian width");
		noExtraArgs=false;
	}
	else if (strcmp(argv[1],"nloc_av")==0) {
		PolRelation=POL_NLOC_AV;
		if (Narg!=2) NargErrorSub(Narg,"pol nloc_av","1");
		ScanDoubleError(argv[2],&polNlocRp);
		TestNonNegative(polNlocRp,"Gaussian width");
		noExtraArgs=false;
	}
	else if (strcmp(argv[1],"rrc")==0) PolRelation=POL_RRC;
	else if (strcmp(argv[1],"so")==0) PolRelation=POL_SO;
	/* TO ADD NEW POLARIZABILITY FORMULATION
	 * add the line to else-if sequence above in the alphabetical order, analogous to the ones already present. The
	 * variable parts of the line are its name used in command line and its descriptor, defined in const.h. If
	 * subarguments are used, test their quantity explicitly, process (scan) subarguments, and set noExtraArgs to false.
	 * See "nloc" for example. You may also need to change the test for Narg in the beginning of this function.
	 */
	else NotSupported("Polarizability relation",argv[1]);
	TestExtraNarg(Narg,noExtraArgs,argv[1]);
}
PARSE_FUNC(prognosis)
{
	prognosis=true;
	run_name="test";
}
PARSE_FUNC(prop)
{
	double tmp;

	ScanDoubleError(argv[1],prop_0);
	ScanDoubleError(argv[2],prop_0+1);
	ScanDoubleError(argv[3],prop_0+2);
	tmp=DotProd(prop_0,prop_0);
	if (tmp==0) PrintErrorHelp("Given propagation vector is null");
	vMultScal(1/sqrt(tmp),prop_0,prop_0);
	prop_used=true;
}
PARSE_FUNC(recalc_resid)
{
	recalc_resid=true;
}
#ifndef SPARSE
PARSE_FUNC(save_geom)
{
	if (Narg>1) NargError(Narg,"0 or 1");
	save_geom=true;
	if (Narg==1) save_geom_fname=ScanStrError(argv[1],MAX_FNAME);
}
#endif
PARSE_FUNC(scat)
{
	if (strcmp(argv[1],"dr")==0) ScatRelation=SQ_DRAINE;
	else if (strcmp(argv[1],"fin")==0) ScatRelation=SQ_FINDIP;
	else if (strcmp(argv[1],"igt_so")==0) ScatRelation=SQ_IGT_SO;
	else if (strcmp(argv[1],"so")==0) ScatRelation=SQ_SO;
	else NotSupported("Scattering quantities relation",argv[1]);
}
PARSE_FUNC(scat_grid_inp)
{
	scat_grid_parms=ScanStrError(argv[1],MAX_FNAME);
}
PARSE_FUNC(scat_matr)
{
	if (strcmp(argv[1],"muel")==0) {
		store_mueller=true;
		store_ampl=false;
	}
	else if (strcmp(argv[1],"ampl")==0) {
		store_mueller=false;
		store_ampl=true;
	}
	else if (strcmp(argv[1],"both")==0) store_mueller=store_ampl=true;
	else if (strcmp(argv[1],"none")==0) store_mueller=store_ampl=false;
	else NotSupported("Scattering matrix specifier",argv[1]);
}
PARSE_FUNC(scat_plane)
{
	scat_plane_used = true;
	scat_plane = true;
}
#ifndef SPARSE
PARSE_FUNC(sg_format)
{
	if (strcmp(argv[1],"text")==0) sg_format=SF_TEXT;
	else if (strcmp(argv[1],"text_ext")==0) sg_format=SF_TEXT_EXT;
	else if (strcmp(argv[1],"ddscat6")==0) sg_format=SF_DDSCAT6;
	else if (strcmp(argv[1],"ddscat7")==0) sg_format=SF_DDSCAT7;
	/* TO ADD NEW FORMAT OF SHAPE FILE
	 * Based on argument of command line option '-sg_format' assign value to variable 'sg_format' (one of handles
	 * defined in const.h).
	 */
	else NotSupported("Geometry format",argv[1]);
}
#endif
PARSE_FUNC(shape)
{
	int i,j,need;
	bool found;

	if (Narg==0) PrintErrorHelp("Suboption (argument) is missing for '-shape' option");
	Narg--;
	found=false;
	i=-1;
	while (shape_opt[++i].name!=NULL) if (strcmp(argv[1],shape_opt[i].name)==0) {
		// set shape and shape option index
		shape=(enum sh)shape_opt[i].type;
		opt.l2=i;
		opt_sh=opt;
		sh_Npars=Narg;
		need=shape_opt[i].narg;
		// check number of arguments
		switch (shape) {
			case SH_COATED: if (Narg!=1 && Narg!=4) NargError(Narg,"1 or 4"); break;
			case SH_BOX: if (Narg!=0 && Narg!=2) NargError(Narg,"0 or 2"); break;
			default: TestNarg(Narg,need); break;
		}
		/* TO ADD NEW SHAPE
		 * If the shape accepts variable number of arguments (UNDEF was used in shape definition above) add a check of
		 * number of received arguments as a new case above. Use NargError function similarly as done in existing tests.
		 */
		// either parse filename or parse all parameters as float; consistency is checked later
		if (!ScanFnamesError(Narg,need,argv+2,&shape_fname,NULL))
			for (j=0;j<Narg;j++) ScanDoubleError(argv[j+2],sh_pars+j);
		// stop search
		found=true;
		break;
	}
	if (!found) NotSupported("Shape type",argv[1]);
	// set shape name; takes place only if shape name was matched above
	shapename=argv[1];
}
PARSE_FUNC(size)
{
	ScanDoubleError(argv[1],&sizeX);
	TestPositive(sizeX,"particle size");
}
PARSE_FUNC(store_beam)
{
	store_beam = true;
}
PARSE_FUNC(store_dip_pol)
{
	store_dip_pol=true;
}
PARSE_FUNC(store_force)
{
	store_force = true;
	calc_mat_force = true;
}
#ifndef SPARSE
PARSE_FUNC(store_grans)
{
	store_grans=true;
}
#endif
PARSE_FUNC(store_int_field)
{
	store_int_field=true;
}
PARSE_FUNC(store_scat_grid)
{
	store_scat_grid = true;
}
PARSE_FUNC(surf)
{
	double mre,mim;
	if (Narg!=2 && Narg!=3) NargError(Narg,"2 or 3");
	ScanDoubleError(argv[1],&hsub);
	TestPositive(hsub,"height above surface");
	if (strcmp(argv[2],"inf")==0) {
		if (Narg>2) PrintErrorHelp("Additional arguments are not allowed for '-surf ... inf'");
		msubInf=true;
	}
	else {
		if (Narg!=3) NargError(Narg,"total 3 or second argument 'inf'");
		ScanDoubleError(argv[2],&mre);
		ScanDoubleError(argv[3],&mim);
		msub = mre + I*mim;
	}
	surface = true;
	symZ=false;
}
PARSE_FUNC(sym)
{
	if (strcmp(argv[1],"auto")==0) sym_type=SYM_AUTO;
	else if (strcmp(argv[1],"no")==0) sym_type=SYM_NO;
	else if (strcmp(argv[1],"enf")==0) sym_type=SYM_ENF;
	else NotSupported("Symmetry option",argv[1]);
}
PARSE_FUNC(test)
{
	run_name="test";
}
PARSE_FUNC(V)
{
	char copyright[]="\n\nCopyright (C) 2006-2014 ADDA contributors\n"
		"This program is free software; you can redistribute it and/or modify it under the terms of the GNU General "
		"Public License as published by the Free Software Foundation; either version 3 of the License, or (at your "
		"option) any later version.\n\n"
		"This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the "
		"implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License "
		"for more details.\n\n"
		"You should have received a copy of the GNU General Public License along with this program. If not, see "
		"<http://www.gnu.org/licenses/>.\n";
	char ccver_str[MAX_LINE];
#if defined(__DECC)
	char cctype;
#elif defined(__BORLANDC__)
	int ccver;
#endif
	size_t num;
	int bits,len;

	if (IFROOT) {
		// compiler & version (works only for selected compilers)
		// Intel
#if defined(__ICC) || defined(__INTEL_COMPILER)
#	define COMPILER "Intel"
#	ifdef __INTEL_COMPILER
#		define CCVERSION __INTEL_COMPILER
#	else
#		define CCVERSION __ICC
#	endif
		sprintf(ccver_str,"%d.%d",CCVERSION/100,CCVERSION%100);
		// DEC (Compaq)
#elif defined(__DECC)
#	define COMPILER "DEC (Compaq)"
		cctype=(__DECC_VER/10000)%10;
		if (cctype==6) cctype='T';
		else if (cctype==8) cctype='S';
		else if (cctype==9) cctype='V';
		else cctype=' ';
		sprintf(ccver_str,"%c%d.%d-%d",cctype,__DECC_VER/10000000,(__DECC_VER/100000)%100,__DECC_VER%1000);
		// Borland
#elif defined(__BORLANDC__)
#	define COMPILER "Borland"
		sprintf(ccver_str,"%x",__BORLANDC__);
		sscanf(ccver_str,"%d",&ccver);
		sprintf(ccver_str,"%d.%d",ccver/100,ccver%100);
		// Microsoft
#elif defined(_MSC_VER)
#	define COMPILER "Microsoft"
		sprintf(ccver_str,"%d.%d",_MSC_VER/100,_MSC_VER%100);
		// GNU
#elif defined(__GNUC__)
#	define COMPILER "GNU"
		sprintf(ccver_str,"%d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
		// IBM
#elif defined(__xlc__)
#	define COMPILER "IBM"
		strncpy(ccver_str,__xlc__,MAX_LINE-1);
		// unknown compiler
#else
#	define COMPILER_UNKNOWN
#	define COMPILER "unknown"
#endif
		// print version, MPI standard, type and compiler information, bit-mode
#ifdef SVNREV // revision number is printed if available, but only here
		printf("ADDA v."ADDA_VERSION" (r"SVNREV")\n");
#else
		printf("ADDA v."ADDA_VERSION"\n");
#endif
#ifdef OPENCL
#	if defined(CL_VERSION_1_2)
#		define OCL_VERSION "1.2"
#	elif defined(CL_VERSION_1_1)
#		define OCL_VERSION "1.1"
#	elif defined(CL_VERSION_1_0)
#		define OCL_VERSION "1.0"
#	else // this should never happen, since minimum OpenCL version is checked in oclcore.h
#		error "OpenCL version not recognized"
#	endif
		printf("GPU-accelerated version conforming to OpenCL standard "OCL_VERSION"\n");
#	ifdef CLFFT_AMD
		printf("Linked to clFFT version %d.%d.%d\n",clfftVersionMajor,clfftVersionMinor,clfftVersionPatch);
#	endif
#	ifdef OCL_BLAS
		printf("Linked to clBLAS version %d.%d.%d\n",clblasVersionMajor,clblasVersionMinor,
			clblasVersionPatch);
#	endif
#elif defined(ADDA_MPI)
		// Version of MPI standard is specified
		printf("Parallel version conforming to MPI standard %d.%d\n",RUN_MPI_VER_REQ,RUN_MPI_SUBVER_REQ);
#	ifdef MPICH2
		printf("Linked to MPICH2 version "MPICH2_VERSION"\n");
#	elif defined(OPEN_MPI)
		printf("Linked to OpenMPI version %d.%d.%d\n",OMPI_MAJOR_VERSION,OMPI_MINOR_VERSION,OMPI_RELEASE_VERSION);
#	endif
		// Additional debug information about MPI implementation
#	ifndef SUPPORT_MPI_BOOL
		D("No full support for MPI_C_BOOL (emulated)");
#	endif
#	ifndef SUPPORT_MPI_COMPLEX
		D("No support for complex MPI datatypes (emulated)");
#	else
#		ifndef SUPPORT_MPI_COMPLEX
		D("Complex MPI datatypes are supported, but their use in reduce operations is not supported (emulated)");
#		endif
#	endif
#else
		printf("Sequential version\n");
#endif
		printf("Built with "COMPILER" compilers");
#ifndef COMPILER_UNKNOWN
		printf(" version %s",ccver_str);
#endif
		// determine number of bits in size_t; not the most efficient way, but should work robustly
		num=SIZE_MAX;
		bits=1;
		while(num>>=1) bits++;
		printf(" (%d-bit)\n",bits);
#ifdef __MINGW64_VERSION_STR
		printf("      using MinGW-64 environment version "__MINGW64_VERSION_STR"\n");
#elif defined(__MINGW32_VERSION)
		printf("      using MinGW-32 environment version %g\n",__MINGW32_VERSION);
#endif
		// extra build flags - can it be shortened by some clever defines?
		const char build_opts[]=
#ifdef DEBUGFULL
		"DEBUGFULL, "
#elif defined(DEBUG)
		"DEBUG, "
#endif
#ifdef FFT_TEMPERTON
		"FFT_TEMPERTON, "
#endif
#ifdef PRECISE_TIMING
		"PRECISE_TIMING, "
#endif
#ifdef NOT_USE_LOCK
		"NOT_USE_LOCK, "
#elif defined(ONLY_LOCKFILE)
		"ONLY_LOCKFILE, "
#endif
#ifdef NO_FORTRAN
		"NO_FORTRAN, "
#endif
#ifdef NO_CPP
		"NO_CPP, "
#endif
#ifdef OVERRIDE_STDC_TEST
		"OVERRIDE_STDC_TEST, "
#endif
#ifdef OCL_READ_SOURCE_RUNTIME
		"OCL_READ_SOURCE_RUNTIME, "
#endif
#ifdef CLFFT_APPLE
		"CLFFT_APPLE, "
#endif
#ifdef SPARSE
		"SPARSE, "
#endif
#ifdef USE_SSE3
		"USE_SSE3, "
#endif
#ifdef OCL_BLAS
		"OCL_BLAS, "
#endif
#ifdef NO_SVNREV
		"NO_SVNREV, "
#endif
		"";
		printf("Extra build options: ");
		len=strlen(build_opts);
		if (len==0) printf("none");
		else printf("%.*s",len-2,build_opts); // trailing comma and space are omitted
		// print copyright information
		WrapLines(copyright);
		printf("%s",copyright);
	}
	// exit
	Stop(EXIT_SUCCESS);
#undef COMPILER
#undef CCVERSION
#undef COMPILER_UNKNOWN
#undef BUILD_OPTS_SIZE
}
PARSE_FUNC(vec)
{
	calc_vec = true;
}
PARSE_FUNC(yz)
{
	yz_used = true;
	yzplane = true;
}
/* TO ADD NEW COMMAND LINE OPTION
 * add a function definition to this list. It should start with 'PARSE_FUNC(option_name)', and contain all necessary
 * logic (code) to implement (or register) command line option. Typically, it should check number of parameters (if
 * UNDEF was used instead of number of parameters in the option definition), scanf input parameters (if any) and check
 * them for consistency and set the values of some internal variables or flags. These variables should be defined in
 * the beginning of this file and initialized in InitVariables() below. It is recommended to use functions from param.h
 * and those defined above since they are designed to be overflow-safe for any input parameters and would automatically
 * produce informative output in case of error. If a new command line option contains suboptions, you may easily parse
 * them using the suboption structure (which should be defined in the beginning of this file). See PARSE_FUNC(beam) and
 * PARSE_FUNC(shape) above for examples. Potential conflicts or interactions between different command line options
 * should be processed in VariablesInterconnect() below. If any output in log files or stdout should be produced based
 * on the values of command line option, this should be implemented in PrintInfo() below.
 */
#undef PAR
#undef PARSE_FUNC
#undef PARSE_NAME

// end of parsing functions
//=============================================================

static FILEHANDLE CreateLockFile(const char * restrict fname ONLY_FOR_LOCK)
// create lock file; works only if USE_LOCK is enabled
{
#ifdef USE_LOCK
	FILEHANDLE fd;
	int i;

#	ifdef WINDOWS
	i=0;
	while ((fd=CreateFile(fname,GENERIC_WRITE,FILE_SHARE_WRITE,NULL,CREATE_NEW,FILE_ATTRIBUTE_NORMAL,NULL))
		    ==INVALID_HANDLE_VALUE) {
		Sleep(LOCK_WAIT*1000);
		if (i++ == MAX_LOCK_WAIT_CYCLES) LogError(ONE_POS,"Lock file %s permanently exists",fname);
	}
#	elif defined(POSIX)
#		ifdef LOCK_FOR_NFS
	struct flock lock;
#		endif
	// open file exclusively
	i=0;
	while ((fd=open(fname,O_WRONLY | O_CREAT | O_EXCL,0666))==-1) {
		sleep(LOCK_WAIT);
		if (i++ == MAX_LOCK_WAIT_CYCLES) LogError(ONE_POS,"Lock file %s permanently exists",fname);
	}
#		ifdef LOCK_FOR_NFS
	// specify lock
	lock.l_type=F_WRLCK;
	lock.l_whence=SEEK_SET;
	lock.l_start=0;
	lock.l_len=0;
	// obtain lock*/
	i=0;
	while (fcntl(fd,F_SETLK,&lock)==-1) {
		// if locked by another process wait and try again
		if (errno==EACCES || errno==EAGAIN) {
			sleep(LOCK_WAIT);
			if (i++ == MAX_LOCK_WAIT_CYCLES) LogError(ONE_POS,"Lock file %s permanently exists",fname);
		}
		else { // otherwise produce a message and continue
			if (errno==EOPNOTSUPP || errno==ENOLCK || errno==ENOSYS)
				LogWarning(EC_WARN,ONE_POS,"Advanced file locking is not supported by the file system");
			else LogWarning(EC_WARN,ONE_POS,"Unknown problem with file locking ('%s').",strerror(errno));
			break;
		}
	}
#		endif
#	endif
	// return file handle
	return fd;
#else
	return 0;
#endif
}

//======================================================================================================================

static void RemoveLockFile(FILEHANDLE fd ONLY_FOR_LOCK,const char * restrict fname ONLY_FOR_LOCK)
// closes and remove lock file; works only if USE_LOCK is enabled
{
#ifdef USE_LOCK
#	ifdef WINDOWS
	// close file
	CloseHandle(fd);
#	elif defined(POSIX)
	// close file; all locks are automatically released
	close(fd);
#	endif
	// remove lock file
	RemoveErr(fname,ONE_POS);
#endif
}

//======================================================================================================================

static void UpdateSymVec(const double a[static 3])
// tests whether vector a satisfies a number of symmetries (with round-off errors) and cancels the failing symmetries
{
	if (fabs(a[0])>ROUND_ERR) symX=false;
	if (fabs(a[1])>ROUND_ERR) symY=false;
	if (fabs(a[2])>ROUND_ERR) symZ=false;
	if (!vAlongZ(a)) symR=false;
}

//======================================================================================================================

void InitVariables(void)
// some defaults are specified also in const.h
{
	prop_used=false;
	orient_used=false;
	directory="";
	lambda=TWO_PI;
	// initialize ref_index of scatterer
	Nmat=Nmat_given=1;
	ref_index[0]=1.5;
	// initialize to null to determine further whether it is initialized
	logfile=NULL;
	boxX=boxY=boxZ=UNDEF;
	sizeX=UNDEF;
	a_eq=UNDEF;
	dpl=UNDEF;
	run_name="run";
	nTheta=UNDEF;
	iter_eps=1E-5;
	shape=SH_SPHERE;
	shapename="sphere";
	store_int_field=false;
	store_dip_pol=false;
	PolRelation=POL_LDR;
	avg_inc_pol=false;
	ScatRelation=SQ_DRAINE;
	IntRelation=G_POINT_DIP;
	IterMethod=IT_QMR_CS;
	sym_type=SYM_AUTO;
	prognosis=false;
	maxiter=UNDEF;
	jagged=1;
	beamtype=B_PLANE;
	alldir_parms=FD_ALLDIR_PARMS;
	avg_parms=FD_AVG_PARMS;
	scat_grid_parms=FD_SCAT_PARMS;
	chp_dir=FD_CHP_DIR;
	chp_time=UNDEF;
	chp_type=CHP_NONE;
	orient_avg=false;
	alph_deg=bet_deg=gam_deg=0.0;
	volcor=true;
	reduced_FFT=true;
	save_geom=false;
	save_geom_fname="";
	yzplane=false;
	yz_used=false;
	scat_plane=false;
	scat_plane_used=false;
	all_dir=false;
	scat_grid=false;
	phi_integr=false;
	store_scat_grid=false;
	calc_Cext=true;
	calc_Cabs=true;
	calc_Csca=false;
	calc_vec=false;
	calc_asym=false;
	calc_mat_force=false;
	store_force=false;
	store_mueller=true;
	store_ampl=false;
	store_grans=false;
	load_chpoint=false;
	sh_granul=false;
	symX=symY=symZ=symR=true;
	anisotropy=false;
	save_memory=false;
	sg_format=SF_TEXT;
	memory=0;
	memPeak=0;
	Ncomp=1;
	igt_lim=UNDEF;
	igt_eps=UNDEF;
	InitField=IF_AUTO;
	recalc_resid=false;
	surface=false;
	msubInf=false;
	int_surf_used=false;
	// sometimes the following two are left uninitialized
	beam_fnameX=NULL;
	infi_fnameX=NULL;
#ifdef OPENCL
	gpuInd=0;
#endif
	/* TO ADD NEW COMMAND LINE OPTION
	 * If you use some new variables, flags, etc. you should specify their default values here. This value will be used
	 * if new option is not specified in the command line.
	 */
}

//======================================================================================================================

void ParseParameters(const int argc,char **argv)
// parses input parameters
{
	int i,j,Narg,tmp;
	bool found;
	char *p1,*p2;

	// try to determine terminal width
	if ((p1=getenv("COLUMNS"))!=NULL && sscanf(p1,"%d",&tmp)==1 && tmp>=MIN_TERM_WIDTH) term_width=tmp;
	// get name of executable; remove all path overhead
	if ((p1=strrchr(argv[0],'\\'))==NULL) p1=argv[0];
	if ((p2=strrchr(argv[0],'/'))==NULL) p2=argv[0];
	exename=MAX(p1,p2)+1;
	// initialize option
	opt.l1=UNDEF;
	// check first argument
	if (argc>1 && !IsOption(argv[1])) PrintErrorHelpSafe("Illegal format of first argument '%s'",argv[1]);
	// read command line
	for (i=1;i<argc;i++) {
		// get number of arguments
		Narg=0;
		while ((i+(++Narg))<argc && !IsOption(argv[i+Narg]));
		Narg--;

		argv[i]++; // shift to remove "-" in the beginning of the string
		found=false;
		opt.l1=opt.l2=UNDEF;
		for (j=0;j<LENGTH(options);j++) if (strcmp(argv[i],options[j].name)==0) {
			opt.l1=j;
			// check consistency, if enabled for this parameter
			TestNarg(Narg,options[j].narg);
			// parse this parameter
			(*options[j].func)(Narg,argv+i);
			// check duplicate options; it is safe since at this point argv[i] is known to be normal
			if (options[j].used) PrintError("Option '-%s' is used more than once",argv[i]);
			else options[j].used=true;
			// stop search
			found=true;
			break;
		}
		if(!found) PrintErrorHelpSafe("Unknown option '-%s'",argv[i]);
		argv[i]--; // shift back
		i+=Narg;
	}
}

//======================================================================================================================

void VariablesInterconnect(void)
// finish parameters initialization based on their interconnections
{
	double temp;

	// initialize WaveNum ASAP
	WaveNum = TWO_PI/lambda;
	// set default incident direction, which is +z for all configurations
	if (!prop_used) {
		prop_0[0]=prop_0[1]=0;
		prop_0[2]=1;
	}
	// parameter interconnections
	if (IntRelation==G_SO) {
		reduced_FFT=false;
		// this limitation is due to assumption of reciprocity in DecayCross()
		if (beamtype==B_DIPOLE) PrintError("'-beam dipole' and '-int so' can not be used together");
	}
	/* TO ADD NEW INTERACTION FORMULATION
	 * If the new Green's tensor is non-symmetric (which is very unlikely) add it to the definition of reduced_FFT
	 */
	if (calc_Csca || calc_vec) all_dir = true;
	// by default, one of the scattering options is activated
	if (store_scat_grid || phi_integr) scat_grid = true;
	else if (!scat_plane_used && !yz_used) { // default values when no scattering plane is explicitly specified
		if (surface) scat_plane = true;
		else yzplane = true;
	}
	// if not initialized before, IGT precision is set to that of the iterative solver
	if (igt_eps==UNDEF) igt_eps=iter_eps;
	// parameter incompatibilities
	if (scat_plane && yzplane) PrintError("Currently '-scat_plane' and '-yz' cannot be used together.");
	if (orient_avg) {
		if (prop_used) PrintError("'-prop' and '-orient avg' can not be used together");
		if (store_int_field) PrintError("'-store_int_field' and '-orient avg' can not be used together");
		if (store_dip_pol) PrintError("'-store_dip_pol' and '-orient avg' can not be used together");
		if (store_beam) PrintError("'-store_beam' and '-orient avg' can not be used together");
		if (beamtype==B_READ) PrintError("'-beam read' and '-orient avg' can not be used together");
		if (scat_grid) PrintError("'-orient avg' can not be used with calculation of scattering for a grid of angles");
		// TODO: this limitation should be removed in the future
		if (all_dir) PrintError("Currently '-orient avg' can not be used with calculation of asym or Csca");
		if (calc_mat_force) PrintError("Currently '-orient avg' can not be used with calculation of Cpr");
		if (store_force) PrintError("'-store_force' and '-orient avg' can not be used together");
		if (!store_mueller && store_ampl) {
			store_ampl=false;
			LogWarning(EC_WARN,ONE_POS,"Amplitude matrix can not be averaged over orientations. So switching off "
				"calculation of Mueller matrix results in no calculation of scattering matrices at all.");
		}
		if (scat_plane && yzplane) PrintError("Only one scattering plane can be used during orientation averaging, so "
			"a combination of '-scat_plane' and '-yz' cannot be used together with '-orient avg'.");
	}
	if (phi_integr && !store_mueller) PrintError("Integration over phi can only be performed for Mueller matrix. "
		"Hence, '-phi_integr' is incompatible with '-scat_matr {ampl|none}'");
	if (!store_mueller && !store_ampl) {
		if (nTheta!=UNDEF) PrintError("'-scat_matr none' turns off calculation of angle-resolved quantities making "
			"'-ntheta' irrelevant. Hence, these two options are incompatible.");
		if (store_scat_grid) PrintError("'-scat_matr none' turns off calculation of angle-resolved quantities. Hence, "
			"it is incompatible with '-store_scat_grid'.");
		if (yz_used) PrintError("'-scat_matr none' turns off calculation of angle-resolved quantities. Hence, it is "
			"incompatible with '-yz'.");
		if (scat_plane_used) PrintError("'-scat_matr none' turns off calculation of angle-resolved quantities. Hence, "
			"it is incompatible with '-scat_plane'.");
		nTheta=0;
		yzplane=false;
		scat_plane=false;
		scat_grid=false;
	}
	if (anisotropy) {
		if (PolRelation==POL_CLDR) PrintError("'-anisotr' is incompatible with '-pol cldr'");
		if (PolRelation==POL_SO) PrintError("'-anisotr' is incompatible with '-pol so'");
		if (ScatRelation==SQ_SO) PrintError("'-anisotr' is incompatible with '-scat so'");
		if (IntRelation==G_SO) PrintError("'-anisotr' is incompatible with '-int so'");
		/* TO ADD NEW POLARIZABILITY FORMULATION
		 * If the new polarizability formulation is incompatible with anisotropic material, add an exception here
		 */
		if (Nmat%3!=0)
			PrintError("When '-anisotr' is used 6 numbers (3 complex values) should be given per each domain");
		else Nmat=Nmat/3;
	}
	if (chp_type!=CHP_NONE) {
		if (chp_time==UNDEF && chp_type!=CHP_ALWAYS) PrintError("You must specify time for this checkpoint type");
		// TODO: this limitation should be removed in the future
		if (orient_avg) PrintError("Currently checkpoint is incompatible with '-orient avg'");
	}
	if (sizeX!=UNDEF && a_eq!=UNDEF) PrintError("'-size' and '-eq_rad' can not be used together");
	if (calc_mat_force && beamtype!=B_PLANE)
		PrintError("Currently radiation forces can not be calculated for non-plane incident wave");
	if (InitField==IF_WKB) {
		if (prop_used) PrintError("Currently '-init_field wkb' and '-prop' can not be used together");
		if (beamtype!=B_PLANE) PrintError("'-init_field wkb' is incompatible with non-plane incident wave");
	}
	if (surface) { // currently a lot of limitations for particles near surface
		if (orient_used) PrintError("Currently '-orient' and '-surf' can not be used together");
		if (calc_mat_force) PrintError("Currently calculation of radiation forces is incompatible with '-surf'");
		if (InitField==IF_WKB) PrintError("'-init_field wkb' and '-surf' can not be used together");
		if (!int_surf_used) ReflRelation = msubInf ? GR_IMG : GR_SOM;
		else if (msubInf && ReflRelation!=GR_IMG) PrintError("For perfectly reflecting surface interaction is always "
			"computed through an image dipole. So this case is incompatible with other options to '-int_surf ...'");
		/* TO ADD NEW REFLECTION FORMULATION
		 * Take a look at the above logic, and revise if the new formulation is not fully consistent with it
		 */
	}
	InteractionRealArgs=(beamtype==B_DIPOLE); // other cases may be added here in the future (e.g. nearfields)
#ifdef SPARSE
	if (shape==SH_SPHERE) PrintError("Sparse mode requires shape to be read from file (-shape read ...)");
#endif
#if defined(PARALLEL) && !defined(SPARSE)
	/* Transpose of the non-symmetric interaction matrix can't be done in MPI mode, due to existing memory distribution
	 * among different processors. This causes problems for iterative solvers, which require a product of Hermitian
	 * transpose (calculated through the standard transpose) of the matrix with vector (currently, only CGNR). To remove
	 * this limitation (issue 174) memory distribution should be changed to have both x and gridX-x on the same
	 * processor.
	 */
	if (!reduced_FFT && IterMethod==IT_CGNR) PrintError("Non-symmetric interaction matrix (e.g., -no_reduced_fft) "
		"can not be used together with CGNR iterative solver in the MPI mode");
	/* TO ADD NEW ITERATIVE SOLVER
	 * add the new iterative solver to the above line, if it requires calculation of product of Hermitian transpose
	 * of the matrix with vector (i.e. calls MatVec function with 'true' as the fourth argument)
	 */
#endif
	// scale boxes by jagged; should be completely robust to overflows
#define JAGGED_BOX(a) { \
	if (a!=UNDEF) { \
		if ((BOX_MAX/(size_t)jagged)<(size_t)a) \
			LogError(ONE_POS,"Derived grid size (" #a ") is too large (>%d)",BOX_MAX); \
		else a*=jagged; \
	} }

	if (jagged!=1) {
		JAGGED_BOX(boxX);
		JAGGED_BOX(boxY);
		JAGGED_BOX(boxZ);
	}
#undef JAGGED_BOX
	/* Determine two incident polarizations. Equivalent to rotation of X,Y,Z basis by angles theta and phi from (0,0,1)
	 * to given propagation vector.
	 */
	if (vAlongZ(prop_0)) {
		propAlongZ=incPolX_0[0]=copysign(1.0,prop_0[2]);
		incPolY_0[1]=1;
		incPolX_0[1]=incPolX_0[2]=incPolY_0[0]=incPolY_0[2]=0.0;
	}
	else {
		propAlongZ=0;
		if (sym_type==SYM_ENF && scat_plane) PrintError("Enforcing symmetry for incident propagation not along the "
			"z-axis in combination with defining scattering directions relative to z-axis is not supported");
		temp=sqrt(1-prop_0[2]*prop_0[2]);
		incPolX_0[0]=prop_0[0]*prop_0[2]/temp;
		incPolX_0[1]=prop_0[1]*prop_0[2]/temp;
		incPolX_0[2]=-temp;
		incPolY_0[0]=-prop_0[1]/temp;
		incPolY_0[1]=prop_0[0]/temp;
		incPolY_0[2]=0.0;
	}
	// initialize beam description
	InitBeam();
	// initialize averaging over orientation
	if (orient_avg) {
		ReadAvgParms(avg_parms);
		symX=symY=symZ=symR=false;
		avg_inc_pol=true;
	}
	else { // initialize rotation stuff and test symmetries of the beam in particle reference frame
		InitRotation();
		UpdateSymVec(prop);
		if (beam_asym) UpdateSymVec(beam_center);
	}
	ipr_required=(IterMethod==IT_BICGSTAB || IterMethod==IT_CGNR);
	/* TO ADD NEW ITERATIVE SOLVER
	 * add the new iterative solver to the above line, if it requires inner product calculation during matrix-vector
	 * multiplication (i.e. calls MatVec function with non-NULL third argument)
	 */

	/* TO ADD NEW COMMAND LINE OPTION
	 * If a new command line option may potentially conflict or interact with other options, add here code to implement
	 * corresponding tests or cross-dependence.
	 */
}

//======================================================================================================================

void FinalizeSymmetry(void) {
	// finalize symmetries
	if (sym_type==SYM_NO) symX=symY=symZ=symR=false;
	else if (sym_type==SYM_ENF) symX=symY=symZ=symR=true;
	// test based on SR^2 = SX*SY; uses handmade XOR
	if (symR && ((symX&&!symY) || (symY&&!symX))) LogError(ONE_POS,"Inconsistency in internally defined symmetries");
	// additional tests in case of two polarization runs
	if (!(symR && !scat_grid)) {
		if (beamtype==B_READ && beam_fnameX==NULL)
			PrintError("Only one beam file is specified, while two incident polarizations need to be considered");
		if (InitField==IF_READ && infi_fnameX==NULL) PrintError("Only one file with initial field is specified, while "
			"two incident polarizations need to be considered");
		// the following limitation should be removed in the future
		if (chp_type!=CHP_NONE) PrintError("Currently checkpoints can be used when internal fields are calculated only "
			"once, i.e. for a single incident polarization");
	}
}

//======================================================================================================================

void DirectoryLog(const int argc,char **argv)
// create input directory and start logfile
{
	int i,Nexp;
	FILE * restrict Nexpfile;
	char *compname;
	FILEHANDLE lockid;
#ifdef PARALLEL
	char *ptmp,*ptmp2;
#endif

	// devise directory name (for output files)
	if (directory[0]==0) {
		// root processor works with ExpCount
		if (IFROOT) {
			// lock file
			lockid=CreateLockFile(F_EXPCOUNT_LCK);
			// read ExpCount
			if ((Nexpfile=fopen(F_EXPCOUNT,"r"))!=NULL) {
				if (fscanf(Nexpfile,"%i",&Nexp)!=1) Nexp=0;
				FCloseErr(Nexpfile,F_EXPCOUNT,ONE_POS);
			}
			else Nexp=0;
			// put new number in Nexpfile
			Nexpfile=FOpenErr(F_EXPCOUNT,"w",ONE_POS);
			fprintf(Nexpfile,"%i",Nexp+1);
			FCloseErr(Nexpfile,F_EXPCOUNT,ONE_POS);
			// unlock
			RemoveLockFile(lockid,F_EXPCOUNT_LCK);
		}
		// cast Nexp to all processors
		MyBcast(&Nexp,int_type,1,NULL);
		/* create automatic directory name; It is stored in the following buffer. MAX_LINE should be enough for
		 * auto-name, however up to MAX_DIRNAME can be obtained from '-dir ...'. So the latter size is considered in all
		 * relevant buffers (for filenames or messages).
		 */
		static char sbuffer[MAX_LINE];
		sprintf(sbuffer,"%s%03i_%s_g%i_m"GFORM_RI_DIRNAME,run_name,Nexp,shapename,boxX,creal(ref_index[0]));
#ifdef PARALLEL
		// add PBS, SGE or SLURM job id to the directory name if available
		if ((ptmp=getenv("PBS_JOBID"))!=NULL || (ptmp=getenv("JOB_ID"))!=NULL || (ptmp=getenv("SLURM_JOBID"))!=NULL) {
			// job ID is truncated at first ".", probably can happen only for PBS
			if ((ptmp2=strchr(ptmp,'.'))!=NULL) *ptmp2=0;
			sprintf(sbuffer+strlen(sbuffer),"_id%s",ptmp);
		}
#endif
		directory=sbuffer;
	}
	// make new directory and print info
	if (IFROOT) {
		MkDirErr(directory,ONE_POS);
		printf("all data is saved in '%s'\n",directory);
	}
	// make logname; do it for all processors to enable additional logging in LogError
	if (IFROOT) SnprintfErr(ONE_POS,logfname,MAX_FNAME,"%s/"F_LOG,directory);
	else SnprintfErr(ALL_POS,logfname,MAX_FNAME,"%s/"F_LOG_ERR,directory,ringid);
	// start logfile
	if (IFROOT) {
		// open logfile
		logfile=FOpenErr(logfname,"w",ONE_POS);
		// log version number
		fprintf(logfile,"Generated by ADDA v."ADDA_VERSION"\n");
		// get computer name
#ifdef WINDOWS
		TCHAR cname[MAX_COMPUTERNAME_LENGTH+1];
		DWORD cname_size=MAX_COMPUTERNAME_LENGTH+1;
		GetComputerName(cname,&cname_size);
		compname=cname;
#else // POSIX and others
		compname=getenv("HOST");
#endif
		// write number of processors and computer name
#ifdef PARALLEL
		// write number of processors
		fprintf(logfile,"The program was run on: %d processors (cores)",nprocs);
		// add PBS or SGE host name if present, otherwise use compname
		if ((ptmp=getenv("PBS_O_HOST"))!=NULL || (ptmp=getenv("SGE_O_HOST"))!=NULL) fprintf(logfile," from %s\n",ptmp);
		else if (compname!=NULL) fprintf(logfile," from %s\n",compname);
		else fprintf(logfile,"\n");
#else // sequential
		if (compname!=NULL) fprintf(logfile,"The program was run on: %s\n",compname);
#endif
		// log command line
		fprintf(logfile,"command: '");
		for(i=0;i<argc;i++) fprintf(logfile,"%s ",argv[i]);
		fprintf(logfile,"'\n");
	}
	Synchronize(); // needed to wait for creation of the output directory
	LogPending();
}

//======================================================================================================================

void PrintInfo(void)
// print info to stdout and logfile
{
	int i;
	char sbuffer[MAX_LINE];

	if (IFROOT) {
		// print basic parameters
		printf("box dimensions: %ix%ix%i\n",boxX,boxY,boxZ);
		printf("lambda: "GFORM"   Dipoles/lambda: "GFORMDEF"\n",lambda,dpl);
		printf("Required relative residual norm: "GFORMDEF"\n",iter_eps);
		printf("Total number of occupied dipoles: %zu\n",nvoid_Ndip);
		// log basic parameters
		fprintf(logfile,"lambda: "GFORM"\n",lambda);
		fprintf(logfile,"shape: ");
		fprintf(logfile,"%s"GFORM"%s\n",sh_form_str1,sizeX,sh_form_str2);
#ifndef SPARSE
		if (sh_granul) fprintf(logfile,"  domain %d is filled with %d granules of diameter "GFORMDEF"\n"
			"    volume fraction: specified - "GFORMDEF", actual - "GFORMDEF"\n",gr_mat+1,gr_N,gr_d,gr_vf,gr_vf_real);
#endif // SPARSE
		fprintf(logfile,"box dimensions: %ix%ix%i\n",boxX,boxY,boxZ);
		if (anisotropy) {
			fprintf(logfile,"refractive index (diagonal elements of the tensor):\n");
			if (Nmat==1) fprintf(logfile,"    "CFORM3V"\n",REIM3V(ref_index));
			else {
				for (i=0;i<Nmat;i++) {
					if (i<Nmat_given) fprintf(logfile,"    %d. "CFORM3V"\n",i+1,REIM3V(ref_index+3*i));
					else fprintf(logfile,"   %d. not specified\n",i+1);
				}
			}
		}
		else {
			fprintf(logfile,"refractive index: ");
			if (Nmat==1) fprintf(logfile,CFORM"\n",REIM(ref_index[0]));
			else {
				fprintf(logfile,"1. "CFORM"\n",REIM(ref_index[0]));
				for (i=1;i<Nmat;i++) {
					if (i<Nmat_given) fprintf(logfile,"                  %d. "CFORM"\n",i+1,REIM(ref_index[i]));
					else fprintf(logfile,"                  %d. not specified\n",i+1);
				}
			}
		}
		if (surface) {
			if (msubInf) fprintf(logfile,"Particle is placed near the perfectly reflecting substrate\n");
			else fprintf(logfile,"Particle is placed near the substrate with refractive index "CFORM",\n",REIM(msub));
			fprintf(logfile,"  height of the particle center: "GFORMDEF"\n",hsub);
		}
		fprintf(logfile,"Dipoles/lambda: "GFORMDEF"\n",dpl);
		if (volcor_used) fprintf(logfile,"\t(Volume correction used)\n");
		fprintf(logfile,"Required relative residual norm: "GFORMDEF"\n",iter_eps);
		fprintf(logfile,"Total number of occupied dipoles: %zu\n",nvoid_Ndip);
		if (Nmat>1) {
			fprintf(logfile,"  per domain: 1. %zu\n",mat_count[0]);
			for (i=1;i<Nmat;i++) fprintf(logfile,"              %d. %zu\n",i+1,mat_count[i]);
		}
		fprintf(logfile,"Volume-equivalent size parameter: "GFORM"\n",ka_eq);
		// log incident beam and polarization
		fprintf(logfile,"\n---In laboratory reference frame:---\nIncident beam: %s\n",beam_descr);
		fprintf(logfile,"Incident propagation vector: "GFORMDEF3V"\n",COMP3V(prop_0));
		if (beamtype==B_DIPOLE) fprintf(logfile,"(dipole orientation)\n");
		else { // polarizations are not shown for dipole incident field
			fprintf(logfile,"Incident polarization Y(par): "GFORMDEF3V"\n",COMP3V(incPolY_0));
			fprintf(logfile,"Incident polarization X(per): "GFORMDEF3V"\n",COMP3V(incPolX_0));
		}
		if (surface && beamtype!=B_DIPOLE) { // include surface-specific vectors
			fprintf(logfile,"Reflected propagation vector: "GFORMDEF3V"\n",COMP3V(prIncRefl));
			fprintf(logfile,"Transmitted propagation vector: ");
			if (msubInf) fprintf (logfile,"none\n");
			else {
				fprintf(logfile,GFORMDEF3V,COMP3V(prIncTran));
				if (prIncTran[2]==0) fprintf(logfile," (evanescent)\n");
				else fprintf(logfile,"\n");
			}
		}
		fprintf(logfile,"\n");
		// log particle orientation
		if (orient_avg) fprintf(logfile,"Particle orientation - averaged\n%s\n",avg_string);
		else {
			// log incident polarization after transformation
			if (alph_deg!=0 || bet_deg!=0 || gam_deg!=0) {
				fprintf(logfile,"Particle orientation (deg): alpha="GFORMDEF", beta="GFORMDEF", gamma="GFORMDEF"\n\n"
					"---In particle reference frame:---\n",alph_deg,bet_deg,gam_deg);
				if (beam_asym) fprintf(logfile,"Incident Beam center position: "GFORMDEF3V"\n",
					beam_center[0],beam_center[1],beam_center[2]);
				fprintf(logfile,"Incident propagation vector: "GFORMDEF3V"\n",COMP3V(prop));
				if (beamtype==B_DIPOLE) fprintf(logfile,"(dipole orientation)\n");
				else { // polarizations are not shown for dipole incident field
					fprintf(logfile,"Incident polarization Y(par): "GFORMDEF3V"\n",COMP3V(incPolY));
					fprintf(logfile,"Incident polarization X(per): "GFORMDEF3V"\n",COMP3V(incPolX));
				}
			}
			else fprintf(logfile,"Particle orientation: default\n");
			fprintf(logfile,"\n");
		}
		// log type of scattering matrices
		if (store_mueller) {
			if (store_ampl) fprintf(logfile,"Calculating both amplitude and Mueller scattering matrices\n");
			else fprintf(logfile,"Calculating only Mueller scattering matrix\n");
		}
		else {
			if (store_ampl) fprintf(logfile,"Calculating only amplitude scattering matrix\n");
			else fprintf(logfile,"Calculating no scattering matrices\n");
		}
		// log polarizability relation
		fprintf(logfile,"Polarizability relation: ");
		switch (PolRelation) {
			case POL_CLDR: fprintf(logfile,"'Corrected Lattice Dispersion Relation'\n"); break;
			case POL_CM: fprintf(logfile,"'Clausius-Mossotti'\n"); break;
			case POL_DGF: fprintf(logfile,"'Digitized Green's Function'\n"); break;
			case POL_FCD: fprintf(logfile,"'Filtered Coupled Dipoles'\n"); break;
			case POL_IGT_SO: fprintf(logfile,"'Integration of Green's Tensor [approximation O(kd^2)]'\n"); break;
			case POL_LAK: fprintf(logfile,"'by Lakhtakia'\n"); break;
			case POL_LDR:
				fprintf(logfile,"'Lattice Dispersion Relation'");
				if (avg_inc_pol) fprintf(logfile," (averaged over incident polarization)");
				fprintf(logfile,"\n");
				break;
			case POL_NLOC:
				fprintf(logfile,"'Non-local' (based on lattice sum, Gaussian width Rp="GFORMDEF")\n",polNlocRp);
				break;
			case POL_NLOC_AV:
				fprintf(logfile,"'Non-local' (averaged, Gaussian width Rp="GFORMDEF")\n",polNlocRp);
				break;
			case POL_RRC: fprintf(logfile,"'Radiative Reaction Correction'\n"); break;
			case POL_SO: fprintf(logfile,"'Second Order'\n"); break;
		}
		/* TO ADD NEW POLARIZABILITY FORMULATION
		 * add a case above in the alphabetical order, analogous to the ones already present. The variable parts of the
		 * case are descriptor, defined in const.h, and its plain-text description (to be shown in log).
		 */
		// log Scattering Quantities formulae
		fprintf(logfile,"Scattering quantities formulae: ");
		switch (ScatRelation) {
			case SQ_DRAINE: fprintf(logfile,"'by Draine'\n"); break;
			case SQ_FINDIP: fprintf(logfile,"'Finite Dipoles'\n"); break;
			case SQ_IGT_SO: fprintf(logfile,"'Integration of Green's Tensor [approximation O(kd^2)]'\n"); break;
			case SQ_SO: fprintf(logfile,"'Second Order'\n"); break;
		}
		// log Interaction term prescription
		fprintf(logfile,"Interaction term prescription: ");
		switch (IntRelation) {
			case G_FCD: fprintf(logfile,"'Filtered Green's tensor'\n"); break;
			case G_FCD_ST: fprintf(logfile,"'Filtered Green's tensor (quasistatic)'\n"); break;
			case G_IGT:
				fprintf(logfile,"'Integrated Green's tensor' (accuracy "GFORMDEF", ",igt_eps);
				if (igt_lim==UNDEF) fprintf(logfile,"no distance limit)\n");
				else fprintf(logfile,"for distance < "GFORMDEF" dipole sizes)\n",igt_lim);
				break;
			case G_IGT_SO: fprintf(logfile,"'Integrated Green's tensor [approximation O(kd^2)]'\n"); break;
			case G_NLOC: fprintf(logfile,"'Non-local' (point-value, Gaussian width Rp="GFORMDEF")\n",nloc_Rp); break;
			case G_NLOC_AV: fprintf(logfile,"'Non-local' (averaged, Gaussian width Rp="GFORMDEF")\n",nloc_Rp); break;
			case G_POINT_DIP: fprintf(logfile,"'as Point dipoles'\n"); break;
			case G_SO: fprintf(logfile,"'Second Order'\n"); break;
		}
		/* TO ADD NEW INTERACTION FORMULATION
		 * add a case above in the alphabetical order, analogous to the ones already present. The variable parts of the
		 * case are descriptor, defined in const.h, and its plain-text description (to be shown in log).
		 */
		// log reflected Green'tensor formulae
		if (surface) {
			fprintf(logfile,"Reflected Green's tensor formulae: ");
			switch (ReflRelation) {
				case GR_IMG: fprintf(logfile,"'Image-dipole approximation'\n"); break;
				case GR_SOM: fprintf(logfile,"'Sommerfeld integrals'\n"); break;
			}
		}
		/* TO ADD NEW REFLECTION FORMULATION
		 * add a case above in the alphabetical order, analogous to the ones already present. The variable parts of the
		 * case are descriptor, defined in const.h, and its plain-text description (to be shown in log).
		 */
		// log FFT and (if needed) clFFT method
		fprintf(logfile,"FFT algorithm: ");
#ifdef FFTW3
		fprintf(logfile,"FFTW3\n");
#elif defined(FFT_TEMPERTON)
		fprintf(logfile,"by C.Temperton\n");
#elif defined(SPARSE)
		fprintf(logfile,"none (sparse mode)\n");
#endif
#if defined(OPENCL) && !defined(SPARSE)
		fprintf(logfile,"OpenCL FFT algorithm: ");
#	ifdef CLFFT_AMD
		fprintf(logfile,"by AMD\n");
#	elif defined(CLFFT_APPLE)
		fprintf(logfile,"by Apple\n");
#	endif
#endif
		// log Iterative Method
		fprintf(logfile,"Iterative Method: ");
		switch (IterMethod) {
			case IT_BCGS2: fprintf(logfile,"Enhanced Bi-CG Stabilized(2)\n"); break;
			case IT_BICG_CS: fprintf(logfile,"Bi-CG (complex symmetric)\n"); break;
			case IT_BICGSTAB: fprintf(logfile,"Bi-CG Stabilized\n"); break;
			case IT_CGNR: fprintf(logfile,"CGNR\n"); break;
			case IT_CSYM: fprintf(logfile,"CSYM\n"); break;
			case IT_QMR_CS: fprintf(logfile,"QMR (complex symmetric)\n"); break;
			case IT_QMR_CS_2: fprintf(logfile,"2-term QMR (complex symmetric)\n"); break;
		}
		/* TO ADD NEW ITERATIVE SOLVER
		 * add a case above in the alphabetical order, analogous to the ones already present. The variable parts of the
		 * case are descriptor, defined in const.h, and its plain-text description (to be shown in log).
		 */
		// log Symmetry options
		switch (sym_type) {
			case SYM_AUTO:
				fprintf(logfile,"Symmetries: ");
				if (symX) fprintf(logfile,"X");
				if (symY) fprintf(logfile,"Y");
				if (symZ) fprintf(logfile,"Z");
				if (symR) fprintf(logfile,"R");
				if (!(symX||symY||symZ||symR)) fprintf(logfile,"none");
				fprintf(logfile,"\n");
				break;
			case SYM_NO: fprintf(logfile,"Symmetries: cancelled by user\n"); break;
			case SYM_ENF: fprintf(logfile,"Symmetries: enforced by user (warning!)\n"); break;
		}
		// log optimization method
		if (save_memory) fprintf(logfile,"Optimization is done for minimum memory usage\n");
		else fprintf(logfile,"Optimization is done for maximum speed\n");
		// log Checkpoint options
		if (load_chpoint) fprintf(logfile,"Simulation is continued from a checkpoint\n");
		if (chp_type!=CHP_NONE) {
			fprintf(logfile,"Checkpoint is turned on:\n");
			switch (chp_type) {
				case CHP_NONE: break; // redundant; to keep a set of cases complete
				case CHP_NORMAL: fprintf(logfile,"    type = normal\n"); break;
				case CHP_REGULAR: fprintf(logfile,"    type = regular\n"); break;
				case CHP_ALWAYS: fprintf(logfile,"    type = always\n"); break;
			}
			if (chp_time==UNDEF) fprintf(logfile,"    time = no limit\n");
			else {
				PrintTime(sbuffer,&chp_time);
				// chp_time is converted to long to avoid problems with definition of time_t (can be either int or long)
				fprintf(logfile,"    time = %s(%ld sec)\n",sbuffer,(long)chp_time);
			}
		}
		if (load_chpoint || chp_type!=CHP_NONE) fprintf(logfile,"    directory = '%s'\n",chp_dir);
		/* TO ADD NEW COMMAND LINE OPTION
		 * If a new command line option requires additional output to log file or stdout, implement this functionality
		 * here, choosing appropriate place for this information among the above lines.
		 */
	}
}
