/* File: const.h
 * $Author$
 * $Date::                            $
 * Descr: all the constants used by ADDA code, including enum constants, also defines some
 *        useful macros
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * Copyright (C) 2009 Institute of Chemical Kinetics and Combustion & University of Amsterdam
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __const_h
#define __const_h

// version number (string)
#define ADDA_VERSION "0.80a2"

/* ADDA uses certain C99 extensions, which are widely supported by GNU and Intel compilers. However,
 * they may be not completely supported by e.g. Microsoft Visual Studio compiler. Therefore, we
 * check the version of the standard here and produce a strong warning, if it is not satisfied.
 * The list of C99 features, used by ADDA, include (but may be not limited to):
 * stdbool.h, snprintf, %z argument in printf, '//' comments, restricted pointers, variadic macros
*/
# if !defined(OVERRIDE_STDC_TEST) && (!defined(__STDC_VERSION__) || (__STDC_VERSION__ < 199901L))
#   error Support for C99 standard (at least many of its parts) is strongly recommended for \
          compilation. Otherwise the compilation will may fail or produce wrong results. If you \
          still want to try, you may enable an override in the Makefile.
#endif


// basic constants
#define UNDEF -1 // should be used only for variables, which are naturally non-negative
	// denotes that shape accepts single filename argument; used in definitions of suboptions
#define FNAME_ARG -2

// simple functions
#define MIN(A,B) (((A) > (B)) ? (B) : (A))
#define MAX(A,B) (((A) < (B)) ? (B) : (A))
#define IS_EVEN(A) (((A)%2) == 0)
#define LENGTH(A) ((int)(sizeof(A)/sizeof(A[0]))) // length of any array (converted to int)

// parallel definitions
#ifdef MPI
#define PARALLEL
#endif

/* ringid of root processor. Using ROOT!=0 should work, however it was not thoroughly tested.
 * Hence do not change without necessity.
 */
#define ROOT 0

// math constants rounded for 32 decimals
#define PI                  3.1415926535897932384626433832795
#define TWO_PI              6.283185307179586476925286766559
#define FOUR_PI             12.566370614359172953850573533118
#define EIGHT_PI            25.132741228718345907701147066236
#define FOUR_PI_OVER_THREE  4.1887902047863909846168578443727
#define PI_OVER_TWO         1.5707963267948966192313216916398
#define PI_OVER_FOUR        0.78539816339744830961566084581988
#define PI_OVER_SIX         0.52359877559829887307710723054658
#define INV_PI              0.31830988618379067153776752674503
#define TWO_OVER_PI         0.63661977236758134307553505349006
#define THREE_OVER_FOUR_PI  0.23873241463784300365332564505877
#define SIX_OVER_PI         1.9098593171027440292266051604702
#define ONE_THIRD           0.33333333333333333333333333333333
#define PI_OVER_180         0.017453292519943295769236907684886
#define INV_PI_180          57.295779513082320876798154814105
#define EULER               0.57721566490153286060651209008241
#define FULL_ANGLE          360.0

/* determines the maximum number representable by size_t. Actually, SIZE_MAX is part of the C99
 * standard, but we leave this code to be.
 */
#ifndef SIZE_MAX
#define SIZE_MAX ((size_t)-1)
#endif

// sets the maximum box size; otherwise 'position' should be changed
#define BOX_MAX USHRT_MAX

// sizes of some arrays
#define MAX_NMAT         15   // maximum number of different refractive indices (<256)
#define MAX_N_SH_PARMS   25   // maximum number of shape parameters
#define MAX_N_BEAM_PARMS 10   // maximum number of beam parameters

// sizes of filenames and other strings
#define MAX_DIRNAME      300 // maximum length of dirname; increase THIS if any errors appear
#define MAX_FNAME_SH     100 // maximum length of filename (used for known names)
#define MAX_TMP_FNAME_SH  15 // maximum length of names of temporary files (short)
#define MAX_SYSTEM_CALL   10 // maximum string length of system call (itself)
#define MAX_WORD          10 // maximum length of a short word
#define MAX_LINE          50 // maximum length of a line
	// size of buffer for reading lines (longer lines are handled robustly)
#define BUF_LINE      150
#define MAX_PARAGRAPH 600 // maximum length of a paragraph (few lines)

// derived sizes
	// maximum string to create directory
#define MAX_DIRSYS  (MAX_DIRNAME + MAX_SYSTEM_CALL)
	// maximum length of filename (including directory name)
#define MAX_FNAME   (MAX_DIRNAME + MAX_FNAME_SH)
	// maximum length of temporary filename (including directory name)
#define MAX_TMP_FNAME   (MAX_DIRNAME + MAX_TMP_FNAME_SH)
	// maximum message that may include a filename (for PrintError)
#define MAX_MESSAGE (MAX_FNAME + MAX_PARAGRAPH)
	// maximum message that may include 2 filenames (for LogError)
#define MAX_MESSAGE2 (2*MAX_FNAME + MAX_PARAGRAPH)

// widths of terminal used for output
#define DEF_TERM_WIDTH 80 // default
#define MIN_TERM_WIDTH 20 // ADDA never takes value less than that from environmental variables

enum sh { // shape types
	SH_AXISYMMETRIC, // axisymmetric
	SH_BOX,          // box (may be rectangular)
	SH_CAPSULE,      // capsule
	SH_COATED,       // coated sphere
	SH_CYLINDER,     // cylinder
	SH_EGG,          // egg
	SH_ELLIPSOID,    // general ellipsoid
	SH_LINE,         // line with width of one dipole
	SH_PRISMA,       // prisma (triangular) -- not operational
	SH_RBC,          // Red Blood Cell
	SH_READ,         // read from file
	SH_SDISK_ROT,    // disc cut of a sphere -- not operational
	SH_SPHERE,       // sphere
	SH_SPHEREBOX     // sphere in a box
	/* TO ADD NEW SHAPE
	 * add an identifier starting with 'SH_' and a descriptive comment to this list in alphabetical
	 * order.
	 */
};

enum pol { // which way to calculate coupleconstant
	POL_CM,   // Clausius-Mossotti
	POL_RRC,  // Radiative Reaction correction
	POL_LDR,  // Lattice Dispersion Relation
	POL_CLDR, // Corrected Lattice Dispersion Relation
	POL_FCD,  // Filtered Coupled Dipoles
	POL_SO    // Second Order formulation
};

enum scat { // how to calculate scattering quantities
	SQ_DRAINE, // classical, as Draine
	SQ_FINDIP, /* Same as Draine, but with correction of radiation energy of a _finite_ dipole when
	            * calculating absorption cross section
	            */
	SQ_SO      // Second Order formulation
};

enum inter { // how to calculate interaction term
	G_POINT_DIP, // as point dipoles
	G_FCD,       // Filtered Green's tensor (Filtered Coupled Dipoles)
	G_FCD_ST,    // quasi-static version of FCD
	G_SO         // Second Order formulation
};

// ldr constants
#define LDR_B1  1.8915316
#define LDR_B2 -0.1648469
#define LDR_B3  1.7700004

// 2nd_order constants
#define SO_B1 1.5867182
#define SO_B2 0.13488017
#define SO_B3 0.11895826

// two boundaries for separation between G_SO 'close', 'median', and 'far'
#define G_BOUND_CLOSE  1 // k*R^2/d < GB_CLOSE => 'close'
#define G_BOUND_MEDIAN 1 // k*R < GB_MEDIAN => 'median'

enum iter { // iterative methods
	IT_CGNR,     // Conjugate Gradient for Normalized equations minimizing Residual norm
	IT_BICGSTAB, // Bi-Conjugate Gradient Stabilized
	IT_BICG_CS,  // Bi-Conjugate Gradient for Complex-Symmetric matrices
	IT_QMR_CS    // Quasi-minimal residual for Complex-Symmetric matrices
};

enum Eftype { // type of E field calculation
	CE_NORMAL, // normal
	CE_PARPER  /* use symmetry to calculate both incident polarizations from one calculation of
	            * internal fields
	            */
};

// path and size of tables
#define TAB_PATH     "tables/"
#define TAB_FNAME(a) "t" #a "f.dat" // a is a number, e.g. TAB_FNAME(2) -> "t2f.dat"
#define TAB_SIZE     142
#define TAB_RMAX     10

enum beam { // beam types
	B_BARTON5, // 5th order description of the Gaussian beam
	B_DAVIS3,  // 3rd order description of the Gaussian beam
	B_LMINUS,  // 1st order description of the Gaussian beam
	B_PLANE    // infinite plane wave
	/* TO ADD NEW BEAM
	 * add an identifier starting with 'B_' and a descriptive comment to this list in alphabetical
	 * order.
	 */
};

enum scatgrid { // types of scattering grid
	SG_GRID, // grid of angles
	SG_PAIRS // set of independent pairs
};
enum angleset {// types of angles set
	AS_RANGE, // range with uniformly spaced points
	AS_VALUES // any set of values
};

// types of phi_integr (should be different one-bit numbers)
#define PHI_UNITY 1 // just integrate
#define PHI_COS2  2 // integrate with cos(2*phi)
#define PHI_SIN2  4 // integrate with sin(2*phi)
#define PHI_COS4  8 // integrate with cos(4*phi)
#define PHI_SIN4 16 // integrate with sin(4*phi)

enum sym { // ways to treat particle symmetries
	SYM_AUTO, // automatic
	SYM_NO,   // do not take into account
	SYM_ENF   // enforce
};

enum chpoint { // types of checkpoint (to save)
	CHP_NONE,    // do not save checkpoint
	CHP_NORMAL,  // save checkpoint if not finished in time and exit
	CHP_REGULAR, // save checkpoints in regular time intervals (until finished or halted)
	CHP_ALWAYS   /* save checkpoint if either simulation is finished or time elapsed and calculate
	              * all scattering quantities
                  */
};

// return values for functions
#define CHP_EXIT -2 // exit after saving checkpoint

// default values; other are specified in InitVariables (param.c)
#define DEF_GRID       (16*jagged)
#define MIN_AUTO_GRID  16 // minimum grid, when set from default dpl

// numbers less than this value (compared to unity) are considered to be zero
#define ROUND_ERR 1E-15

// output and input file and directory names (can only be changed at compile time)
#define F_EXPCOUNT      "ExpCount"
#define F_EXPCOUNT_LCK  F_EXPCOUNT ".lck"
#define F_CS            "CrossSec"
#define F_FRP           "VisFrp"
#define F_INTFLD        "IntField"
#define F_DIPPOL        "DipPol"
#define F_BEAM          "IncBeam"
#define F_GRANS         "granules"
	// suffixes
#define F_XSUF          "-X"
#define F_YSUF          "-Y"
	// logs
#define F_LOG           "log"
#define F_LOG_ERR       "logerr.%d"    // ringid as argument
#define F_LOG_ORAVG     "log_orient_avg"
#define F_LOG_INT_CSCA  "log_int_Csca"
#define F_LOG_INT_ASYM  "log_int_asym"
	// log suffixes
#define F_LOG_X         "_x"
#define F_LOG_Y         "_y"
#define F_LOG_Z         "_z"
	// Mueller files
#define F_MUEL          "mueller"
#define F_MUEL_SG       "mueller_scatgrid"
#define F_MUEL_INT      "mueller_integr"
#define F_MUEL_C2       "mueller_integr_c2"
#define F_MUEL_S2       "mueller_integr_s2"
#define F_MUEL_C4       "mueller_integr_c4"
#define F_MUEL_S4       "mueller_integr_s4"
	// temporary files; used in printf with ringid as argument
#define F_BEAM_TMP      "b%d.tmp"
#define F_INTFLD_TMP    "f%d.tmp"
#define F_DIPPOL_TMP    "p%d.tmp"
#define F_GEOM_TMP      "g%d.tmp"
	// checkpoint files
#define F_CHP_LOG       "chp.log"
#define F_CHP           "chp.%d"   // ringid as argument

// default file and directory names; can be changed by command line options
#define FD_ALLDIR_PARMS "alldir_params.dat"
#define FD_AVG_PARMS    "avg_params.dat"
#define FD_SCAT_PARMS   "scat_params.dat"
#define FD_CHP_DIR      "chpoint"

/* number of components of D. Really, it can't be easily changed, but using constant instead of 6
 * adds more understanding for reader
 */
#define NDCOMP 6

// shape formats; numbers should be nonnegative
#define SF_TEXT     0 // ADDA text format for one-domain particles
#define SF_TEXT_EXT 1 // ADDA text format for multi-domain particles
#define SF_DDSCAT   2 // DDSCAT 6.1 format (FRMFIL), produced by calltarget

#define POSIT __FILE__,__LINE__ // position of the error in source code

enum enwho { // who is calling
	ALL, // each processor may report this error
	ONE  // only root processor reports an error
};

// derived; for simplicity
#define ALL_POS ALL,POSIT
#define ONE_POS ONE,POSIT

enum ec { // error codes
	EC_OK,    // no error, corresponds to zero exit code
	EC_ERROR, // error
	EC_WARN,  // warning
	EC_INFO   // slight warning, that does not interfere at all with normal execution
};

#endif // __const_h
