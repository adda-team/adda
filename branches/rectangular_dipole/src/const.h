/* File: const.h
 * $Date::                            $
 * Descr: all the constants used by ADDA code, including enum constants, also defines some useful macros
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
#ifndef __const_h
#define __const_h

// version number (string)
#define ADDA_VERSION "1.3b6"

/* ADDA uses certain C99 extensions, which are widely supported by GNU and Intel compilers. However, they may be not
 * completely supported by e.g. Microsoft Visual Studio compiler. Therefore, we check the version of the standard here
 * and produce a strong warning, if it is not satisfied. The list of C99 features, used by ADDA, include (but may be not
 * limited to): stdbool.h, snprintf, %z argument in printf, '//' comments, restricted pointers, variadic macros
*/
# if !defined(OVERRIDE_STDC_TEST) && (!defined(__STDC_VERSION__) || (__STDC_VERSION__ < 199901L))
#   error "Support for C99 standard (at least many of its parts) is strongly recommended for compilation. Otherwise \
	the compilation will may fail or produce wrong results. If you still want to try, you may enable an override in \
	the Makefile."
#endif

/* The following is to ensure that mingw64 with "-std=c99" will use c99-compliant printf-family functions. For some
 * (philosophical) reasons mingw64 developers have not implemented this behavior as the default one. So we need to set
 * it manually. This macro should be defined before any system includes, hence inclusion of "const.h" should be the
 * first one in all sources. This is also convenient for testing c99 standard above. However, there is no simple way
 * then to test for MinGW64 at this point, since, e.g., __MINGW64_VERSION_STR is defined by the system header (not by
 * the compiler itself). So the code executes always. Not the most reliable way, but seems the only way to keep the
 * following definition in one place.
 */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L)
#	define __USE_MINGW_ANSI_STDIO 1
#endif

// basic constants
#define UNDEF -1 // should be used only for variables, which are naturally non-negative
	// denotes that shape accepts filename arguments; used in definitions of options
#define FNAME_ARG     -2  // single filename
#define FNAME_ARG_2   -3  // two filenames
#define FNAME_ARG_1_2 -4  // 1 or 2 filenames
	// macro to test for occurrence of one of FNAME_ARG
#define IS_FNAME_ARG(A) (((A)==FNAME_ARG) || ((A)==FNAME_ARG_2) || ((A)==FNAME_ARG_1_2))

// simple functions
#define MIN(A,B) (((A) > (B)) ? (B) : (A))
#define MAX(A,B) (((A) < (B)) ? (B) : (A))
#define MAXIMIZE(A,B) {if ((A)<(B)) (A)=(B);}
#define IS_EVEN(A) (((A)%2) == 0)
#define DIV_CEILING(A,B) (((A)%(B)==0) ? (A)/(B) : ((A)/(B))+1 ) // valid only for nonnegative A and B
#define LENGTH(A) ((int)(sizeof(A)/sizeof(A[0]))) // length of any array (converted to int)
#define STRINGIFY(A) #A
#define GREATER_EQ2(a1,a2,b1,b2) ( (a1)>(b1) || ( (a1)==(b1) && (a2)>=(b2) )) // a1.a2>=b1.b2

// parallel definitions
#ifdef ADDA_MPI
#	define PARALLEL
#endif

/* ringid of root processor. Using ADDA_ROOT!=0 should work, however it was not thoroughly tested. Hence do not change
 * without necessity.
 */
#define ADDA_ROOT 0

// math constants rounded for 32 decimals; C99 standard specifies that they are encoded as double
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
#define SQRT_PI             1.7724538509055160272981674833411
#define TWO_OVER_SQRT_PI    1.1283791670955125738961589031215
#define SQRT2               1.4142135623730950488016887242097
#define SQRT3               1.7320508075688772935274463415059
#define SQRT1_2             0.70710678118654752440084436210485
#define SQRT1_2PI           0.39894228040143267793994605993438
#define SQRT2_9PI           0.26596152026762178529329737328959
#define EULER               0.57721566490153286060651209008241
#define FULL_ANGLE          360.0
#define MICRO               1E-6

// sets the maximum box size; otherwise 'position' should be changed
#define BOX_MAX USHRT_MAX

// sizes of some arrays
#define MAX_NMAT         15   // maximum number of different refractive indices (<256)
#define MAX_N_SH_PARMS   25   // maximum number of shape parameters
#define MAX_N_BEAM_PARMS 10   // maximum number of beam parameters

// sizes of filenames and other strings
/* There is MAX_PATH constant that equals 260 on Windows. However, even this OS allows ways to override this limit. On
 * POSIX this constant may be much larger and have even less reliability. So the values below are conservatively high,
 * but are not guaranteed to suffice. However, the functions in the code are designed to survive if this buffer won't
 * suffice.
 */
#define MAX_DIRNAME      300 // maximum length of dirname; increase THIS if any errors appear
#define MAX_FNAME_SH     100 // maximum length of filename (used for known names)
#define MAX_TMP_FNAME_SH  15 // maximum length of names of temporary files (short)
#define MAX_SYSTEM_CALL   10 // maximum string length of system call (itself)
#define MAX_WORD          10 // maximum length of a short word
#define MAX_LINE         100 // maximum length of a line
#define BUF_LINE         300 // size of buffer for reading lines (longer lines are handled robustly)
#define MAX_PARAGRAPH    600 // maximum length of a paragraph (few lines)

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

// formats for outputs of float values
#define EFORM "%.10E"        // fixed width
#define GFORM "%.10g"        // variable width (showing significant digits)
#define GFORMDEF "%g"        // default output for non-precise values
#define GFORM_DEBUG "%.2g"   // for debug and error output
#define CFORM "%.10g%+.10gi" // for complex numbers; may be defined in terms of GFORM
	// derived formats; starting "" is to avoid redundant syntax errors in Eclipse
#define GFORM3V "("GFORM","GFORM","GFORM")"
#define GFORM3L ""GFORM" "GFORM" "GFORM
#define GFORM6L ""GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM
#define GFORM7L ""GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM
#define GFORM10L ""GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM" "GFORM
#define GFORMDEF3V "("GFORMDEF","GFORMDEF","GFORMDEF")"
#define CFORM3V "("CFORM","CFORM","CFORM")"
	// macros to shorten printing of all vector components
#define COMP3V(a) (a)[0],(a)[1],(a)[2]
#define COMP16V(a) (a)[0],(a)[1],(a)[2],(a)[3],(a)[4],(a)[5],(a)[6],(a)[7],(a)[8],(a)[9],(a)[10],(a)[11],(a)[12],\
	(a)[13],(a)[14],(a)[15]

enum sh { // shape types
	SH_AXISYMMETRIC, // axisymmetric
	SH_BICOATED,     // two coated spheres
	SH_BIELLIPSOID,  // two general ellipsoids
	SH_BISPHERE,     // two spheres
	SH_BOX,          // box (may be rectangular)
	SH_CAPSULE,      // capsule
	SH_CHEBYSHEV,    // Chebyshev particle (axisymmetric)
	SH_COATED,       // coated sphere
	SH_CYLINDER,     // cylinder
	SH_EGG,          // egg
	SH_ELLIPSOID,    // general ellipsoid
	SH_LINE,         // line with width of one dipole
	SH_PLATE,        // plate
	SH_PRISM,        // right rectangular prism
	SH_RBC,          // Red Blood Cell
	SH_READ,         // read from file
	//SH_SDISK_ROT,  // disc cut of a sphere -- not operational
	SH_SPHERE,       // sphere
	SH_SPHEREBOX     // sphere in a box
	/* TO ADD NEW SHAPE
	 * add an identifier starting with 'SH_' and a descriptive comment to this list in alphabetical order.
	 */
};

enum pol { // which way to calculate coupleconstant
	POL_CLDR,    // Corrected Lattice Dispersion Relation
	POL_CM,      // Clausius-Mossotti
	POL_DGF,     // Digitized Green's Function (second order approximation of LAK)
	POL_FCD,     // Filtered Coupled Dipoles
	POL_IGT_SO,  // Second order approximation to Green's tensor integrated over a cube
	POL_LAK,     // Exact result of IGT for sphere
	POL_LDR,     // Lattice Dispersion Relation
	POL_NLOC,    // non-local extension (Gaussian dipole-density, formula based on lattice sums)
	POL_NLOC_AV, // same as NLOC, but based on averaging of Gaussian over the dipole volume
	POL_RRC,     // Radiative Reaction correction
	POL_SO       // Second Order formulation
	/* TO ADD NEW POLARIZABILITY FORMULATION
	 * add an identifier starting with 'POL_' and a descriptive comment to this list in the alphabetical order.
	 */
};
// in alphabetical order

enum scat { // how to calculate scattering quantities
	SQ_DRAINE, // classical, as Draine
	SQ_FINDIP, /* Same as Draine, but with correction of radiation energy of a _finite_ dipole when calculating
	              absorption cross section */
	SQ_IGT_SO, // Integration of Green's tensor (second order in kd approximation)
	SQ_SO      // Second Order formulation
};
// in alphabetical order

enum inter { // how to calculate interaction term
	G_FCD,       // Filtered Green's tensor (Filtered Coupled Dipoles)
	G_FCD_ST,    // quasi-static version of FCD
	G_IGT,       // (direct) integration of Green's tensor
	G_IGT_SO,    // approximate integration of Green's tensor (based on ideas of SO)
	G_NLOC,      // non-local extension (interaction of Gaussian dipole-densities)
	G_NLOC_AV,   // same as NLOC, but based on averaging of Gaussian over the dipole volume
	G_POINT_DIP, // as point dipoles
	G_SO         // Second Order formulation
	/* TO ADD NEW INTERACTION FORMULATION
	 * add an identifier starting with 'G_' and a descriptive comment to this list in the alphabetical order.
	 */
};
enum refl { // how to calculate interaction of dipoles through the nearby surface (reflected G)
	GR_IMG,       // approximate expression based on a single image dipole
	GR_SOM        // direct evaluation of Sommerfeld integrals
	/* TO ADD NEW REFLECTION FORMULATION
	 * add an identifier starting with 'GR_' and a descriptive comment to this list in the alphabetical order.
	 */
};// in alphabetical order

// ldr constants
#define LDR_B1  1.8915316
#define LDR_B2 -0.1648469
#define LDR_B3  1.7700004

// 2nd_order constants; derived from c1=ln(5+3^(3/2))-ln(2)/2-pi/4 and c2=pi/6
#define SO_B1 1.5867182426530356710958782335228 // 4c1/3
#define SO_B2 0.13488017286410948123541594310740 // c1/3 - c2/2
#define SO_B3 0.11895825700597042937085940122438 // (5/2)c1 - c2

// other constants for polarizability
#define DGF_B1 1.6119919540164696407169668466385  // (4pi/3)^(1/3)
#define LAK_C  0.62035049089940001666800681204778 // (4pi/3)^(-1/3)

// two boundaries for separation between G_SO 'close', 'median', and 'far'
#define G_BOUND_CLOSE  1 // k*R^2/d < GB_CLOSE => 'close'
#define G_BOUND_MEDIAN 1 // k*R < GB_MEDIAN => 'median'

enum iter { // iterative methods
	IT_BCGS2,    // Enhanced Bi-Conjugate Gradient Stabilized (2)
	IT_BICG_CS,  // Bi-Conjugate Gradient for Complex-Symmetric matrices
	IT_BICGSTAB, // Bi-Conjugate Gradient Stabilized
	IT_CGNR,     // Conjugate Gradient for Normalized equations minimizing Residual norm
	IT_CSYM,     // Algorithm CSYM
	IT_QMR_CS,   // Quasi-minimal residual for Complex-Symmetric matrices
	IT_QMR_CS_2  // 2-term QMR (better roundoff properties)
	/* TO ADD NEW ITERATIVE SOLVER
	 * add an identifier starting with 'IT_' and a descriptive comment to this list in the alphabetical order.
	 */
};

enum Eftype { // type of E field calculation
	CE_NORMAL, // normal
	CE_PARPER  // use symmetry to calculate both incident polarizations from one calculation of internal fields
};

enum incpol {
	INCPOL_Y, // y-polarization, assumed to be used first
	INCPOL_X  // x-polarization
};
// path and size of tables
#define TAB_PATH     "tables/"
#define TAB_FNAME(a) "t"#a"f.dat" // a is a number, e.g. TAB_FNAME(2) -> "t2f.dat"
#define TAB_SIZE     142
#define TAB_RMAX     10

enum beam { // beam types
	B_BARTON5, // 5th order description of the Gaussian beam
	B_DAVIS3,  // 3rd order description of the Gaussian beam
	B_DIPOLE,  // field of a point dipole
	B_LMINUS,  // 1st order description of the Gaussian beam
	B_PLANE,   // infinite plane wave
	B_READ     // read from file
	/* TO ADD NEW BEAM
	 * add an identifier starting with 'B_' and a descriptive comment to this list in alphabetical order.
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
	CHP_ALWAYS   /* save checkpoint if either simulation is finished or time elapsed and calculate all scattering
	                quantities */
};

enum init_field { // how to calculate initial field to be used in the iterative solver
	IF_AUTO, // automatically choose from ZERO or INC (based on lower residual value)
	IF_ZERO, // zero
	IF_INC,  // equal to incident field
	IF_READ, // read from file
	IF_WKB   // from WKB approximation (incident field corrected for phase shift in the particle)
};

// return values for functions
#define CHP_EXIT -2 // exit after saving checkpoint

// default values; other are specified in InitVariables (param.c)
#define DEF_GRID       (16*jagged)
#define MIN_AUTO_GRID  16 // minimum grid, when set from default dpl

// numbers less than this value (compared to unity) are considered to be zero (approximately 10*DBL_EPSILON)
#define ROUND_ERR 1E-15
#define SQRT_RND_ERR 3E-8 // sqrt(ROUND_ERR)

// output and input file and directory names (can only be changed at compile time)
#define F_EXPCOUNT      "ExpCount"
#define F_EXPCOUNT_LCK  F_EXPCOUNT ".lck"
#define F_CS            "CrossSec"
#define F_FRP           "RadForce"
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
	// files for amplitude matrix
#define F_AMPL          "ampl"
#define F_AMPL_SG       "ampl_scatgrid"
	// temporary files; used in printf with ringid as argument
#define F_FRP_TMP       "rf%d.tmp"
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

/* number of components of D. Really, it can't be easily changed, but using constant instead of 6 adds more
 * understanding for reader
 */
#define NDCOMP 6

// shape formats; numbers should be nonnegative
enum shform {
	SF_DDSCAT6,  // DDSCAT 6 format (FRMFIL), produced by calltarget
	SF_DDSCAT7,  // DDSCAT 7 format (FRMFIL), produced by calltarget
	SF_TEXT,     // ADDA text format for one-domain particles
	SF_TEXT_EXT  // ADDA text format for multi-domain particles
	/* TO ADD NEW FORMAT OF SHAPE FILE
	 * add an identifier starting with 'SF_' and a descriptive comment to this list in alphabetical order.
	 */
};

#define POSIT __FILE__,__LINE__ // position of the error in source code

enum enwho { // who is calling
	ALL, // each processor may report this error
	ONE  // only root processor reports an error
};

// derived; for simplicity
#define ALL_POS ALL,POSIT
#define ONE_POS ONE,POSIT
#define ALL_POS_FUNC ALL_POS,__func__
#define ONE_POS_FUNC ONE_POS,__func__

enum ec { // error codes
	EC_OK,    // no error, corresponds to zero exit code
	EC_ERROR, // error
	EC_WARN,  // warning
	EC_INFO   // slight warning, that does not interfere at all with normal execution
};

#endif // __const_h
