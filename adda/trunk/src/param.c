/* FILE: param.c
 * AUTH: Maxim Yurkin
 * DESCR: Initialization, parsing and handling of input parameters.
 *        Also printout general information. Contains file locking routines.
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <limits.h>
#include "os.h"
#include "io.h"
#include "const.h"
#include "comm.h"
#include "vars.h"
#include "crosssec.h"
#include "fft.h"
#include "param.h"
#include "cmplx.h"
#include "function.h"
#include "parbas.h"

/* definitions for file locking */
#ifdef USE_LOCK
# ifdef WINDOWS
#  define FILEHANDLE HANDLE
# elif defined(POSIX)
#  include <unistd.h>
#  include <fcntl.h>
#  ifdef LOCK_FOR_NFS
#   include <errno.h>    /* for error handling of fcntl call */
#  endif
#  define FILEHANDLE int
# else
#  error *** Unknown operation system. Creation of lock files is not supported. ***
# endif
# define LOCK_WAIT 1              /* in seconds */
# define MAX_LOCK_WAIT_CYCLES 60
#else
# define FILEHANDLE int
#endif

/* GLOBAL VARIABLES */

opt_index opt;   /* main option index */

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in crosssec.c */
extern const char avg_string[];
/* defined and initialized in GenerateB.c */
extern const char beam_descr[];
/* defined and initialized in make_particle.c */
extern const int volcor_used;
extern const char sh_form_str[];
extern const int gr_N;
extern const double gr_vf_real;
extern double mat_count[];

/* used in CalculateE.c */
int store_int_field; /* save full internal fields to text file */
int store_dip_pol;   /* save dipole polarizations to text file */
int store_beam;      /* save incident beam to file */
int store_scat_grid; /* Store the scattered field for grid of angles */
int calc_Cext;       /* Calculate the extinction cross-section - allways do */
int calc_Cabs;       /* Calculate the absorption cross-section - allways do */
int calc_Csca;       /* Calculate the scattering cross-section by integration */
int calc_vec;        /* Calculate the unnormalized asymmetry-parameter */
int calc_asym;       /* Calculate the asymmetry-parameter */
int calc_mat_force;  /* Calculate the scattering force by matrix-evaluation */
int store_force;     /* Write radiation pressure per dipole to file */
int phi_int_type;    /* type of phi integration (each bit determines
                          whether to calculate with different multipliers) */
/* used in calculator.c */
int avg_inc_pol;                 /* whether to average CC over incident polarization */
char alldir_parms[MAX_FNAME];    /* name of file with alldir parameters */
char scat_grid_parms[MAX_FNAME]; /* name of file with parameters of scattering grid */
/* used in crosssec.c */
double prop_0[3];  /* initial incident direction (in laboratory reference frame) */
double incPolX_0[3],incPolY_0[3]; /* initial incident polarizations (in lab RF)*/
int ScatRelation;                 /* type of formulae for scattering quantities */
/* used in GenerateB.c */
int beam_Npars;
double beam_pars[MAX_N_BEAM_PARMS]; /* beam parameters */
/* used in io.c */
char logname[MAX_FNAME]="";      /* name of logfile */
/* used in iterative.c */
double eps;                      /* relative error to reach */
/* used in make_particle.c */
int shape;                       /* particle shape definition */
int sh_Npars;                    /* number of shape parameters */
double sh_pars[MAX_N_SH_PARMS];  /* storage for shape parameters */
int symmetry_enforced;           /* enforce use of all symmetries; suppresses NoSymmetry */
double sizeX;                    /* size of particle along x-axis */
double dpl;                      /* number of dipoles per lambda (wavelength) */
double lambda;                   /* incident wavelength (in vacuum) */
int jagged;                      /* size of big dipoles, used to construct a particle */
char aggregate_file[MAX_FNAME];  /* name of aggregate file */
char save_geom_fname[MAX_FNAME]; /* geometry file name to save dipole configuration */
char shapename[MAX_LINE];        /* name of the shape used */
int volcor;                      /* whether to use volume correction */
int save_geom;                   /* whether to save dipole configuration in .geom file */
opt_index opt_sh;                /* option index of shape option used */
double gr_vf;                    /* granules volume fraction */
double gr_d;                     /* granules diameter */
int gr_mat;                      /* domain number to granulate */

/* LOCAL VARIABLES */

static char run_name[MAX_WORD];   /* first part of the dir name ('run' or 'test') */
static char avg_parms[MAX_FNAME]; /* name of file with orientation averaging parameters */
static char *exename;             /* name of executable (adda or adda.exe) */
static int Nmat_given;            /* number of refractive indices given in the command line */
  /* structure definitions */
struct subopt_struct {
  const char *name;   /* name of option */
  const char *usage;  /* how to use (argument list) */
  const char *help;   /* help string */
  const int narg;     /* possible number of argumetns ; UNDEF -> should not be checked */
  const int type;     /* type of suboption */
};
struct opt_struct {
  const char *name;         /* name of option */
  void (*func)(int Narg,char **argv);  /* pointer to a function, that parse this parameter */
  int used;                 /* flag to indicate, if the option was allready used */
  const char *usage;        /* how to use (argument list) */
  const char *help;         /* help string */
  const int narg;           /* possible number of argumetns ; UNDEF -> should not be checked */
  const struct subopt_struct *sub;  /* suboptions */
};
  /* const string for usage of ADDA */
static const char exeusage[]="[-<opt1> [<args1>] [-<opt2> [<args2>]...]]";
  /* initializations of suboptions; should be 'NULL terminated'
     each row contains: suboption name, usage string, help string, number of arguments
       (UNDEF = not checked automatically, identifier (number)  */
static const struct subopt_struct beam_opt[]={
  {"plane","","Infinite plane wave",0,B_PLANE},
  {"lminus","<width> [<x> <y> <z>]",
     "Simplest approximation of the Gaussian beam. The beam width is obligatory\n"\
     "and x, y, z coordinates of the center of the beam are optional parameters (all in um).\n"\
     "By default beam center coincides with the center of the computational box.",
     UNDEF,B_LMINUS},
  {"davis3","<width> [<x> <y> <z>]",
     "3rd order approximation of the Gaussian beam (by Davis). The beam width is obligatory\n"\
     "and x, y, z coordinates of the center of the beam are optional parameters (all in um).\n"\
     "By default beam center coincides with the center of the computational box.",
     UNDEF,B_DAVIS3},
  {"barton5","<width> [<x> <y> <z>]",
     "5th order approximation of the Gaussian beam (by Barton). The beam width is obligatory\n"\
     "and x, y, z coordinates of the center of the beam are optional parameters (all in um).\n"\
     "By default beam center coincides with the center of the computational box.\n"\
     "This option is recommended for the description of the Gaussian beam.",
     UNDEF,B_BARTON5},
  {NULL,NULL,NULL,0,0}
};
static const struct subopt_struct shape_opt[]={
  {"box","","Homogenous cube (edges along the axes)",0,SH_BOX},
  {"coated","<d_in/d> [<x/d> <y/d> <z/d>]",
     "Sphere with a spherical inclusion; outer sphere has a diameter d (first domain).\n"\
     "The included sphere has a diameter d_in (optional position of the center: x,y,z).",
     UNDEF,SH_COATED},
  {"cylinder","<h/d>",
     "Homogenous cylinder with height (length) h and diameter d (its axis of symmetry\n"\
     "coincides with the z-axis).",1,SH_CYLINDER},
  {"ellipsoid","<y/x> <z/x>","Homogenous general ellipsoid with semi-axes x,y,z",2,SH_ELLIPSOID},
  {"line","","Line along the x-axis with the width of one dipole",0,SH_LINE},
  {"rbc","<h/d> <b/d> <c/d>",
     "Red Blood Cell, an axisymmetric (over z-axis) biconcave homogenous particle,\n"\
     "which is characterized by diameter d, maximum and minimum width h, b, and\n"\
     "diameter at the position of the maximum width c.",3,SH_RBC},
  {"read","<filename>","Read a particle geometry from file <filename>",1,SH_READ},
  {"sphere","","Homogenous sphere",0,SH_SPHERE},
  {"spherebox","<d_sph/Dx>",
     "Sphere (diameter d_sph) in a cube (size Dx, first domain)",1,SH_SPHEREBOX},
/* TO ADD NEW SHAPE
   add a row here, before null-terminating element. It contains:
   shape name (used in command line), usage string (what command line parameters can be used
   for this shape), help string (shown when -h option is used), possible number of parameters
   (use UNDEF if shape can accept different number of parameters, then check it explicitly in
   function InitShape), shape identifier (constant defined in const.h). Number of parameters
   should not be greater than MAX_N_SH_PARMS (defined in const.h). It is recommended to use
   dimensionless shape parameters, e.g. aspect ratios. */

  {NULL,NULL,NULL,0,0}
};

/* EXTERNAL FUNCTIONS */

/* GenerateB.c */
void InitBeam(void);

/*========================================================================*/
   /* declarations of parsing functions; definitions are given below.
      defines are for conciseness */
#define PARSE_NAME(a) parse_##a
#define PARSE_FUNC(a) void PARSE_NAME(a)(int Narg,char **argv)
#define PAR(a) #a,PARSE_NAME(a),FALSE
PARSE_FUNC(alldir_inp);
PARSE_FUNC(anisotr);
PARSE_FUNC(asym);
PARSE_FUNC(beam);
PARSE_FUNC(chp_dir);
PARSE_FUNC(chp_load);
PARSE_FUNC(chp_type);
PARSE_FUNC(chpoint);
PARSE_FUNC(Cpr_mat);
PARSE_FUNC(Csca);
PARSE_FUNC(dir);
PARSE_FUNC(dpl);
PARSE_FUNC(eps);
PARSE_FUNC(granul);
PARSE_FUNC(grid);
PARSE_FUNC(h) ATT_NORETURN;
PARSE_FUNC(int);
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
PARSE_FUNC(prognose);
PARSE_FUNC(prop);
PARSE_FUNC(save_geom);
PARSE_FUNC(scat);
PARSE_FUNC(scat_grid_inp);
PARSE_FUNC(shape);
PARSE_FUNC(size);
PARSE_FUNC(store_beam);
PARSE_FUNC(store_dip_pol);
PARSE_FUNC(store_force);
PARSE_FUNC(store_int_field);
PARSE_FUNC(store_scat_grid);
PARSE_FUNC(sym);
PARSE_FUNC(test);
PARSE_FUNC(V) ATT_NORETURN;
PARSE_FUNC(vec);
PARSE_FUNC(yz);
   /* initialization of options, their usage and help;
      each row contains: PAR(option name),usage string, help string, number of arguments
      (UNDEF = not checked automatically),pointer to suboption (if exist) */
static struct opt_struct options[]={
  {PAR(alldir_inp),"<filename>",
     "Specifies a file with parameters of the grid of scattering angles\n"\
     "for calculating integral scattering quantities.\n"\
     "Default: " FD_ALLDIR_PARMS,1,NULL},
  {PAR(anisotr),"",
     "Specifies that refractive index is anisotropic (its tensor is limited to be diagonal\n"\
     "in particle reference frame). '-m' then accepts 6 arguments per each domain.\n"\
     "Can not be used with CLDR polarizability and all SO formulations.",0,NULL},
  {PAR(asym),"","Calculate the asymmetry vector. Implies '-Csca' and '-vec'",0,NULL},
  {PAR(beam),"<type> [<arg1>...]",
     "Sets a type of the incident beam. Four other float arguments must be specified\n"\
     "for all beam types except 'plane'. These are the width and x, y, z coordinates\n"\
     "of the center of the beam respectively (all in um).\n"\
     "Default: plane",UNDEF,beam_opt},
  {PAR(chp_dir),"<dirname>",
     "Sets directory for the checkpoint (both for saving and loading).\n"\
     "Default: " FD_CHP_DIR,1,NULL},
  {PAR(chp_load),"","Restart a simulation from a checkpoint",0,NULL},
  {PAR(chp_type),"{normal|regular|always}",
     "Sets type of the checkpoint. All types, except 'always', require '-chpoint'.\n"\
     "Default: normal",1,NULL},
  {PAR(chpoint),"<time>",
     "Specifies the time for checkpoints in format '#d#h#m#s'. All fields are optional,\n"\
     "numbers are integers, 's' can be omitted, the format is not case sensitive.\n"\
     "Examples: 12h30M, 1D10s, 3600",1,NULL},
  {PAR(Cpr_mat),"","Calculate the total radiation force",0,NULL},
  {PAR(Csca),"","Calculate scattering cross section (by integrating the scattered field)",0,NULL},
  {PAR(dir),"<dirname>",
     "Sets directory for output files.\n"\
     "Default: constructed automatically",1,NULL},
  {PAR(dpl),"<arg>",
     "Sets parameter 'dipoles per lambda', float.\n"\
     "Default: 10|m|, where 'm' is the maximum (by absolute value) refractive index\n"\
     "         specified by the '-m' option.",1,NULL},
  {PAR(eps),"<arg>",
     "Specifies the stopping criterion for the iterative solver by setting the\n"\
     "relative norm of the residual 'epsilon' to reach. <arg> is an exponent\n"\
     "of base 10 (float), i.e. epsilon=10^(-<arg>).\n"\
     "Default: 5 (epsilon=1E-5)",1,NULL},
  {PAR(granul),"<vol_frac> <diam> [<dom_number>]",
     "Specifies that one particle domain should be randomly filled with spherical granules\n"\
     "with specified diameter <diam> and volume fraction <vol_frac>. Domain number to fill\n"\
     "is given by the last optional argument. Algorithm may fail for volume fractions > 30-50%.\n"\
     "Default <dom_number>: 1",UNDEF,NULL},
  {PAR(grid),"<nx> [<ny> <nz>]",
     "Sets dimensions of the computation grid. Arguments should be even integers.\n"\
     "In most cases <ny> and <nz> can be omitted (they are automatically determined\n"\
     "by <nx> based on the proportions of the scatterer). This command line option\n"\
     "is not relevant when particle geometry is read from a file ('-shape read').\n"\
     "If '-jagged' option is used the grid dimension is effectively multiplied\n"\
     "by the specified number.\n"\
     "Default: 16 (if  size is not specified) or defined by\n"\
     "         '-size', '-lambda', and '-dpl'.",UNDEF,NULL},
  {PAR(h),"[<opt> [<subopt>]]",
     "Shows help. If used without arguments, ADDA shows a list of all available\n"\
     "command line options. If first argument is specified, help on specific command\n"\
     "line option <opt> is shown (only the name of the option should be given\n"\
     "without preceding dash). For some options (e.g. '-beam' or '-shape') specific\n"\
     "help on a particular suboption <subopt> may be shown.\n"\
     "Example: shape coated",UNDEF,NULL},
  {PAR(int),"{poi|fcd|fcd_st|so}",
     "Sets prescription to calculate interaction term. 'so' is under development\n"\
     "and incompatible with '-anisotr'. 'fcd' requires dpl to be larger than 2.\n"\
     "Default: poi",1,NULL},
  {PAR(iter),"{cgnr|bicg|bicgstab|qmr}",
     "Sets the iterative solver.\n"\
     "Default: qmr",1,NULL},
  {PAR(jagged),"<arg>",
     "Sets a size of a big dipole in units of small dipoles, integer. It is used\n"\
     "to improve the discretization of the particle without changing the shape.\n"\
     "Default: 1",1,NULL},
  {PAR(lambda),"<arg>",
     "Sets incident wavelength in um, float.\n"\
     "Default: 2*pi",1,NULL},
  {PAR(m),"{<m1Re> <m1Im> [...]|<m1xxRe> <m1xxIm> <m1yyRe> <m1yyIm> <m1zzRe> <m1zzIm> [...]}",
     "Sets refractive indices, float. Each pair of arguments specifies real and\n"\
     "imaginary part of the refractive index of one of the domains. If '-anisotr' is\n"\
     "specified, three refractive indices correspond to one domain (diagonal elements of\n"
     "refractive index tensor in particle reference frame). Maximum number of different\n"\
     "refractive indices is defined at compilation time by the parameter MAX_NMAT.\n"\
     "None of the refractive indices can be equal to 1+0i.\n"\
     "in file const.h (by default, 15).\n"\
     "Default: 1.5 0",UNDEF,NULL},
  {PAR(maxiter),"<arg>",
     "Sets the maximum number of iterations of the iterative solver, integer.\n"\
     "Default: very large, not realistic value",1,NULL},
  {PAR(no_reduced_fft),"",
     "Do not use symmetry of the interaction matrix to reduce the storage space\n"\
     "for the Fourier-transformed matrix.",0,NULL},
  {PAR(no_vol_cor),"",
     "Do not use 'dpl correction', which ensures (if used) that the volume of\n"\
     "the dipole representation of the particle is exactly correct.",0,NULL},
  {PAR(ntheta),"<arg>",
     "Sets the number of intervals, into which the range of scattering angles [0,180]\n"\
     "(degrees) is equally divided, integer. This is used for scattering angles in\n"\
     "yz-plane. If particle is not symmetric and orientation averaging is not used,\n"\
     "the range is extended to 360 degrees (with the same length of elementary interval,\n"\
     "i.e. number of intervals is doubled).\n"\
     "Default: from 90 to 720 depending on the size of the computational grid.",1,NULL},
  {PAR(opt),"{speed|mem}",
     "Sets whether ADDA should optimize itself for maximum speed or for minimum memory usage.\n"\
     "Default: speed",1,NULL},
  {PAR(orient),"{<alpha> <beta> <gamma>|avg [<filename>]}",
     "Either sets an orientation of the particle by three Euler angles 'alpha',\n"\
     "'beta','gamma' (in degrees) or specifies that orientation averaging should be\n"\
     "performed. <filename> sets a file with parameters for orientation averaging. Here\n"\
     "zyz-notation (or y-convention) is used for the Euler angles.\n"\
     "Default orientation: 0 0 0\n"\
     "Default <filename>: " FD_AVG_PARMS,UNDEF,NULL},
  {PAR(phi_integr),"<arg>",
     "Turns on and specifies the type of Mueller matrix integration over azimuthal\n"\
     "angle 'phi'. <arg> is an integer from 1 to 31, each bit of which, from lowest\n"\
     "to highest, indicates whether the integration should be performed with\n"\
     "multipliers 1, cos(2*phi), sin(2*phi), cos(4*phi), and sin(4*phi)\n"\
     "respectively.\n"\
     "Examples: 1 (one integration with no multipliers),\n"\
     "          6 (two integration with cos(2*phi) and sin(2*phi) multipliers).",1,NULL},
  {PAR(pol),"{cm|rrc|ldr [avgpol]|cldr|so|fcd}",
     "Type of polarization prescription. An optional flag 'avg' can be added for LDR\n"\
     "- it specifies that LDR polarizability should be averaged over incident\n"\
     "polarizations. 'so' is under development. 'cldr' and 'so' are incompatible\n"\
     "with '-anisotr'. 'fcd' requires dpl to be larger than 2.\n"\
     "Default: ldr (without averaging).",UNDEF,NULL},
  {PAR(prognose),"",
     "Do not actually perform simulation (not even memory allocation) but only\n"\
     "estimate the required RAM. Implies '-test'.",0,NULL},
  {PAR(prop),"<x> <y> <z>",
     "Sets propagation direction of incident radiation, float. Normalization\n"\
     "(to the unity vector) is performed automatically.\n"\
     "Default: 0 0 1",3,NULL},
  {PAR(save_geom),"[<filename>]",
     "Saves dipole configuration to a file <filename> (a path relative to the\n"\
     "output directory). Can be used with '-prognose'.\n"\
     "Default: <type>.geom (<type> is a first argument to the '-shape' option; '_gran' is added\n"\
     "                      if '-granul' option is used).",UNDEF,NULL},
  {PAR(scat),"{dr|so}",
     "Sets prescription to calculate scattering quantities.\n"\
     "'so' is under development and incompatible with '-anisotr'.\n"\
     "Default: dr",1,NULL},
  {PAR(scat_grid_inp),"<filename>",
     "Specifies a file with parameters of the grid of scattering angles for\n"\
     "calculating Mueller matrix (possibly integrated over 'phi').\n"\
     "Default: " FD_SCAT_PARMS,1,NULL},
  {PAR(shape),"<type> [<arg1>...]",
     "Sets shape of the particle, either predefined or 'read' from file.\n"\
     "All the parameters of predefined shapes are floats.\n"\
     "Default: sphere",UNDEF,shape_opt},
  {PAR(size),"<arg>",
     "Sets the size of the computational grid along the x-axis in um, float.\n"\
     "Default: determined by the values of '-grid', '-dpl', and '-lambda'.",1,NULL},
  {PAR(store_beam),"","Save incident beam to a file",0,NULL},
  {PAR(store_dip_pol),"","Save dipole polarizations to a file",0,NULL},
  {PAR(store_force),"","Calculate the radiation force on each dipole. Requires '-Cpr_mat'",0,NULL},
  {PAR(store_int_field),"","Save internal fields to a file",0,NULL},
  {PAR(store_scat_grid),"",
     "Calculate Mueller matrix for a grid of scattering angles and save it to a file.",0,NULL},
  {PAR(sym),"{no|enf}",
     "Do not take into account ('no') or enforce ('enf') all particle symmetries",1,NULL},
  {PAR(test),"","Begin name of the output directory with 'test' instead of 'run'",0,NULL},
  {PAR(V),"",
     "Show ADDA version, compiler used to build this executable,\n"\
     "and copyright information",0,NULL},
  {PAR(vec),"","Calculate the not-normalized asymmetry vector",0,NULL},
  {PAR(yz),"",
     "Calculate the Mueller matrix in yz-plane even if it is calculated for a\n"\
     "scattering grid. If the latter option is not enabled, scattering in yz-plane\n"\
     "is always calculated.",0,NULL}
};
      /* auxiliary functions */
/*============================================================*/

static const char *OptionName(void)
   /* produces full option name for error messages */
{
   static char buf[MAX_LINE];

   if (opt.l2==UNDEF) return options[opt.l1].name;
   else {
     sprintf(buf,"%s %s",options[opt.l1].name,options[opt.l1].sub[opt.l2].name);
     return buf;
   }
}

/*============================================================*/

void PrintErrorHelp(const char *fmt, ... )
   /* print anything to stderr (on root processor), then help on the arguments used, and stop;
      assumes that all processors call it */
{
  va_list args;
  const char *optname,*use;

  if (ringid==ROOT) {
    /* produce error message */
    va_start(args,fmt);
    fprintf(stderr,"ERROR: ");
    vfprintf(stderr,fmt,args);
    fprintf(stderr,"\n");
    va_end(args);
    /* add help message */
    if (opt.l1==UNDEF)     /* no option is found */
      fprintf(stderr,"Usage: %s %s\n"\
                     "Type '%s -h' for help\n",exename,exeusage,exename);
    else {  /* at least option is found */
      if (opt.l2==UNDEF) use=options[opt.l1].usage;
      else use=options[opt.l1].sub[opt.l2].usage;
      optname=OptionName();
      fprintf(stderr,"Usage: -%s %s\n"\
                     "Type '%s -h %s' for details\n",optname,use,exename,optname);
    }
    fflush(stderr);
  }
  /* wait for root to generate an error message */
  Synchronize();
  Stop(1);
}

/*============================================================*/

static void NargError(const int Narg,const char *expec)
      /* Print error of illegal number of arguments to an option (suboption);
         and display correct usage information */
{
  char buf[MAX_WORD]; /* not to allocate memory if needed */

  if (expec==NULL) {
    if (opt.l2==UNDEF) sprintf(buf,"%d",options[opt.l1].narg);
    else sprintf(buf,"%d",options[opt.l1].sub[opt.l2].narg);
    expec=buf;
  }
  PrintErrorHelp("Illegal number of arguments (%d) to '-%s' option (%s expected)",
                 Narg,OptionName(),expec);
}

/*============================================================*/
                /* following two functions are interfaces to NargError */
INLINE void TestNarg(const int Narg)
     /* check if Narg given to an option is correct */
{
  if (options[opt.l1].narg!=UNDEF && Narg!=options[opt.l1].narg)
    NargError(Narg,NULL);
}

/*============================================================*/

INLINE void TestNarg_sub(const int Narg)
     /* check if Narg given to a suboption is correct */
{
  if (options[opt.l1].sub[opt.l2].narg!=UNDEF && Narg!=options[opt.l1].sub[opt.l2].narg)
    NargError(Narg,NULL);
}

/*============================================================*/

static void NotSupported(const char *type,const char *given)
     /* print error message that "type 'given' is not supported"
        type should start with a capital letter */
{
  PrintErrorHelp("%s '%s' is not supported",type,given);
}

/*============================================================*/

INLINE void TestStrLength(const char *str,const unsigned int size)
    /* check if string fits in buffer of size 'size', otherwise produces error message
       'opt' is command line option that checks its argument */
{
  if (strlen(str)>=size)
    PrintErrorHelp("Too long argument to '-%s' option (only %ud chars allowed).\n"\
                   "If you really need it you may increase MAX_DIRNAME in const.h and recompile",
                   OptionName(),size-1);
}

/*============================================================*/

INLINE void ScanfDoubleError(const char *str,double *res)
    /* scanf an option argument and checks for errors */
{
  if (sscanf(str,"%lf",res)!=1)
    PrintErrorHelp("Non-numeric argument (%s) is given to the option '-%s'",str,OptionName());
}

/*============================================================*/

INLINE void ScanfIntError(const char *str,int *res)
    /* scanf an option argument and checks for errors */
{
  double tmp;

  if (sscanf(str,"%lf",&tmp)!=1)
    PrintErrorHelp("Non-numeric argument (%s) is given to the option '-%s'",str,OptionName());
  if (tmp <INT_MIN || tmp>INT_MAX) PrintErrorHelp(
    "Argumenent value (%s) of the option '-%s' is out of integer bounds",str,OptionName());
  if (sscanf(str,"%d",res)!=1)
    PrintErrorHelp("Error reading argument (%s) of the option '-%s'",str,OptionName());
}

/*============================================================*/

INLINE int IsOption(const char *str)
   /* checks if string is an option. First should be '-' and then letter (any case);
      it enables use of negative numbers as subparameters */
{
  /* conversion to int is needed to remove warnings caused by the fact
     that str[1] is _signed_ char */
  return (str[0]=='-' && isalpha((int)(str[1])));
}
/*============================================================*/

static int TimeField(const char c)
   /* analyze one time multiplier */
{
  if (c=='d' || c=='D') return 86400;
  else if (c=='h' || c=='H') return 3600;
  else if (c=='m' || c=='M') return 60;
  else if (c=='s' || c=='S' || c==0) return 1;
  else PrintErrorHelp("Illegal time format specifier (%c)",c);
  /* never reached */
  return 0;
}

/*============================================================*/

static int ScanTime(const char *str)
   /* scans time in seconds from a string "%d[d,D[%d]][h,H[%d]][m,M[%d]][s,S] */
{
#define TIME_N_TYPES 4   /* not so easy to change */
  int tim,t[TIME_N_TYPES],n,i;
  char c[TIME_N_TYPES];

  for (i=0;i<TIME_N_TYPES;i++) c[i]=0;
  n=sscanf(str,"%d%c%d%c%d%c%d%c",t,c,t+1,c+1,t+2,c+2,t+3,c+3);
  if (n<1) PrintErrorHelp("Wrong time format '%s'",str);
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

/*============================================================*/

static void PrintTime(char *s,const time_t *time_ptr)
{
   struct tm *t;

   t=gmtime(time_ptr);
   s[0]=0; /* initialize string */
   if (t->tm_yday>0) sprintf(s,"%dd ",t->tm_yday);
   if (t->tm_hour>0) sprintf(s+strlen(s),"%dh ",t->tm_hour);
   if (t->tm_min>0) sprintf(s+strlen(s),"%dm ",t->tm_min);
   if (t->tm_sec>0) sprintf(s+strlen(s),"%ds ",t->tm_sec);
}

/*========================================================================*/
   /* parsing functions definitions*/
PARSE_FUNC(alldir_inp)
{
  TestStrLength(argv[1],MAX_FNAME);
  strcpy(alldir_parms,argv[1]);
}
PARSE_FUNC(anisotr)
{
  anisotropy = TRUE;
  Ncomp=3;
}
PARSE_FUNC(asym)
{
  calc_asym = TRUE;
  calc_vec = TRUE;
  calc_Csca = TRUE;
}
PARSE_FUNC(beam)
{
  int i,j,found;

  Narg--;
  found=FALSE;
  i=-1;
  while (beam_opt[++i].name!=NULL) if (strcmp(argv[1],beam_opt[i].name)==0) {
    /* set suboption and beamtype */
    opt.l2=i;
    beamtype=beam_opt[i].type;
    beam_Npars=Narg;
    /* check number of arguments */
    TestNarg_sub(Narg);
    if (beamtype!=B_PLANE) {
      if (Narg!=1 && Narg!=4) NargError(Narg,"1 or 4");
    }
    /* parse and check consistency */
    for (j=0;j<Narg;j++) ScanfDoubleError(argv[j+2],beam_pars+j);
    if (Narg>0) TestPositive(beam_pars[0],"beam width");
    /* stop search */
    found=TRUE;
    break;
  }
  if(!found) NotSupported("Beam type",argv[1]);
}
PARSE_FUNC(chp_dir)
{
  TestStrLength(argv[1],MAX_DIRNAME);
  strcpy(chp_dir,argv[1]);
}
PARSE_FUNC(chp_load)
{
  load_chpoint = TRUE;
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
PARSE_FUNC(Cpr_mat)
{
  calc_mat_force = TRUE;
}
PARSE_FUNC(Csca)
{
  calc_Csca = TRUE;
}
PARSE_FUNC(dir)
{
  TestStrLength(argv[1],MAX_DIRNAME);
  strcpy(directory,argv[1]);
}
PARSE_FUNC(dpl)
{
  ScanfDoubleError(argv[1],&dpl);
  TestPositive(dpl,"dpl");
}
PARSE_FUNC(eps)
{
  double tmp;

  ScanfDoubleError(argv[1],&tmp);
  TestPositive(tmp,"eps exponent");
  eps=pow(10,-tmp);
}
PARSE_FUNC(granul)
{
  if (Narg!=2 && Narg!=3) NargError(Narg,"2 or 3");
  ScanfDoubleError(argv[1],&gr_vf);
  TestRange(gr_vf,"volume fraction",0,PI_OVER_SIX);
  ScanfDoubleError(argv[2],&gr_d);
  TestPositive(gr_d,"diameter");
  if (Narg==3) {
    ScanfIntError(argv[3],&gr_mat);
    TestPositive_i(gr_mat,"domain number");
  }
  else gr_mat=1;
  gr_mat--;  /* converted to usual indexing starting from 0 */
  sh_granul=TRUE;
}
PARSE_FUNC(grid)
{
  if (Narg!=1 && Narg!=3) NargError(Narg,"1 or 3");
  ScanfIntError(argv[1],&boxX);  /* boxes are further multiplied by jagged if needed */
  TestRange_i(boxX,"gridX",1,BOX_MAX);
  if (Narg==3) {
    ScanfIntError(argv[2],&boxY);
    TestRange_i(boxY,"gridY",1,BOX_MAX);
    ScanfIntError(argv[3],&boxZ);
    TestRange_i(boxY,"gridY",1,BOX_MAX);
  }
}
PARSE_FUNC(h)
{
  int i,j,found;

  if (Narg>2) NargError(Narg,"not more than 2");
  /* do all output on root processor */
  if (ringid==ROOT) {
    found=FALSE;
    if (Narg>=1) {
      for(i=0;i<LENGTH(options);i++) if (strcmp(argv[1],options[i].name)==0) {
        if (Narg==2) {
          if (options[i].sub==NULL)
            printf("No specific help is available for suboptions of this option\n\n");
          else {
            j=-1;
            while (options[i].sub[++j].name!=NULL) if (strcmp(argv[2],options[i].sub[j].name)==0) {
              printf("  -%s %s %s\n%s\n",options[i].name,options[i].sub[j].name,
                     options[i].sub[j].usage,options[i].sub[j].help);
              found=TRUE;
              break;
            }
            if (!found) printf("Unknown suboption '%s'\n\n",argv[2]);
          }
        }
        if (!found) {
          printf("  -%s %s\n%s\n",options[i].name,options[i].usage,options[i].help);
          if (options[i].sub!=NULL) {
            printf("Available suboptions:\n");
            j=-1;
            while (options[i].sub[++j].name!=NULL)
              printf("  %s %s\n",options[i].sub[j].name,options[i].sub[j].usage);
            printf("Type '%s -h %s <subopt>' for details\n",exename,options[i].name);
          }
        }
        found=TRUE;
        break;
      }
      if (!found) printf("Unknown option '%s'\n\n",argv[1]);
    }
    if (!found) {
      printf("Usage: '%s %s'\n"\
             "Available options:\n",exename,exeusage);
      for (i=0;i<LENGTH(options);i++) printf("  -%s %s\n",options[i].name,options[i].usage);
      printf("Type '%s -h <opt>' for details\n",exename);
    }
  }
  /* exit */
  Stop(0);
}
PARSE_FUNC(int)
{
  if (strcmp(argv[1],"poi")==0) IntRelation=G_POINT_DIP;
  else if (strcmp(argv[1],"fcd")==0) IntRelation=G_FCD;
  else if (strcmp(argv[1],"fcd_st")==0) IntRelation=G_FCD_ST;
  else if (strcmp(argv[1],"so")==0) IntRelation=G_SO;
  else NotSupported("Interaction term prescription",argv[1]);
}
PARSE_FUNC(iter)
{
  if (strcmp(argv[1],"cgnr")==0) IterMethod=IT_CGNR;
  else if (strcmp(argv[1],"bicgstab")==0) IterMethod=IT_BICGSTAB;
  else if (strcmp(argv[1],"bicg")==0) IterMethod=IT_BICG_CS;
  else if (strcmp(argv[1],"qmr")==0) IterMethod=IT_QMR_CS;
  else NotSupported("Iterative method",argv[1]);
}
PARSE_FUNC(jagged)
{
  ScanfIntError(argv[1],&jagged);
  TestRange_i(jagged,"jagged",1,BOX_MAX);
}
PARSE_FUNC(lambda)
{
  ScanfDoubleError(argv[1],&lambda);
  TestPositive(lambda,"wavelength");
}
PARSE_FUNC(m)
{
  int i;

  if (Narg%2!=0 || Narg==0) NargError(Narg,"even");
  Nmat=Nmat_given=Narg/2;
  if (Nmat>MAX_NMAT)
    PrintErrorHelp("Too many materials (%d), maximum %d are supported.\n"\
                   "You may increase parameter MAX_NMAT in const.h and recompile",Nmat,MAX_NMAT);
  for (i=0;i<Nmat;i++) {
    ScanfDoubleError(argv[2*i+1],&ref_index[i][RE]);
    ScanfDoubleError(argv[2*i+2],&ref_index[i][IM]);
    if (ref_index[i][RE]==1 && ref_index[i][IM]==0)
      PrintErrorHelp("Given refractive index #%d is that of vacuum, which is not supported.\n"\
                     "Consider using, for instance, 1.0001 instead.",i+1);
  }
}
PARSE_FUNC(maxiter)
{
  ScanfIntError(argv[1],&maxiter);
  TestPositive_i(maxiter,"maximum number of iterations");
}
PARSE_FUNC(no_reduced_fft)
{
  reduced_FFT=FALSE;
}
PARSE_FUNC(no_vol_cor)
{
  volcor=FALSE;
}
PARSE_FUNC(ntheta)
{
  ScanfIntError(argv[1],&nTheta);
  TestPositive_i(nTheta,"number of theta intervals");
  nTheta++;
}
PARSE_FUNC(opt)
{
  if (strcmp(argv[1],"speed")==0) save_memory=FALSE;
  else if (strcmp(argv[1],"mem")==0) save_memory=TRUE;
  else NotSupported("Optimization method",argv[1]);
}
PARSE_FUNC(orient)
{
  if (Narg==0) NargError(Narg,"at least 1");
  if (strcmp(argv[1],"avg")==0) {
    if (Narg>2) PrintErrorHelp(
      "Illegal number of arguments (%d) to '-orient avg' option (0 or 1 expected)",Narg-1);
    orient_avg=TRUE;
    if (Narg==2) {
      TestStrLength(argv[2],MAX_FNAME);
      strcpy(avg_parms,argv[2]);
    }
  }
  else {
    if (Narg!=3) NargError(Narg,"3");
    ScanfDoubleError(argv[1],&alph_deg);
    ScanfDoubleError(argv[2],&bet_deg);
    ScanfDoubleError(argv[3],&gam_deg);
  }
}
PARSE_FUNC(phi_integr)
{
  phi_integr = TRUE;
  ScanfIntError(argv[1],&phi_int_type);
  TestRange_i(phi_int_type,"type of integration over phi",1,31);
}
PARSE_FUNC(pol)
{
  if (Narg!=1 && Narg!=2) NargError(Narg,"1 or 2");
  if (strcmp(argv[1],"cm")==0) PolRelation=POL_CM;
  else if (strcmp(argv[1],"rrc")==0) PolRelation=POL_RR;
  else if (strcmp(argv[1],"ldr")==0) PolRelation=POL_LDR;
  else if (strcmp(argv[1],"cldr")==0) PolRelation=POL_CLDR;
  else if (strcmp(argv[1],"fcd")==0) PolRelation=POL_FCD;
  else if (strcmp(argv[1],"so")==0) PolRelation=POL_SO;
  else NotSupported("Polarization relation",argv[1]);
  if (Narg==2) {
    if (strcmp(argv[2],"avgpol")==0) avg_inc_pol=TRUE;
    else PrintErrorHelp("Unknown argument '%s' to '-pol %s' option",argv[2],argv[1]);
  }
}
PARSE_FUNC(prognose)
{
  prognose=TRUE;
  strcpy(run_name,"test");
}
PARSE_FUNC(prop)
{
  double tmp;

  ScanfDoubleError(argv[1],prop_0);
  ScanfDoubleError(argv[2],prop_0+1);
  ScanfDoubleError(argv[3],prop_0+2);
  tmp=DotProd(prop_0,prop_0);
  if (tmp==0) PrintErrorHelp("Given propagation vector is null");
  tmp=1/sqrt(tmp);
  prop_0[0]*=tmp;
  prop_0[1]*=tmp;
  prop_0[2]*=tmp;
}
PARSE_FUNC(save_geom)
{
  if (Narg>1) NargError(Narg,"0 or 1");
  save_geom=TRUE;
  if (Narg==1) {
    TestStrLength(argv[1],MAX_FNAME);
    strcpy(save_geom_fname,argv[1]);
  }
}
PARSE_FUNC(scat)
{
  if (strcmp(argv[1],"dr")==0) ScatRelation=SQ_DRAINE;
  else if (strcmp(argv[1],"so")==0) ScatRelation=SQ_SO;
  else NotSupported("Scattering quantities relation",argv[1]);
}
PARSE_FUNC(scat_grid_inp)
{
  TestStrLength(argv[1],MAX_FNAME);
  strcpy(scat_grid_parms,argv[1]);
}
PARSE_FUNC(shape)
{
  int i,j,found;

  Narg--;
  found=FALSE;
  i=-1;
  while (shape_opt[++i].name!=NULL) if (strcmp(argv[1],shape_opt[i].name)==0) {
    /* set shape and shape option index */
    shape=shape_opt[i].type;
    opt.l2=i;
    opt_sh=opt;
    sh_Npars=Narg;
    /* check number of arguments */
    TestNarg_sub(Narg);
    if (shape==SH_COATED) {
      if (Narg!=1 && Narg!=4) NargError(Narg,"1 or 4");
    }
    /* parse; consistency of shape arguments is checked in InitShape() */
    if (shape==SH_READ) {
      TestStrLength(argv[2],MAX_FNAME);
      strcpy(aggregate_file,argv[2]);
    }
    else for (j=0;j<Narg;j++) ScanfDoubleError(argv[j+2],sh_pars+j);
    /* stop search */
    found=TRUE;
    break;
  }
  if(!found) NotSupported("Shape type",argv[1]);
  /* set shapename; takes place only if shapename was matched above */
  strcpy(shapename,argv[1]);
}
PARSE_FUNC(size)
{
  ScanfDoubleError(argv[1],&sizeX);
  TestPositive(sizeX,"particle size");
}
PARSE_FUNC(store_beam)
{
  store_beam = TRUE;
}
PARSE_FUNC(store_dip_pol)
{
  store_dip_pol=TRUE;
}
PARSE_FUNC(store_force)
{
  store_force = TRUE;
}
PARSE_FUNC(store_int_field)
{
  store_int_field=TRUE;
}
PARSE_FUNC(store_scat_grid)
{
  store_scat_grid = TRUE;
}
PARSE_FUNC(sym)
{
  if (strcmp(argv[1],"no")==0) NoSymmetry=TRUE;
  else if (strcmp(argv[1],"enf")==0) symmetry_enforced=TRUE;
  else NotSupported("Symmetry option",argv[1]);
}
PARSE_FUNC(test)
{
  strcpy(run_name,"test");
}
PARSE_FUNC(V)
{
  char ccver_str[MAX_LINE];
#if defined(__DECC)
  char cctype;
#elif defined(__BORLANDC__)
  int ccver;
#endif

  if (ringid==ROOT) {
    /* compiler & version (works only for selected compilers) */
    /* Intel */
#if defined(__ICC) || defined(__INTEL_COMPILER)
# define COMPILER "Intel"
# ifdef __INTEL_COMPILER
#  define CCVERSION __INTEL_COMPILER
# else
#  define CCVERSION __ICC
# endif
    sprintf(ccver_str,"%d.%d",CCVERSION/100,CCVERSION%100);
    /* DEC (Compaq) */
#elif defined(__DECC)
# define COMPILER "DEC (Compaq)"
    cctype=(__DECC_VER/10000)%10;
    if (cctype==6) cctype='T';
    else if (cctype==8) cctype='S';
    else if (cctype==9) cctype='V';
    else cctype=' ';
    sprintf(ccver_str,"%c%d.%d-%d",cctype,__DECC_VER/10000000,(__DECC_VER/100000)%100,
                                   __DECC_VER%1000);
    /* Borland */
#elif defined(__BORLANDC__)
# define COMPILER "Borland"
    sprintf(ccver_str,"%x",__BORLANDC__);
    sscanf(ccver_str,"%d",&ccver);
    sprintf(ccver_str,"%d.%d",ccver/100,ccver%100);
    /* Microsoft */
#elif defined(_MSC_VER)
# define COMPILER "Microsoft"
    sprintf(ccver_str,"%d.%d",_MSC_VER/100,_MSC_VER%100);
    /* GNU */
#elif defined(__GNUC__)
# define COMPILER "GNU"
    sprintf(ccver_str,"%d.%d.%d",__GNUC__,__GNUC_MINOR__,__GNUC_PATCHLEVEL__);
    /* unknown compiler */
#else
# define COMPILER_UNKNOWN
# define COMPILER "unknown"
#endif
    /* print version, type and compiler information */
    printf("'Amsterdam DDA' v." ADDA_VERSION "\n");
#ifdef MPI
   /* Version of MPI standard is specified, requires MPI 1.2 */
    printf("Parallel version conforming to MPI standard %d.%d\n",
           MPI_VERSION,MPI_SUBVERSION);
#else
    printf("Sequential version\n");
#endif
    printf("Built with " COMPILER " C compiler");
#ifndef COMPILER_UNKNOWN
    printf(", version %s",ccver_str);
#endif
    /* print copyright information; split in two to eliminate warning of a long string */
    printf("\n\nCopyright (C) 2006-2008 University of Amsterdam\n"\
           "This program is free software; you can redistribute it and/or modify\n"\
           "it under the terms of the GNU General Public License as published by\n"\
           "the Free Software Foundation; either version 2 of the License, or\n"\
           "(at your option) any later version.\n\n");
    printf("This program is distributed in the hope that it will be useful,\n"\
           "but WITHOUT ANY WARRANTY; without even the implied warranty of\n"\
           "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n"\
           "GNU General Public License for more details.\n\n"\
           "You should have received a copy of the GNU General Public License\n"\
           "along with this program; if not, write to the Free Software\n"\
           "Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n");
  }
  /* exit */
  Stop(0);
#undef COMPILER
#undef CCVERSION
#undef COMPILER_UNKNOWN
}
PARSE_FUNC(vec)
{
  calc_vec = TRUE;
}
PARSE_FUNC(yz)
{
  yzplane = TRUE;
}
#undef PAR
#undef PARSE_FUNC
#undef PARSE_NAME

       /* end of parsing functions */
/*=============================================================*/

static FILEHANDLE CreateLockFile(const char *fname)
   /* create locks file
      works only if USE_LOCK is enabled */
{
#ifdef USE_LOCK
  FILEHANDLE fd;
  int i;

# ifdef WINDOWS
  i=0;
  while ((fd=CreateFile(fname,GENERIC_WRITE,FILE_SHARE_WRITE,NULL,CREATE_NEW,
                        FILE_ATTRIBUTE_NORMAL,NULL))==INVALID_HANDLE_VALUE) {
    Sleep(LOCK_WAIT*1000);
    if (i++ == MAX_LOCK_WAIT_CYCLES)
      LogError(EC_ERROR,ONE_POS,"Lock file %s permanently exists",fname);
  }
# elif defined(POSIX)
#  ifdef LOCK_FOR_NFS
  struct flock lock;
#  endif
  /* open file exclusively */
  i=0;
  while ((fd=open(fname,O_WRONLY | O_CREAT | O_EXCL,0666))==-1) {
    sleep(LOCK_WAIT);
    if (i++ == MAX_LOCK_WAIT_CYCLES)
      LogError(EC_ERROR,ONE_POS,"Lock file %s permanently exists",fname);
  }
#  ifdef LOCK_FOR_NFS
  /* specify lock */
  lock.l_type=F_WRLCK;
  lock.l_whence=SEEK_SET;
  lock.l_start=0;
  lock.l_len=0;
  /* obtain lock*/
  i=0;
  while (fcntl(fd,F_SETLK,&lock)==-1) {
    /* if locked by another process wait and try again */
    if (errno==EACCES || errno==EAGAIN) {
      sleep(LOCK_WAIT);
      if (i++ == MAX_LOCK_WAIT_CYCLES)
        LogError(EC_ERROR,ONE_POS,"Lock file %s permanently exists",fname);
    }
    else { /* otherwise produce a message and continue */
      if (errno==EOPNOTSUPP || errno==ENOLCK) LogError(EC_WARN,ONE_POS,
        "Advanced file locking is not supported by the filesystem");
      else LogError(EC_WARN,ONE_POS,"Unknown problem with file locking ('%s').",strerror(errno));
      break;
    }
  }
#  endif
# endif
  /* return file handle */
  return fd;
#else
  return 0;
#endif
}

/*============================================================*/

static void RemoveLockFile(FILEHANDLE fd,const char *fname)
   /* closes and remove lock file
      works only if USE_LOCK is enabled */
{
#ifdef USE_LOCK
# ifdef WINDOWS
  /* close file */
  CloseHandle(fd);
# elif defined(POSIX)
  /* close file; all locks are automatically released */
  close(fd);
# endif
  /* remove lock file */
  RemoveErr(fname,ONE_POS);
#endif
}

/*============================================================*/

void InitVariables(void)
   /* some defaults are specified also in const.h */
{
  /* defaults */
  prop_0[0]=0;        /* by default beam propagates along z-axis */
  prop_0[1]=0;
  prop_0[2]=1;
  directory[0]=0;
  lambda=TWO_PI;
    /* initialize ref_index of scatterer */
  Nmat=Nmat_given=1;
  ref_index[0][RE]=1.5;
  ref_index[0][IM]=0.0;
    /* initialize to null to determine further whether it is initialized */
  logfile=NULL;
  logname[0]=0;

  boxX=boxY=boxZ=UNDEF;
  sizeX=UNDEF;
  dpl=UNDEF;
  strcpy(run_name,"run");
  nTheta=UNDEF;
  eps=1E-5;
  shape=SH_SPHERE;
  strcpy(shapename,"sphere");
  store_int_field=FALSE;
  store_dip_pol=FALSE;
  PolRelation=POL_LDR;
  ScatRelation=SQ_DRAINE;
  IntRelation=G_POINT_DIP;
  IterMethod=IT_QMR_CS;
  NoSymmetry=FALSE;
  symmetry_enforced=FALSE;
  prognose=FALSE;
  maxiter=UNDEF;
  jagged=1;
  beamtype=B_PLANE;
  strcpy(alldir_parms,FD_ALLDIR_PARMS);
  strcpy(avg_parms,FD_AVG_PARMS);
  strcpy(scat_grid_parms,FD_SCAT_PARMS);
  strcpy(chp_dir,FD_CHP_DIR);
  chp_time=UNDEF;
  chp_type=CHP_NONE;
  orient_avg=FALSE;
  alph_deg=bet_deg=gam_deg=0.0;
  volcor=TRUE;
  reduced_FFT=TRUE;
  save_geom=FALSE;
  save_geom_fname[0]=0;
  yzplane=UNDEF;
  all_dir=FALSE;
  scat_grid=FALSE;
  phi_integr=FALSE;
  store_scat_grid=FALSE;
  calc_Cext=TRUE;
  calc_Cabs=TRUE;
  calc_Csca=FALSE;
  calc_vec=FALSE;
  calc_asym=FALSE;
  calc_mat_force=FALSE;
  store_force=FALSE;
  load_chpoint=FALSE;
  sh_granul=FALSE;
  symX=symY=symZ=symR=TRUE;
  anisotropy=FALSE;
  save_memory=FALSE;
  memory=0;
  Ncomp=1;
}

/*============================================================*/

void ParseParameters(const int argc,char **argv)
  /* parses input parameters */
{
  int i,j,Narg;
  int found;
  char *p1,*p2;

  /* get name of executable; remove all path overhead */
  if ((p1=strrchr(argv[0],'\\'))==NULL) p1=argv[0];
  if ((p2=strrchr(argv[0],'/'))==NULL) p2=argv[0];
  exename=MAX(p1,p2)+1;
  /* initialize option */
  opt.l1=UNDEF;
  /* check first argument */
  if (argc>1 && !IsOption(argv[1]))
    PrintErrorHelp("Illegal format of first argument '%s'",argv[1]);
  /* read command line */
  for (i=1;i<argc;i++) {
    /* get number of arguments */
    Narg=0;
    while ((i+(++Narg))<argc && !IsOption(argv[i+Narg]));
    Narg--;

    argv[i]++; /* shift to remove "-" in the beginning of the string */
    found=FALSE;
    opt.l1=opt.l2=UNDEF;
    for (j=0;j<LENGTH(options);j++) if (strcmp(argv[i],options[j].name)==0) {
      opt.l1=j;
      /* check consistency, if enabled for this parameter */
      TestNarg(Narg);
      /* parse this parameter */
      (*options[j].func)(Narg,argv+i);
      /* check duplicate options */
      if (options[j].used) PrintError("Option '-%s' is used more than once",argv[i]);
      else options[j].used=TRUE;
      /* stop search */
      found=TRUE;
      break;
    }
    if(!found) PrintErrorHelp("Unknown option '-%s'",argv[i]);
    argv[i]--; /* shift back */
    i+=Narg;
  }  /* end of reading command line arguments */
}

/*============================================================*/

void VariablesInterconnect(void)
  /* finish parameters initialization based on their interconnections */
{
  double temp;

  /* initialize WaveNum ASAP */
  WaveNum  = TWO_PI/lambda;
  /* parameter interconnections */
  if (IntRelation==G_SO) reduced_FFT=FALSE;
  if (calc_Csca || calc_vec) all_dir = TRUE;
  if (store_scat_grid || phi_integr) {
    scat_grid = TRUE;
    if (yzplane==UNDEF) yzplane = FALSE;
  }
  else if (yzplane==UNDEF) yzplane = TRUE;
  /* parameter incompatibilities */
  if (orient_avg) {
    if (prop_0[2]!=1) PrintError("'-prop' and '-orient avg' can not be used together");
    if (store_int_field)
      PrintError("'-store_int_field' and '-orient avg' can not be used together");
    if (store_dip_pol)
      PrintError("'-store_dip_pol' and '-orient avg' can not be used together");
    if (store_beam) PrintError("'-store_beam' and '-orient avg' can not be used together");
    if (scat_grid) PrintError(
      "'-orient avg' can not be used with calculation of scattering for a grid of angles");
    /* this limitation should be removed in the future */
    if (all_dir)
      PrintError("Currently '-orient avg' can not be used with calculation of asym or Csca");
  }
  if (anisotropy) {
    if (PolRelation==POL_CLDR) PrintError("'-anisotr' is incompatible with '-pol cldr'");
    if (PolRelation==POL_SO) PrintError("'-anisotr' is incompatible with '-pol so'");
    if (ScatRelation==SQ_SO) PrintError("'-anisotr' is incompatible with '-scat so'");
    if (IntRelation==G_SO) PrintError("'-anisotr' is incompatible with '-int so'");
    if (Nmat%3!=0) PrintError(
      "When '-anisotr' is used 6 numbers (3 complex values) should be given per each domain");
    else Nmat=Nmat/3;
  }
  if (chp_type!=CHP_NONE) {
    if (chp_time==UNDEF && chp_type!=CHP_ALWAYS)
      PrintError("You must specify time for this checkpoint type");
    /* this limitation should be removed in the future */
    if (orient_avg) PrintError("Currently checkpoint is incompatible with '-orient avg'");
  }
  /* scale boxes by jagged; should be completely robust to overflows */
#define JAGGED_BOX(a) if (a!=UNDEF) { if ((BOX_MAX/(size_t)jagged)<(size_t)a) \
       LogError(EC_ERROR,ONE_POS,"Derived grid size (" #a ") is too large (>%d)",BOX_MAX); \
     else a*=jagged; }
  if (jagged!=1) {
    JAGGED_BOX(boxX);
    JAGGED_BOX(boxY);
    JAGGED_BOX(boxZ);
  }
#undef JAGGED_BOX
  /*determine two incident polarizations. Equivalent to rotation of X,Y,Z basis by
  angles Theta and Phi from (0,0,1) to given propagation vector */
  if (fabs(prop_0[2])>=1) {     /* can not be >1 except for machine precision */
    incPolX_0[0]=prop_0[2];
    incPolY_0[1]=1;
    incPolX_0[1]=incPolX_0[2]=incPolY_0[0]=incPolY_0[2]=0.0;
  }
  else {
    temp=sqrt(1-prop_0[2]*prop_0[2]);
    incPolX_0[0]=prop_0[0]*prop_0[2]/temp;
    incPolX_0[1]=prop_0[1]*prop_0[2]/temp;
    incPolX_0[2]=-temp;
    incPolY_0[0]=-prop_0[1]/temp;
    incPolY_0[1]=prop_0[0]/temp;
    incPolY_0[2]=0.0;
  }
  /* initialize beam description */
  InitBeam();
  /* initialize averaging over orientation */
  if (orient_avg) {
    ReadAvgParms(avg_parms);
    NoSymmetry=TRUE;
    avg_inc_pol=TRUE;
  }
  else {
    /* else - initialize rotation stuff */
    InitRotation();
    /* if not default incidence, break the symmetry completely. This can be improved to
       account for some special cases, however, then symmetry of Gaussian beam should be
       treated more thoroughly than now. */
    if (prop[2]!=1) NoSymmetry=TRUE;
  }
}

/*============================================================*/

void DirectoryLog(const int argc,char **argv)
   /* create input directory and start logfile */
{
  int  i,Nexp;
  FILE *Nexpfile;
  char sbuffer[MAX_LINE];
  char *ptmp,*compname;
  FILEHANDLE lockid;
#ifdef PARALLEL
  char *ptmp2;
#endif
#ifdef WINDOWS  /* for obtaining computer name */
  TCHAR cname[MAX_COMPUTERNAME_LENGTH+1];
  DWORD cname_size=MAX_COMPUTERNAME_LENGTH+1;
#endif

  /* devise directory name (for output files) */
  if (directory[0]==0) {
    /* ROOT processor works with ExpCount */
    if (ringid==ROOT) {
      /* lock file */
      lockid=CreateLockFile(F_EXPCOUNT_LCK);
      /* read ExpCount */
      if ((Nexpfile=fopen(F_EXPCOUNT,"r"))!=NULL) {
        if (fscanf(Nexpfile,"%i",&Nexp)!=1) Nexp=0;
        FCloseErr(Nexpfile,F_EXPCOUNT,ONE_POS);
      }
      else Nexp=0;
      /* put new number in Nexpfile */
      Nexpfile=FOpenErr(F_EXPCOUNT,"w",ONE_POS);
      fprintf(Nexpfile,"%i",Nexp+1);
      FCloseErr(Nexpfile,F_EXPCOUNT,ONE_POS);
      /* unlock */
      RemoveLockFile(lockid,F_EXPCOUNT_LCK);
    }
    /* cast Nexp to all processors */
    MyBcast(&Nexp,int_type,1,NULL);
    /* create directory name */
    sprintf(sbuffer,"m%.4g",ref_index[0][RE]);
    ptmp=strchr(sbuffer,'.');
    if (ptmp!=NULL) *ptmp='_';
    sprintf(directory,"%s%03i_%s_g%i%s",run_name,Nexp,shapename,boxX,sbuffer);
#ifdef PARALLEL
    /* add PBS or SGE job id to the directory name if available */
    if ((ptmp=getenv("PBS_JOBID"))!=NULL) {
      /* jobid is truncated at first "." */
      if ((ptmp2=strchr(ptmp,'.'))!=NULL) *ptmp2=0;
    }
    else ptmp=getenv("JOB_ID");
    if (ptmp!=NULL) sprintf(directory+strlen(directory),"id%s",ptmp);
#endif
  }
  /* make new directory and print info */
  if (ringid==ROOT) {
    MkDirErr(directory,ONE_POS);
    printf("all data is saved in '%s'\n",directory);
  }
  /* make logname; do it for all processors to enable additional logging in LogError */
  if (ringid==ROOT) sprintf(logname,"%s/" F_LOG,directory);
  else sprintf(logname,"%s/" F_LOG_ERR,directory,ringid);
  /* start logfile */
  if (ringid==ROOT) {
    /* open logfille */
    logfile=FOpenErr(logname,"w",ONE_POS);
    /* log version number */
    fprintf(logfile,"Generated by ADDA v." ADDA_VERSION "\n");
    /* get computer name */
#ifdef WINDOWS
    GetComputerName(cname,&cname_size);
    compname=cname;
#else /* POSIX and others */
    compname=getenv("HOST");
#endif
    /* write number of processors and computer name */
#ifdef PARALLEL
    /* write number of processors */
    fprintf(logfile,"The program was run on: %d processors (cores)",nprocs);
    /* add PBS or SGE host name if present, otherwise use compname */
    if ((ptmp=getenv("PBS_O_HOST"))!=NULL || (ptmp=getenv("SGE_O_HOST"))!=NULL)
      fprintf(logfile," from %s\n",ptmp);
    else if (compname!=NULL) fprintf(logfile," from %s\n",compname);
    else fprintf(logfile,"\n");
#else  /* sequential */
    if (compname!=NULL) fprintf(logfile,"The program was run on: %s\n",compname);
#endif
    /* log command line */
    fprintf(logfile,"command: '");
    for(i=0;i<argc;i++) fprintf(logfile,"%s ",argv[i]);
    fprintf(logfile,"'\n");
  }
  Synchronize(); /* needed to wait for creation of the output directory */
  LogPending();
}

/*============================================================*/

void PrintInfo(void)
  /* print info to stdout and logfile */
{
   int  i;
   char sbuffer[MAX_LINE];

   if (ringid==ROOT) {
    /* print basic parameters */
    printf("lambda: %.10g   Dipoles/lambda: %g\n",lambda,dpl);
    printf("Required relative residual norm: %g\n",eps);
    printf("Total number of occupied dipoles: %.0f\n",nvoid_Ndip);
    /* log basic parameters */
    fprintf(logfile,"lambda: %.10g\n",lambda);
    fprintf(logfile,"shape: ");
    fprintf(logfile,sh_form_str,sizeX);
    if (sh_granul) fprintf(logfile,
      "\n  domain %d is filled with %d granules of diameter %g\n"\
      "    volume fraction: specified - %g, actual - %g",gr_mat+1,gr_N,gr_d,gr_vf,gr_vf_real);
    fprintf(logfile,"\nbox dimensions: %ix%ix%i\n",boxX,boxY,boxZ);
    if (anisotropy) {
      fprintf(logfile,"refractive index (diagonal elements of the tensor):\n");
      if (Nmat==1) fprintf(logfile,"    (%.10g%+.10gi,%.10g%+.10gi,%.10g%+.10gi)\n",
                           ref_index[0][RE],ref_index[0][IM],ref_index[1][RE],ref_index[1][IM],
                           ref_index[2][RE],ref_index[2][IM]);
      else {
        for (i=0;i<Nmat;i++) {
          if (i<Nmat_given) fprintf(logfile,"    %d. (%.10g%+.10gi,%.10g%+.10gi,%.10g%+.10gi)\n",
                                    i+1,ref_index[3*i][RE],ref_index[3*i][IM],ref_index[3*i+1][RE],
                                    ref_index[3*i+1][IM],ref_index[3*i+2][RE],ref_index[3*i+2][IM]);
          else fprintf(logfile,"   %d. not specified\n",i+1);
        }
      }
    }
    else {
      fprintf(logfile,"refractive index: ");
      if (Nmat==1) fprintf(logfile,"%.10g%+.10gi\n",ref_index[0][RE],ref_index[0][IM]);
      else {
        fprintf(logfile,"1. %.10g%+.10gi\n",ref_index[0][RE],ref_index[0][IM]);
        for (i=1;i<Nmat;i++) {
          if (i<Nmat_given) fprintf(logfile,
            "                  %d. %.10g%+.10gi\n",i+1,ref_index[i][RE],ref_index[i][IM]);
          else fprintf(logfile,"                  %d. not specified\n",i+1);
        }
      }
    }
    fprintf(logfile,"Dipoles/lambda: %g\n",dpl);
    if (volcor_used) fprintf(logfile,"\t(Volume correction used)\n");
    fprintf(logfile,"Required relative residual norm: %g\n",eps);
    fprintf(logfile,"Total number of occupied dipoles: %.0f\n",nvoid_Ndip);
    if (Nmat>1) {
      fprintf(logfile,"  per domain: 1. %.0f\n",mat_count[0]);
      for (i=1;i<Nmat;i++) fprintf(logfile,"              %d. %.0f\n",i+1,mat_count[i]);
    }
    fprintf(logfile,"Volume-equivalent size parameter: %.10g\n",ka_eq);
    /* log incident beam and polarization polarization */
    fprintf(logfile,"\n---In laboratory reference frame:---\nIncident beam: %s\n",beam_descr);
    fprintf(logfile,"Incident propagation vector: (%g,%g,%g)\n",
            prop_0[0],prop_0[1],prop_0[2]);
    fprintf(logfile,"Incident polarization Y(par): (%g,%g,%g)\n",
            incPolY_0[0],incPolY_0[1],incPolY_0[2]);
    fprintf(logfile,"Incident polarization X(per): (%g,%g,%g)\n\n",
            incPolX_0[0],incPolX_0[1],incPolX_0[2]);
    /* log particle orientation */
    if (orient_avg) fprintf(logfile,"Particle orientation - averaged\n%s\n",avg_string);
    else {
      /* log incident polarization after transformation */
      if (alph_deg!=0 || bet_deg!=0 || gam_deg!=0) {
        fprintf(logfile,"Particle orientation (deg): alpha=%g, beta=%g, gamma=%g\n\n"\
                        "---In particle reference frame:---\n",alph_deg,bet_deg,gam_deg);
        if (beam_asym) fprintf(logfile,"Incident Beam center position: (%g,%g,%g)\n",
                               beam_center[0],beam_center[1],beam_center[2]);
        fprintf(logfile,"Incident propagation vector: (%g,%g,%g)\n",
                prop[0],prop[1],prop[2]);
        fprintf(logfile,"Incident polarization Y(par): (%g,%g,%g)\n",
                incPolY[0],incPolY[1],incPolY[2]);
        fprintf(logfile,"Incident polarization X(per): (%g,%g,%g)\n\n",
                incPolX[0],incPolX[1],incPolX[2]);
      }
      else fprintf(logfile,"Particle orientation: default\n\n");
    }
    /* log Polarization relation */
    if (PolRelation==POL_CM)
      fprintf(logfile,"Polarization relation: 'Clausius-Mossotti'\n");
    else if (PolRelation==POL_RR)
      fprintf(logfile,"Polarization relation: 'Radiative Reaction Correction'\n");
    else if (PolRelation==POL_LDR) {
      fprintf(logfile,"Polarization relation: 'Lattice Dispersion Relation'");
      if (avg_inc_pol) fprintf(logfile," (averaged over incident polarization)");
      fprintf(logfile,"\n");
    }
    else if (PolRelation==POL_CLDR)
      fprintf(logfile,"Polarization relation: 'Corrected Lattice Dispersion Relation'\n");
    else if (PolRelation==POL_FCD)
      fprintf(logfile,"Polarization relation: 'Filtered Coupled Dipoles'\n");
    else if (PolRelation==POL_SO)
      fprintf(logfile,"Polarization relation: 'Second Order'\n");
    /* log Scattering Quantities formulae */
    if (ScatRelation==SQ_DRAINE)
      fprintf(logfile,"Scattering quantities formulae: 'by Draine'\n");
    else if (ScatRelation==SQ_SO)
      fprintf(logfile,"Scattering quantities formulae: 'Second Order'\n");
    /* log Interaction term prescription */
    if (IntRelation==G_POINT_DIP)
      fprintf(logfile,"Interaction term prescription: 'as Point dipoles'\n");
    else if (IntRelation==G_FCD)
      fprintf(logfile,"Interaction term prescription: 'Filtered Green's tensor'\n");
    else if (IntRelation==G_FCD_ST)
      fprintf(logfile,"Interaction term prescription: 'Filtered Green's tensor (quasistatic)'\n");
    else if (IntRelation==G_SO)
      fprintf(logfile,"Interaction term prescription: 'Second Order'\n");
    /* log FFT method */
#ifdef FFTW3
    fprintf(logfile,"FFT algorithm: FFTW3\n");
#elif defined(FFT_TEMPERTON)
    fprintf(logfile,"FFT algorithm: by C.Temperton\n");
#endif
    /* log Iterative Method */
    if (IterMethod==IT_CGNR)
      fprintf(logfile,"Iterative Method: CGNR\n");
    else if (IterMethod==IT_BICGSTAB)
      fprintf(logfile,"Iterative Method: Bi-CG Stabilized\n");
    else if (IterMethod==IT_BICG_CS)
      fprintf(logfile,"Iterative Method: Bi-CG (complex symmetric)\n");
    else if (IterMethod==IT_QMR_CS)
      fprintf(logfile,"Iterative Method: QMR (complex symmetric)\n");
    /* log Symmetry options */
    if (symmetry_enforced) fprintf(logfile,"Symmetry is enforced by user (warning!)\n");
    else if (NoSymmetry) fprintf(logfile,"No symmetries are used\n");
    /* log optimization method */
    if (save_memory) fprintf(logfile,"Optimization is done for minimum memory usage\n");
    else fprintf(logfile,"Optimization is done for maximum speed\n");
    /* log Checkpoint options */
    if (load_chpoint) fprintf(logfile,"Simulation is continued from a checkpoint\n");
    if (chp_type!=CHP_NONE) {
      fprintf(logfile,"Checkpoint is turned on:\n");
      if (chp_type==CHP_NORMAL) fprintf(logfile,"    type = normal\n");
      else if (chp_type==CHP_REGULAR) fprintf(logfile,"    type = regular\n");
      else if (chp_type==CHP_ALWAYS) fprintf(logfile,"    type = always\n");
      if (chp_time==UNDEF) fprintf(logfile,"    time = no limit\n");
      else {
        PrintTime(sbuffer,&chp_time);
        /* chp_time is converted to long to avoid problems with definition of time_t
           (can be either int or long) */
        fprintf(logfile,"    time = %s(%ld sec)\n",sbuffer,(long)chp_time);
      }
    }
    if (load_chpoint || chp_type!=CHP_NONE)
      fprintf(logfile,"    directory = '%s'\n",chp_dir);
  }
}
