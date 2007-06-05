/* FILE : Romberg.c
 * AUTH : Maxim Yurkin
 * DESCR: 2D-integration routine, based on the qromb-routine from
 *        Numerical Recipes. The main differences from the original
 *        code are :
 *        1) Romberg1D works with precalculated values
 *           Romberg2D is adaptive, works through pointer to a function,
 *           can handle also precalculated values
 *        2) From 1D-integration to 2D-integration. The function-values for
 *           OuterQromb are themselves the result of an integration
 *           performed by InnerQromb.
 *        The types declared in types.h deserve some attention :
 *        - Parms_1D : this structure holds the parameters for integration
 *             in one dimension. They must be set outside of the
 *             Romberg routine.
 *        All the routines normalize the result on the interval width, i.e.
 *           actually averaging takes place.
 *
 *        Previous version by Martin Frijlink
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include "vars.h"
#include "cmplx.h"
#include "types.h"
#include "comm.h"
#include "const.h"
#include "Romberg.h"
#include "memory.h"
#include "io.h"

/* LOCAL VARIABLES */

static int dim;                /* dimension of the data (integrated simultaneously) */
static int N_eval;             /* number of function evaluations in inner cycle */
static int N_tot_eval;         /* total number of function evaluation */
static int no_convergence;     /* number of inner integrals that did not converge */
static FILE *file;             /* file to print info */
/* used in InnerQromb */
static double **s_in,           /* Intermediate results */
              *dss_in,          /* extrapolation error in final result */
              **st_in,          /* Dummy carrying transposed s[][] */
              *h_in;            /* squared stepsize */
/* used in OuterQromb */
static double **s_out,*dss_out,**st_out,*h_out;  /* analogous to the above */
/* used in InnerTrapzd */
static double *dummy_in,        /* save function values */
              *sum_in;          /* sums of function values */
/* used in OuterTrapzd */
static double *dummy_out,*sum_out;  /* analogous to the above */
/* pointer to the function that is integrated */
static void (*func)(int theta,int phi,double *res);
static const Parms_1D *input;  /* parameters of integration */
static double dtheta,dphi;     /* lengths of intervals */

/*============================================================*/

static void AllocateAll(void)
   /* allocates all needed memory */
{
  s_in = dMatrix(1,input[PHI].Jmax+1,0,dim-1);         /* InnerQromb */
  st_in = dMatrix(0,dim-1,1,input[PHI].Jmax+1);
  dss_in = (double *) malloc(dim*sizeof(double));
  h_in = dVector(1,input[PHI].Jmax+1);

  s_out = dMatrix(1,input[THETA].Jmax+1,0,dim-1);      /* OuterQromb */
  st_out = dMatrix(0,dim-1,1,input[THETA].Jmax+1);
  dss_out = (double *) malloc(dim*sizeof(double));
  h_out = dVector(1,input[THETA].Jmax+1);

  dummy_in = (double *) malloc(dim*sizeof(double));    /* InnerTrapzd */
  sum_in = (double *) malloc(dim*sizeof(double));

  dummy_out = (double *) malloc(dim*sizeof(double));   /* OuterTrapzd */
  sum_out = (double *) malloc(dim*sizeof(double));
}

/*============================================================*/

static void FreeAll(void)
   /* frees all memory */
{
  Free_dMatrix(s_in,1,input[PHI].Jmax+1,0);      /* InnerQromb */
  Free_dMatrix(st_in,0,dim-1,1);
  free(dss_in);
  Free_dVector(h_in,1);

  Free_dMatrix(s_out,1,input[THETA].Jmax+1,0);   /* OuterQromb */
  Free_dMatrix(st_out,0,dim-1,1);
  free(dss_out);
  Free_dVector(h_out,1);

  free(dummy_in);                                /* InnerTrapzd */
  free(sum_in);

  free(dummy_out);                               /* OuterTrapzd */
  free(sum_out);
}

/*============================================================*/

static void PolInterp(const double xa[],const double ya[],const int n,
                      const double x,double *y,double *dy)
  /* Inter- or extrapolate the n function-values ya[], at
   * corresponding gridpoints xa[], to y, at gridpoint x
   * with Neville's algorithm. The estimated error in y
   * will be stored in dy. */
{
  int i,m,     /* Loop indices */
      ns;
  double
    dif,     /* distance(xa[i],x) */
    ho,hp,   /* dummy variables */
    den,w;   /* idem */
  double *c,*d;

  c = dVector(1,n);
  d = dVector(1,n);

  /* Find the index of the gridpoint closest to x */
  for (i=1, dif=fabs(x-xa[1]); i<=n; ++i) {
    double dift; /* temporary dif */
    if ( (dift=fabs(x-xa[i])) < dif ) {
      ns=i;
      dif=dift;
    }
    c[i] = ya[i];
    d[i] = ya[i];
  }
  /* Search-algorithm */
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;++i) { /* Calculate corrections for polynomial of one degree higher */
      ho = xa[i]-x;
      hp = xa[i+m]-x;
      w = c[i+1] - d[i];
      if ((den=ho-hp)==0) /* Catch for identical (to within roundoff) xa's. */
         LogError(EC_ERROR,ONE_POS,"Error in routine POLINT"); /* Preventing division by zero */
      den  = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    /* update *y and *dy */
    *dy = (2*ns < (n-m)) ? c[ns+1] : d[ns--];
    *y += *dy;
  }
  Free_dVector(c,1);
  Free_dVector(d,1);
}

/*============================================================*/

static void Trapzd1D(/* Calculate the n'th refinement stage of the integration */
		 const Parms_1D param,  /* parameters of integration */
                 const int size,        /* size of block of data */
                 const double *data,    /* written as sequential blocks */
                 double *s,             /* previous result (replaced by new) */
		 double *sum,           /* allocated memory for buffer */
		 const int n)
{
  int j,comp,step,index;
  double tnm,delta;

  delta=param.max-param.min;

  if (n == 1) {    /* Initialise the algorithm */    
    if (param.equival)
      for (comp=0;comp<size;++comp) s[comp] = delta*data[comp]; 
    else {
      index=(param.Grid_size-1)*size;
      for (comp=0;comp<size;++comp)
        s[comp] = 0.5*delta*(data[comp]+data[index+comp]);
    }
  }
  else {  /* refinement stage n */
    for (comp=0;comp<size;++comp) sum[comp]=0;
    step = (param.Grid_size-1) >> (n-2);

    for (j=step>>1;j<param.Grid_size;j+=step) {   /* initialise refinement */
      index=j*size;
      for (comp=0;comp<size;++comp) sum[comp] += data[index+comp];
    }
    /* add refinement */
    tnm = pow(2,2-n);
    for (comp=0;comp<size;++comp)
      s[comp] = 0.5*(s[comp] + delta*sum[comp]*tnm);
  }
}

/*============================================================*/

double Romberg1D(/* Performs integration of data */
                    const Parms_1D param,/* parameters of integration */
                    const int size,      /* size of block of data */
                    const double *data,  /* written as sequential blocks */
                    double *ss)          /* where to put result */

   /* Since all values are already calculated, no adaptation is used (all data is used)
      Result is normalized on the interval width, i.e. actually averaging takes place.
      returns mean-square error (mean) */
{
  double **s,           /* Intermediate results */
         *dss,          /* extrapolation error in final result */
         **st,          /* Dummy carrying transposed s[][] */
         *h,            /* squared stepsize */
         *sum;          /* workspace for Trapzd1D */
  int j,comp;
  double abs_ss,abs_dss; /* ...-norms(squared) of ss and dss */
  double scale;

  /* allocate memory */
  s = dMatrix(1,param.Jmax+1,0,size-1);
  st = dMatrix(0,size-1,1,param.Jmax+1);
  dss = (double *) malloc(size*sizeof(double));
  h = dVector(1,param.Jmax+1);
  sum = (double *) malloc(size*sizeof(double));

  if (param.min==param.max) {
    memcpy(ss,data,size*sizeof(double));
    return 0;
  }
  h[1]=1.0;
  for (j=1;j<=param.Jmax;++j) {
    Trapzd1D(param,size,data,s[j],sum,j);
    for (comp=0;comp<size;++comp) st[comp][j] = s[j][comp];
    memcpy(s[j+1],s[j],size*sizeof(double));
    h[j+1] = 0.25*h[j];
  }
  abs_ss=abs_dss=0;
  for (comp=0; comp<size; ++comp) {
    PolInterp(&h[param.Jmax-param.K],&st[comp][param.Jmax-param.K],
              param.K,0.0,&ss[comp],&dss[comp]);
    abs_ss+=ss[comp]*ss[comp];
    abs_dss+=dss[comp]*dss[comp];
  }
  /* free memory */
  Free_dMatrix(s,1,param.Jmax+1,0);
  Free_dMatrix(st,0,size-1,1);
  free(dss);
  Free_dVector(h,1);
  free(sum);

  scale= 1/(param.max - param.min);
  for(comp=0;comp<size;comp++) ss[comp]*=scale;

  return (sqrt(abs_dss/abs_ss));
}

/*============================================================*/

static void InnerTrapzd(const int fixed,double *s,const int n)
  /* Calculate the n'th refinement stage of the integration
     of func over phi_min < phi < phi_max for fixed theta = th_f */
{
  int j,comp,step;
  double tnm;

  if (n == 1) {    /* Initialise the algorithm */
    (*func)(fixed,0,dummy_in);
    N_eval ++;

    if (input[PHI].equival)
      for (comp=0;comp<dim;++comp) s[comp] = dphi*dummy_in[comp];
    else {
      (*func)(fixed,input[PHI].Grid_size-1,sum_in);
      N_eval ++;
      for (comp=0;comp<dim;++comp)
      s[comp] = 0.5*dphi*(dummy_in[comp]+sum_in[comp]);
    }
  }
  else {  /* refinement stage n */
    for (comp=0;comp<dim;++comp) sum_in[comp]=0;
    step = (input[PHI].Grid_size-1) >> (n-2);

    for (j=step>>1;j<input[PHI].Grid_size;j+=step) {   /* initialise refinement */
      (*func)(fixed,j,dummy_in);
      N_eval++;
      for (comp=0;comp<dim;++comp) sum_in[comp] += dummy_in[comp];
    }
    /* add refinement */
    tnm = pow(2,2-n);
    for (comp=0;comp<dim;++comp)
      s[comp] = 0.5*(s[comp] + dphi*sum_in[comp]*tnm);
  }
}

/*============================================================*/

static void InnerQromb(const int fixed,double *ss)
  /* Integrate func for fixed theta=fixed */
{
  int j,k,comp;
  double abs_ss,abs_dss; /* ...-norms of ss and dss */

  if (input[PHI].min==input[PHI].max) {
    (*func)(fixed,0,ss);
    N_eval++;
    return;
  }
  h_in[1]=1.0;
  for (j=1;j<=input[PHI].Jmax;++j) {
    InnerTrapzd(fixed,s_in[j],j);
    for (comp=0;comp<dim;++comp) st_in[comp][j] = s_in[j][comp];

    if (j >= ROMB_KMIN) {
      k = MIN(j,input[PHI].K);
      for (comp=0; comp<dim; ++comp) {
	PolInterp(&h_in[j-k],&st_in[comp][j-k],k,0.0,&ss[comp],&dss_in[comp]);
      }
      abs_ss = fabs(ss[0]);
      abs_dss = fabs(dss_in[0]);
      if ( abs_dss < input[PHI].eps*abs_ss || abs_ss==0) return;
    }

    memcpy(s_in[j+1],s_in[j],dim*sizeof(double));
    h_in[j+1] = 0.25*h_in[j];
  }

  fprintf(file,"Inner_qromb converged only to d=%g for cosine value #%d\n",abs_dss/abs_ss,fixed);
  fflush(file);
  no_convergence++;
}

/*============================================================*/

static void OuterTrapzd(double *s,const int n)
  /* Calculate n'th refinement for the outer integration of func */
{
  int j,comp,step;
  double tnm;

  if (n == 1) {    /* Initialise the algorithm */
    InnerQromb(0,dummy_out);

    if (input[THETA].equival)
      for (comp=0;comp<dim;++comp) s[comp] = dtheta*dummy_out[comp];
    else {
      InnerQromb(input[THETA].Grid_size-1,sum_out);
      for (comp=0;comp<dim;++comp)
        s[comp] = 0.5*dtheta*(dummy_out[comp]+sum_out[comp]);
    }
  }
  else {  /* refinement stage n */
    for (comp=0;comp<dim;++comp) sum_out[comp]=0;
    step = (input[THETA].Grid_size-1) >> (n-2);

    for (j=step>>1;j<input[THETA].Grid_size;j+=step) {   /* initialise refinement */
      InnerQromb(j,dummy_out);
      for (comp=0;comp<dim;++comp) sum_out[comp] += dummy_out[comp];
    }
    /* add refinement */
    tnm = pow(2,2-n);
    for (comp=0;comp<dim;++comp)
      s[comp] = 0.5*(s[comp] + dtheta*sum_out[comp]*tnm);
  }
}

/*============================================================*/

static double OuterQromb(double *ss)
  /* Performs outer integration */
{
  int j,k,comp;
  double abs_ss,abs_dss; /* ...-norms of ss and dss */

  if (input[THETA].min==input[THETA].max) {
    InnerQromb(0,ss);
    N_tot_eval = N_eval;
    return 0;
  }
  h_out[1]=1.0;
  for (j=1;j<=input[THETA].Jmax;++j) {
    N_eval = 0;
    OuterTrapzd(s_out[j],j);
    for (comp=0;comp<dim;++comp) st_out[comp][j] = s_out[j][comp];
    fprintf(file,"%d\t\t%d integrand-values were used.\n",j,N_eval);
    fflush(file);
    N_tot_eval += N_eval;

    if (j >= ROMB_KMIN) {
      k = MIN(j,input[THETA].K);
      for (comp=0; comp<dim; ++comp) {
	PolInterp(&h_out[j-k],&st_out[comp][j-k],k,0.0,&ss[comp],&dss_out[comp]);
      }
      abs_ss = fabs(ss[0]);
      abs_dss = fabs(dss_out[0]);
      if ( abs_dss < input[THETA].eps*abs_ss || abs_ss==0) return 0;
    }

    memcpy(s_out[j+1],s_out[j],dim*sizeof(double));
    h_out[j+1] = 0.25*h_out[j];
  }
  return (abs_dss/abs_ss);
}

/*============================================================*/

void Romberg2D(const Parms_1D parms_input[2],
	       void (*func_input)(int theta,int phi,double *res),
	       const int dim_input,double *ss,const char *fname)
/* Integrate 2D func with Romberg's method according to input's parameters.
   Argument dim_input gives the number of components of (double *).
   Consistency between func and dim_input is the user's responsibilty.
   Result is normalized on the interval widths, i.e. actually averaging takes place. */
{
  double error, scale;
  int comp;
  char buf1[MAX_WORD],buf2[MAX_WORD];

  /* initialize global values */
  dim = dim_input;
  func = func_input;
  input = parms_input;  /* conversion from constant to non-constant pointer */
  file=FOpenErr(fname,"w",ONE_POS);
  no_convergence = 0;
  N_tot_eval=0;
  dphi=input[PHI].max-input[PHI].min;
  dtheta=input[THETA].max-input[THETA].min;

  AllocateAll();      /* allocate memory */

  if (orient_avg) {
    strcpy(buf1,"BETA");
    strcpy(buf2,"GAMMA");
  }
  else {
    strcpy(buf1,"THETA");
    strcpy(buf2,"PHI");
  }
  /* print info */
  fprintf(file,"                   %4s(rad)   cos(%s)\n"\
               "EPS                    %-7g   %g\n"\
               "Maximum number of\n"\
               "refinement-stages      %-7d   %d\n"\
               "Number of evaluations\n"\
               "for an extrapolation   %-7d   %d\n"\
               "lower boundary         %-7g   %g\n"\
               "upper boundary         %-7g   %g\n",
               buf2,buf1,
               input[PHI].eps,input[THETA].eps,
               input[PHI].Jmax,input[THETA].Jmax,
               input[PHI].K,input[THETA].K,
	       input[PHI].min,input[THETA].min,
	       input[PHI].max,input[THETA].max);
  fprintf(file,"\n\nOuter-Loop\tInner Loop\n");
  fflush(file);

  error=OuterQromb(ss);    /* main calculation */

  /* scale the result */
  if (dtheta==0) scale=1;
  else scale= dtheta;       /* theta is really cos(theta) */
  if (dphi!=0) scale *= dphi;
  scale=1/scale;
  for(comp=0;comp<dim;comp++) ss[comp]*=scale;

  /* finalize log */
  if (error==0) {
    if (no_convergence==0) PrintBoth(file,"All inner integrations converged\n"\
    		                          "The outer integration converged\n");
    else PrintBoth(file,"%d inner integrations did not converge.\n"\
                        "The outer integration converged\n",no_convergence);
  }
  else {
    if (no_convergence==0) PrintBoth(file,"Only the outer integration did not converge \n"\
 	                                  "It reached d=%g\n",error);
    else PrintBoth(file,
	           "%d inner integrations did not converge.\n"\
                   "The outer integration did not converge\n"\
	           "The outer integration reached d=%g\n",no_convergence,error);
  }
  PrintBoth(file,"In total %d evaluations were used\n",N_tot_eval);
  FCloseErr(file,fname,ONE_POS);

  FreeAll();  /* free all memory */
}

