/* FILE : Romberg.c
 * AUTH : Maxim Yurkin (based on old code by Martijn Frijlink)
 * DESCR: 2D-integration routine, based on the qromb-routine from 
 *        Numerical Recipes. The main differences from the original
 *        code are :
 *        1) Romberg1D works with precalculated values
 *           Romberg2D is adaptive, works through pointer to a function,
 *           can handle also precalculated values
 *        2) From 1D-integration to 2D-integration. The function-values for
 *           outer_qromb are themselves the result of an integration
 *           performed by inner_qromb.
 *        The types declared in types.h deserve some attention :
 *        - Parms_1D : this structure holds the parameters for integration
 *             in one dimension. They must be set outside of the 
 *             Romberg routine.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "cmplx.h"
#include "types.h"
#include "comm.h"
#include "const.h"
#include "Romberg.h"
#include "memory.h"

int dim;
int N_eval;
int N_tot_eval;
int no_convergence;
FILE *file;

double **s_in,           /* Intermediate results */
       *dss_in,          /* extrapolation error in final result */
       **st_in,          /* Dummy carrying transposed s[][] */
       *h_in;            /* squared stepsize */

double **s_out,*dss_out,**st_out,*h_out;        
double *dummy_in,*sum_in;
double *dummy_out,*sum_out;

void (*func)(int theta,int phi,double *res);
Parms_1D *input;

double dtheta,dphi; /* lengths of intervals */

/*============================================================*/

void allocate_all(void)
   /* allocates all needed memory */
{
  s_in = dmatrix(1,input[PHI].JMAX+1,0,dim-1);         /* inner_qromb */
  st_in = dmatrix(0,dim-1,1,input[PHI].JMAX+1);
  dss_in = (double *) malloc(dim*sizeof(double));
  h_in = dvector(1,input[PHI].JMAX+1);

  s_out = dmatrix(1,input[THETA].JMAX+1,0,dim-1);      /* outer qromb */
  st_out = dmatrix(0,dim-1,1,input[THETA].JMAX+1);
  dss_out = (double *) malloc(dim*sizeof(double));
  h_out = dvector(1,input[THETA].JMAX+1);

  dummy_in = (double *) malloc(dim*sizeof(double));    /* inner_trapzd */ 
  sum_in = (double *) malloc(dim*sizeof(double));                         
                                                       
  dummy_out = (double *) malloc(dim*sizeof(double));   /* outrer_trapzd */ 
  sum_out = (double *) malloc(dim*sizeof(double));
}

/*============================================================*/

void free_all(void)
   /* frees all memory */
{
  free_dmatrix(s_in,1,input[PHI].JMAX+1,0);      /* inner_qromb */
  free_dmatrix(st_in,0,dim-1,1);
  free(dss_in);
  free_dvector(h_in,1);
  
  free_dmatrix(s_out,1,input[THETA].JMAX+1,0);   /* outer qromb */
  free_dmatrix(st_out,0,dim-1,1);
  free(dss_out);
  free_dvector(h_out,1);
  
  free(dummy_in);                                /* inner_trapzd */
  free(sum_in);                                  
  
  free(dummy_out);                               /* outrer_trapzd */
  free(sum_out);
}

/*============================================================*/

void polint(/* Inter- or extrapolate the n function-values ya[], at
	     * corresponding gridpoints xa[], to y, at gridpoint x 
	     * with Neville's algorithm. The estimated error in y
	     * will be stored in dy. */
	    double xa[],double ya[],int n,
	    double x,double *y,double *dy)
{
  int i,m,     /* Loop indices */
      ns;    
  double 
    dif,     /* distance(xa[i],x) */
    ho,hp,   /* dummy variables */
    den,w;   /* idem */
  double *c,*d;

  c = dvector(1,n);
  d = dvector(1,n);

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
      if ( (den=ho-hp) == 0.0)                                   /* Catch for identical (to within roundoff) xa's. */
         LogError(EC_ERROR,ONE,POSIT,"Error in routine POLINT"); /* Preventing division by zero */
      den  = w/den;
      d[i] = hp*den;
      c[i] = ho*den;
    }
    /* update *y and *dy */ 
    *dy = (2*ns < (n-m)) ? c[ns+1] : d[ns--];
    *y += *dy;
  }
  free_dvector(c,1);
  free_dvector(d,1);
}

/*============================================================*/

void trapzd1D(/* Calculate the n'th refinement stage of the integration */
		 Parms_1D param,  /* parameters of integration */
                 int size,        /* size of block of data */
                 double *data,    /* written as sequential blocks */
                 double *s,       /* previous result (replaced by new) */
		 double *sum,     /* allocated memory for buffer */
		 int n)
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
                    Parms_1D param,  /* parameters of integration */
                    int size,      /* size of block of data */
                    double *data,  /* written as sequential blocks */
                    double *ss)    /* where to put result */

   /* Since all values are already calculated, no adaptation is used (all data is used)
      returns mean-square error (mean) */
{
  double **s,           /* Intermediate results */
         *dss,          /* extrapolation error in final result */
         **st,          /* Dummy carrying transposed s[][] */
         *h,            /* squared stepsize */
         *sum;          /* workspace for trapzd1D */
  int j,comp;
  double abs_ss,abs_dss; /* ...-norms(squared) of ss and dss */
  double scale;

  /* allocate memory */
  s = dmatrix(1,param.JMAX+1,0,size-1);      
  st = dmatrix(0,size-1,1,param.JMAX+1);
  dss = (double *) malloc(size*sizeof(double));
  h = dvector(1,param.JMAX+1);
  sum = (double *) malloc(size*sizeof(double));  

  if (param.min==param.max) {
    memcpy(ss,data,size*sizeof(double)); 
    return 0;
  }
  h[1]=1.0;
  for (j=1;j<=param.JMAX;++j) { 
    trapzd1D(param,size,data,s[j],sum,j);
    for (comp=0;comp<size;++comp) st[comp][j] = s[j][comp];
    memcpy(s[j+1],s[j],size*sizeof(double));
    h[j+1] = 0.25*h[j];
  }
  abs_ss=abs_dss=0;
  for (comp=0; comp<size; ++comp) {
    polint(&h[param.JMAX-param.K],&st[comp][param.JMAX-param.K],param.K,0.0,&ss[comp],&dss[comp]);
    abs_ss+=ss[comp]*ss[comp];
    abs_dss+=dss[comp]*dss[comp];
  }
  /* free memory */
  free_dmatrix(s,1,param.JMAX+1,0);   
  free_dmatrix(st,0,size-1,1);
  free(dss);
  free_dvector(h,1);
  free(sum);

  scale= 1/(param.max - param.min);
  for(comp=0;comp<size;comp++) ss[comp]*=scale;
    
  return (sqrt(abs_dss/abs_ss));
}

/*============================================================*/

void inner_trapzd(/* Calculate the n'th refinement stage of the integration
		   * of func over phi_min < phi < phi_max for 
		   * fixed theta = th_f */
		  int fixed,
		  double *s,
		  int n)
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

void inner_qromb(/* Integrate func for fixed theta=fixed */
		    int fixed,
		    double *ss)
{
  int j,k,comp;
  double abs_ss,abs_dss; /* ...-norms of ss and dss */

  if (input[PHI].min==input[PHI].max) {
    (*func)(fixed,0,ss);
    N_eval++;
    return;
  }
  h_in[1]=1.0;
  for (j=1;j<=input[PHI].JMAX;++j) {
    inner_trapzd(fixed,s_in[j],j);
    for (comp=0;comp<dim;++comp) st_in[comp][j] = s_in[j][comp];

    if (j >= KMIN) {
      k = MIN(j,input[PHI].K);
      for (comp=0; comp<dim; ++comp) {
	polint(&h_in[j-k],&st_in[comp][j-k],k,0.0,&ss[comp],&dss_in[comp]);
      }
      abs_ss = fabs(ss[0]);
      abs_dss = fabs(dss_in[0]);
      if ( abs_dss < input[PHI].INT_EPS*abs_ss || abs_ss==0) return;
    }
    
    memcpy(s_in[j+1],s_in[j],dim*sizeof(double));
    h_in[j+1] = 0.25*h_in[j];
  }
  
  fprintf(file,"Inner_qromb converged only to d=%lg for cosine value #%d\n",
	  (double) abs_dss/abs_ss,fixed);
  fflush(file);
  no_convergence++;
}

/*============================================================*/

void outer_trapzd(/* Calculate n'th refinement for the outer integration
		   * of func */
		  double *s,
		  int n)
{
  int j,comp,step;
  double tnm;
  
  if (n == 1) {    /* Initialise the algorithm */
    inner_qromb(0,dummy_out);

    if (input[THETA].equival)
      for (comp=0;comp<dim;++comp) s[comp] = dtheta*dummy_out[comp];
    else {
      inner_qromb(input[THETA].Grid_size-1,sum_out);
      for (comp=0;comp<dim;++comp)
        s[comp] = 0.5*dtheta*(dummy_out[comp]+sum_out[comp]);
    }
  }
  else {  /* refinement stage n */
    for (comp=0;comp<dim;++comp) sum_out[comp]=0;
    step = (input[THETA].Grid_size-1) >> (n-2);

    for (j=step>>1;j<input[THETA].Grid_size;j+=step) {   /* initialise refinement */
      inner_qromb(j,dummy_out);
      for (comp=0;comp<dim;++comp) sum_out[comp] += dummy_out[comp];
    }
    /* add refinement */
    tnm = pow(2,2-n);
    for (comp=0;comp<dim;++comp)
      s[comp] = 0.5*(s[comp] + dtheta*sum_out[comp]*tnm);
  }
}

/*============================================================*/

double outer_qromb(/* Performs outer integration */ 
                    double *ss)
{
  int j,k,comp;
  double abs_ss,abs_dss; /* ...-norms of ss and dss */

  if (input[THETA].min==input[THETA].max) {
    inner_qromb(0,ss);
    N_tot_eval = N_eval;
    return 0;
  }
  h_out[1]=1.0;
  for (j=1;j<=input[THETA].JMAX;++j) { 
    N_eval = 0;
    outer_trapzd(s_out[j],j);
    for (comp=0;comp<dim;++comp) st_out[comp][j] = s_out[j][comp];
    fprintf(file,"%d\t\t%d integrand-values were used.\n",j,N_eval);
    fflush(file);
    N_tot_eval += N_eval;

    if (j >= KMIN) {
      k = MIN(j,input[THETA].K);
      for (comp=0; comp<dim; ++comp) {
	polint(&h_out[j-k],&st_out[comp][j-k],k,0.0,&ss[comp],&dss_out[comp]);
      }
      abs_ss = fabs(ss[0]);
      abs_dss = fabs(dss_out[0]);
      if ( abs_dss < input[THETA].INT_EPS*abs_ss || abs_ss==0) return 0;
    }

    memcpy(s_out[j+1],s_out[j],dim*sizeof(double));
    h_out[j+1] = 0.25*h_out[j];
  }
  return (abs_dss/abs_ss);
}

/*============================================================*/

void Romberg2D(/* Integrate 2D func with Romberg's method
		* according to input's parameters. Argument
		* dim_input gives the number of components of double *.
		* Consistency between func and dim_input is the user's
		* responsibilty */
		Parms_1D parms_input[Narg],
		void (*func_input)(int theta,int phi,double *res),
		int dim_input,
                double *ss,
                char *fname)
{
  double error, scale;
  int comp;
  char buf1[10],buf2[10];

  extern int orient_avg;

  /* initialize global values */
  dim = dim_input;
  func = func_input;
  input = parms_input;
  if ((file=fopen(fname,"w"))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"Failed to open file '%s'",fname);
  no_convergence = 0;
  N_tot_eval=0;
  dphi=input[PHI].max-input[PHI].min;
  dtheta=input[THETA].max-input[THETA].min;

  /* check for errors */
  if (input[THETA].JMAX<KMIN && input[THETA].min != input[THETA].max)
    LogError(EC_ERROR,ONE,POSIT,"Too small JMAX (%d) for THETA(BETA)",input[THETA].JMAX);
  if (input[PHI].JMAX<KMIN && input[PHI].min != input[PHI].max)
    LogError(EC_ERROR,ONE,POSIT,"Too small JMAX (%d) for PHI(GAMMA)",input[PHI].JMAX);

  allocate_all();      /* allocate memory */

  if (orient_avg) {
    strcpy(buf1,"GAMMA");
    strcpy(buf2,"BETA");
  }
  else {
    strcpy(buf1,"THETA");
    strcpy(buf2,"PHI");
  }
  /* print info */
  fprintf(file,"\t\t\t%s\tcos(%s)\n"\
                  "EPS\t\t\t%.8lg\t%.8lg\n"\
                  "Maximal number of\nrefinement-stages\t%d\t%d\n"\
                  "Number of evaluations\nfor an extrapolation\t%d\t%d\n"\
                  "lower boundary\t\t%.8lg\t%.8lg\n"\
                  "upper boundary\t\t%.8lg\t%.8lg\n",
                  buf2,buf1,
                  (double) input[PHI].INT_EPS,(double) input[THETA].INT_EPS,
                  input[PHI].JMAX,input[THETA].JMAX,
                  input[PHI].K,input[THETA].K,
	          (double) input[PHI].min,(double) input[THETA].min,
	          (double) input[PHI].max,(double) input[THETA].max);
  fprintf(file,"\n\nOuter-Loop\tInner Loop\n");
  fflush(file);

  error=outer_qromb(ss);    /* main calculation */

  /* scale the result */
  if (dtheta==0) scale=1;
  else scale= dtheta;       /* theta is really cos(theta) */
  if (dphi!=0) scale *= dphi;
  scale=1/scale;
  for(comp=0;comp<dim;comp++) ss[comp]*=scale;

  /* finalize log */
  if (error==0) {
    if (no_convergence==0) {
      printf("All inner integrations converged\n"\
    	     "The outer integration converged\n");
      fprintf(file,"All inner integrations converged\n"\
		      "The outer integration converged\n");
    }
    else {
      printf("%d inner integrations did not converge.\n"\
             "The outer integration converged\n",no_convergence);
      fprintf(file,"%d inner integrations did not converge.\n"\
              "The outer integration converged\n",no_convergence);
    }
  }
  else {
    if (no_convergence==0) {
      printf("Only the outer integration did not converge \n"\
	     "It reached d=%lg\n",error);
      fprintf(file,"Only the outer integration did not converge \n"\
	      "It reached d=%lg\n",error);
    }
    else {
      printf("%d inner integrations did not converge.\n"\
	     "The outer integration did not converge\n"\
	     "The outer integration reached d=%lg\n",no_convergence,error);
      fprintf(file,
	      "%d inner integrations did not converge.\n"\
	      "The outer integration did not converge\n"\
	      "The outer integration reached d=%lg\n",no_convergence,error);

    }
  }
  printf("In total %d evaluations were used\n",N_tot_eval);
  fprintf(file,"In total %d evaluations were used\n",N_tot_eval);
  fclose(file);

  free_all();  /* free all memory */
}


