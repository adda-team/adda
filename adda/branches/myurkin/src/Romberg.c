/* FILE : Romberg.c
 * Auth : Martijn Frijlink
 * Hist : 1) Created 1998
 * Descr: 2D-integration routine, based on the qromb-routine from 
 *        Numerical Recipes. The main differences from the original
 *        code are :
 *        1)The function-values are not calculated on demand, but before
 *          the outer_qromb is actually callled. So instead of determining
 *          the proper argument-values, the proper array-index must be 
 *          determined. The array is not imported directly into to the 
 *          routine, but accessed through a function, supplied as an 
 *          argument to outer_qromb.
 *        2)From 1D-integration to 2D-integration. The function-values for
 *          outer_qromb are themselves the result of an integration
 *          performed by inner_qromb.
 *        The types declared in types.h deserve some attention :
 *        - Rvector : originally this routine was designed to integrate
 *             a 3D function of two variables AND a scalar function of 
 *             two variables. In order to hold both scalar and 3D vector
 *             values, Rvector was made a pointer.
 *        - Parms_1D : this structure holds the parameters for integration
 *             in one dimension. They must be set outside of the 
 *             outer_qromb-routine.
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "cmplx.h"
#include "types.h"
#include "Romberg.h"

int DIM;
int N_eval;
int no_convergence;
FILE *log_int;

void nrerror(char error_text[])
{
  printf("Runtime error: %s.\n",error_text);
  printf("Now aborting the program\n");
  fprintf(log_int,error_text);
  fclose(log_int);
  exit(1);
}

REAL *alloc_Rvector(int nl,int nh)
{
  REAL *v;

  v = (REAL *)malloc((unsigned)(nh-nl+1)*sizeof(REAL));
  if (!v) nrerror("memory allocation failure in Rvector()");

  return v-nl;
}

void freeRv(REAL *v,int nl,int nh)
{
  free((v+nl));
}





void polint(/* Inter- or extrapolate the n function-values ya[], at
	     * corresponding gridpoints xa[], to y, at gridpoint x 
	     * with Neville's algorithm. The estimated error in y
	     * will be stored in dy. */
	    REAL xa[],REAL ya[],int n,
	    REAL x,REAL *y,REAL *dy)
{
  int 
    i,m,     /* Loop indices */
    ns;    
  REAL
    dif,     /* distance(xa[i],x) */
    *c,*d,   /* polynomial corrections */
    ho,hp,   /* dummy variables */
    den,w;   /* idem */

  c = alloc_Rvector(1,n);
  d = alloc_Rvector(1,n);

  /* Find the index of the gridpoint closest to x */
  for (i=1, dif=fabs(x-xa[1]);
       i<=n; ++i)
    {
      REAL dift; /* temporary dif */

      if ( (dift=fabs(x-xa[i])) < dif )
	{
	  ns=i;
	  dif=dift;
	}

      c[i] = ya[i];
      d[i] = ya[i];
    }

  /* Search-algorithm */
  *y=ya[ns--];
  for (m=1;m<n;m++)
    {
      for (i=1;i<=n-m;++i)
	/* Calculate corrections for polynomial of one degree higher */
	{
	  ho = xa[i]-x;
	  hp = xa[i+m]-x;
	  w = c[i+1] - d[i];

	  if ( (den=ho-hp) == 0.0)
	    /* Catch for identical (to within roundoff) xa's.
	     * Preventing division by zero */
	    nrerror("Error in routine POLINT");
	  den  = w/den;
	  d[i] = hp*den;
	  c[i] = ho*den;
	}

      /* update *y and *dy */ 
      *dy = (2*ns < (n-m)) ? c[ns+1] : d[ns--];
      *y += *dy;
    }

  freeRv(c,1,n);
  freeRv(d,1,n);
}




Rvector outer_qromb(/* Integrate 2D FUNC with Romberg's method
		     * according to input's parameters. Argument 
		     * dim gives the number of components of Rvector.
		     * Consistency between FUNC and dim is the user's
		     * responsibilty */
		    Parms_1D input[Narg],
		    Rvector (*FUNC)(int theta,int phi),
		    int dim,
                    char logfile[])
{
  int 
    j,
    comp,
    N_tot_eval=0;   /* total number of integrand-evaluations */
  char
    *io_control,
    buffer[220];
  Rvector
    *s,             /* Intermediate results */
    ss,             /* final result */
    dss;            /* extrapolation error in final result */
  REAL 
    abs_ss,abs_dss, /* ...-norms of ss and dss */
    **s_pi,         /* Dummy carrying transposed ss[] */
    *h;             /* squared stepsize */
  double
    error;

  DIM = dim;
  log_int=fopen(logfile,"w");

  fprintf(log_int,"\t\t\tPHI\tTHETA\n"\
                  "EPS\t\t\t%lg\t%lg\n"\
                  "Maximal number of\nrefinement-stages\t%d\t%d\n"\
                  "Number of evaluations\nfor an extrapolation\t%d\t%d\n"\
                  "lower boundary\t\t%lg\t%lg\n"\
                  "upper boundary\t\t%lg\t%lg\n",
                  (double) input[PHI].INT_EPS,(double) input[THETA].INT_EPS,
                  input[PHI].JMAX,input[THETA].JMAX,
                  input[PHI].K,input[THETA].K,
	          (double) input[PHI].min,(double) input[THETA].min,
	          (double) input[PHI].max,(double) input[THETA].max);

  /* memory-allocation */
  s = (Rvector *) malloc((input[THETA].JMAX+2)*sizeof(Rvector));
  for (j=0;j<input[THETA].JMAX+2;++j)
    s[j].x = (REAL *) malloc(DIM*sizeof(REAL));
  ss.x = (REAL *) malloc(DIM*sizeof(REAL));
  dss.x = (REAL *) malloc(DIM*sizeof(REAL));
  s_pi = (REAL **) malloc(DIM*sizeof(REAL *));
  for (comp=0;comp<DIM;++comp)
    s_pi[comp] = (REAL *) malloc((input[THETA].JMAX+2)*sizeof(REAL));
  h = (REAL *) malloc((input[THETA].JMAX+2)*sizeof(REAL));

  fprintf(log_int,"\n\nOuter-Loop\tInner Loop\n");  
  no_convergence = 0;
  h[1] = 1.0;
  for (j=1;j<=input[THETA].JMAX;++j)
    {
      N_eval = 0;
      outer_trapzd(input,FUNC,&s[j],j);
      for (comp=0;comp<DIM;++comp)
	s_pi[comp][j] = s[j].x[comp];
      fprintf(log_int,"%d\t\t%d integrand-values were used.\n",j,N_eval);
      N_tot_eval += N_eval;

      if (j >= input[THETA].K)
	{
	  for (comp=0, abs_ss=0, abs_dss=0;
	       comp<DIM;
	       ++comp)
	    {
	      polint(&h[j-input[THETA].K],&s_pi[comp][j-input[THETA].K],
		     input[THETA].K,0.0,&(ss.x[comp]),&(dss.x[comp]));
	      abs_ss += fabs(ss.x[comp]);
	      abs_dss += fabs(dss.x[comp]);
	    }

	  if ( abs_dss < input[THETA].INT_EPS*abs_ss)
	    {
	      if (!no_convergence)
		{
		  printf("All inner integrations converged\n"\
			 "The outer integration converged\n"\
			 "In total %d evaluations were used\n",
			 N_tot_eval);
		  fprintf(log_int,
			  "All inner integrations converged\n"\
			  "The outer integration converged\n"\
			  "In total %d evaluations were used\n",
			  N_tot_eval);
		}
	      else
		{
		  printf("%d inner integrations did not converge.\n"\
			 "In total %d evaluations were used\n",
			 no_convergence,
			 N_tot_eval);
		  fprintf(log_int,
			  "%d inner integrations did not converge.\n"\
			  "In total %d evaluations were used\n",
			  no_convergence,
			  N_tot_eval);
		}

	      fclose(log_int);
	      return ss;
	    }
	}

      for (comp=0;comp<DIM;++comp)
	s[j+1].x[comp] = s[j].x[comp];
      h[j+1]   = 0.25*h[j];
    }
  
  /* Calculate percentual error */
  error = (double) abs_dss/abs_ss;

  if (no_convergence == 0)
    {
      printf("Only the outer integration did not converge \n"\
	     "It reached d=%lg\n",
	     error);
      fprintf(log_int,
	      "Only the outer integration did not converge \n"\
	      "It reached d=%lg\n",
	      error);
    }
  else
    {
      printf("%d inner integrations did not converge.\n"\
	     "The outer integration did not converge\n"\
	     "The outer integration reached d=%lg\n",
	     no_convergence,error);
      fprintf(log_int,
	      "%d inner integrations did not converge.\n"\
	      "The outer integration did not converge\n"\
	      "The outer integration reached d=%lg\n",
	      no_convergence,error);
    }
  printf("In total %d evaluations were used\n",
	  N_tot_eval);
  fprintf(log_int,"In total %d evaluations were used\n",
	  N_tot_eval);

  fclose(log_int);
  return ss;
}

void outer_trapzd(/* Calculate n'th refinement for the outer integration
		   * of FUNC */
		  Parms_1D input[Narg],
		  Rvector (*FUNC)(int theta,int phi),
		  Rvector *s,int n)
{
  if (n == 1)
    /* Initialise the algorithm */
    {
      int
	comp;
      Rvector
	dummy[2];

      /* memory allocation */
      dummy[0].x = (REAL *) malloc(DIM*sizeof(REAL));
      dummy[1].x = (REAL *) malloc(DIM*sizeof(REAL));

      dummy[0] = inner_qromb(input,FUNC,0);
      dummy[1] = inner_qromb(input,FUNC,input[THETA].Grid_size-1);

      for (comp=0;comp<DIM;++comp)
	s->x[comp] = 0.5*(input[THETA].max-input[THETA].min)*
	  (dummy[0].x[comp] + dummy[1].x[comp]);
   }
  else
    /* refinement stage n */
    {
      REAL
	tnm;       /* number of intervals per refinement stage */
      int
	j,
	step,
	comp;
      Rvector
	sum,       /* refinement */
	dummy;

      /* memory allocation */
      dummy.x = (REAL *) malloc(DIM*sizeof(REAL));
      sum.x = (REAL *) calloc(DIM,sizeof(REAL));

      step = (input[THETA].Grid_size-1) >> n-1;

      for (j = step >> 1;  /* initialise refinement */
	   j < input[THETA].Grid_size;
	   j += step)
	{
	  dummy = inner_qromb(input,FUNC,j);
	  for (comp=0;comp<DIM;++comp)
	    sum.x[comp] += dummy.x[comp];
	}

      /* add refinement */
      tnm = pow(2,n-1);
      for (comp=0;comp<DIM;++comp)
	s->x[comp] = 0.5*(s->x[comp] + 
			  (input[THETA].max-input[THETA].min)*sum.x[comp]/tnm);
    }
}

Rvector inner_qromb(/* Integrate FUNC for fixed theta=fixed */
		    Parms_1D input[Narg],
		    Rvector (*FUNC)(int theta,int phi),
		    int fixed)
{
  int
    j,
    comp;
  Rvector
    *s,             /* Intermediate results */
    ss,             /* final result */
    dss;            /* extrapolation error in final result */
  REAL 
    abs_ss,abs_dss, /* ...-norms of ss and dss */
    **s_pi,         /* Dummy carrying transposed ss[] */
    *h;             /* squared stepsize */

  /* memory-allocation */
  s = (Rvector *) malloc((input[PHI].JMAX+2)*sizeof(Rvector));
  for (j=0;j<input[PHI].JMAX+2;++j)
    s[j].x = (REAL *) malloc(DIM*sizeof(REAL));
  ss.x = (REAL *) malloc(DIM*sizeof(REAL));
  dss.x = (REAL *) malloc(DIM*sizeof(REAL));
  
  s_pi = (REAL **) malloc(DIM*sizeof(REAL *));
  for (comp=0;comp<DIM;++comp)
    s_pi[comp] = (REAL *) malloc((input[PHI].JMAX+2)*sizeof(REAL));
  h = (REAL *) malloc((input[PHI].JMAX+2)*sizeof(REAL));

  h[1]=1.0;
  for (j=1;j<=input[PHI].JMAX;++j)
    {
      inner_trapzd(input,FUNC,fixed,&s[j],j);
      for (comp=0;comp<DIM;++comp)
	s_pi[comp][j] = s[j].x[comp];

      if (j >= input[PHI].K)
	{
	  for (comp=0, abs_ss=0, abs_dss=0;
	       comp<DIM;
	       ++comp)
	    {
	      polint(&h[j-input[PHI].K],&s_pi[comp][j-input[PHI].K],
		     input[PHI].K,0.0,&(ss.x[comp]),&(dss.x[comp]));
	      abs_ss += fabs(ss.x[comp]);
	      abs_dss += fabs(dss.x[comp]);
	    }

	  if ( abs_dss < input[PHI].INT_EPS*abs_ss) 
	    {
	      return ss;
	    }
	}

      for (comp=0;comp<DIM;++comp)
	s[j+1].x[comp] = s[j].x[comp];
      h[j+1] = 0.25*h[j];
    }

  fprintf(log_int,"Inner_qromb converged to d=%lg for theta = %d\n",
	  (double) abs_dss/abs_ss,fixed);
  no_convergence++;
  
  return ss;
}

void inner_trapzd(/* Calculate the n'th refinement stage of the integration
		   * of FUNC over phi_min < phi < phi_max for 
		   * fixed theta = th_f */
		  Parms_1D input[Narg],
		  Rvector (*FUNC)(int theta,int phi),
		  int fixed,
		  Rvector *s,int n)
{
  if (n == 1)
    /* Initialise the algorithm */
    {
      int
	comp;
      Rvector
	dummy[2];

      /* memory allocation */
      dummy[0].x = (REAL *) malloc(DIM*sizeof(REAL));
      dummy[1].x = (REAL *) malloc(DIM*sizeof(REAL));

      dummy[0] = FUNC(fixed,0);
      dummy[1] = FUNC(fixed,input[PHI].Grid_size-1);
      N_eval += 2;

      for (comp=0;comp<DIM;++comp)
      s->x[comp] = 0.5*(input[PHI].max-input[PHI].min)*
	(dummy[0].x[comp]+dummy[1].x[comp]);
    }
  else
    /* refinement stage n */
    {
      int
	j,
	comp,
	step;
      REAL
	tnm;
      Rvector
	dummy,
	sum;     /* refinement */

 
      /* memory allocation */
      dummy.x = (REAL *) malloc(DIM*sizeof(REAL));
      sum.x = (REAL *) calloc(DIM,sizeof(REAL));

      step = (input[PHI].Grid_size-1) >> n-1;

      for (j = step >> 1;  /* initialise refinement */
	   j < input[PHI].Grid_size;
	   j += step)
	{
	  dummy = FUNC(fixed,j);
	  ++N_eval;
	  for (comp=0;comp<DIM;++comp)
	    sum.x[comp] += dummy.x[comp];
	}

      /* add refinement */
      tnm = pow(2,n-1);
      for (comp=0;comp<DIM;++comp)
	s->x[comp] = 0.5*(s->x[comp] + 
			  (input[PHI].max-input[PHI].min)*sum.x[comp]/tnm);
    }
}
