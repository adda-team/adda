/* FILE: CGNR.c
 * AUTH: Alfons Hoekstra
 * DATE: -
 */

/* the crude harness for the CGNR implementation, again, a complete
   mess as far as the variables are concerned
*/

#include <stdio.h>
#include "cmplx.h"
#include "types.h"
#include "comm.h"
#include "debug.h"

extern dcomplex *x;		/* defined inin */
extern dcomplex *r;		/* defined in main */
extern dcomplex *p;		/* defined in main */
extern dcomplex *buffer;	/* defined in main */
extern FILE *logfile;           /* defined in main */

extern int nDip;		/* defined in calculator */
extern int nRing;		/* defined in calculator */
extern int RingId;		/* defined in calculator */
extern int nlocalDip;		/* defined in calculator */
extern dcomplex *Avecbuffer;	/* defined in calculator */
extern double **DipoleCoord;	/* defined in calculator */  
extern double WaveNum;		/* defined in calcularor */

extern double inprodB;       	/* defined in CalculateE */
extern double epsB;          	/* defined in CalculateE */
extern double inprodR;       	/* defined in CalculateE */

extern int Nmat;
extern int *material;

extern void MatVec_nim ( );
extern void MatVecHer_nim ( );

int
CGNR (int max_iter,	/* maximum number of allowed iterations */
      int max_count)    /* maximum number of allowed residual increases */

{
  int i;
  int count = 0, counter = 0;	/* two counters to control the convergence */
  double inprodRplus1;		/* inner product of rk+1 */
  double alpha, numeratorAlpha, denumeratorAlpha;
  double beta, numeratorBeta;
  FILE *teststop;
  char sbuffer[200];
  extern char directory[200];

  extern int maxiter;
  unsigned int tstart, tstop, tstart2, tstop2;
  extern unsigned long Timing_OneIter,
                       Timing_OneIterCalc,
                       Timing_OneIterComm;

  if (inprodR > epsB) {    	/* x0 is not good enough */
    
    /* calculate numerator of first alpha */
    numeratorAlpha = 0.0;
    for (i = 0; i < local_Ndip; ++i) 
      if(material[i]<Nmat-1) {
	numeratorAlpha += p[3*i].r * p[3*i].r + p[3*i].i * p[3*i].i;
	numeratorAlpha += p[3*i+1].r * p[3*i+1].r + p[3*i+1].i * p[3*i+1].i;
	numeratorAlpha += p[3*i+2].r * p[3*i+2].r + p[3*i+2].i * p[3*i+2].i;
      }
    
    my_inner_product(&numeratorAlpha); /* ACCUMULATE NUMERATOR OF ALPHA */
    
    /* calculate denominator of alpha, and alpha itself */
    MatVecAndInp_nim (p, Avecbuffer, &denumeratorAlpha,
		      nlocalDip, nDip, RingId, nRing, 
		      WaveNum);
    
    alpha = numeratorAlpha / denumeratorAlpha;
    
    /* calculate xk+1 */
    for (i = 0; i < local_Ndip; ++i)  {
      if(material[i]<Nmat-1) {
	x[3*i].r += alpha * p[3*i].r;
	x[3*i].i += alpha * p[3*i].i;
	x[3*i+1].r += alpha * p[3*i+1].r;
	x[3*i+1].i += alpha * p[3*i+1].i;
	x[3*i+2].r += alpha * p[3*i+2].r;
	x[3*i+2].i += alpha * p[3*i+2].i;
      }
    }
    
    
    /* calculate rk+1 and |rk+1| */
    inprodRplus1 = 0.0;
    for (i = 0; i < local_Ndip; ++i) {
      if(material[i]<Nmat-1) {
	r[3*i].r -= alpha * Avecbuffer[3*i].r;
	r[3*i].i -= alpha * Avecbuffer[3*i].i;
	inprodRplus1 += r[3*i].r * r[3*i].r + r[3*i].i * r[3*i].i;
	r[3*i+1].r -= alpha * Avecbuffer[3*i+1].r;
	r[3*i+1].i -= alpha * Avecbuffer[3*i+1].i;
	inprodRplus1 += r[3*i+1].r * r[3*i+1].r + r[3*i+1].i * r[3*i+1].i;
	r[3*i+2].r -= alpha * Avecbuffer[3*i+2].r;
	r[3*i+2].i -= alpha * Avecbuffer[3*i+2].i;
	inprodRplus1 += r[3*i+2].r * r[3*i+2].r + r[3*i+2].i * r[3*i+2].i;
      }
    }
    
    tstart2=extime();
    my_inner_product(&inprodRplus1);	/* ACCUMULATE INPRODUCT_PLUS_1 */
    tstop2=extime();
    Timing_OneIterComm += tstop2-tstart2;
    
    /* cheque the innerproduct for ongoing convergence */
    if (inprodRplus1 < inprodR) 
      inprodR = inprodRplus1;
    else
      counter = 1;
    
    count = 1;
    
    /* PRINT ALPHA AND |R| */
    if (RingId == 0) {
      fprintz(logfile,"\nalpha01   = %1.12e\n",alpha);
      fprintz(logfile,"r01       = %1.7e\n", inprodRplus1);
      printz("r01       = %1.7e\n", inprodRplus1);
      fflush(stdout); fflush(logfile);
    }

    /* start the iteration */
    while ((inprodR >= epsB) && (count <= max_iter) && count<maxiter && 
	   (counter <= max_count) ) {
      strcpy(sbuffer,directory);
      strcat(sbuffer,"/stop");
      Timing_OneIterComm  = 0;	/* time spent communicating */
      
      tstart = extime();
      
      /* start of the iterations */
      
      /* calculation of beta */
      MatVecAndInpHer_nim (r, Avecbuffer, &numeratorBeta,
			   nlocalDip, nDip, RingId, nRing,
			   WaveNum);
      
      
      beta = numeratorBeta / numeratorAlpha;
      
      /* calculate pk+1 */
      for (i = 0; i < local_Ndip; ++i)  {
        if(material[i]<Nmat-1){
	  p[3*i].r = beta * p[3*i].r + Avecbuffer[3*i].r;
	  p[3*i].i = beta * p[3*i].i + Avecbuffer[3*i].i;
	  p[3*i+1].r = beta * p[3*i+1].r + Avecbuffer[3*i+1].r;
	  p[3*i+1].i = beta * p[3*i+1].i + Avecbuffer[3*i+1].i;
	  p[3*i+2].r = beta * p[3*i+2].r + Avecbuffer[3*i+2].r;
	  p[3*i+2].i = beta * p[3*i+2].i + Avecbuffer[3*i+2].i;
	}
      }
      
      /* calculate alpha */
      MatVecAndInp_nim (p, Avecbuffer, &denumeratorAlpha,
			nlocalDip, nDip, RingId, nRing,
			WaveNum);
      
      numeratorAlpha = numeratorBeta;
      alpha = numeratorAlpha / denumeratorAlpha;
      
      /* calculate xk+1 */
      for (i = 0; i < local_Ndip; ++i)  {
        if(material[i]<Nmat-1) {
	  x[3*i].r += alpha * p[3*i].r;
	  x[3*i].i += alpha * p[3*i].i;
	  x[3*i+1].r += alpha * p[3*i+1].r;
	  x[3*i+1].i += alpha * p[3*i+1].i;
	  x[3*i+2].r += alpha * p[3*i+2].r;
	  x[3*i+2].i += alpha * p[3*i+2].i;
	}
      }
      
      
      /* calculate rk+1 and |rk+1| */
      inprodRplus1 = 0.0;
      for (i = 0; i < local_Ndip; ++i) {
        if(material[i]<Nmat-1) {
	  r[3*i].r -= alpha * Avecbuffer[3*i].r;
	  r[3*i].i -= alpha * Avecbuffer[3*i].i;
	  inprodRplus1 += r[3*i].r * r[3*i].r + r[3*i].i * r[3*i].i;
	  r[3*i+1].r -= alpha * Avecbuffer[3*i+1].r;
	  r[3*i+1].i -= alpha * Avecbuffer[3*i+1].i;
	  inprodRplus1 += r[3*i+1].r * r[3*i+1].r + r[3*i+1].i * r[3*i+1].i;
	  r[3*i+2].r -= alpha * Avecbuffer[3*i+2].r;
	  r[3*i+2].i -= alpha * Avecbuffer[3*i+2].i;
	  inprodRplus1 += r[3*i+2].r * r[3*i+2].r + r[3*i+2].i * r[3*i+2].i;
	}
      }
      
      tstart2=extime();
      my_inner_product(&inprodRplus1); /* ACCUMULATE INPRODUKT RK+1 */
      tstop2=extime();
      Timing_OneIterComm += tstop2-tstart2;
      
      /* cheque the innerproduct for ongoing convergence */
      ++count;
      if (inprodRplus1 < inprodR) {
	inprodR = inprodRplus1;
	counter = 0;
      }
      else
	++counter;
      
      
      tstop = extime ();
      if (RingId == 0) {	/* PRINT BETA, ALPHA AND |R| TO FILE */
        fprintz(logfile,"beta%02d    = %1.12e\n\n",count-1, beta);
        fprintz(logfile,"r%02d       = %1.12e\n",count,inprodRplus1);
	printz("r%02d       = %1.12e\n",count,inprodRplus1);
        fprintz(logfile,"alpha%02d   = %1.12e\n",count, alpha);
        fflush(stdout); fflush(logfile);
      }
      
    } /* end of the big while loop */
  } /* end of the big if construct */
  
  /* so, let us check for convergence */
  
  Timing_OneIter = tstop - tstart;
  Timing_OneIterCalc = Timing_OneIter - Timing_OneIterComm;
  
  if (count > max_iter) {
    /* no convergence, too many iterations */
    return (-1);
  } else if (counter > max_count) {
    /* no convergence, residuals increase too much */
    return (-2);
  } else {
    /* convergence */
    return (0);
  }
}
