#include <stdio.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "types.h"
#include "comm.h"

void
calc_field (dcomplex *x,
	    dcomplex *ebuff,
            REAL **rdip,
	    double *n,
            double kk,
            dcomplex **f,
	    int nlocalDip
	    )
     
{
  /*  Near-optimal routine to compute the scattered fields at one specific
   *  angle.
   *
   *  Michel Grimminck 1995
   */
  
  double kr,rt00, rt01, rt02, rt11, rt12, rt22;
  dcomplex a;
  double sumrr[3],sumri[3],sumir[3],sumii[3];
  int k,l,m,j,jjj;
  double xp,yp,zp,temp;
  extern short int *position;
  extern double gridspaceX,gridspaceY,gridspaceZ,centreX,centreY,centreZ;
  extern int Nmat;
  extern int *material;
  extern doublecomplex cc[10];
  doublecomplex tbuff[3];
  
  /* calculate projection matrix */
  rt00=n[1]*n[1]+n[2]*n[2];
  rt11=n[0]*n[0]+n[2]*n[2];
  rt22=n[0]*n[0]+n[1]*n[1];
  rt01=-n[1]*n[0];
  rt02=-n[2]*n[0];
  rt12=-n[2]*n[1];
  
  for(k=0;k<3;k++) n[k]=-kk*n[k];
  
  for(m=0;m<Nmat;m++) { /* for each refraction index */
    for(k=0;k<3;k++) sumrr[k]=sumri[k]=sumir[k]=sumii[k]=0.0;
    
    for (j = local_d0; j < local_d1; ++j) 
      if (material[j-local_d0]==m) { /* for each dipole */
	/* calculate field in cartesian coordinates */
	jjj = 3 * (j-local_d0);
	
	/* calculate fase of field */
	
	kr=rdip[j-local_d0][0]*n[0]+rdip[j-local_d0][1]*n[1]+rdip[j-local_d0][2]*n[2];
	while (kr > 2 * PI)
	  kr -= 2 * PI;
      
	a.r = cos (kr);
	a.i = sin (kr);
	
	for(k=0;k<3;k++) {
	  sumrr[k]+=x[jjj+k].r*a.r; sumir[k]+=x[jjj+k].i*a.r;
	  sumri[k]+=x[jjj+k].r*a.i; sumii[k]+=x[jjj+k].i*a.i;
	}
      } /* end for j */
    f[0][0].r = rt00;  f[0][0].i = rt00;
    f[0][1].r = f[1][0].r = rt01;  f[0][1].i = f[1][0].i = rt01;
    f[0][2].r = f[2][0].r = rt02;  f[0][2].i = f[2][0].i = rt02;
    f[1][1].r = rt11;  f[1][1].i = rt11;
    f[1][2].r = f[2][1].r = rt12;  f[1][2].i = f[2][1].i = rt12;
    f[2][2].r = rt22;  f[2][2].i = rt22;
    
    for (k = 0; k < 3; ++k) {
      tbuff[k].r=tbuff[k].i=0.0;
      for (l = 0; l < 3; ++l) {
	tbuff[k].r += f[k][l].r * sumrr[l] -
	  f[k][l].i * sumii[l];
	tbuff[k].i += f[k][l].r * sumri[l] +
	  f[k][l].i * sumir[l];
      }
    }
    /* multiply with coupleconstant */
    for (k = 0; k < 3; ++k) {
      temp = tbuff[k].r;
      tbuff[k].r = temp * cc[m].r -
	tbuff[k].i * cc[m].i;
      tbuff[k].i = temp * cc[m].i +
	tbuff[k].i * cc[m].r;
    }
    /* add it to the E field */
    for(k = 0; k < 3; ++k) {
      ebuff[k].r+=tbuff[k].r*kk*kk;
      ebuff[k].i+=tbuff[k].i*kk*kk;
    }
  } /* end m */
}


