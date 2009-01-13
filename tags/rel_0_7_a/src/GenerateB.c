/* FILE: GenerateB.c
 * AUTH: Alfons Hoekstra
 * DESCR: generate a incident beam
 *        original by A. Hoekstra, rewritten by Grimminck.
 *        plane wave and buggy type by Alfons Hoekstra, others by
 *        Michel Grimminck.
 *  
 *        Currently is developed by Maxim Yurkin
 *
 *        included incident propagation vector,
 *        but only for plane wave. Others work correctly only for default incidence
 */
#include <math.h>
#include <stdio.h>
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "debug.h"

void GenerateB (char which,		/* x - or y polarized incident light   */
                doublecomplex *b)	/* the b vector for the incident field */

{ 
  int i;
  
  extern double prop[3],incPolX[3],incPolY[3]; /*incident propagation and 
						basic vectors for polarization*/
  extern int shape;
  extern double *DipoleCoord;
  extern double WaveNum;
  extern double beam_w0,beam_x0,beam_y0,beam_z0;
  extern int beamtype;

  int p1,p2,p3;
  doublecomplex psi0, Q;
  double l = WaveNum * beam_w0 * beam_w0;       /* spreading length */
  double hplus, Qdenom, xsidiff;
  double x, y, z;
  double co, si, ex, tmpArg;
  doublecomplex ctemp;  
  double *incPol;	  	 /*used incident polarization*/
  
  /* polarisation is in the p1 direction - old version; new - choose from one of the basic 
  (latter only for plane wave, for now) */
  if (which=='Y') {p1=1; p2=0; p3=2; incPol=incPolY;}
  if (which=='X') {p1=0; p2=1; p3=2; incPol=incPolX;}
  if (beamtype==PLANE)
    for (i=0;i<local_nvoid_Ndip;i++) {
      tmpArg=WaveNum*DotProd(DipoleCoord+3*i,prop);
      ctemp[re]=cos(tmpArg);
      ctemp[im]=sin(tmpArg);
      cScalMultRVec(incPol,ctemp,b+3*i); /* b[i]=ctemp*incPol */
    } 
  else if (beamtype==LMINUS) {
    double s;

    s=beam_w0/l;
    printz("beam confinement factor s:%f\n",s);
    for (i=0;i<local_nvoid_Ndip;i++) {
      x = DipoleCoord[3*i];
      y = DipoleCoord[3*i+1];
      z = DipoleCoord[3*i+2];

      hplus = ((x-beam_x0)*(x-beam_x0)+(y-beam_y0)*(y-beam_y0))/(beam_w0*beam_w0);

      /* calculate Q */
      xsidiff = -(z - beam_z0)/l;
      Qdenom = 1.0 + 4 * xsidiff * xsidiff;
      Q[re] = 2.0 * xsidiff / Qdenom;
      Q[im] = -1.0 / Qdenom;

      /* calculate psi0 */
      co = cos(Q[re]*hplus);
      si = sin(Q[re]*hplus);
      ex = exp(Q[im]*hplus);
      psi0[re] = (-Q[im]*co + Q[re]*si)*ex;
      psi0[im] = (Q[re]*co + Q[im]*si)*ex; 
      
      /* calculate the fields */
      co = cos(WaveNum*(z-beam_z0));
      si = sin(WaveNum*(z-beam_z0));
      b[3*i+p1][re] = psi0[re]*co - psi0[im]*si;
      b[3*i+p1][im] = psi0[im]*co + psi0[re]*si;
      b[3*i+p2][re] = 0.0;
      b[3*i+p2][im] = 0.0;
      b[3*i+p3][re] = 0.0;
      b[3*i+p3][im] = 0.0;
    }
  }
  else if (beamtype==BUGGY) {
    double s;
    s=beam_w0/l;
    printz("beam confinement factor s:%f\n",s);
    for (i=0;i<local_nvoid_Ndip;i++) {
      x = DipoleCoord[3*i];
      y = DipoleCoord[3*i+1];
      z = DipoleCoord[3*i+2];

      hplus = ((x-beam_x0)*(x-beam_x0)+(y-beam_y0)*(y-beam_y0))/(beam_w0*beam_w0);
      
      /* calculate Q */
      xsidiff = (z - beam_z0)/l;
      Qdenom = 1.0 + 4 * xsidiff * xsidiff;
      Q[re] = 2.0 * xsidiff / Qdenom;
      Q[im] = -1.0 / Qdenom;
      
      /* calculate psi0 */
      co = cos(Q[re]*hplus);
      si = sin(Q[re]*hplus);
      ex = exp(Q[im]*hplus);
      psi0[re] = (-Q[im]*co + Q[re]*si)*ex;
      psi0[im] = (Q[re]*co + Q[im]*si)*ex; 
      
      /* calculate the fields */
      co = cos(WaveNum*(z-beam_z0));
      si = sin(WaveNum*(z-beam_z0));
      b[3*i+p1][re] = psi0[re]*co - psi0[im]*si;
      b[3*i+p1][im] = psi0[im]*co + psi0[re]*si;
      b[3*i+p2][re] = 0.0;
      b[3*i+p2][im] = 0.0;
      if (which=='Y') {
	b[3*i+p3][re] = -2.0*(y-beam_y0)*(Q[re]*b[3*i+p1][re] - Q[im]*b[3*i+p1][im])/l;
	b[3*i+p3][im] = -2.0*(y-beam_y0)*(Q[im]*b[3*i+p1][re] + Q[re]*b[3*i+p1][im])/l;
      }
      else {
	b[3*i+p3][re] = -2.0*(x-beam_x0)*(Q[re]*b[3*i+p1][re] - Q[im]*b[3*i+p1][im])/l;
	b[3*i+p3][im] = -2.0*(x-beam_x0)*(Q[im]*b[3*i+p1][re] + Q[re]*b[3*i+p1][im])/l;
      }
    }
  }
  else if (beamtype==DAVIS1 || beamtype==BARTON1 || beamtype==BARTON3
	   || beamtype==DAVIS3 || beamtype==BARTON5) {
    double s;
    doublecomplex ex,ey,ez;
    double zeta,xsi,eta,nu;
    doublecomplex d,d2,d3,d4,d5,d6;
    doublecomplex temp,temp2,temp3;

    s=beam_w0/l;
    printz("beam confinement factor s:%f\n",s);
    for (i=0;i<local_nvoid_Ndip;i++) {
      if (which=='X') {
        x = DipoleCoord[3*i]-beam_x0;
        y = DipoleCoord[3*i+1]-beam_y0;
      } else {
        x = DipoleCoord[3*i+1]-beam_y0;
        y = DipoleCoord[3*i]-beam_x0;
      }

      z = -(DipoleCoord[3*i+2]-beam_z0); /* minus to give beam the right direction */
      xsi=x/beam_w0;
      eta=y/beam_w0;
      zeta=s*z/beam_w0;
      nu=xsi*xsi+eta*eta;
      d[re]=1.0/(1.0+4*zeta*zeta); d[im]=2.0*zeta/(1.0+4*zeta*zeta);
      
      if (beamtype==DAVIS1 || beamtype==BARTON1) {
	ex[re]=1; ex[im]=0;
	ey[re]=0; ey[im]=0;
	ez[re]=s; ez[im]=0;
      }
      else if (beamtype==BARTON3) {
	cSquare(d,d2);
	cMult(d2,d,d3);
	ex[re]=1+s*s*(2*xsi*xsi*d2[re]+nu*d2[re]-nu*nu*d3[re]);
	ex[im]=s*s*(2*xsi*xsi*d2[im]+nu*d2[im]-nu*nu*d3[im]);
	ey[re]=s*s*2*xsi*eta*d2[re];
	ey[im]=s*s*2*xsi*eta*d2[im];
	ez[re]=s+s*s*s*(3*nu*d2[re]-nu*nu*d3[re]);
	ez[im]=s*s*s*(3*nu*d2[im]-nu*nu*d3[im]);
      }
      else if (beamtype==DAVIS3) {
	cSquare(d,d2);
	cMult(d2,d,d3);
	ex[re]=1+s*s*(4*xsi*xsi*d2[re]-nu*nu*d3[re]);
	ex[im]=s*s*(4*xsi*xsi*d2[im]-nu*nu*d3[im]);
	ey[re]=s*s*4*xsi*eta*d2[re];
	ey[im]=s*s*4*xsi*eta*d2[im];
	ez[re]=s+s*s*s*(-2*d[re]+4*nu*d2[re]-nu*nu*d3[re]);
	ez[im]=s*s*s*(-2*d[im]+4*nu*d2[im]-nu*nu*d3[im]);
      }
      else if (beamtype==BARTON5) {
	cSquare(d,d2);
	cMult(d2,d,d3);
	cMult(d3,d,d4);
	cMult(d4,d,d5);
	cMult(d5,d,d6);
	ex[re]=1+s*s*(2*xsi*xsi*d2[re]+nu*d2[re]-nu*nu*d3[re]) +
	           s*s*s*s*(2*nu*nu*d4[re]+8*xsi*xsi*nu*d4[re]-3*nu*nu*nu*d5[re] -
		   2*xsi*xsi*nu*nu*d5[re]+.5*nu*nu*nu*nu*d6[re]);
	ex[im]=s*s*(2*xsi*xsi*d2[im]+nu*d2[im]-nu*nu*d3[im]) +
	           s*s*s*s*(2*nu*nu*d4[im]+8*xsi*xsi*nu*d4[im]-3*nu*nu*nu*d5[im]-
		   2*xsi*xsi*nu*nu*d5[im]+.5*nu*nu*nu*nu*d6[im]);
	ey[re]=s*s*2*xsi*eta*d2[re] +
	           s*s*s*s*(2*xsi*eta*d2[re])*(4*nu*d2[re]-nu*nu*d3[re]);
	ey[im]=s*s*2*xsi*eta*d2[im] +
	           s*s*s*s*(2*xsi*eta*d2[im])*(4*nu*d2[im]-nu*nu*d3[im]);
	ez[re]=s+s*s*s*(3*nu*d2[re]-nu*nu*d3[re]) +
	           s*s*s*s*s*(10*nu*nu*d4[re]-5*nu*nu*nu*d5[re]+.5*nu*nu*nu*nu*d6[re]);
	ez[im]=s*s*s*(3*nu*d2[im]-nu*nu*d3[im]) +
	           s*s*s*s*s*(10*nu*nu*d4[im]-5*nu*nu*nu*d5[im]+.5*nu*nu*nu*nu*d6[im]);
      }
      temp[re]=cos(-WaveNum*z);     /* temp=exp(-i*k*z) */
      temp[im]=sin(-WaveNum*z);     
      temp2[re]=temp[re]*d[re]-temp[im]*d[im];      /* temp2=exp(-i*k*z)*D */
      temp2[im]=temp[re]*d[im]+temp[im]*d[re];
      temp[re]=exp(-nu*d[re])*cos(-nu*d[im]);      /* temp=exp(-nu*d) */
      temp[im]=exp(-nu*d[re])*sin(-nu*d[im]);
      temp3[re]=temp[re]*temp2[re]-temp[im]*temp2[im]; /* temp3=(i*k*z)*D*exp(-nu*d) */
      temp3[im]=temp[re]*temp2[im]+temp[im]*temp2[re];
      
      temp[re]=d[re]*ez[re]-d[im]*ez[im];              /* temp=d*ez */
      temp[im]=d[re]*ez[im]+d[im]*ez[re];

      temp[re]*=2*xsi;                         /* temp=2*xsi*d*ez */
      temp[im]*=2*xsi;

      temp2[re]=-temp[im]*temp3[re]-temp[re]*temp3[im];  /* temp2=2*i*xsi*d*ez*temp3 */
      temp2[im]=temp[re]*temp3[re]-temp[im]*temp3[im];

      /* fill in E field */
      b[3*i+p1][re] =temp3[re]*ex[re]-temp3[im]*ex[im];
      b[3*i+p1][im] =temp3[re]*ex[im]+temp3[im]*ex[re];
      b[3*i+p2][re] =temp3[re]*ey[re]-temp3[im]*ey[im];
      b[3*i+p2][im] =temp3[re]*ey[im]+temp3[im]*ey[re];
      b[3*i+p3][re] =-temp2[re];
      b[3*i+p3][im] =-temp2[im];
      
      /* correction for carvature of the field */
      /*corr=1+.63*.63*(xsi*xsi/12-1/24)/beam_w0/beam_w0;
      b[3*i+p1][re]*=corr; b[3*i+p1][im]*=corr;
      b[3*i+p2][re]*=corr; b[3*i+p2][im]*=corr;
      b[3*i+p3][re]*=corr; b[3*i+p3][im]*=corr;*/

      /* show intensity of incident field along x-axis */
      if (shape==LINE && x>0) printf("%f %.8g\n",x,
			  b[3*i+p1][re]*b[3*i+p1][re]+b[3*i+p1][im]*b[3*i+p1][im]+
			  b[3*i+p2][re]*b[3*i+p2][re]+b[3*i+p2][im]*b[3*i+p2][im]+
			  b[3*i+p3][re]*b[3*i+p3][re]+b[3*i+p3][im]*b[3*i+p3][im]);
    }
  }
}
