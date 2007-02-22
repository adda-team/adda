/* FILE: GenerateB.c
 * AUTH: Alfons Hoekstra
 * DESCR: generate a incident beam
 *        original by A. Hoekstra, rewritten by Grimminck.
 *        plane wave by Alfons Hoekstra, others by Michel Grimminck.
 *
 *        Currently is developed by Maxim Yurkin
 *
 *        included incident propagation vector,
 *        but only for plane wave. Others work correctly only for default incidence
 *
 *        Davis beams are based on:
 *        L. W. Davis, "Theory of electromagnetic beams," Phys. Rev. A 19, 1177-1179 (1979).
 *
 *	  Barton beams are based on:
 *        J. P. Barton and D. R. Alexander, "Fifth-order corrected electromagnetic-field
 *        components for a fundamental gaussian-beam," J.Appl.Phys. 66, 2800-2802 (1989).
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include "vars.h"
#include "cmplx.h"
#include "const.h"
#include "comm.h"

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in param.c */
extern const int beamtype;
extern const double beam_w0,beam_x0,beam_y0,beam_z0;

/*============================================================*/

void GenerateB (const char which,	/* x - or y polarized incident light   */
                doublecomplex *b)	/* the b vector for the incident field */
   /* generates incident beam at every dipole */
{ 
  int i;
  
  int p1,p2,p3;
  doublecomplex psi0, Q;
  double l = WaveNum * beam_w0 * beam_w0;       /* spreading length */
  double hplus, Qdenom, xsidiff;
  double x, y, z;
  double co,si;
  doublecomplex ctemp;
  double const *incPol;	  	 /*used incident polarization*/

  /* polarisation is in the p1 direction - old version; new - choose from one of the basic
  (latter only for plane wave, for now) */
  if (which=='Y') {p1=1; p2=0; p3=2; incPol=incPolY;}
  if (which=='X') {p1=0; p2=1; p3=2; incPol=incPolX;}
  if (beamtype==B_PLANE)
    for (i=0;i<local_nvoid_Ndip;i++) {
      cExp(WaveNum*DotProd(DipoleCoord+3*i,prop),ctemp); /* ctemp=exp(ik*r.a) */
      cScalMultRVec(incPol,ctemp,b+3*i); /* b[i]=ctemp*incPol */
    }
  else if (beamtype==B_LMINUS) {
    double s,ex;

    s=beam_w0/l;
    PRINTZ("beam confinement factor s:%g\n",s);
    for (i=0;i<local_nvoid_Ndip;i++) {
      x = DipoleCoord[3*i];
      y = DipoleCoord[3*i+1];
      z = DipoleCoord[3*i+2];

      hplus = ((x-beam_x0)*(x-beam_x0)+(y-beam_y0)*(y-beam_y0))/(beam_w0*beam_w0);

      /* calculate Q */
      xsidiff = -(z - beam_z0)/l;
      Qdenom = 1.0 + 4 * xsidiff * xsidiff;
      Q[RE] = 2.0 * xsidiff / Qdenom;
      Q[IM] = -1.0 / Qdenom;

      /* calculate psi0 */
      co = cos(Q[RE]*hplus);
      si = sin(Q[RE]*hplus);
      ex = exp(Q[IM]*hplus);
      psi0[RE] = (-Q[IM]*co + Q[RE]*si)*ex;
      psi0[IM] = (Q[RE]*co + Q[IM]*si)*ex; 
      
      /* calculate the fields */
      co = cos(WaveNum*(z-beam_z0));
      si = sin(WaveNum*(z-beam_z0));
      b[3*i+p1][RE] = psi0[RE]*co - psi0[IM]*si;
      b[3*i+p1][IM] = psi0[IM]*co + psi0[RE]*si;
      b[3*i+p2][RE] = 0.0;
      b[3*i+p2][IM] = 0.0;
      b[3*i+p3][RE] = 0.0;
      b[3*i+p3][IM] = 0.0;
    }
  }
  else if (beamtype==B_DAVIS1 || beamtype==B_BARTON1 || beamtype==B_BARTON3
	   || beamtype==B_DAVIS3 || beamtype==B_BARTON5) {
    double s;
    doublecomplex ex,ey,ez;
    double zeta,xsi,eta,nu;
    doublecomplex d,d2,d3,d4,d5,d6;
    doublecomplex temp,temp2,temp3;

    s=beam_w0/l;
    PRINTZ("beam confinement factor s:%g\n",s);
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
      d[RE]=1.0/(1.0+4*zeta*zeta); d[IM]=2.0*zeta/(1.0+4*zeta*zeta);
      
      if (beamtype==B_DAVIS1 || beamtype==B_BARTON1) {
	ex[RE]=1; ex[IM]=0;
	ey[RE]=0; ey[IM]=0;
	ez[RE]=s; ez[IM]=0;
      }
      else if (beamtype==B_BARTON3) {
	cSquare(d,d2);
	cMult(d2,d,d3);
	ex[RE]=1+s*s*(2*xsi*xsi*d2[RE]+nu*d2[RE]-nu*nu*d3[RE]);
	ex[IM]=s*s*(2*xsi*xsi*d2[IM]+nu*d2[IM]-nu*nu*d3[IM]);
	ey[RE]=s*s*2*xsi*eta*d2[RE];
	ey[IM]=s*s*2*xsi*eta*d2[IM];
	ez[RE]=s+s*s*s*(3*nu*d2[RE]-nu*nu*d3[RE]);
	ez[IM]=s*s*s*(3*nu*d2[IM]-nu*nu*d3[IM]);
      }
      else if (beamtype==B_DAVIS3) {
	cSquare(d,d2);
	cMult(d2,d,d3);
	ex[RE]=1+s*s*(4*xsi*xsi*d2[RE]-nu*nu*d3[RE]);
	ex[IM]=s*s*(4*xsi*xsi*d2[IM]-nu*nu*d3[IM]);
	ey[RE]=s*s*4*xsi*eta*d2[RE];
	ey[IM]=s*s*4*xsi*eta*d2[IM];
	ez[RE]=s+s*s*s*(-2*d[RE]+4*nu*d2[RE]-nu*nu*d3[RE]);
	ez[IM]=s*s*s*(-2*d[IM]+4*nu*d2[IM]-nu*nu*d3[IM]);
      }
      else if (beamtype==B_BARTON5) {
	cSquare(d,d2);
	cMult(d2,d,d3);
	cMult(d3,d,d4);
	cMult(d4,d,d5);
	cMult(d5,d,d6);
	ex[RE]=1+s*s*(2*xsi*xsi*d2[RE]+nu*d2[RE]-nu*nu*d3[RE]) +
	           s*s*s*s*(2*nu*nu*d4[RE]+8*xsi*xsi*nu*d4[RE]-3*nu*nu*nu*d5[RE] -
		   2*xsi*xsi*nu*nu*d5[RE]+.5*nu*nu*nu*nu*d6[RE]);
	ex[IM]=s*s*(2*xsi*xsi*d2[IM]+nu*d2[IM]-nu*nu*d3[IM]) +
	           s*s*s*s*(2*nu*nu*d4[IM]+8*xsi*xsi*nu*d4[IM]-3*nu*nu*nu*d5[IM]-
		   2*xsi*xsi*nu*nu*d5[IM]+.5*nu*nu*nu*nu*d6[IM]);
	ey[RE]=s*s*2*xsi*eta*d2[RE] +
	           s*s*s*s*(2*xsi*eta*d2[RE])*(4*nu*d2[RE]-nu*nu*d3[RE]);
	ey[IM]=s*s*2*xsi*eta*d2[IM] +
	           s*s*s*s*(2*xsi*eta*d2[IM])*(4*nu*d2[IM]-nu*nu*d3[IM]);
	ez[RE]=s+s*s*s*(3*nu*d2[RE]-nu*nu*d3[RE]) +
	           s*s*s*s*s*(10*nu*nu*d4[RE]-5*nu*nu*nu*d5[RE]+.5*nu*nu*nu*nu*d6[RE]);
	ez[IM]=s*s*s*(3*nu*d2[IM]-nu*nu*d3[IM]) +
	           s*s*s*s*s*(10*nu*nu*d4[IM]-5*nu*nu*nu*d5[IM]+.5*nu*nu*nu*nu*d6[IM]);
      }
      temp[RE]=cos(-WaveNum*z);     /* temp=exp(-i*k*z) */
      temp[IM]=sin(-WaveNum*z);     
      temp2[RE]=temp[RE]*d[RE]-temp[IM]*d[IM];      /* temp2=exp(-i*k*z)*D */
      temp2[IM]=temp[RE]*d[IM]+temp[IM]*d[RE];
      temp[RE]=exp(-nu*d[RE])*cos(-nu*d[IM]);      /* temp=exp(-nu*d) */
      temp[IM]=exp(-nu*d[RE])*sin(-nu*d[IM]);
      temp3[RE]=temp[RE]*temp2[RE]-temp[IM]*temp2[IM]; /* temp3=(i*k*z)*D*exp(-nu*d) */
      temp3[IM]=temp[RE]*temp2[IM]+temp[IM]*temp2[RE];
      
      temp[RE]=d[RE]*ez[RE]-d[IM]*ez[IM];              /* temp=d*ez */
      temp[IM]=d[RE]*ez[IM]+d[IM]*ez[RE];

      temp[RE]*=2*xsi;                         /* temp=2*xsi*d*ez */
      temp[IM]*=2*xsi;

      temp2[RE]=-temp[IM]*temp3[RE]-temp[RE]*temp3[IM];  /* temp2=2*i*xsi*d*ez*temp3 */
      temp2[IM]=temp[RE]*temp3[RE]-temp[IM]*temp3[IM];

      /* fill in E field */
      b[3*i+p1][RE] =temp3[RE]*ex[RE]-temp3[IM]*ex[IM];
      b[3*i+p1][IM] =temp3[RE]*ex[IM]+temp3[IM]*ex[RE];
      b[3*i+p2][RE] =temp3[RE]*ey[RE]-temp3[IM]*ey[IM];
      b[3*i+p2][IM] =temp3[RE]*ey[IM]+temp3[IM]*ey[RE];
      b[3*i+p3][RE] =-temp2[RE];
      b[3*i+p3][IM] =-temp2[IM];

      /* correction for carvature of the field */
      /*corr=1+.63*.63*(xsi*xsi/12-1/24)/beam_w0/beam_w0;
      b[3*i+p1][RE]*=corr; b[3*i+p1][IM]*=corr;
      b[3*i+p2][RE]*=corr; b[3*i+p2][IM]*=corr;
      b[3*i+p3][RE]*=corr; b[3*i+p3][IM]*=corr;*/
    }
  }
}

