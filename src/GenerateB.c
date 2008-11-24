#include "cmplx.h"
#include <math.h>
#include "const.h"
#include <stdio.h>
#include "types.h"
#include "comm.h"

void cmult(doublecomplex a,doublecomplex b,doublecomplex *c)
     /* complex multiplycation; c=ab */
{
  c->r=a.r*b.r - a.i*b.i;
  c->i=a.i*b.r + a.r*b.i;
}

/* generate a incident beam
 * original by A. Hoekstra, rewritten by Grimminck.
 * plane wave and buggy type by Alfons Hoekstra, others by
 * Michel Grimminck.
 */
void
GenerateB (char which,		/* x - or y polarized incident light */
           dcomplex *b,		/* the b vector for the incident field */
           REAL **dip_coord,	/* matrix containing dipole coordinates */
           int nldip,		/* number of local dipoles */
           int ndip,		/* total number of dipoles */
           int procid,		/* processor ring number */
           int nproc,		/* total number of processors in the ring */
           double k,		/* wavenumber */
	   int beamtype,
           double w0,           /* beam waist radius */
           double x0,           /* x coordinate of beam center */
           double y0,           /* y coordinate of beam center */
           double z0)           /* z coordinate of beam center */

{ 
  int i;
  
  int p1,p2,p3;
  dcomplex psi0, Q;
  double l = k * w0 * w0;       /* spreading length */
  double hplus, Qdenom, xsidiff;
  double x, y, z;
  double co, si, ex;
  extern int shape;
  FILE *debug;
  double corr;
  
  /* polarisation is in the p1 direction */
  if (which=='Y') {p1=1-3*local_d0; p2=0-3*local_d0; p3=2-3*local_d0;}
  if (which=='X') {p1=0-3*local_d0; p2=1-3*local_d0; p3=2-3*local_d0;}
  
  if (beamtype==PLANE)
    for (i = local_d0; i < local_d1; ++i) {
      b[3*i+p1].r = cos(k * dip_coord[i-local_d0][2]);
      b[3*i+p1].i = sin(k * dip_coord[i-local_d0][2]);
      b[3*i+p2].r = 0.0;
      b[3*i+p2].i = 0.0;
      b[3*i+p3].r = 0.0;
      b[3*i+p3].i = 0.0;
    } 
  else if (beamtype==LMINUS) {
    double s;

    s=w0/l;
    printf("beam confinement factor s:%f\n",s);
    for (i = local_d0; i < local_d1; ++i) {
      x = dip_coord[i-local_d0][0];
      y = dip_coord[i-local_d0][1];
      z = dip_coord[i-local_d0][2];
      
      hplus = ((x-x0)*(x-x0)+(y-y0)*(y-y0))/(w0*w0);
      
      /* calculate Q */
      xsidiff = -(z - z0)/l;
      Qdenom = 1.0 + 4 * xsidiff * xsidiff;
      Q.r = 2.0 * xsidiff / Qdenom;
      Q.i = -1.0 / Qdenom;
      
      /* calculate psi0 */
      co = cos(Q.r*hplus);
      si = sin(Q.r*hplus);
      ex = exp(Q.i*hplus);
      psi0.r = (-Q.i*co + Q.r*si)*ex;
      psi0.i = (Q.r*co + Q.i*si)*ex; 
      
      /* calculate the fields */
      co = cos(k*(z-z0));
      si = sin(k*(z-z0));
      b[3*i+p1].r = psi0.r*co - psi0.i*si;
      b[3*i+p1].i = psi0.i*co + psi0.r*si;
      b[3*i+p2].r = 0.0;
      b[3*i+p2].i = 0.0;
      b[3*i+p3].r = 0.0;
      b[3*i+p3].i = 0.0;
      
      /*if (which=='Y') {
	b[3*i+p3].r = 2.0*(y-y0)*(Q.r*b[3*i+p1].r - Q.i*b[3*i+p1].i)/l;
	b[3*i+p3].i = 2.0*(y-y0)*(Q.i*b[3*i+p1].r + Q.r*b[3*i+p1].i)/l;
      }
      else {
	b[3*i+p3].r = 2.0*(x-x0)*(Q.r*b[3*i+p1].r - Q.i*b[3*i+p1].i)/l;
	b[3*i+p3].i = 2.0*(x-x0)*(Q.i*b[3*i+p1].r + Q.r*b[3*i+p1].i)/l;
      }*/
    }
  }
  else if (beamtype==BUGGY) {
    double s;
    s=w0/l;
    printf("beam confinement factor s:%f\n",s);
    for (i = local_d0; i < local_d1; ++i) {
      x = dip_coord[i-local_d0][0];
      y = dip_coord[i-local_d0][1];
      z = dip_coord[i-local_d0][2];
      
      hplus = ((x-x0)*(x-x0)+(y-y0)*(y-y0))/(w0*w0);
      
      /* calculate Q */
      xsidiff = (z - z0)/l;
      Qdenom = 1.0 + 4 * xsidiff * xsidiff;
      Q.r = 2.0 * xsidiff / Qdenom;
      Q.i = -1.0 / Qdenom;
      
      /* calculate psi0 */
      co = cos(Q.r*hplus);
      si = sin(Q.r*hplus);
      ex = exp(Q.i*hplus);
      psi0.r = (-Q.i*co + Q.r*si)*ex;
      psi0.i = (Q.r*co + Q.i*si)*ex; 
      
      /* calculate the fields */
      co = cos(k*(z-z0));
      si = sin(k*(z-z0));
      b[3*i+p1].r = psi0.r*co - psi0.i*si;
      b[3*i+p1].i = psi0.i*co + psi0.r*si;
      b[3*i+p2].r = 0.0;
      b[3*i+p2].i = 0.0;
      if (which=='Y') {
	b[3*i+p3].r = -2.0*(y-y0)*(Q.r*b[3*i+p1].r - Q.i*b[3*i+p1].i)/l;
	b[3*i+p3].i = -2.0*(y-y0)*(Q.i*b[3*i+p1].r + Q.r*b[3*i+p1].i)/l;
      }
      else {
	b[3*i+p3].r = -2.0*(x-x0)*(Q.r*b[3*i+p1].r - Q.i*b[3*i+p1].i)/l;
	b[3*i+p3].i = -2.0*(x-x0)*(Q.i*b[3*i+p1].r + Q.r*b[3*i+p1].i)/l;
      }
    }
  }
  else if (beamtype==DAVIS1 || beamtype==BARTON1 || beamtype==BARTON3
	   || beamtype==DAVIS3 || beamtype==BARTON5) {
    double s;
    doublecomplex ex,ey,ez,bx,by,bz,he,hb;
    double zeta,xsi,eta,nu;
    doublecomplex d,d2,d3,d4,d5,d6;
    double theta,phi,r;
    doublecomplex temp,temp2,temp3;

    s=w0/l;
    printf("beam confinement factor s:%f\n",s);
    for (i = local_d0; i < local_d1; ++i) {
      if (which=='X') {
        x = dip_coord[i-local_d0][0]-x0;
        y = dip_coord[i-local_d0][1]-y0;
      } else {
        x = dip_coord[i-local_d0][1]-y0;
        y = dip_coord[i-local_d0][0]-x0;
      }

      z = -(dip_coord[i-local_d0][2]-z0); /* minus to give beam the right direction */
      xsi=x/w0;
      eta=y/w0;
      zeta=s*z/w0;
      nu=xsi*xsi+eta*eta;
      d.r=1.0/(1.0+4*zeta*zeta); d.i=2.0*zeta/(1.0+4*zeta*zeta);
      
      if (beamtype==DAVIS1 || beamtype==BARTON1) {
	ex.r=1; ex.i=0;
	ey.r=0; ey.i=0;
	ez.r=s; ez.i=0;
      }
      else if (beamtype==BARTON3) {
	cmult(d,d,&d2);
	cmult(d2,d,&d3);
	ex.r=1+s*s*(2*xsi*xsi*d2.r+nu*d2.r-nu*nu*d3.r);
	ex.i=s*s*(2*xsi*xsi*d2.i+nu*d2.i-nu*nu*d3.i);
	ey.r=s*s*2*xsi*eta*d2.r;
	ey.i=s*s*2*xsi*eta*d2.i;
	ez.r=s+s*s*s*(3*nu*d2.r-nu*nu*d3.r);
	ez.i=s*s*s*(3*nu*d2.i-nu*nu*d3.i);
      }
      else if (beamtype==DAVIS3) {
	cmult(d,d,&d2);
	cmult(d2,d,&d3);
	ex.r=1+s*s*(4*xsi*xsi*d2.r-nu*nu*d3.r);
	ex.i=s*s*(4*xsi*xsi*d2.i-nu*nu*d3.i);
	ey.r=s*s*4*xsi*eta*d2.r;
	ey.i=s*s*4*xsi*eta*d2.i;
	ez.r=s+s*s*s*(-2*d.r+4*nu*d2.r-nu*nu*d3.r);
	ez.i=s*s*s*(-2*d.i+4*nu*d2.i-nu*nu*d3.i);
      }
      else if (beamtype==BARTON5) {
	cmult(d,d,&d2);
	cmult(d2,d,&d3);
	cmult(d3,d,&d4);
	cmult(d4,d,&d5);
	cmult(d5,d,&d6);
	ex.r=1+s*s*(2*xsi*xsi*d2.r+nu*d2.r-nu*nu*d3.r) +
	           s*s*s*s*(2*nu*nu*d4.r+8*xsi*xsi*nu*d4.r-3*nu*nu*nu*d5.r -
		   2*xsi*xsi*nu*nu*d5.r+.5*nu*nu*nu*nu*d6.r);
	ex.i=s*s*(2*xsi*xsi*d2.i+nu*d2.i-nu*nu*d3.i) +
	           s*s*s*s*(2*nu*nu*d4.i+8*xsi*xsi*nu*d4.i-3*nu*nu*nu*d5.i-
		   2*xsi*xsi*nu*nu*d5.i+.5*nu*nu*nu*nu*d6.i);
	ey.r=s*s*2*xsi*eta*d2.r +
	           s*s*s*s*(2*xsi*eta*d2.r)*(4*nu*d2.r-nu*nu*d3.r);
	ey.i=s*s*2*xsi*eta*d2.i +
	           s*s*s*s*(2*xsi*eta*d2.i)*(4*nu*d2.i-nu*nu*d3.i);
	ez.r=s+s*s*s*(3*nu*d2.r-nu*nu*d3.r) +
	           s*s*s*s*s*(10*nu*nu*d4.r-5*nu*nu*nu*d5.r+.5*nu*nu*nu*nu*d6.r);
	ez.i=s*s*s*(3*nu*d2.i-nu*nu*d3.i) +
	           s*s*s*s*s*(10*nu*nu*d4.i-5*nu*nu*nu*d5.i+.5*nu*nu*nu*nu*d6.i);
      }
      temp.r=cos(-k*z); temp.i=sin(-k*z);     /* temp=exp(-i*k*z) */
      temp2.r=temp.r*d.r-temp.i*d.i;      /* temp2=exp(-i*k*z)*D */
      temp2.i=temp.r*d.i+temp.i*d.r;
      temp.r=exp(-nu*d.r)*cos(-nu*d.i);      /* temp=exp(-nu*d) */
      temp.i=exp(-nu*d.r)*sin(-nu*d.i);
      temp3.r=temp.r*temp2.r-temp.i*temp2.i; /* temp3=(i*k*z)*D*exp(-nu*d) */
      temp3.i=temp.r*temp2.i+temp.i*temp2.r;
      
      temp.r=d.r*ez.r-d.i*ez.i;              /* temp=d*ez */
      temp.i=d.r*ez.i+d.i*ez.r;

      temp.r*=2*xsi;                         /* temp=2*xsi*d*ez */
      temp.i*=2*xsi;

      temp2.r=-temp.i*temp3.r-temp.r*temp3.i;/* temp2=2*i*xsi*d*ez*temp3 */
      temp2.i=temp.r*temp3.r-temp.i*temp3.i;

      /* fill in E field */
      b[3*i+p1].r =temp3.r*ex.r-temp3.i*ex.i;
      b[3*i+p1].i =temp3.r*ex.i+temp3.i*ex.r;
      b[3*i+p2].r =temp3.r*ey.r-temp3.i*ey.i;
      b[3*i+p2].i =temp3.r*ey.i+temp3.i*ey.r;
      b[3*i+p3].r =-temp2.r;
      b[3*i+p3].i =-temp2.i;
      
      /* correction for carvature of the field */
      /*corr=1+.63*.63*(xsi*xsi/12-1/24)/w0/w0;
      b[3*i+p1].r*=corr; b[3*i+p1].i*=corr;
      b[3*i+p2].r*=corr; b[3*i+p2].i*=corr;
      b[3*i+p3].r*=corr; b[3*i+p3].i*=corr;*/

      /* show intensity of incident field along x-axis */
      if (shape==LINE && x>0) printf("%f %g\n",x,
			  b[3*i+p1].r*b[3*i+p1].r+b[3*i+p1].i*b[3*i+p1].i+
			  b[3*i+p2].r*b[3*i+p2].r+b[3*i+p2].i*b[3*i+p2].i+
			  b[3*i+p3].r*b[3*i+p3].r+b[3*i+p3].i*b[3*i+p3].i);
    }
  }
  else {printf("requested beam type is not implemented\n"); stop(1);}
}







