/* FILE : chebyshev.c
 * $Date::                            $
 * Descr: routines for determining parameters of the chebyshev particles
 *
 * Copyright (C) 2011-2013 ADDA contributors
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
#include "const.h" // keep this first
// project headers
#include "io.h"
// system headers 
#include <math.h>

// LOCAL VARIABLES

#define MAX_NEWTON_ITERATIONS 20
// The following corresponds to sqrt(DBL_EPSILON)/10. It should lead to DBL_EPSILON precision in the function value
#define NEWTON_EPS 1e-9
enum newt_func { // type of function to calculate in Newton's maximization
	NF_VALUE, // f(x)
	NF_DRATIO // f'(x)/f''(x)
};
static double eps;
static int n;

//======================================================================================================================

static double Zfunc(double x,enum newt_func mode)
// computes either f(x) or f'(x)/f''(x), where f(x)=[1+e*cos(nx)]cos(x)
{
	double e,xn,s,c,sn,cn;

	e=-fabs(eps); // assumes that e is negative; solution is trivial otherwise
	xn=n*x;
	switch (mode) {
		case NF_VALUE:
			c=cos(x);
			cn=cos(xn);
			return (1+e*cn)*c;
		case NF_DRATIO:
			s=sin(x);
			c=cos(x);
			sn=sin(xn);
			cn=cos(xn);
			return (-e*n*sn*c-(1+e*cn)*s)/(2*e*n*sn*s-(1+e*(1+n*n)*cn)*c);
	}
	LogError(ONE_POS,"Unknown mode %d for calling Zfunc",(int)mode);
}

//======================================================================================================================

static double Xfunc(double x,enum newt_func mode)
// computes either f(x) or f'(x)/f''(x), where f(x)=[1+e*sin(nx)]cos(x)
{
	double e,xn,s,c,sn,cn;

	e=fabs(eps); // assumes that e is positive (maximum value is invariant to sign change)
	xn=n*x;
	switch (mode) {
		case NF_VALUE:
			c=cos(x);
			sn=sin(xn);
			return (1+e*sn)*c;
		case NF_DRATIO:
			s=sin(x);
			c=cos(x);
			sn=sin(xn);
			cn=cos(xn);
			return (e*n*cn*c-(1+e*sn)*s)/(-2*e*n*cn*s-(1+e*(1+n*n)*sn)*c);
	}
	LogError(ONE_POS,"Unknown mode %d for calling Zfunc",(int)mode);
}

//======================================================================================================================

static double NewtonMaximum(double x0,double x1,double x2,double (*func)(double x,enum newt_func))
/* finds a maximum (minimum) of a function f(x) using the Newton's method; x0 - should be a good guess, and [x1,x2] is
 * a bounding interval (preliminary analysis is usually required for this); func - is a reference to function, which
 * implements both f(x) and f'(x)/f''(x) (depending on the second argument)
 */
{
	int i;
	double x,dx;

	x=x0;
	for (i=0;i<MAX_NEWTON_ITERATIONS;i++) {
		dx=(*func)(x,NF_DRATIO);
		x-=dx;
		if (x<x1 || x>x2) LogError(ONE_POS,"Newton's method fell out of bounds");
		if (fabs(dx)<=NEWTON_EPS*x) return (*func)(x,NF_VALUE);
	}
	LogError(ONE_POS,"Newton's method failed to converge in %d iterations",i);
}

//======================================================================================================================

static double Zmax(double e)
// finds maximum of function [1+e*cos(n*x)]cos(x)
{
	double res,n2,tmp,arg,x01,x02,x0;

	n2=n*n;
	tmp=e*(n2+1);
	if (tmp>=-1) res=1+e;
	/* two special cases for small n are based on direct solution of f'(x)=0; the formulae are self-derived but agree
	 * (for n=2) with A. Mugnai and W.J. Wiscombe, “Scattering of radiation by moderately nonspherical particles,”
	 * J. Atmos. Sci. 37, 1291-1307 (1980).
	 */
	else if (n==1) res=-1/(4*e);
	else if (n==2) res=2*(1-e)*sqrt((e-1)/(6*e))/3;
	else { // a few other values of n allow for a closed-form solution, but they are cumbersome
		arg=PI/n;
		x01=arg - tan(arg)*(1-e)/(1-tmp);
		x02=sqrt(6*(1+tmp)/(1+e*(1+6*n2+n2*n2)));
		// this is not the most efficient empirics, but should work fine
		x0=Zfunc(x01,NF_VALUE)>Zfunc(x02,NF_VALUE) ? x01 : x02;
		res=NewtonMaximum(x0,0,arg,Zfunc);
	}
	return res;
}

//======================================================================================================================

static double Xmax(void)
// finds maximum of function [1+e*sin(n*x)]cos(x); in principle should work for any n, but is called only for odd n
{
	double res,e,s,arg;
	e=fabs(eps); // result does not depend on sign of eps anyway

	if (n==1) {
		s=2*e/(1+sqrt(1+8*e*e)); // solution of quadratic equation on sin(x) from f'(x)=0
		res=(1+e*s)*sqrt(1-s*s);
	} // a few other values of n allow for a closed-form solution, but they are cumbersome
	else {
		arg=PI/(2*n);
		res=NewtonMaximum(arg-tan(arg)*(1+e)/(1+e*(1+n*n)),0,arg,Xfunc);
	}
	return res;
}

//======================================================================================================================

void ChebyshevParams(double eps_in,int n_in,double *dx,double *dz,double *sz,double *vr)
// finds geometric parameters of the Chebyshev particles: r(th)=1+eps*cos(n*th)
{
	double n2,tmp1,tmp2,zmax,zmin,xmax;
	eps=eps_in;
	n=n_in;
	n2=n*n;

	// In the following, Zmax is called a few times, but never with the same arguments
	zmax=Zmax(eps);
	if (n&1) { // n=2k+1
		zmin=-Zmax(-eps);
		xmax=Xmax();
	}
	else { // n=2k
		zmin=-zmax; // symmetric case
		if (n&2) xmax=Zmax(-eps); // n=4k+2
		else xmax=zmax;
	}
	*dz=zmax-zmin;
	*sz=(zmax+zmin)/2;
	*dx=2*xmax;
	/* determine volume fraction; the formula is self-derived but agrees with A. Mugnai and W.J. Wiscombe, “Scattering
	 * of radiation by moderately nonspherical particles,” J. Atmos. Sci. 37, 1291-1307 (1980).
	 */
	tmp1=eps*eps/4;
	tmp2=1+6*tmp1*(4*n2-2)/(4*n2-1);
	if (!(n&1)) tmp2-=eps*(3*(1+tmp1)/(n2-1)+tmp1/(9*n2-1)); // n=2k
	*vr = FOUR_PI_OVER_THREE*tmp2/((*dx)*(*dx)*(*dx));
}
