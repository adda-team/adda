/* File: cmplx.h
 * $Date::                            $
 * Descr: inline complex functions plus few auxiliary functions
 *
 *        'const' can be used for many more function variables, however it doesn't work in
 *        combination with 'doublecomplex *' or more nested lists. That seems to be a principal
 *        limitation of C standard (some compilers may work, some produce warnings). A few changes
 *        for reliability and stability were made according to the ideas of section 5.5 of the
 *        Numerical Recipes, 3rd edition.
 *
 * Copyright (C) 2006-2008,2010,2012 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __cmplx_h
#define __cmplx_h

// project headers
#include "const.h"    // for math constants
#include "function.h" // for static inline
#include "types.h"    // for doublecomplex
// system headers
#include <math.h>     // for cos, sin
#include <string.h>   // for memcpy

#ifdef USE_SSE3
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#endif

//============================================================
// operations on complex numbers

static inline void cEqual(const doublecomplex a,doublecomplex b)
// performs b=a; !!! pointers b and a must be different !!!
{
	memcpy(b,a,sizeof(doublecomplex));
}

//============================================================

static inline double cAbs2(const doublecomplex a)
// square of absolute value of complex number; |a|^2
{
	return (a[RE]*a[RE] + a[IM]*a[IM]);
}

//============================================================

static inline double cAbs(const doublecomplex a)
// absolute value of complex number |a|, specially designed to avoid overflow
{
	double u,v,w;
	u=fabs(a[RE]);
	v=fabs(a[IM]);

	if (u==0 && v==0) return 0;
	else {
		if (u>=v) {
			w=v/u;
			return (u*sqrt(1+w*w));
		}
		else {
			w=u/v;
			return (v*sqrt(1+w*w));
		}
	}
}

//============================================================

static inline void cConj(const doublecomplex a,doublecomplex b)
// complex conjugate; b=a*; b and a may coincide
{
	b[RE] = a[RE];
	b[IM] = - a[IM];
}

//============================================================

static inline void cAdd(const doublecomplex a,const doublecomplex b,doublecomplex c)
// add two complex numbers; c=a+b; any two or three pointers may coincide
{
	c[RE] = a[RE] + b[RE];
	c[IM] = a[IM] + b[IM];
}

//============================================================

static inline void cSubtr(const doublecomplex a,const doublecomplex b,doublecomplex c)
// subtract two complex numbers; c=a-b; any two or three pointers may coincide
{
	c[RE] = a[RE] - b[RE];
	c[IM] = a[IM] - b[IM];
}

//============================================================

static inline void cSquare(const doublecomplex a,doublecomplex b)
// square of complex number; b=a^2; !!! pointers b and a must me different !!!
{
	b[RE]=a[RE]*a[RE] - a[IM]*a[IM];
	b[IM]=2*a[IM]*a[RE];
}

//============================================================

static inline void cMultReal(const double a,const doublecomplex b,doublecomplex c)
// complex multiplication by real; c=ab; pointers b and c may coincide
{
	c[RE]=a*b[RE];
	c[IM]=a*b[IM];
}
//============================================================

static inline void cMultRealConj(const double a,const doublecomplex b,doublecomplex c)
// conjugation with complex multiplication by real; c=ab(*); pointers b and c may coincide
{
	c[RE]=a*b[RE];
	c[IM]=-a*b[IM];
}

//============================================================

static inline void cMult_i(doublecomplex c)
// complex multiplication by i; c=i*c
{
	double tmp;
	tmp=c[RE];
	c[RE]=-c[IM];
	c[IM]=tmp;
}

//============================================================

static inline void cMult_i2(doublecomplex a,doublecomplex b)
// complex multiplication by i; b=i*a; !!! pointers b and a must be different !!!
{
	b[RE]=-a[IM];
	b[IM]=a[RE];
}
//============================================================

static inline void cMult(const doublecomplex a,const doublecomplex b,doublecomplex c)
// complex multiplication; c=ab; !!! pointer c must be different from a and b !!!
{
	c[RE]=a[RE]*b[RE] - a[IM]*b[IM];
	c[IM]=a[IM]*b[RE] + a[RE]*b[IM];
}

//============================================================

static inline void cMultSelf(doublecomplex a,const doublecomplex b)
// complex multiplication; a*=b; !!! pointers a and b must be different !!!
{
	double tmp;

	tmp=a[RE];
	a[RE]=a[RE]*b[RE] - a[IM]*b[IM];
	a[IM]=a[IM]*b[RE] + tmp*b[IM];
}

//============================================================

static inline void cMultConj(const doublecomplex a,const doublecomplex b,doublecomplex c)
// complex multiplication by conjugate; c=a*b(*); !!! pointer c must be different from a and b !!!
{
	c[RE]=a[RE]*b[RE] + a[IM]*b[IM];
	c[IM]=a[IM]*b[RE] - a[RE]*b[IM];
}

//============================================================

static inline double cMultConRe(const doublecomplex a,const doublecomplex b)
// complex multiplication; returns real(a*b_conjugated); pointers a and b may coincide
{
	return (a[RE]*b[RE] + a[IM]*b[IM]);
}

//============================================================

static inline double cMultConIm(const doublecomplex a,const doublecomplex b)
// complex multiplication; returns imag(a*b_conjugated); pointers a and b may coincide
{
	return (a[IM]*b[RE] - a[RE]*b[IM]);
}

//============================================================

static inline void cLinComb(const doublecomplex a,const doublecomplex b,
                     const double c1,const double c2,doublecomplex c)
// linear combination of two complex numbers; c=c1*a+c2*b; any two or three pointers may coincide
{
	c[RE]=c1*a[RE]+c2*b[RE];
	c[IM]=c1*a[IM]+c2*b[IM];
}

//============================================================

static inline void cInvSign(doublecomplex a)
// change sign of complex number; a*=-1;
{
	a[RE] = - a[RE];
	a[IM] = - a[IM];
}

//============================================================

static inline void cInvSign2(const doublecomplex a,doublecomplex b)
// change sign of complex number and store to different address; b=-a; pointers a and b may coincide
{
	b[RE] = - a[RE];
	b[IM] = - a[IM];
}

//============================================================

static inline void cInv(const doublecomplex a,doublecomplex b)
// complex inversion; b=1/a; designed to avoid under and overflows; pointers a and b may coincide
{
	double tmp;

	if (fabs(a[RE])>=fabs(a[IM])) {
		tmp=a[IM]/a[RE];
		b[RE]=1/(a[RE]+a[IM]*tmp);
		b[IM]=-b[RE]*tmp;
	}
	else {
		tmp=a[RE]/a[IM];
		b[IM]=-1/(a[RE]*tmp+a[IM]);
		b[RE]=-b[IM]*tmp;
	}
}

//============================================================

static inline double cInvIm(const doublecomplex a)
// returns Im of inverse of a; designed to avoid under and overflows
{
	double tmp;

	if (fabs(a[RE])>=fabs(a[IM])) {
		tmp=a[IM]/a[RE];
		return (-tmp/(a[RE]+a[IM]*tmp));
	}
	else {
		tmp=a[RE]/a[IM];
		return (-1/(a[RE]*tmp+a[IM]));
	}
}

//============================================================

static inline void cDiv(const doublecomplex a,const doublecomplex b,doublecomplex c)
/* complex division; c=a/b; designed to avoid under and overflows
 * !!! pointer c must be different from a !!!
 */
{
	double u,v;

	if (fabs(b[RE])>=fabs(b[IM])) {
		u=b[IM]/b[RE];
		v=1/(b[RE]+b[IM]*u);
		c[RE]=(a[RE]+a[IM]*u)*v;
		c[IM]=(a[IM]-a[RE]*u)*v;
	}
	else {
		u=b[RE]/b[IM];
		v=1/(b[RE]*u+b[IM]);
		c[RE]=(a[RE]*u+a[IM])*v;
		c[IM]=(a[IM]*u-a[RE])*v;
	}
}

//============================================================

static inline void cDivSelf(doublecomplex a,const doublecomplex b)
// complex division; a/=b; designed to avoid under and overflows; pointers a and b may coincide
{
	double u,v,w;

	w=a[RE];
	if (fabs(b[RE])>=fabs(b[IM])) {
		u=b[IM]/b[RE];
		v=1/(b[RE]+b[IM]*u);
		a[RE]=(w+a[IM]*u)*v;
		a[IM]=(a[IM]-w*u)*v;
	}
	else {
		u=b[RE]/b[IM];
		v=1/(b[RE]*u+b[IM]);
		a[RE]=(w*u+a[IM])*v;
		a[IM]=(a[IM]*u-w)*v;
	}
}

//============================================================

static inline void cSqrt(const doublecomplex a,doublecomplex b)
/* complex square root; b=sqrt(a); designed to avoid under and overflows;
 * branch cut discontinuity is (-inf,0) - b[RE]>=0; pointers a and b may coincide
 */
{
	double u,v,w,r;

	u=fabs(a[RE]);
	v=fabs(a[IM]);
	if (u==0 && v==0) b[RE]=b[IM]=0;
	else {
		// first determine w
		if (u>=v) {
			r=v/u;
			w=sqrt(u)*sqrt((1+sqrt(1+r*r))/2);
		}
		else {
			r=u/v;
			w=sqrt(v)*sqrt((r+sqrt(1+r*r))/2);
		}
		// compute the result
		if (a[RE]>=0) {
			b[RE]=w;
			b[IM]=a[IM]/(2*w);
		}
		else {
			b[RE]=v/(2*w);
			if (a[IM]>=0) b[IM]=w;
			else b[IM]=-w;
		}
	}
}

//============================================================

static inline void imExp(const double arg,doublecomplex c)
// exponent of imaginary argument c=Exp(i*arg); optimization is performed by compiler
{
	c[RE]=cos(arg);
	c[IM]=sin(arg);
}

//============================================================

static inline void imExp_arr(const double arg,const int size,doublecomplex *c)
/* construct an array of exponent of imaginary argument c=Exp(i*k*arg)
 * where k=0,1,...,size-1. Uses stable recurrence from Numerical Recipes.
 * Optimization of the initial simultaneous calculation of sin and cos is performed
 * by compiler; It is assumed that size is at least 1
 */
{
	int k;
	double a,b;

	c[0][RE]=1;
	c[0][IM]=0;
	if (size>1) {
		// set a=2*sin^2(arg/2), b=sin(arg)
		a=sin(arg/2);
		b=cos(arg/2);
		b*=2*a;
		a*=2*a;
		// this a bit faster than in the main cycle
		c[1][RE]=1-a;
		c[1][IM]=b;
		// main cycle
		for (k=2;k<size;k++) {
			/* potentially compiler may open brackets to accelerate calculation but lose significant
			 * digits. We hope it doesn't happen, but it should not be a big problem anyway
			 */
			c[k][RE]=c[k-1][RE]-(a*c[k-1][RE]+b*c[k-1][IM]);
			c[k][IM]=c[k-1][IM]-(a*c[k-1][IM]-b*c[k-1][RE]);
		}
	}
}

//============================================================

static inline void cExp(const doublecomplex arg,doublecomplex c)
/* complex exponent of complex argument c=Exp(arg); optimization is performed by compiler;
 * !!! pointer c must be different from arg !!!
 */
{
	c[RE]=c[IM]=exp(arg[RE]);
	c[RE]*=cos(arg[IM]);
	c[IM]*=sin(arg[IM]);
}

//============================================================

static inline void cExpSelf(doublecomplex arg)
/* complex exponent of complex argument arg=Exp(arg); result is stored in the argument itself
 * Optimization is performed by compiler
 */
{
	double tmp;

	tmp=arg[IM];
	arg[RE]=arg[IM]=exp(arg[RE]);
	arg[RE]*=cos(tmp);
	arg[IM]*=sin(tmp);
}

//============================================================
// operations on complex vectors

static inline void cvMultScal(const double a,doublecomplex b[static restrict 3],
	doublecomplex c[static restrict 3])
// multiplication of vector[3] by real scalar; c=ab; !!! pointers b and c must not alias !!!
{
	c[0][RE] = a*b[0][RE];
	c[0][IM] = a*b[0][IM];
	c[1][RE] = a*b[1][RE];
	c[1][IM] = a*b[1][IM];
	c[2][RE] = a*b[2][RE];
	c[2][IM] = a*b[2][IM];
}

//============================================================
/* Here and further we leave the possibility for the compiler to consider double and doublecomplex
 * types compatible, not to require strict aliasing rule upon them. Therefore we declare double
 * pointers as restrict wherever possible. This doesn't add any problems to the rest of the code,
 * since no conversion between double and doublecomplex should be performed.
 */

static inline void cScalMultRVec(const double a[static restrict 3],const doublecomplex b,
	doublecomplex c[static restrict 3])
// complex scalar- real vector[3] multiplication; c=b*a; !!! pointers a,b,c must not alias !!!
{
	c[0][RE] = b[RE]*a[0];
	c[0][IM] = b[IM]*a[0];
	c[1][RE] = b[RE]*a[1];
	c[1][IM] = b[IM]*a[1];
	c[2][RE] = b[RE]*a[2];
	c[2][IM] = b[IM]*a[2];
}

//============================================================

static inline void cvMultScal_cmplx(const doublecomplex a,doublecomplex b[static restrict 3],
	doublecomplex c[static restrict 3])
// multiplication of vector[3] by complex scalar; c=ab; !!! pointers a,b,c must not alias !!!
{
	c[0][RE] = a[RE]*b[0][RE] - a[IM]*b[0][IM];
	c[0][IM] = a[RE]*b[0][IM] + a[IM]*b[0][RE];
	c[1][RE] = a[RE]*b[1][RE] - a[IM]*b[1][IM];
	c[1][IM] = a[RE]*b[1][IM] + a[IM]*b[1][RE];
	c[2][RE] = a[RE]*b[2][RE] - a[IM]*b[2][IM];
	c[2][IM] = a[RE]*b[2][IM] + a[IM]*b[2][RE];
}

//============================================================

static inline double cvNorm2(doublecomplex a[static 3])
// square of the norm of a complex vector[3]
{
	return ( a[0][RE]*a[0][RE] + a[0][IM]*a[0][IM]
	       + a[1][RE]*a[1][RE] + a[1][IM]*a[1][IM]
	       + a[2][RE]*a[2][RE] + a[2][IM]*a[2][IM] );
}


//============================================================

static inline void cDotProd(doublecomplex a[static restrict 3],doublecomplex b[static restrict 3],
	doublecomplex c)
/* conjugate dot product of two complex vector[3]; c=a.b = a[0]*b*[0]+...+a[2]*b*[2]
 * !!! pointers a,b,c must not alias !!!)
 */
{
	c[RE] = a[0][RE]*b[0][RE] + a[0][IM]*b[0][IM]
	      + a[1][RE]*b[1][RE] + a[1][IM]*b[1][IM]
	      + a[2][RE]*b[2][RE] + a[2][IM]*b[2][IM];
	c[IM] = a[0][IM]*b[0][RE] - a[0][RE]*b[0][IM]
	      + a[1][IM]*b[1][RE] - a[1][RE]*b[1][IM]
	      + a[2][IM]*b[2][RE] - a[2][RE]*b[2][IM];
}

//============================================================

static inline double cDotProd_Re(doublecomplex a[static 3],doublecomplex b[static 3])
// real part of dot product of two complex vector[3]; c=Re(a.b); a and b may alias
{
	return ( a[0][RE]*b[0][RE] + a[0][IM]*b[0][IM]
	       + a[1][RE]*b[1][RE] + a[1][IM]*b[1][IM]
	       + a[2][RE]*b[2][RE] + a[2][IM]*b[2][IM] );
}

//============================================================

static inline double cDotProd_Im(doublecomplex a[static 3],doublecomplex b[static 3])
// imaginary part of dot product of two complex vector[3]; c=Im(a.b); a and b may alias
{
	return ( a[0][IM]*b[0][RE] - a[0][RE]*b[0][IM]
	       + a[1][IM]*b[1][RE] - a[1][RE]*b[1][IM]
	       + a[2][IM]*b[2][RE] - a[2][RE]*b[2][IM] );
}

//============================================================

static inline void cDotProd_conj(doublecomplex a[static restrict 3],doublecomplex b[static restrict 3],
	doublecomplex c)
/* dot product of two complex vector[3]; c=a.b* = a[0]*b[0]+...+a[2]*b[2]
 * !!! pointers a,b,c must not alias !!!
 */
{
	c[RE] = a[0][RE]*b[0][RE] - a[0][IM]*b[0][IM]
	      + a[1][RE]*b[1][RE] - a[1][IM]*b[1][IM]
	      + a[2][RE]*b[2][RE] - a[2][IM]*b[2][IM];
	c[IM] = a[0][IM]*b[0][RE] + a[0][RE]*b[0][IM]
	      + a[1][IM]*b[1][RE] + a[1][RE]*b[1][IM]
	      + a[2][IM]*b[2][RE] + a[2][RE]*b[2][IM];
}

//============================================================

static inline double cDotProd_conj_Re(doublecomplex a[static 3],doublecomplex b[static 3])
// real part of dot product of two complex vector[3]; c=Re(a.b*); a and b may alias
{
	return ( a[0][RE]*b[0][RE] - a[0][IM]*b[0][IM]
	       + a[1][RE]*b[1][RE] - a[1][IM]*b[1][IM]
	       + a[2][RE]*b[2][RE] - a[2][IM]*b[2][IM] );
}

//============================================================

static inline double cDotProd_conj_Im(doublecomplex a[static 3],doublecomplex b[static 3])
// imaginary part of dot product of two complex vector[3]; c=Im(a.b*); a and b may alias
{
	return ( a[0][IM]*b[0][RE] + a[0][RE]*b[0][IM]
	       + a[1][IM]*b[1][RE] + a[1][RE]*b[1][IM]
	       + a[2][IM]*b[2][RE] + a[2][RE]*b[2][IM] );
}

//============================================================

static inline void cvAdd(doublecomplex a[static restrict 3],doublecomplex b[static restrict 3],
	doublecomplex c[static restrict 3])
/* add two complex vector[3]; c=a+b; !!! pointers a,b,c must not alias !!!
 * (coincidence of, e.g., a and c is logically possible but disallowed for 'restrict' optimization)
 */
{
	c[0][RE] = a[0][RE] + b[0][RE];
	c[0][IM] = a[0][IM] + b[0][IM];
	c[1][RE] = a[1][RE] + b[1][RE];
	c[1][IM] = a[1][IM] + b[1][IM];
	c[2][RE] = a[2][RE] + b[2][RE];
	c[2][IM] = a[2][IM] + b[2][IM];
}

//============================================================

static inline void cvAdd2Self(doublecomplex a[static restrict 3],doublecomplex b[static restrict 3],
	doublecomplex c[static restrict 3])
/* increment one complex vector[3] by sum of other two; a+=b+c;!!! pointers a,b,c must not alias !!!
 * (coincidence of, e.g., a and b is logically possible but disallowed for 'restrict' optimization)
 */
{
	a[0][RE] += b[0][RE] + c[0][RE];
	a[0][IM] += b[0][IM] + c[0][IM];
	a[1][RE] += b[1][RE] + c[1][RE];
	a[1][IM] += b[1][IM] + c[1][IM];
	a[2][RE] += b[2][RE] + c[2][RE];
	a[2][IM] += b[2][IM] + c[2][IM];
}

//============================================================

static inline void cvSubtrSelf(doublecomplex a[static restrict 3],doublecomplex b[static restrict 3])
/* inverse sign of the complex vector[3] b and then increment by complex vector a; b=a-b
 * rather uncommon operation that appears only in one place in the code and is subject to compiler
 * optimization. !!! b and a must not alias !!! (and their coincidence has little sense)
 */
{
	b[0][RE] = a[0][RE] - b[0][RE];
	b[0][IM] = a[0][IM] - b[0][IM];
	b[1][RE] = a[1][RE] - b[1][RE];
	b[1][IM] = a[1][IM] - b[1][IM];
	b[2][RE] = a[2][RE] - b[2][RE];
	b[2][IM] = a[2][IM] - b[2][IM];
}

//============================================================

static inline void crDotProd(doublecomplex a[static 3],const double b[static restrict 3],
	doublecomplex c)
// dot product of complex and real vectors[3]; c=a.b; pointers a and c may alias, but that is weird
{
	c[RE] = a[0][RE]*b[0] + a[1][RE]*b[1] + a[2][RE]*b[2];
	c[IM] = a[0][IM]*b[0] + a[1][IM]*b[1] + a[2][IM]*b[2];
}

//============================================================

static inline double crDotProd_Re(doublecomplex a[static 3],const double b[static 3])
// real part of dot product of complex and real vectors[3]; c=Re(a.b); a and b may alias
{
	return (a[0][RE]*b[0] + a[1][RE]*b[1] + a[2][RE]*b[2]);
}

//============================================================

static inline double crDotProd_Im(doublecomplex a[static restrict 3],const double b[static restrict 3])
// imaginary part of dot product of complex and real vectors[3]; c=Im(a.b); a and b may alias
{
	return (a[0][IM]*b[0] + a[1][IM]*b[1] + a[2][IM]*b[2]);
}

//============================================================

static inline void cvIncremScaled_cmplx(doublecomplex a[static restrict 3],const doublecomplex b,
	doublecomplex c[static restrict 3])
/* increment of complex vectors[3] by complex-scaled other vector; c+=b*a;
 * !!! pointers a,b,c must not alias !!!
 */
{
	c[0][RE] += b[RE]*a[0][RE] - b[IM]*a[0][IM];
	c[0][IM] += b[RE]*a[0][IM] + b[IM]*a[0][RE];
	c[1][RE] += b[RE]*a[1][RE] - b[IM]*a[1][IM];
	c[1][IM] += b[RE]*a[1][IM] + b[IM]*a[1][RE];
	c[2][RE] += b[RE]*a[2][RE] - b[IM]*a[2][IM];
	c[2][IM] += b[RE]*a[2][IM] + b[IM]*a[2][RE];
}
//============================================================

static inline void cvMultAdd(doublecomplex a[static restrict 3],const doublecomplex b,
	doublecomplex c[static restrict 3])
/* multiply complex vectors[3] with complex constant and add another vector; c=b*c+a
 * !!! pointers a,b,c must not alias !!!
 */
{
	double tmp;
	tmp=c[0][RE];
	c[0][RE] = c[0][RE]*b[RE] - c[0][IM]*b[IM] + a[0][RE];
	c[0][IM] = tmp*b[IM] + c[0][IM]*b[RE] + a[0][IM];
	tmp=c[1][RE];
	c[1][RE] = c[1][RE]*b[RE] - c[1][IM]*b[IM] + a[1][RE];
	c[1][IM] = tmp*b[IM] + c[1][IM]*b[RE] + a[1][IM];
	tmp=c[2][RE];
	c[2][RE] = c[2][RE]*b[RE] - c[2][IM]*b[IM] + a[2][RE];
	c[2][IM] = tmp*b[IM] + c[2][IM]*b[RE] + a[2][IM];
}

//============================================================

static inline void cvLinComb1(doublecomplex a[static restrict 3],doublecomplex b[static restrict 3],
	const double c1,doublecomplex c[static restrict 3])
/* linear combination of complex vectors[3]; second coefficient is unity; c=c1*a+b
 * !!! pointers a,b,c must not alias !!!
 * (coincidence of, e.g., a and c is logically possible but disallowed for 'restrict' optimization)
 */
{
	c[0][RE] = c1*a[0][RE] + b[0][RE];
	c[0][IM] = c1*a[0][IM] + b[0][IM];
	c[1][RE] = c1*a[1][RE] + b[1][RE];
	c[1][IM] = c1*a[1][IM] + b[1][IM];
	c[2][RE] = c1*a[2][RE] + b[2][RE];
	c[2][IM] = c1*a[2][IM] + b[2][IM];
}

//============================================================

static inline void cvLinComb1_cmplx(doublecomplex a[static restrict 3],doublecomplex b[static restrict 3],
	const doublecomplex c1,doublecomplex c[static restrict 3])
/* linear combination of complex vectors[3] with complex coefficients;
 * second coefficient is unity; c=c1*a+b; !!! pointers a,b,c,c1 must not alias !!!
 */
{
	c[0][RE] = a[0][RE]*c1[RE] - a[0][IM]*c1[IM] + b[0][RE];
	c[0][IM] = a[0][RE]*c1[IM] + a[0][IM]*c1[RE] + b[0][IM];
	c[1][RE] = a[1][RE]*c1[RE] - a[1][IM]*c1[IM] + b[1][RE];
	c[1][IM] = a[1][RE]*c1[IM] + a[1][IM]*c1[RE] + b[1][IM];
	c[2][RE] = a[2][RE]*c1[RE] - a[2][IM]*c1[IM] + b[2][RE];
	c[2][IM] = a[2][RE]*c1[IM] + a[2][IM]*c1[RE] + b[2][IM];
}

//============================================================

static inline void cSymMatrVec(doublecomplex matr[static restrict 6],doublecomplex vec[static restrict 3],
	doublecomplex res[static restrict 3])
/* multiplication of complex symmetric matrix[6] by complex vec[3]; res=matr.vec
 * !!! pointers matr,vec,res must not alias !!!
 */
{
	res[0][RE] = matr[0][RE]*vec[0][RE] - matr[0][IM]*vec[0][IM]
	           + matr[1][RE]*vec[1][RE] - matr[1][IM]*vec[1][IM]
	           + matr[2][RE]*vec[2][RE] - matr[2][IM]*vec[2][IM];
	res[0][IM] = matr[0][RE]*vec[0][IM] + matr[0][IM]*vec[0][RE]
	           + matr[1][RE]*vec[1][IM] + matr[1][IM]*vec[1][RE]
	           + matr[2][RE]*vec[2][IM] + matr[2][IM]*vec[2][RE];

	res[1][RE] = matr[1][RE]*vec[0][RE] - matr[1][IM]*vec[0][IM]
	           + matr[3][RE]*vec[1][RE] - matr[3][IM]*vec[1][IM]
	           + matr[4][RE]*vec[2][RE] - matr[4][IM]*vec[2][IM];
	res[1][IM] = matr[1][RE]*vec[0][IM] + matr[1][IM]*vec[0][RE]
	           + matr[3][RE]*vec[1][IM] + matr[3][IM]*vec[1][RE]
	           + matr[4][RE]*vec[2][IM] + matr[4][IM]*vec[2][RE];

	res[2][RE] = matr[2][RE]*vec[0][RE] - matr[2][IM]*vec[0][IM]
	           + matr[4][RE]*vec[1][RE] - matr[4][IM]*vec[1][IM]
	           + matr[5][RE]*vec[2][RE] - matr[5][IM]*vec[2][IM];
	res[2][IM] = matr[2][RE]*vec[0][IM] + matr[2][IM]*vec[0][RE]
	           + matr[4][RE]*vec[1][IM] + matr[4][IM]*vec[1][RE]
	           + matr[5][RE]*vec[2][IM] + matr[5][IM]*vec[2][RE];
}

//============================================================
// operations on real vectors

static inline void vAdd(const double a[static 3],const double b[static 3],double c[static 3])
// adds two real vectors; c=a+b; vectors may alias
{
	c[0]=a[0]+b[0];
	c[1]=a[1]+b[1];
	c[2]=a[2]+b[2];
}

//============================================================
static inline void vInvSign(double a[static restrict 3])
// inverts the sign in the double vector[3]
{
	a[0]=-a[0];
	a[1]=-a[1];
	a[2]=-a[2];
}

//============================================================

static inline void vMultScal(const double a,const double b[static restrict 3],double c[static restrict 3])
/* multiplication of real vector by scalar; c=a*b;
 * !!! pointers b and c must not alias !!! (for 'restrict' optimization)
 */
{
	c[0]=a*b[0];
	c[1]=a*b[1];
	c[2]=a*b[2];
}

//============================================================
static inline void vMultScalSelf(const double a,double b[static restrict 3])
// multiplication of real vector by scalar; b*=a;
{
	b[0]*=a;
	b[1]*=a;
	b[2]*=a;
}

//============================================================

static inline void vMult(const double a[static restrict 3],const double b[static restrict 3],
	double c[static restrict 3])
/* multiplication of two vectors (by elements); c[i]=a[i]*b[i];
 * !!! pointers a,b,c must not alias !!!
 * (coincidence of, e.g., a and c is logically possible but disallowed for 'restrict' optimization)
 */
{
	c[0]=a[0]*b[0];
	c[1]=a[1]*b[1];
	c[2]=a[2]*b[2];
}

//============================================================

static inline void vSquare(const double a[static restrict 3],double b[static restrict 3])
/* element-wise square of a vector; b[i]=a[i]*a[i]; !!! pointers a and b must not alias !!!
 * (coincidence of a and b is logically possible but disallowed for 'restrict' optimization)
 */
{
	b[0]=a[0]*a[0];
	b[1]=a[1]*a[1];
	b[2]=a[2]*a[2];
}

//============================================================

static inline double DotProd(const double a[static 3],const double b[static 3])
// dot product of two real vectors[3]; a and b may alias
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

//============================================================

static inline void LinComb(const double a[static restrict 3],const double b[static restrict 3],
	const double c1,const double c2, double c[static restrict 3])
/* linear combination of real vectors[3]; c=c1*a+c2*b
 * !!! pointers a,b,c must not alias !!!
 * (coincidence of, e.g., a and c is logically possible but disallowed for 'restrict' optimization)
 *
 */
{
	c[0]=c1*a[0]+c2*b[0];
	c[1]=c1*a[1]+c2*b[1];
	c[2]=c1*a[2]+c2*b[2];
}

//============================================================

static inline double TrSym(const double a[static restrict 6])
// trace of a symmetric matrix stored as a vector of size 6
{
	return (a[0]+a[2]+a[5]);
}

//============================================================

static inline double QuadForm(const double matr[static 6],const double vec[static 3])
/* value of a quadratic form matr (symmetric matrix stored as a vector of size 6) over a vector vec;
 * matr and vec may alias
 */
{
	return ( vec[0]*vec[0]*matr[0] + vec[1]*vec[1]*matr[2] + vec[2]*vec[2]*matr[5]
	       + 2*(vec[0]*vec[1]*matr[1] + vec[0]*vec[2]*matr[3] + vec[1]*vec[2]*matr[4]) );
}

//============================================================
static inline void MatrVec(double matr[static restrict 3][3],const double vec[static restrict 3],
	double res[static restrict 3])
/* multiplication of matrix[3][3] by vec[3] (all real); res=matr.vec;
 * !!! none of inputs may alias !!!
 */
{
	res[0]=matr[0][0]*vec[0]+matr[0][1]*vec[1]+matr[0][2]*vec[2];
	res[1]=matr[1][0]*vec[0]+matr[1][1]*vec[1]+matr[1][2]*vec[2];
	res[2]=matr[2][0]*vec[0]+matr[2][1]*vec[1]+matr[2][2]*vec[2];
}

//============================================================

static inline void Permutate(double vec[static restrict 3],const int ord[static restrict 3])
/* permutate double vector vec using permutation ord; !!! vec and ord must not alias !!!
 * (coincidence of them is logically possible but disallowed for 'restrict' optimization)
 */
{
	double buf[3];

	memcpy(buf,vec,3*sizeof(double));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//============================================================

static inline void Permutate_i(int vec[static restrict 3],const int ord[static restrict 3])
/* permutate int vector vec using permutation ord; !!! vec and ord must not alias !!!
 * (coincidence of them is logically possible but disallowed for 'restrict' optimization)
 */
{
	int buf[3];

	memcpy(buf,vec,3*sizeof(int));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//============================================================
// Auxiliary functions

static inline double Deg2Rad(const double deg)
// transforms angle in degrees to radians
{
	return (PI_OVER_180*deg);
}

//============================================================

static inline double Rad2Deg(const double rad)
// transforms angle in radians to degrees
{
	return (INV_PI_180*rad);
}

#ifdef USE_SSE3

//============================================================

static inline __m128d cmul(__m128d a,__m128d b)
// complex multiplication
{
	__m128d t;
	t = _mm_movedup_pd(a);
	a = _mm_shuffle_pd(a,a,3);
	t = _mm_mul_pd(b,t);
	b = _mm_shuffle_pd(b,b,1);
	b = _mm_mul_pd(a,b);
	return _mm_addsub_pd(t,b);
}

//============================================================

static inline __m128d cadd(__m128d a,__m128d b)
// complex addition
{
	return _mm_add_pd(a,b);
}

#endif // USE_SSE3

#endif // __cmplx_h
