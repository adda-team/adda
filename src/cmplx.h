/* File: cmplx.h
 * $Author$
 * $Date::                            $
 * Descr: inline complex functions
 *        plus few auxiliary functions
 *
 *        'const' can be used for many more function variables, however it doesn't
 *        work in combination with 'doublecomplex *' or more nested lists. That seems
 *        to be a principal limitation of C standard (some compilers may work, some produce
 *        warnings)
 *        A few changes for reliability and stability were made according to the ideas of
 *        section 5.5 of the Numerical Recipes, 3rd edition
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __cmplx_h
#define __cmplx_h

#include <string.h>   // for memcpy
#include <math.h>     // for cos, sin
#include "const.h"    // for math constants
#include "types.h"    // for doublecomplex
#include "function.h" // for INLINE

//============================================================
// operations on complex numbers

INLINE void cEqual(const doublecomplex a,doublecomplex b)
// performs b=a
{
	memcpy(b,a,sizeof(doublecomplex));
}

//============================================================

INLINE double cAbs2(const doublecomplex a)
// square of absolute value of complex number; |a|^2
{
	return (a[RE]*a[RE] + a[IM]*a[IM]);
}

//============================================================

INLINE double cAbs(const doublecomplex a)
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

INLINE void cConj(const doublecomplex a,doublecomplex b)
// complex conjugate; b=a*
{
	b[RE] = a[RE];
	b[IM] = - a[IM];
}

//============================================================

INLINE void cAdd(const doublecomplex a,const doublecomplex b,doublecomplex c)
// add two complex numbers; c=a+b
{
	c[RE] = a[RE] + b[RE];
	c[IM] = a[IM] + b[IM];
}

//============================================================

INLINE void cSubtr(const doublecomplex a,const doublecomplex b,doublecomplex c)
// subtract two complex numbers; c=a-b
{
	c[RE] = a[RE] - b[RE];
	c[IM] = a[IM] - b[IM];
}

//============================================================

INLINE void cSquare(const doublecomplex a,doublecomplex b)
// square of complex number; b=a^2
{
	b[RE]=a[RE]*a[RE] - a[IM]*a[IM];
	b[IM]=2*a[IM]*a[RE];
}

//============================================================

INLINE void cMultReal(const double a,const doublecomplex b,doublecomplex c)
// complex multiplication by real; c=ab
{
	c[RE]=a*b[RE];
	c[IM]=a*b[IM];
}

//============================================================

INLINE void cMult_i(doublecomplex c)
// complex multiplication by i; c=i*c
{
	double tmp;
	tmp=c[RE];
	c[RE]=-c[IM];
	c[IM]=tmp;
}

//============================================================

INLINE void cMult_i2(doublecomplex a,doublecomplex b)
// complex multiplication by i; b=i*a; !!! b and c should be different !!!
{
	b[RE]=-a[IM];
	b[IM]=a[RE];
}
//============================================================

INLINE void cMult(const doublecomplex a,const doublecomplex b,doublecomplex c)
// complex multiplication; c=ab; !!! c should be different from a and b !!!
{
	c[RE]=a[RE]*b[RE] - a[IM]*b[IM];
	c[IM]=a[IM]*b[RE] + a[RE]*b[IM];
}

//============================================================

INLINE void cMultSelf(doublecomplex a,const doublecomplex b)
// complex multiplication; a*=b
{
	double tmp;

	tmp=a[RE];
	a[RE]=a[RE]*b[RE] - a[IM]*b[IM];
	a[IM]=a[IM]*b[RE] + tmp*b[IM];
}

//============================================================

INLINE double cMultConRe(const doublecomplex a,const doublecomplex b)
// complex multiplication; returns real(a*b_conjugated)
{
	return (a[RE]*b[RE] + a[IM]*b[IM]);
}

//============================================================

INLINE double cMultConIm(const doublecomplex a,const doublecomplex b)
// complex multiplication; returns imag(a*b_conjugated)
{
	return (a[IM]*b[RE] - a[RE]*b[IM]);
}

//============================================================

INLINE void cLinComb(const doublecomplex a,const doublecomplex b,
                     const double c1,const double c2,doublecomplex c)
// linear combination of two complex numbers; c=c1*a+c2*b
{
	c[RE]=c1*a[RE]+c2*b[RE];
	c[IM]=c1*a[IM]+c2*b[IM];
}

//============================================================

INLINE void cInvSign(doublecomplex a)
// change sign of complex number; a*=-1;
{
	a[RE] = - a[RE];
	a[IM] = - a[IM];
}

//============================================================

INLINE void cInvSign2(const doublecomplex a,doublecomplex b)
// change sign of complex number and store to different address; b=-a;
{
	b[RE] = - a[RE];
	b[IM] = - a[IM];
}

//============================================================

INLINE void cInv(const doublecomplex a,doublecomplex b)
// complex inversion; b=1/a; designed to avoid under and overflows
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

INLINE double cInvIm(const doublecomplex a)
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

INLINE void cDiv(const doublecomplex a,const doublecomplex b,doublecomplex c)
/* complex division; c=a/b; designed to avoid under and overflows
 * !!! c should be different from a !!!
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

INLINE void cDivSelf(doublecomplex a,const doublecomplex b)
// complex division; a/=b; designed to avoid under and overflows
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

INLINE void cSqrt(const doublecomplex a,doublecomplex b)
/* complex square root; b=sqrt(a); designed to avoid under and overflows;
 * branch cut discontinuity is (-inf,0) - b[RE]>=0
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

INLINE void imExp(const double arg,doublecomplex c)
// exponent of imaginary argument c=Exp(i*arg); optimization is performed by compiler
{
	c[RE]=cos(arg);
	c[IM]=sin(arg);
}

//============================================================

INLINE void imExp_arr(const double arg,const int size,doublecomplex *c)
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

INLINE void cExp(const doublecomplex arg,doublecomplex c)
/* complex exponent of complex argument c=Exp(arg); optimization is performed by compiler;
 * !!! c should be different from arg !!!
 */
{
	c[RE]=c[IM]=exp(arg[RE]);
	c[RE]*=cos(arg[IM]);
	c[IM]*=sin(arg[IM]);
}

//============================================================

INLINE void cExpSelf(doublecomplex arg)
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

INLINE void cvMultScal(const double a,doublecomplex *b,doublecomplex *c)
// multiplication of vector by real scalar; c=ab
{
	c[0][RE] = a*b[0][RE];
	c[0][IM] = a*b[0][IM];
	c[1][RE] = a*b[1][RE];
	c[1][IM] = a*b[1][IM];
	c[2][RE] = a*b[2][RE];
	c[2][IM] = a*b[2][IM];
}

//============================================================

INLINE void cScalMultRVec(const double *a,const doublecomplex b,doublecomplex *c)
// complex scalar- real vector[3] multiplication; c=b*a
{
	c[0][RE] = b[RE]*a[0];
	c[0][IM] = b[IM]*a[0];
	c[1][RE] = b[RE]*a[1];
	c[1][IM] = b[IM]*a[1];
	c[2][RE] = b[RE]*a[2];
	c[2][IM] = b[IM]*a[2];
}

//============================================================

INLINE void cvMultScal_cmplx(const doublecomplex a,doublecomplex *b,doublecomplex *c)
// multiplication of vector[3] by complex scalar; c=ab
{
	c[0][RE] = a[RE]*b[0][RE] - a[IM]*b[0][IM];
	c[0][IM] = a[RE]*b[0][IM] + a[IM]*b[0][RE];
	c[1][RE] = a[RE]*b[1][RE] - a[IM]*b[1][IM];
	c[1][IM] = a[RE]*b[1][IM] + a[IM]*b[1][RE];
	c[2][RE] = a[RE]*b[2][RE] - a[IM]*b[2][IM];
	c[2][IM] = a[RE]*b[2][IM] + a[IM]*b[2][RE];
}

//============================================================

INLINE double cvNorm2(doublecomplex *a)
// square of the norm of a complex vector[3]
{
	return ( a[0][RE]*a[0][RE] + a[0][IM]*a[0][IM]
	       + a[1][RE]*a[1][RE] + a[1][IM]*a[1][IM]
	       + a[2][RE]*a[2][RE] + a[2][IM]*a[2][IM] );
}


//============================================================

INLINE void cDotProd(doublecomplex *a,doublecomplex *b,doublecomplex c)
// conjugate dot product of two complex vector[3]; c=a.b = a[0]*b*[0]+...+a[2]*b*[2]*/
{
	c[RE] = a[0][RE]*b[0][RE] + a[0][IM]*b[0][IM]
	      + a[1][RE]*b[1][RE] + a[1][IM]*b[1][IM]
	      + a[2][RE]*b[2][RE] + a[2][IM]*b[2][IM];
	c[IM] = a[0][IM]*b[0][RE] - a[0][RE]*b[0][IM]
	      + a[1][IM]*b[1][RE] - a[1][RE]*b[1][IM]
	      + a[2][IM]*b[2][RE] - a[2][RE]*b[2][IM];
}

//============================================================

INLINE double cDotProd_Re(doublecomplex *a,doublecomplex *b)
// real part of dot product of two complex vector[3]; c=Re(a.b)
{
	return ( a[0][RE]*b[0][RE] + a[0][IM]*b[0][IM]
	       + a[1][RE]*b[1][RE] + a[1][IM]*b[1][IM]
	       + a[2][RE]*b[2][RE] + a[2][IM]*b[2][IM] );
}

//============================================================

INLINE double cDotProd_Im(doublecomplex *a,doublecomplex *b)
// imaginary part of dot product of two complex vector[3]; c=Im(a.b)
{
	return ( a[0][IM]*b[0][RE] - a[0][RE]*b[0][IM]
	       + a[1][IM]*b[1][RE] - a[1][RE]*b[1][IM]
	       + a[2][IM]*b[2][RE] - a[2][RE]*b[2][IM] );
}

//============================================================

INLINE void cDotProd_conj(doublecomplex *a,doublecomplex *b,doublecomplex c)
// dot product of two complex vector[3]; c=a.b* = a[0]*b[0]+...+a[2]*b[2]
{
	c[RE] = a[0][RE]*b[0][RE] - a[0][IM]*b[0][IM]
	      + a[1][RE]*b[1][RE] - a[1][IM]*b[1][IM]
	      + a[2][RE]*b[2][RE] - a[2][IM]*b[2][IM];
	c[IM] = a[0][IM]*b[0][RE] + a[0][RE]*b[0][IM]
	      + a[1][IM]*b[1][RE] + a[1][RE]*b[1][IM]
	      + a[2][IM]*b[2][RE] + a[2][RE]*b[2][IM];
}

//============================================================

INLINE double cDotProd_conj_Re(doublecomplex *a,doublecomplex *b)
// real part of dot product of two complex vector[3]; c=Re(a.b*)
{
	return ( a[0][RE]*b[0][RE] - a[0][IM]*b[0][IM]
	       + a[1][RE]*b[1][RE] - a[1][IM]*b[1][IM]
	       + a[2][RE]*b[2][RE] - a[2][IM]*b[2][IM] );
}

//============================================================

INLINE double cDotProd_conj_Im(doublecomplex *a,doublecomplex *b)
// imaginary part of dot product of two complex vector[3]; c=Im(a.b*)
{
	return ( a[0][IM]*b[0][RE] + a[0][RE]*b[0][IM]
	       + a[1][IM]*b[1][RE] + a[1][RE]*b[1][IM]
	       + a[2][IM]*b[2][RE] + a[2][RE]*b[2][IM] );
}

//============================================================

INLINE void cvAdd(doublecomplex *a,doublecomplex *b,doublecomplex *c)
// add two complex vector[3]; c=a+b
{
	c[0][RE] = a[0][RE] + b[0][RE];
	c[0][IM] = a[0][IM] + b[0][IM];
	c[1][RE] = a[1][RE] + b[1][RE];
	c[1][IM] = a[1][IM] + b[1][IM];
	c[2][RE] = a[2][RE] + b[2][RE];
	c[2][IM] = a[2][IM] + b[2][IM];
}

//============================================================

INLINE void cvSubtr(doublecomplex *a,doublecomplex *b,doublecomplex *c)
// subtraction of two complex vector[3]; c=a-b
{
	c[0][RE] = a[0][RE] - b[0][RE];
	c[0][IM] = a[0][IM] - b[0][IM];
	c[1][RE] = a[1][RE] - b[1][RE];
	c[1][IM] = a[1][IM] - b[1][IM];
	c[2][RE] = a[2][RE] - b[2][RE];
	c[2][IM] = a[2][IM] - b[2][IM];
}

//============================================================

INLINE void crDotProd(doublecomplex *a,const double *b,doublecomplex c)
// dot product of complex and real vectors[3]; c=a.b
{
	c[RE] = a[0][RE]*b[0] + a[1][RE]*b[1] + a[2][RE]*b[2];
	c[IM] = a[0][IM]*b[0] + a[1][IM]*b[1] + a[2][IM]*b[2];
}

//============================================================

INLINE double crDotProd_Re(doublecomplex *a,const double *b)
// real part of dot product of complex and real vectors[3]; c=Re(a.b)
{
	return (a[0][RE]*b[0] + a[1][RE]*b[1] + a[2][RE]*b[2]);
}

//============================================================

INLINE double crDotProd_Im(doublecomplex *a,const double *b)
// imaginary part of dot product of complex and real vectors[3]; c=Im(a.b)
{
	return (a[0][IM]*b[0] + a[1][IM]*b[1] + a[2][IM]*b[2]);
}

//============================================================

INLINE void cvIncremScaled_cmplx(doublecomplex *a,const doublecomplex b,doublecomplex *c)
// increment of complex vectors[3] by complex-scaled other vector; c+=b*a
{
	c[0][RE] += b[RE]*a[0][RE] - b[IM]*a[0][IM];
	c[0][IM] += b[RE]*a[0][IM] + b[IM]*a[0][RE];
	c[1][RE] += b[RE]*a[1][RE] - b[IM]*a[1][IM];
	c[1][IM] += b[RE]*a[1][IM] + b[IM]*a[1][RE];
	c[2][RE] += b[RE]*a[2][RE] - b[IM]*a[2][IM];
	c[2][IM] += b[RE]*a[2][IM] + b[IM]*a[2][RE];
}
//============================================================

INLINE void cvMultAdd(doublecomplex *a,const doublecomplex b,doublecomplex *c)
// multiply complex vectors[3] with complex constant and add another vector; c=b*c+a
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

INLINE void cvLinComb1(doublecomplex *a,doublecomplex *b,
                       const double c1,doublecomplex *c)
// linear combination of complex vectors[3]; second coefficient is unity; c=c1*a+b
{
	c[0][RE] = c1*a[0][RE] + b[0][RE];
	c[0][IM] = c1*a[0][IM] + b[0][IM];
	c[1][RE] = c1*a[1][RE] + b[1][RE];
	c[1][IM] = c1*a[1][IM] + b[1][IM];
	c[2][RE] = c1*a[2][RE] + b[2][RE];
	c[2][IM] = c1*a[2][IM] + b[2][IM];
}

//============================================================

INLINE void cvLinComb1_cmplx(doublecomplex *a,doublecomplex *b,
                             const doublecomplex c1,doublecomplex *c)
/* linear combination of complex vectors[3] with complex coefficients;
 * second coefficient is unity; c=c1*a+b   !!! c!=a
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

INLINE void cSymMatrVec(doublecomplex *matr,doublecomplex *vec,doublecomplex *res)
// multiplication of complex symmetric matrix[6] by complex vec[3]; res=matr.vec
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

INLINE void MultScal(const double a,const double *b,double *c)
// multiplication of vector by scalar; c=a*b
{
	c[0]=a*b[0];
	c[1]=a*b[1];
	c[2]=a*b[2];
}

//============================================================

INLINE void vMult(const double *a,const double *b,double *c)
// multiplication of two vectors (by elements); c[i]=a[i]*b[i]
{
	c[0]=a[0]*b[0];
	c[1]=a[1]*b[1];
	c[2]=a[2]*b[2];
}

//============================================================

INLINE double DotProd(const double *a,const double *b)
// dot product of two real vectors[3]
{
	return (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]);
}

//============================================================

INLINE void LinComb(const double *a,const double *b,
	const double c1,const double c2, double *c)
// linear combination of real vectors[3]; c=c1*a+c2*b
{
	c[0]=c1*a[0]+c2*b[0];
	c[1]=c1*a[1]+c2*b[1];
	c[2]=c1*a[2]+c2*b[2];
}

//============================================================

INLINE double TrSym(const double *a)
// trace of a symmetric matrix stored as a vector of size 6
{
	return (a[0]+a[2]+a[5]);
}

//============================================================

INLINE double QuadForm(const double *matr,const double *vec)
// value of a quadratic form matr (symmetric matrix stored as a vector of size 6) over a vector vec
{
	return ( vec[0]*vec[0]*matr[0] + vec[1]*vec[1]*matr[2] + vec[2]*vec[2]*matr[5]
	       + 2*(vec[0]*vec[1]*matr[1] + vec[0]*vec[2]*matr[3] + vec[1]*vec[2]*matr[4]) );
}

//============================================================
INLINE void MatrVec(double matr[3][3],const double *vec, double *res)
// multiplication of matrix[3][3] by vec[3] (all real); res=matr.vec
{
	res[0]=matr[0][0]*vec[0]+matr[0][1]*vec[1]+matr[0][2]*vec[2];
	res[1]=matr[1][0]*vec[0]+matr[1][1]*vec[1]+matr[1][2]*vec[2];
	res[2]=matr[2][0]*vec[0]+matr[2][1]*vec[1]+matr[2][2]*vec[2];
}

//============================================================

INLINE void Permutate(double *vec,const int *ord)
// permutate double vector vec using permutation ord
{
	double buf[3];

	memcpy(buf,vec,3*sizeof(double));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//============================================================

INLINE void Permutate_i(int *vec,const int *ord)
// permutate int vector vec using permutation ord
{
	int buf[3];

	memcpy(buf,vec,3*sizeof(int));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//============================================================
// Auxiliary functions

INLINE double Deg2Rad(const double deg)
// transforms angle in degrees to radians
{
	return (PI_OVER_180*deg);
}

//============================================================

INLINE double Rad2Deg(const double rad)
// transforms angle in radians to degrees
{
	return (INV_PI_180*rad);
}

#endif // __cmplx_h
