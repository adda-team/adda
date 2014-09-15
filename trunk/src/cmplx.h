/* File: cmplx.h
 * $Date::                            $
 * Descr: inline complex functions, functions on length-3 real and complex vectors, and several auxiliary functions
 *
 * Copyright (C) 2006-2008,2010,2012-2014 ADDA contributors
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
#ifndef __cmplx_h
#define __cmplx_h

// project headers
#include "const.h"    // for math constants
#include "types.h"    // for doublecomplex
// system headers
#include <math.h>     // for cos, sin
#include <string.h>   // for memcpy

#ifdef USE_SSE3
#include <xmmintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#endif

// useful macro for printing complex numbers and matrices
#define REIM(a) creal(a),cimag(a)
#define REIM3V(a) REIM((a)[0]),REIM((a)[1]),REIM((a)[2])

#ifdef ACCIMEXP
doublecomplex * restrict imexptable;
double xp;
int ixp;
double x;
doublecomplex imexp;
doublecomplex res;
double par[2];
#endif

//======================================================================================================================
// operations on complex numbers

static inline double cAbs2(const doublecomplex a)
// square of absolute value of complex number; |a|^2
{
	return creal(a)*creal(a) + cimag(a)*cimag(a);
}

//======================================================================================================================

static inline doublecomplex cSqrtCut(const doublecomplex a)
// square root of complex number, with explicit handling of branch cut (not to depend on sign of zero of imaginary part)
/* It is designed for calculating normal component of the transmitted wavevector when passing through the plane
 * interface. However, such choice of branch cut (while physically correct) leads to all kind of weird consequences.
 *
 * For instance, the electric field above the interface for plane wave propagating from a slightly absorbing substrate
 * at large incident angle (larger than critical angle for purely real refractive index) is unexpectedly large. This
 * happens because the wave in the vacuum is inhomogeneous and the real part of wavevector is almost parallel to the
 * surface. So the field above the surface actually comes from distant points on the surface, which has much larger
 * amplitude of the incident wave from below (compared to that under the observation point). Since the distance along
 * the surface (or the corresponding slope) is inversely proportional to the imaginary part of the substrate refractive
 * index, the effect remains finite even in the limit of absorption going to zero. Therefore, in this case there exist
 * a discontinuity when switching from non-absorbing to absorbing substrate. Physically, this fact is a consequence of
 * the infinite lateral extent of the plane wave.
 *
 * Exactly the same issue exist when scattering into the absorbing medium is calculated. At large scattering angles the
 * amplitude becomes very large, which also amplifies a lot the calculated Csca.
 */
{
	if (cimag(a)==0) {
		if (creal(a)>=0) return sqrt(a);
		else return I*sqrt(-a);
	}
	else return csqrt(a);
}

//======================================================================================================================

static inline doublecomplex imExp(const double arg)
// exponent of imaginary argument Exp(i*arg); optimization is performed by compiler
// this may be faster than using generic cexp, since imaginary type is not supported by all compilers
{
#ifndef ACCIMEXP
	return cos(arg) + I*sin(arg);
#else
	xp=arg*0.159154943091895;
	ixp=floor(xp);
	xp=ixp*FULL_ANGLE;
	x=arg*INV_PI_180;
	x=x-xp;	//x° mod 360°
	ixp=round(x);
	imexp=imexptable[ixp];	//first guess
	x=x-ixp;		// residual
	xp=x*x;
	
	par[0]=-3.925831985743094e-14*xp;
	par[1]=1.349601623163255e-11*xp;
	par[0]+=3.866323851562993e-09;
	par[1]+=-8.860961557012979e-07;
	par[0]*=xp;
	par[1]*=xp;
	par[0]+=-1.523087098933543e-04;
	par[1]+=1.74532925199432957214e-2;
	par[0]*=xp;
	par[1]*=x;
	par[0]+=1.0;
	res=par[0]+I*par[1];
	res=imexp*res; 
	return res;
#endif
}

//======================================================================================================================

static inline void imExp_arr(const doublecomplex arg,const int size,doublecomplex *c)
/* construct an array of exponent of imaginary argument c=Exp(i*k*arg), where k=0,1,...,size-1. arg can be complex.
 * Uses stable recurrence from Numerical Recipes. Optimization of the initial simultaneous calculation of sin and cos is
 * performed by compiler; It is assumed that size is at least 1
 */
{
	int k;
	double a,b;
	doublecomplex d,tmp;
	double re,im;

	re=creal(arg);
	im=cimag(arg);
	// handles real part, no special case for re=0
	c[0]=1;
	if (size>1) {
		// set a=2*sin^2(arg/2), b=sin(arg), d = 1 - exp(i*arg)
		a=sin(re/2);
		b=cos(re/2);
		b*=2*a;
		a*=2*a;
		d= a - I*b;
		// this a bit faster than in the main cycle
		c[1]=1-d;
		// main cycle
		for (k=2;k<size;k++) {
			/* potentially compiler may group terms to accelerate calculation but lose significant digits. We hope it
			 * doesn't happen, but it should not be a big problem anyway
			 */
			tmp=c[k-1]*d;
			c[k]=c[k-1]-tmp;
		}
	}
	// handles imaginary part
	if (im!=0) {
		a=exp(-fabs(im));
		if (im>0) for (k=1,b=a;k<size;k++) {
			c[k]*=b;
			b*=a;
		}
		else for (k=size-1,b=exp(-(size-1)*im);k>0;k--) {
			c[k]*=b;
			b*=a;
		}
	}
}

//======================================================================================================================
// operations on complex vectors

static inline void cvInit(doublecomplex a[static 3])
// set complex vector[3] to zero; a=0
{
	a[0] = 0;
	a[1] = 0;
	a[2] = 0;
}

//======================================================================================================================

static inline void vConj(const doublecomplex a[static 3],doublecomplex b[static 3])
// complex conjugate of vector[3]; b=a*
{
	b[0] = conj(a[0]);
	b[1] = conj(a[1]);
	b[2] = conj(a[2]);
}

//======================================================================================================================

static inline void vReal(const doublecomplex a[static 3],double b[static 3])
// takes real part of the complex vector; b=Re(a)
{
	b[0]=creal(a[0]);
	b[1]=creal(a[1]);
	b[2]=creal(a[2]);
}

//======================================================================================================================

static inline void cvBuildRe(const double a[static 3],doublecomplex b[static 3])
// builds complex vector from real part; b=a + i*0
{
	b[0]=a[0];
	b[1]=a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void vInvRefl_cr(const double a[static 3],doublecomplex b[static 3])
/* reflects real vector with respect to the xy-plane and then inverts it, equivalent to reflection about the z-axis
 * result is stored into the complex vector
 */
{
	b[0]=-a[0];
	b[1]=-a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void crCrossProd(const double a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// cross product of real and complex vector; c = a x b; !!! vectors must not alias
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

//======================================================================================================================

static inline void cvMultScal(const double a,const doublecomplex b[static 3],doublecomplex c[static 3])
// multiplication of vector[3] by real scalar; c=ab;
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline void cvMultScal_RVec(const doublecomplex a,const double b[static 3],doublecomplex c[static 3])
// complex scalar - real vector[3] multiplication; c=ab
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline void cvMultScal_cmplx(const doublecomplex a,const doublecomplex b[static 3],doublecomplex c[static 3])
// multiplication of vector[3] by complex scalar; c=ab
{
	c[0] = a*b[0];
	c[1] = a*b[1];
	c[2] = a*b[2];
}

//======================================================================================================================

static inline double cvNorm2(const doublecomplex a[static 3])
// square of the norm of a complex vector[3]
{
	return cAbs2(a[0]) + cAbs2(a[1]) + cAbs2(a[2]);
}


//======================================================================================================================

static inline doublecomplex cDotProd(const doublecomplex a[static 3],const doublecomplex b[static 3])
// conjugate dot product of two complex vector[3]; c=a.b = a[0]*b*[0]+...+a[2]*b*[2]; for one vector use cvNorm2
{
	return a[0]*conj(b[0]) + a[1]*conj(b[1]) + a[2]*conj(b[2]);
}

//======================================================================================================================

static inline double cDotProd_Im(const doublecomplex a[static 3],const doublecomplex b[static 3])
/* imaginary part of dot product of two complex vector[3]; c=Im(a.b)
 * It is not clear whether this way is faster than cimag(cDotProd)
 */
{
	return ( cimag(a[0])*creal(b[0]) - creal(a[0])*cimag(b[0])
	       + cimag(a[1])*creal(b[1]) - creal(a[1])*cimag(b[1])
	       + cimag(a[2])*creal(b[2]) - creal(a[2])*cimag(b[2]) );
}

//======================================================================================================================

static inline doublecomplex cDotProd_conj(const doublecomplex a[static 3],const doublecomplex b[static 3])
// dot product of two complex vector[3] without conjugation; a.(b*) = a[0]*b[0]+...+a[2]*b[2]
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

//======================================================================================================================

static inline void cvAdd(const doublecomplex a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// add two complex vector[3]; c=a+b;
{
	c[0] = a[0] + b[0];
	c[1] = a[1] + b[1];
	c[2] = a[2] + b[2];
}

//======================================================================================================================

static inline void cvSubtr(const doublecomplex a[static 3],const doublecomplex b[static 3],doublecomplex c[static 3])
// add two complex vector[3]; c=a-b;
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

//======================================================================================================================

static inline void cvAdd2Self(doublecomplex a[static 3],const doublecomplex b[static 3],const doublecomplex c[static 3])
// increment one complex vector[3] by sum of other two; a+=b+c
{
	a[0] += b[0] + c[0];
	a[1] += b[1] + c[1];
	a[2] += b[2] + c[2];
}

//======================================================================================================================

static inline doublecomplex crDotProd(const doublecomplex a[static 3],const double b[static 3])
// dot product of complex and real vectors[3]; a.b
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//======================================================================================================================

static inline double crDotProd_Re(const doublecomplex a[static 3],const double b[static 3])
// real part of dot product of complex and real vectors[3]; Re(a.b)
{
	return creal(a[0])*b[0] + creal(a[1])*b[1] + creal(a[2])*b[2];
}

//======================================================================================================================

static inline void cvLinComb1(const doublecomplex a[static 3],const doublecomplex b[static 3],const double c1,
	doublecomplex c[static 3])
// linear combination of complex vectors[3] with real coefficient; second coefficient is unity; c=c1*a+b
{
	c[0] = c1*a[0] + b[0];
	c[1] = c1*a[1] + b[1];
	c[2] = c1*a[2] + b[2];
}

//======================================================================================================================

static inline void cvLinComb1_cmplx(doublecomplex a[static 3],doublecomplex b[static 3],const doublecomplex c1,
	doublecomplex c[static 3])
// linear combination of complex vectors[3] with complex coefficients; second coefficient is unity; c=c1*a+b
{
	c[0] = c1*a[0] + b[0];
	c[1] = c1*a[1] + b[1];
	c[2] = c1*a[2] + b[2];
}

//======================================================================================================================

static inline void cSymMatrVec(const doublecomplex matr[static 6],const doublecomplex vec[static 3],
	doublecomplex res[static 3])
// multiplication of complex symmetric matrix[6] by complex vec[3]; res=matr.vec; !!! vec and res must not alias
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = matr[2]*vec[0] + matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cSymMatrVecReal(const doublecomplex matr[static 6],const double vec[static 3],
	doublecomplex res[static 3])
// multiplication of complex symmetric matrix[6] by  vec[3]; res=matr.vec; !!! vec and res must not alias
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = matr[2]*vec[0] + matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cReflMatrVec(const doublecomplex matr[static 6],const doublecomplex vec[static 3],
	doublecomplex res[static 3])
/* multiplication of matrix[6] by complex vec[3]; res=matr.vec; passed components are the same as for symmetric matrix:
 * 11,12,13,22,23,33, but the matrix has the following symmetry - M21=M12, M31=-M13, M32=-M23
 * !!! vec and res must not alias
 */
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = - matr[2]*vec[0] - matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================

static inline void cReflMatrVecReal(const doublecomplex matr[static 6],const double vec[static 3],
	doublecomplex res[static 3])
/* multiplication of matrix[6] by complex vec[3]; res=matr.vec; passed components are the same as for symmetric matrix:
 * 11,12,13,22,23,33, but the matrix has the following symmetry - M21=M12, M31=-M13, M32=-M23
 * !!! vec and res must not alias
 */
{
	res[0] = matr[0]*vec[0] + matr[1]*vec[1] + matr[2]*vec[2];
	res[1] = matr[1]*vec[0] + matr[3]*vec[1] + matr[4]*vec[2];
	res[2] = - matr[2]*vec[0] - matr[4]*vec[1] + matr[5]*vec[2];
}

//======================================================================================================================
// operations on real vectors

static inline void vInit(double a[static 3])
// set real vector[3] to zero; a=0
{
	a[0] = 0;
	a[1] = 0;
	a[2] = 0;
}

//======================================================================================================================

static inline void vCopy(const double a[static 3],double b[static 3])
// copies one vector into another; b=a
{
	// can be rewritten through memcpy, but compiler should be able to do it itself if needed
	b[0]=a[0];
	b[1]=a[1];
	b[2]=a[2];
}

//======================================================================================================================

static inline void vAdd(const double a[static 3],const double b[static 3],double c[static 3])
// adds two real vectors; c=a+b
{
	c[0]=a[0]+b[0];
	c[1]=a[1]+b[1];
	c[2]=a[2]+b[2];
}

//======================================================================================================================

static inline void vSubtr(const double a[static 3],const double b[static 3],double c[static 3])
// subtracts two real vectors; c=a-b
{
	c[0]=a[0]-b[0];
	c[1]=a[1]-b[1];
	c[2]=a[2]-b[2];
}

//======================================================================================================================

static inline void vInvSign(double a[static 3])
// inverts the sign in the double vector[3]
{
	a[0]=-a[0];
	a[1]=-a[1];
	a[2]=-a[2];
}

//======================================================================================================================

static inline void vRefl(const double inc[static 3],double ref[static 3])
// reflects the incident vector 'inc' with respect to the xy-plane (inverts z-component)
{
	ref[0]=inc[0];
	ref[1]=inc[1];
	ref[2]=-inc[2];
}

//======================================================================================================================

static inline void vMultScal(const double a,const double b[static 3],double c[static 3])
// multiplication of real vector by scalar; c=a*b;
{
	c[0]=a*b[0];
	c[1]=a*b[1];
	c[2]=a*b[2];
}

//======================================================================================================================

static inline void vMult(const double a[static 3],const double b[static 3],double c[static 3])
// multiplication of two vectors (by elements); c[i]=a[i]*b[i];
{
	c[0]=a[0]*b[0];
	c[1]=a[1]*b[1];
	c[2]=a[2]*b[2];
}

//======================================================================================================================

static inline bool vAlongZ(const double a[static 3])
// a robust (with respect to round-off errors) way to test that vector is along the z-axis (+ or -)
{
	return fabs(a[0])<ROUND_ERR && fabs(a[1])<ROUND_ERR;
}

//======================================================================================================================

static inline double DotProd(const double a[static 3],const double b[static 3])
// dot product of two real vectors[3]; use DotProd(x,x) to get squared norm
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//======================================================================================================================

static inline double DotProdSquare(const double a[static 3],const double b[static 3])
// dot product of element-wise squares of two real vectors[3]
{
	return a[0]*a[0]*b[0]*b[0] + a[1]*a[1]*b[1]*b[1] + a[2]*a[2]*b[2]*b[2];
}

//======================================================================================================================

static inline double vNorm(const double a[static 3])
// norm of a real vector[3]
{
	return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
}

//======================================================================================================================

static inline void CrossProd(const double a[static 3],const double b[static 3],double c[static 3])
// cross product of two real vectors; c = a x b; !!! vectors must not alias
{
	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];
}

//======================================================================================================================

static inline void vNormalize(double a[static 3])
// normalize real vector to have unit norm
{
	double c;
	c=1/sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
	a[0]*=c;
	a[1]*=c;
	a[2]*=c;
}

//======================================================================================================================

static inline void LinComb(const double a[static 3],const double b[static 3],const double c1,const double c2,
	double c[static 3])
// linear combination of real vectors[3]; c=c1*a+c2*b
{
	c[0]=c1*a[0]+c2*b[0];
	c[1]=c1*a[1]+c2*b[1];
	c[2]=c1*a[2]+c2*b[2];
}

//======================================================================================================================

static inline void OuterSym(const double a[static 3],double matr[static 6])
// outer product of real vector a with itself, stored in symmetric matrix matr
{
	matr[0] = a[0]*a[0];
	matr[1] = a[0]*a[1];
	matr[2] = a[0]*a[2];
	matr[3] = a[1]*a[1];
	matr[4] = a[1]*a[2];
	matr[5] = a[2]*a[2];
}

//======================================================================================================================

static inline double TrSym(const double a[static 6])
// trace of a symmetric matrix stored as a vector of size 6
{
	return (a[0]+a[2]+a[5]);
}

//======================================================================================================================

static inline double QuadForm(const double matr[static 6],const double vec[static 3])
// value of a quadratic form matr (symmetric matrix stored as a vector of size 6) over a vector vec;
{
	return ( vec[0]*vec[0]*matr[0] + vec[1]*vec[1]*matr[2] + vec[2]*vec[2]*matr[5]
	       + 2*(vec[0]*vec[1]*matr[1] + vec[0]*vec[2]*matr[3] + vec[1]*vec[2]*matr[4]) );
}

//======================================================================================================================

static inline void MatrVec(double matr[static 3][3],const double vec[static 3],double res[static 3])
// multiplication of matrix[3][3] by vec[3] (all real); res=matr.vec;
{
	res[0]=matr[0][0]*vec[0]+matr[0][1]*vec[1]+matr[0][2]*vec[2];
	res[1]=matr[1][0]*vec[0]+matr[1][1]*vec[1]+matr[1][2]*vec[2];
	res[2]=matr[2][0]*vec[0]+matr[2][1]*vec[1]+matr[2][2]*vec[2];
}

//======================================================================================================================

static inline void MatrColumn(double matr[static 3][3],const int ind,double vec[static 3])
// get ind's column of matrix[3][3] and store it into vec[3] (all real, ind starts from zero); vec=matr[.][ind];
{
	vec[0]=matr[0][ind];
	vec[1]=matr[1][ind];
	vec[2]=matr[2][ind];
}

//======================================================================================================================

static inline void Permutate(double vec[static 3],const int ord[static 3])
// permutate double vector vec using permutation ord
{
	double buf[3];

	memcpy(buf,vec,3*sizeof(double));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//======================================================================================================================

static inline void Permutate_i(int vec[static 3],const int ord[static 3])
// permutate int vector vec using permutation ord
{
	int buf[3];

	memcpy(buf,vec,3*sizeof(int));
	vec[0]=buf[ord[0]];
	vec[1]=buf[ord[1]];
	vec[2]=buf[ord[2]];
}

//======================================================================================================================
// Auxiliary functions

static inline double Deg2Rad(const double deg)
// transforms angle in degrees to radians
{
	return (PI_OVER_180*deg);
}

//======================================================================================================================

static inline double Rad2Deg(const double rad)
// transforms angle in radians to degrees
{
	return (INV_PI_180*rad);
}

//======================================================================================================================

static inline bool TestBelowDeg(const double deg)
/* tests if the direction is below the substrate using the degree theta in degrees;
 * if unsure (within rounded error) returns false (above)
 */
{
	return fabs(fmod(fabs(deg),360)-180) < 90*(1-ROUND_ERR);
}

//======================================================================================================================
// functions used for substrate

//======================================================================================================================

static inline doublecomplex FresnelRS(const doublecomplex ki,const doublecomplex kt)
/* reflection coefficient for s-polarized wave (E perpendicular to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 */
{
	return (ki-kt)/(ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelTS(const doublecomplex ki,const doublecomplex kt)
/* transmission coefficient for s-polarized wave (E perpendicular to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 */
{
	return 2*ki/(ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelRP(const doublecomplex ki,const doublecomplex kt,const doublecomplex mr)
/* reflection coefficient for p-polarized wave (E parallel to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 * mr is the ratio of refractive indices (mt/mi)
 */
{
	return (mr*mr*ki-kt)/(mr*mr*ki+kt);
}

//======================================================================================================================

static inline doublecomplex FresnelTP(const doublecomplex ki,const doublecomplex kt,const doublecomplex mr)
/* transmission coefficient for p-polarized wave (E parallel to the main plane),
 * ki,kt are normal (positive) components of incident and transmitted wavevector (with arbitrary mutual scaling)
 * mr is the ratio of refractive indices (mt/mi)
 */
{
	return 2*mr*ki/(mr*mr*ki+kt);
}

#ifdef USE_SSE3

//======================================================================================================================

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

//======================================================================================================================

static inline __m128d cadd(__m128d a,__m128d b)
// complex addition
{
	return _mm_add_pd(a,b);
}

#endif // USE_SSE3

#endif // __cmplx_h
