/* Functions used to calculate the interaction term
 *
 * Copyright (C) ADDA contributors
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
#include "interaction.h" // corresponding headers
// project headers
#include "cmplx.h"
#include "comm.h"
#include "igt_so.h"
#include "io.h"
#include "memory.h"
#include "somnec.h"
#include "vars.h"
// system headers
#include <float.h> // for DBL_EPSILON
#include <stdlib.h>

// GLOBAL VARIABLES

// function pointers, defined as extern in interaction.h
/* Calculates interaction term between two dipoles; given integer distance vector {i,j,k} (in units of d). The acting
 * dipole is placed in the origin, while the field is calculated at position given as argument. All six components of
 * the symmetric matrix are computed at once. The elements in result are: [G11, G12, G13, G22, G23, G33]
 */
void (*InterTerm_int)(const int i,const int j,const int k,doublecomplex result[static restrict 6]);
// same as above, but distance is passed as a double vector (in um)
void (*InterTerm_real)(const double qvec[static restrict 3],doublecomplex result[static restrict 6]);
/* Calculates reflection term between two dipoles; given integer distance vector {i,j,k} (in units of d). k is the _sum_
 * of dipole indices along z with respect to the center of bottom dipoles of the particle. Bottom is considered for the
 * current processor (position) and the whole particle (position_full) in FFT and SPARSE modes respectively. The latter
 * behavior is determined by ZsumShift.
 * The acting dipole is placed in the origin, while the field is calculated at position given as argument.
 * Six components of the matrix are computed at once: [GR11, GR12, GR13, GR22, GR23, GR33].
 * The matrix is not symmetric, but satisfies: GR21=GR12, GR31=-GR13, GR32=-GR23. However the large matrix GR, which
 * acts on the total vector of dipole polarizations is still complex-symmetric, since GR[i,j]=GR^T[j,i] (interchange of
 * i and j changes only the sign of x and y components, but not z, which leads to sign change of 13,31,23,32 components)
 */
void (*ReflTerm_int)(const int i,const int j,const int k,doublecomplex result[static restrict 6]);
/* same as above, but distance is passed as a double vector (in um) and its z-component is the sum of heights of
 * source and probe points above the surface. So qvec_in is the actual distance between probe point and the image of the
 * source point.
 */
void (*ReflTerm_real)(const double qvec[static restrict 3],doublecomplex result[static restrict 6]);

// SEMI-GLOBAL VARIABLES

// defined and initialized in make_particle.c
extern const double ZsumShift;
// defined and initialized in param.c
extern const double igt_lim,igt_eps,nloc_Rp;
extern const bool InteractionRealArgs;

// used in fft.c
int local_Nz_Rm; // number of local layers in Rmatrix, not greater than 2*boxZ-1 (also used in SPARSE)

// LOCAL VARIABLES

/* Performance of the functions in this file can be improved by using static variable declarations (either inside
 * functions or here for the whole source file. However, such optimizations are not thread-safe without additional'
 * efforts. So we prefer not to do it, since parallel (multi-threaded) execution of these functions is probably the
 * first step for the OpenMP implementation.
 */

// KroneckerDelta[mu,nu] - can serve both as multiplier, and as bool
static const double dmunu[6] = {1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
static doublecomplex surfRCn; // reflection coefficient for normal incidence
static bool XlessY; // whether boxX is not larger than boxY (used for SomTable)
static doublecomplex * restrict somTable; // table of Sommerfeld integrals
static size_t * restrict somIndex; // array for indexing somTable (in the xy-plane)

#ifdef USE_SSE3
static __m128d c1, c2, c3, zo, inv_2pi, p360, prad_to_deg;
static __m128d exptbl[361];
#endif //USE_SSE3

// EXTERNAL FUNCTIONS

#ifndef NO_FORTRAN
static double igtLimR2; // IGT distance threshold squared
// fort/propaesplibreintadda.f
void propaespacelibreintadda_(const double *Rij,const double *ka,const double *gridspacex,const double *gridspacey,
	const double *gridspacez,const double *relreq,double *result,int *ifail);
#endif
// sinint.c
void cisi(double x,double *ci,double *si);

// this is used for debugging, should be empty define, when not required
#define PRINT_GVAL /*printf("%s: %d,%d,%d: %g%+gi, %g%+gi, %g%+gi,\n%g%+gi, %g%+gi, %g%+gi\n",__func__,i,j,k,\
	REIM(result[0]),REIM(result[1]),REIM(result[2]),REIM(result[3]),REIM(result[4]),REIM(result[5]));*/

/* the following wrappers incur a lot of redundant compilation (since each of them inlines the same functional part),
 * however this way all the code is visible to compiler in one place, avoiding extra function calls and allowing full
 * optimizations. Overall, this looks like trying to implement C++ features in C.
 */
// wrapper for <name> (direct interaction), based on integer input; arguments are described in .h file
# define INT_WRAPPER_INTER(name) \
void name##_int(const int i,const int j,const int k,doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	UnitsGridToCoord(i,j,k,qvec); \
	name(qvec,result); }

// same as above, but calling function is different from name and accepts additional argument
# define INT_WRAPPER_INTER_3(name,func,arg) \
void name##_int(const int i,const int j,const int k,doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	UnitsGridToCoord(i,j,k,qvec); \
	func(qvec,result,arg); }

// wrapper for <name> (reflected interaction), based on integer input; arguments are described in .h file
# define INT_WRAPPER_REFL(name) \
void name##_int(const int i,const int j,const int k,doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	UnitsGridToCoordShift(i,j,k,qvec); \
	name(qvec,result); }

/* wrapper for <name>, based on real input; arguments are described in .h file; fine both for direct and reflected parts
 * to keep the input argument const, it has to be duplicated since some of the formulations of InterTerms may change it
 */
# define REAL_WRAPPER(name) \
void name##_real(const double qvec_in[static restrict 3],doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	vCopy(qvec_in,qvec); \
	name(qvec,result); }

// same as above, but calling function is different from name and accepts additional argument
# define REAL_WRAPPER_3(name,func,arg) \
void name##_real(const double qvec_in[static restrict 3],doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	vCopy(qvec_in,qvec); \
	func(qvec,result,arg); }

// aggregate defines
#define WRAPPERS_INTER(name) INT_WRAPPER_INTER(name) REAL_WRAPPER(name)
#define WRAPPERS_INTER_3(name,func,arg) INT_WRAPPER_INTER_3(name,func,arg) REAL_WRAPPER_3(name,func,arg)
#define WRAPPERS_REFL(name) INT_WRAPPER_REFL(name) REAL_WRAPPER(name)

/* this macro defines a void (error generating) real-input wrapper for Green's tensor formulations, which are, for
 * instance, based on tables. Doesn't generate 'int' wrapper, since the function itself is used for it.
 */
# define NO_REAL_WRAPPER(name) \
void name##_real(const double qvec_in[static restrict 3] ATT_UNUSED ,doublecomplex result[static restrict 6] ATT_UNUSED ) { \
	LogError(ALL_POS,"Function "#name" to compute dipole interaction does not support real input"); }

//=====================================================================================================================
/* The following two functions are used to do common calculation parts. For simplicity it is best to call them with the
 * same set of parameters - in this respect they can be replaced by macros. But inline functions are probably easier to
 * maintain.
 */

//=====================================================================================================================

static inline void UnitsGridToCoord(const int i,const int j,const int k,double qvec[static 3])
// convert integer coordinates (in dipole sizes) to real distance vector
{
	qvec[0]=i*dsX;
	qvec[1]=j*dsY;
	qvec[2]=k*dsZ;
}

//=====================================================================================================================

static inline void InterParams(double qvec[static 3],double qmunu[static 6],double *rr,double *invr3,double *kr,
	double *kr2)
/* some common variables needed by the interaction functions - needed for all except IGT
 * makes |qvec|=1
 * It will probably break down for qvec=0, so if any interaction term is required for zero argument, a new function need
 * to be created, which will avoid this singularity
 */
{
	double invr;

	*rr=vNorm(qvec);
	invr = 1.0/(*rr);
	vMultScal(invr,qvec,qvec); // finalize qvec (to have unity amplitude)
	*invr3=invr*invr*invr;
	*kr=WaveNum*(*rr);
	*kr2=(*kr)*(*kr);
	OuterSym(qvec,qmunu);
}

#ifdef USE_SSE3

//=====================================================================================================================

static inline __m128d accImExp_pd(const double x)
/* Accelerated sin-cos (or imaginary exp) routine for use in the calculation of the interaction tensor. Returns
 * c*exp(ix) in the resultant vector. The code is adapted from the CEPHES library. The idea is that we have
 * precalculated exp(iy) for some discrete values of y. Then we take the y that is nearest to our x and write
 * exp(ix)=exp(iy+(ix-iy))=exp(iy)exp(i(x-y)). We take exp(y) from the table and exp(i(x-y)) from the Taylor series.
 * This converges very fast since |x-y| is small.
 */
{
	__m128d px = _mm_set_sd(x);
	__m128d ipx = _mm_mul_pd(px,inv_2pi);
	int ix = _mm_cvttsd_si32(ipx); // truncation so the function only works for 0 <= x < 2*pi*2^-31
	ipx = _mm_cvtsi32_sd(ipx,ix);
	ipx = _mm_mul_pd(p360,ipx);
	px = _mm_mul_pd(prad_to_deg,px);
	px = _mm_sub_pd(px,ipx); // this is x (deg) mod 360

	ix = _mm_cvtsd_si32(px); // the closest integer (rounded, not truncated)
	__m128d scx = exptbl[ix]; // the tabulated value

	ipx = _mm_cvtsi32_sd(ipx,ix);
	__m128d pz = _mm_sub_pd(ipx,px); // the residual -z=(ix-x)
	__m128d py = _mm_mul_pd(pz,pz);
	__m128d yy = _mm_shuffle_pd(py,py,0); // now (y,y)
	__m128d zy = _mm_shuffle_pd(yy,pz,1);

	__m128d scz = _mm_mul_pd(c1,yy);	// Taylor series approximation
	scz = _mm_add_pd(c2,scz);
	scz = _mm_mul_pd(yy,scz);
	scz = _mm_add_pd(c3,scz);
	scz = _mm_mul_pd(zy,scz);
	scz = _mm_sub_pd(zo,scz);
	return cmul(scz,scx); // multiply lookup and approximation
}

//=====================================================================================================================

static inline doublecomplex accImExp(const double x)
{
	/* Here and further in the SSE3 part it is assumed that doublecomplex is equivalent to two doubles (that is
	 * specified by the C99 standard). Explicit pointer casts have been put in place, and pragmas to ignore remaining
	 * warnings from strict aliasing.
	 *
	 * !!! TODO: SSE3 code is a nice hack. But it should be considered carefully - is it worth it? After implementation
	 * of the tabulated imExp for the whole code, the further improvement from SSE3 is 10% (1.52 - 1.70 for matvec in
	 * test sparse runs). Moreover, if the following is replaced by call to standard imExp(), the timing is almost the
	 * same with slight improvement if the check for int overflow and negative numbers is turned off (as is the case in
	 *  SSE3 code).
	 *
	 * There seems to be some space for improvement in InterTerm_core (in comparison with SSE3) code, but otherwise we
	 * should move to remove SSE3 code for better maintainability. Interestingly, adding -msse3 to standard code
	 * compilation (without -DSSE3) doesn't help
	 */
	doublecomplex c;
	IGNORE_WARNING(-Wstrict-aliasing);
	_mm_store_pd((double *)(&c),accImExp_pd(x));
	STOP_IGNORE;
	return c;
}

//=====================================================================================================================

static inline void InterTerm_core(const double kr,const double kr2,const double invr3,const double qmunu[static 6],
	doublecomplex *expval,doublecomplex result[static 6])
// Core routine that calculates the point interaction term between two dipoles
{
	const __m128d ie = accImExp_pd(kr);
	const __m128d sc = _mm_mul_pd(_mm_set1_pd(invr3),ie);
	const double t1=(3-kr2), t2=-3*kr, t3=(kr2-1);
	const __m128d v1 = _mm_set_pd(kr,t3);
	const __m128d v2 = _mm_set_pd(t2,t1);
	__m128d qff,im_re;
	IGNORE_WARNING(-Wstrict-aliasing);
	_mm_store_pd((double *)expval,sc);
	STOP_IGNORE;

#undef INTERACT_MUL
#define INTERACT_DIAG(ind) { \
	qff = _mm_set1_pd(qmunu[ind]); \
	im_re = _mm_add_pd(v1,_mm_mul_pd(v2,qff)); \
	im_re = cmul(sc,im_re); \
	_mm_store_pd((double *)(result+ind),im_re); }
#define INTERACT_NONDIAG(ind) { \
	qff = _mm_set1_pd(qmunu[ind]); \
	im_re = _mm_mul_pd(v2,qff); \
	im_re = cmul(sc,im_re); \
	_mm_store_pd((double *)(result+ind),im_re); }

	IGNORE_WARNING(-Wstrict-aliasing);
	INTERACT_DIAG(0);    // xx
	INTERACT_NONDIAG(1); // xy
	INTERACT_NONDIAG(2); // xz
	INTERACT_DIAG(3);    // yy
	INTERACT_NONDIAG(4); // yz
	INTERACT_DIAG(5);    // zz
	STOP_IGNORE;

#undef INTERACT_DIAG
#undef INTERACT_NONDIAG
}

#else //not using SSE3

//=====================================================================================================================

static inline doublecomplex accImExp(const double x)
// Without SSE3, this is just an alias for imExp
{
	return imExp(x);
}

//=====================================================================================================================

static inline void InterTerm_core(const double kr,const double kr2,const double invr3,const double qmunu[static 6],
	doublecomplex *expval,doublecomplex result[static 6])
// Core routine that calculates the point interaction term between two dipoles
{
	const double t1=(3-kr2), t2=-3*kr, t3=(kr2-1);
	*expval=invr3*imExp(kr);

#define INTERACT_DIAG(ind) { result[ind] = ((t1*qmunu[ind]+t3) + I*(kr+t2*qmunu[ind]))*(*expval); }
#define INTERACT_NONDIAG(ind) { result[ind] = (t1+I*t2)*qmunu[ind]*(*expval); }

	INTERACT_DIAG(0);    // xx
	INTERACT_NONDIAG(1); // xy
	INTERACT_NONDIAG(2); // xz
	INTERACT_DIAG(3);    // yy
	INTERACT_NONDIAG(4); // yz
	INTERACT_DIAG(5);    // zz

#undef INTERACT_DIAG
#undef INTERACT_NONDIAG
}

#endif //USE_SSE3

//=====================================================================================================================

static inline void InterTerm_poi(double qvec[static 3],doublecomplex result[static 6])
/* Interaction term between two dipoles using the point-dipoles formulation;
 * qvec is the real distance, result is for produced output
 */
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr3,kr,kr2; // |R|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	InterParams(qvec,qmunu,&rr,&invr3,&kr,&kr2);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);
	PRINT_GVAL;
}

WRAPPERS_INTER(InterTerm_poi)

//=====================================================================================================================

static inline void InterTerm_fcd(double qvec[static 3],doublecomplex result[static 6])
/* Interaction term between two dipoles for FCD. See InterTerm_poi for more details.
 * !!! Works only for cubical dipoles, otherwise careful reconsideration of all formulae is required
 *
 * FCD is based on Gay-Balmaz P., Martin O.J.F. "A library for computing the filtered and non-filtered 3D Green's tensor
 * associated with infinite homogeneous space and surfaces", Comp. Phys. Comm. 144:111-120 (2002), and
 * Piller N.B. "Increasing the performance of the coupled-dipole approximation: A spectral approach",
 * IEEE Trans.Ant.Propag. 46(8): 1126-1137. Here it differs by a factor of 4*pi*k^2.
 *
 * speed of FCD can be improved by using faster version of sici routine, using predefined tables, etc (e.g. as is
 * done in GSL library). But currently extra time for this computation is already smaller than one main iteration.
 */
// If needed, it can be updated to work fine for qvec==0
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr3,kr,kr2; // |R|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double temp,kfr,ci,si,ci1,si1,ci2,si2,brd,g0,g2;
	int comp;
	doublecomplex eikfr; // exp(i*k_F*R)
	// next line should never happen
	if (rectDip) LogError(ONE_POS,"Incompatibility error in InterTerm_fcd");

	InterParams(qvec,qmunu,&rr,&invr3,&kr,&kr2);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	kfr=PI*rr/gridspace; // k_F*r, for FCD
	eikfr=accImExp(kfr);

	// ci,si_1,2 = ci,si_+,- = Ci,Si((k_F +,- k)r)
	cisi(kfr+kr,&ci1,&si1);
	cisi(kfr-kr,&ci2,&si2);
	// ci=ci1-ci2; si=pi-si1-si2
	ci=ci1-ci2;
	si=PI-si1-si2;
	g0=INV_PI*(cimag(expval)*ci+creal(expval)*si);
	g2=INV_PI*(kr*(creal(expval)*ci-cimag(expval)*si)+2*ONE_THIRD*invr3*(kfr*creal(eikfr)-4*cimag(eikfr)))-g0;
	temp=g0*kr2;
	for (comp=0;comp<NDCOMP;comp++) {
		// brd=(delta[mu,nu]*(-g0*kr^2-g2)+qmunu*(g0*kr^2+3g2))/r^3
		brd=qmunu[comp]*(temp+3*g2);
		if (dmunu[comp]) brd-=temp+g2;
		// result=Gp+brd; only the real part of the Green's tensor is corrected
		result[comp]+=brd;
	}
	PRINT_GVAL;
}

WRAPPERS_INTER(InterTerm_fcd)

//=====================================================================================================================

static inline void InterTerm_fcd_st(double qvec[static 3],doublecomplex result[static 6])
/* Interaction term between two dipoles for static FCD (in the limit of k->inf). See InterTerm_fcd for more details.
 * !!! Works only for cubical dipoles, otherwise careful reconsideration of all formulae is required
 */
// If needed, it can be updated to work fine for qvec==0
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr3,kr,kr2; // |R|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double kfr,ci,si,brd;
	int comp;
	doublecomplex eikfr;
	// next line should never happen
	if (rectDip) LogError(ONE_POS,"Incompatibility error in InterTerm_fcd_st");

	InterParams(qvec,qmunu,&rr,&invr3,&kr,&kr2);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	kfr=PI*rr/gridspace; // k_F*r, for FCD
	eikfr=accImExp(kfr);
	// result = Gp*[3*Si(k_F*r)+k_F*r*cos(k_F*r)-4*sin(k_F*r)]*2/(3*pi)
	cisi(kfr,&ci,&si);
	brd=TWO_OVER_PI*ONE_THIRD*(3*si+kfr*creal(eikfr)-4*cimag(eikfr));
	for (comp=0;comp<NDCOMP;comp++) result[comp]*=brd;
	PRINT_GVAL;
}

WRAPPERS_INTER(InterTerm_fcd_st)

//=====================================================================================================================

void InterTerm_igt_so(double qvec[static 3],doublecomplex result[static restrict 6])
{
	CalcIGTso(qvec,WaveNum,dsX,dsY,dsZ,result);
}

WRAPPERS_INTER(InterTerm_igt_so)

//=====================================================================================================================

static inline double lower_gamma52(const double sx,const double expMx)
/* computes the values of lower incomplete gamma function of order 5/2 at sx^2 by upward recursion from erf(x),
 * i.e. sx = sqrt(x) is given together with exp(-sx^2)=exp(-x)
 * it is numerically stable (up to all digits) for x>1
 * if extension to complex input is required, cerf can be obtained from http://apps.jcns.fz-juelich.de/doku/sc/libcerf
 */
{
	return 0.75*SQRT_PI*erf(sx) - sx*(1.5+sx*sx)*expMx;
}

//=====================================================================================================================

static inline double gamma_scaled(const double s,const double x,const double expMx)
/* computes the values of scaled lower incomplete gamma function, using given value of exp(-x)
 * f[s-1/2,x] = g[s,x]/x^s = M(s,s+1,-x)/s = exp(-x)M(1,s+1,x)/s, where M is Kummer's confluent hypergeometric function
 *
 * The code is based on Numerical Recipes, 3rd ed. There it is mentioned that such series are more efficient than
 * continuous fraction for x < s. Series is f[s-1/2,x] = exp(-x)*Sum{k=0->inf,Gamma(s)*x^k/Gamma(s+k+1)}.
 * For s=5/2 to reach double precision it requires from 9 to 16 iterations when 0.1<=x<=1 (and even smaller -
 * for smaller x)
 */
{
	double ap,del,sum;
	ap=s;
	del=sum=1.0/s;
	do {
		ap++;
		del*=x/ap;
		sum+=del;
	} while (fabs(del) > fabs(sum)*DBL_EPSILON);
	return sum*expMx;
}

//=====================================================================================================================

static inline double AverageGaussCube(double x)
/* computes one component of Gaussian averaging over a cube: (multiplied by 2)
 * [2/(sqrt(2pi)*d*Rp)]*Integral[exp(-(x+t)^2/(2Rp^2)),{t,-d/2,d/2}]
 * !!! works only for cubical dipoles (with single size gridspace)
 */
{
	double dif;
	if (x<0) x=-x; // for convenience work only with non-negative input
	const double y=x/(SQRT2*nloc_Rp);
	const double z=gridspace/(2*SQRT2*nloc_Rp);
	// two branches not to lose precision
	if (x<1) dif = erf(y+z) - erf(y-z);
	else dif = erfc(y-z) - erfc(y+z);
	return dif/gridspace;
}

//=====================================================================================================================

static inline void InterTerm_nloc_both(double qvec[static 3],doublecomplex result[static 6],const bool averageH)
/* Interaction term between two dipoles using the non-local interaction;
 * qvec is the real distance, result is for produced output
 * averageH specifies if the h function should be averaged over the dipole (cube) volume
 *
 * !!! Currently only static version is implemented; and only for cubical dipoles
 * !!! Mind the difference in sign with term in quantum-mechanical simulations, which defines the interaction energy
 * G = 4/[3sqrt(PI)R^3]g(5/2,x)[3(RR/R^2)-I] - (4pi/3)h(R), where x=R^2/(2Rp^2), g is lower incomplete  gamma-function.
 * For moderate x, those gamma functions can be easily expressed through erf by upward recursion, since
 * g(1/2,y^2) = sqrt(pi)erf(y) and g(s+1,x) = s*g(s,x) - exp(-x)x^s
 * For small x to save significant digits, we define f(m,x) = g(m+1/2,x)/x^(m+1/2) [no additional coefficient 1/2]
 * (computed by series representation), then g(5/2,x)/R^3 = f(2,x)/[sqrt(2)*Rp]^3
 *
 * (4pi/3)h(R)=exp(-x)*sqrt(2/pi)/(3*Rp^3) - is the non-locality function. If averageH, it is replaced by its integral
 * over cube, which can be easily expressed through erf.
 *
 * Currently this function is faster than FCD, so we should not worry about speed. But if needed, calculation of both
 * h(R) and its integrals can be optimized by tabulating 1D functions (either exp, or combinations of erf)
 */
// If needed, it can be updated to work fine for qvec==0
{
	// standard variable definitions used for functions InterParams
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,rn,invr3,kr,kr2; // |R|, |R/d|, |R|^-3, kR, (kR)^2

	double sx,x,t1,t2,expMx,invRp3;
	// next line should never happen
	if (rectDip) LogError(ONE_POS,"Incompatibility error in InterTerm_nloc_both");

	InterParams(qvec,qmunu,&rr,&invr3,&kr,&kr2);
	rn=rr/gridspace;
	if (nloc_Rp==0) {
		if (rr==0) LogError(ALL_POS,"Non-local interaction is not defined for both R and Rp equal to 0");
		t1=invr3;
		// when averaging, we check if the point r is inside the cube around r0. Conforms with general formula below
		if (averageH && rn<=SQRT3*0.5 && fabs(qvec[0])*rn<=0.5 && fabs(qvec[1])*rn<=0.5 && fabs(qvec[2])*rn<=0.5)
			t2=FOUR_PI_OVER_THREE;
		else t2=0;
	}
	else {
		invRp3=1/(nloc_Rp*nloc_Rp*nloc_Rp);
		sx=SQRT1_2*rr/nloc_Rp;
		x=sx*sx;
		expMx=exp(-sx*sx);
		// the threshold for x is somewhat arbitrary
		if (x>1) t1=(4/(3*SQRT_PI))*lower_gamma52(sx,expMx)*invr3;
		else t1=SQRT2_9PI*x*gamma_scaled(2.5,x,expMx)*invRp3;
		if (averageH)
			t2=PI_OVER_SIX*AverageGaussCube(qvec[0]*rr)*AverageGaussCube(qvec[1]*rr)*AverageGaussCube(qvec[2]*rr);
		else t2=expMx*invRp3*SQRT2_9PI;
	}

//#define INTERACT_DIAG(ind) { result[ind] = ((t1*qmunu[ind]+t3) + I*(kr+t2*qmunu[ind]))*(*expval); }
//#define INTERACT_NONDIAG(ind) { result[ind] = (t1+I*t2)*qmunu[ind]*(*expval); }
#define INTERACT_DIAG(ind) { result[ind] = 3*t1*qmunu[ind] - (t1+t2); }
#define INTERACT_NONDIAG(ind) { result[ind] = 3*t1*qmunu[ind]; }

	INTERACT_DIAG(0);    // xx
	INTERACT_NONDIAG(1); // xy
	INTERACT_NONDIAG(2); // xz
	INTERACT_DIAG(3);    // yy
	INTERACT_NONDIAG(4); // yz
	INTERACT_DIAG(5);    // zz

#undef INTERACT_DIAG
#undef INTERACT_NONDIAG

	PRINT_GVAL;
}

// wrappers both for nloc and nloc_av
WRAPPERS_INTER_3(InterTerm_nloc,InterTerm_nloc_both,false)
WRAPPERS_INTER_3(InterTerm_nloc_av,InterTerm_nloc_both,true)

//=====================================================================================================================
#ifndef NO_FORTRAN

static inline void InterTerm_igt(double qvec[static 3],doublecomplex result[static 6])
// Interaction term between two dipoles with integration of Green's tensor. See InterTerm_poi for more details.
{
	double tmp[12];
	int comp,ifail;

	if (igt_lim==UNDEF || DotProd(qvec,qvec)<=igtLimR2 ) {
		/* passing complex vectors from Fortran to C is not necessarily portable (at least requires extra effort in
		 * the Fortran code. So we do it through double. This is not bad for performance, since double is anyway used
		 * internally for integration in this Fortran routine.
		 */
		propaespacelibreintadda_(qvec,&WaveNum,&dsX,&dsY,&dsZ,&igt_eps,tmp,&ifail);
		if (ifail!=0) {
			if (ifail==1) LogWarning(EC_WARN,ALL_POS,"Failed to reach relative accuracy of %g for Green's tensor "
				"integration for distance "GFORMDEF3V,igt_eps,COMP3V(qvec));
			else LogWarning(EC_WARN,ALL_POS,"Error in Green's tensor integration (code %d - see dcuhre.f) for distance "
				GFORMDEF3V,ifail,COMP3V(qvec));
		}
		for (comp=0;comp<6;comp++) result[comp] = tmp[comp] + I*tmp[comp+6];
		PRINT_GVAL;

// test IGT_SO
//		doublecomplex result1[6];
//		InterTerm_igt_so(qvec, result1);
//		for (comp=0;comp<6;comp++) {
//			double norm_IGT = cAbs2(result[comp]);
//			double norm_SO = cAbs2(result1[comp]);
//			if (norm_IGT == 0 && fabs(norm_SO) > 1e-10) int rrrr = 1;
//			if (norm_IGT != 0) {
//				double errrr  = fabs(cAbs2(result[comp]-result1[comp]))/norm_IGT*100;
//				if (errrr > 0.0001) int rrrr = 1;
//			}
//		}

	}
	else InterTerm_poi(qvec,result);
}

WRAPPERS_INTER(InterTerm_igt)

#endif
/* TO ADD NEW INTERACTION FORMULATION
 * Add above functions that actually perform the calculation of the interaction term, according to the new formulae. At
 * the end you need to have two functions according to the declarations InterTerm_int and InterTerm_real in
 * interaction.h. Their input and output arguments are described there as well. Precise definition of Green's tensor is
 * given in the manual, in particular, it is based on CGS system of units. There are two ways to proceed:
 *
 * 1) Recommended. You create one main function with the following declaration
 * static inline void InterTerm_<name>(double qvec[static 3],doublecomplex result[static 6])
 * which works for double input vector. After the function definition you specify "WRAPPERS_INTER(InterTerm_<name>)",
 * which automatically produces declarations compatible to that in interaction.h. See InterTerm_poi() for example.
 *
 * 2) Create two separate functions named InterTerm_<name>_int and InterTerm_<name>_real with declarations described in
 * interaction.h. If the new formulation does not support arbitrary real input vector (e.g. it is based on tables), then
 * use "NO_REAL_WRAPPER(InterTerm_<name>)" instead of the real-input declaration.
 *
 * In any case you may benefit from existing utility functions InterParams() and InterTerm_core(). Adhering to the
 * naming conventions is important to be able to use macros here and below in InitInteraction().
 */

//=====================================================================================================================

static inline void UnitsGridToCoordShift(const int i,const int j,const int k,double qvec[static 3])
// convert integer coordinates (in dipole sizes) to real distance vector with extra shift along the z-axis
{
	qvec[0]=i*dsX;
	qvec[1]=j*dsY;
	qvec[2]=k*dsZ+ZsumShift;
}

//=====================================================================================================================

static inline void ReflTerm_core(const double kr,const double kr2,const double invr3,const double qmunu[static 6],
	doublecomplex *expval,doublecomplex result[static 6])
// Core routine that calculates the reflection interaction between probe dipole and image of source dipole
{
	// this is a modification of InterTerm_core, multiplying by refl. coef. and inverting z-components
	const double t1=(3-kr2), t2=-3*kr, t3=(kr2-1);
	*expval=invr3*imExp(kr);
	doublecomplex scale=surfRCn*(*expval);

#define INTERACT_DIAG(ind) { result[ind] = ((t1*qmunu[ind]+t3) + I*(kr+t2*qmunu[ind]))*scale; }
#define INTERACT_NONDIAG(ind) { result[ind] = (t1+I*t2)*qmunu[ind]*scale; }
	INTERACT_DIAG(0);    // xx
	INTERACT_NONDIAG(1); // xy
	INTERACT_NONDIAG(2); // xz
	INTERACT_DIAG(3);    // yy
	INTERACT_NONDIAG(4); // yz
	INTERACT_DIAG(5);    // zz
#undef INTERACT_DIAG
#undef INTERACT_NONDIAG
	// invert sign of *z components
	result[2]*=-1;
	result[4]*=-1;
	result[5]*=-1;
	PRINT_GVAL;
}
//=====================================================================================================================

static inline void ReflTerm_img(double qvec[static 3],doublecomplex result[static 6])
/* Reflection term using the image-dipole approximation;
 * qvec is the real distance between probe and image of source, so its z-component is the sum of heights of source and
 * probe points above the surface. result is for produced output
 */
{
	// standard variable definitions used for functions InterParams and ReflTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr3,kr,kr2; // |R|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	InterParams(qvec,qmunu,&rr,&invr3,&kr,&kr2);
	ReflTerm_core(kr,kr2,invr3,qmunu,&expval,result);
}

WRAPPERS_REFL(ReflTerm_img)

//=====================================================================================================================

static inline void SingleSomIntegral(double rho,const double z,doublecomplex vals[static 4])
/* computes a single Sommerfeld integral (4-element array); arguments are in real units
 * currently it is a wrapper around evlua, but in the future it may be replaced by use of interpolation table
 */
{
	// TODO: these scales can be removed by changes in som_init to use proper wavenumber instead of 2*pi
	const double scale=WaveNum/TWO_PI;
	const double isc=pow(scale,3); // this is subject to under/overflow

	evlua(z*scale,rho*scale,vals,vals+1,vals+2,vals+3,0);
	vals[0]*=isc;
	vals[1]*=isc;
	vals[2]*=isc;
	vals[3]*=isc;
}

//=====================================================================================================================

static void CalcSomTable(void)
/* calculates a table of (essential Sommerfeld integrals), which are further combined into reflected Green's tensor
 * For z values - all local grid; for x- and y-values only positive values are considered and additionally y<=x.
 *
 * That is good for FFT code, since all these values are required anyway. A minor improvement can be achieved by
 * locating different pairs of i,j that lead to the same rho (like 3,4 and 5,0), but the fraction of such matching pairs
 * is very small (also see below). Another way for improvement is to set a (coarser) 2D grid in plane of z-rho and
 * perform interpolation on it (as done in Schmehl's thesis). But that will add another free parameter affecting the
 * final accuracy.
 *
 * However, in sparse mode this procedure is inefficient, since incurs (potentially) a lot of unnecessary evaluations
 * of Sommerfeld integrals. For really sparse aggregates the better way is to buildup a lookup table, using only
 * actually used pairs of (z,rho), as is done in DDA-SI code. A hash table can be used for that.
 * TODO: Implement such improvement for the sparse mode (issue 175)
 *
 * Using only actually used value of (z,rho) can be also relevant for FFT mode (consider, e.g. a sphere and z close to 0
 * and to 2*boxZ-1). However, searching through such pairs seems to be O(N^2) operation, which is unacceptable in FFT
 * mode.
 *
 * Another problem in SPARSE mode is that currently all values of z are computed on each processor (in MPI mode). This
 * is related to the current parallelization mode (issue 160), but it can be improved by computing all values in chunks
 * and then gathering them on each processor. Better to combine it with the hash table above.
 */
{
	int i,j,k;
	double z;
	size_t ind;

	/* The logic below is heavily based on dsX=dsY. In principle, it can be extended to integer ratios, but doesn't seem
	 * worth the effort. First, it is hard to find interesting practical cases of particles on substrate that would
	 * benefit from dsX!=dsY. Second, the issue will be solved automatically once optimized routines are implemented
	 * (as discussed above).
	 */
	if (dsX!=dsY) LogError(ONE_POS,"Incompatibility error in CalcSomTable");
	XlessY=(boxX<=boxY);
	// create index for plane x,y; if boxX<=boxY the space above the main diagonal is indexed (so x<=y) and vice versa
	MALLOC_VECTOR(somIndex,sizet,boxY+1,ALL);
	memory+=(boxY+1)*sizeof(size_t);
	somIndex[0]=0;
	for (j=0;j<boxY;j++) somIndex[j+1]=somIndex[j] + (XlessY ? MIN(j+1,boxX) : (boxX-j));
	// allocate and fill the table
	const size_t tmp=4*local_Nz_Rm*somIndex[boxY];
	memory+=tmp*sizeof(doublecomplex);
	if (!prognosis) {
		MALLOC_VECTOR(somTable,complex,tmp,ALL);
		if (IFROOT) PRINTFB("Calculating table of Sommerfeld integrals\n");
		ind=0;
		for (k=0;k<local_Nz_Rm;k++) {
			z=k*dsZ+ZsumShift;
			for (j=0;j<boxY;j++) {
				if (XlessY) for (i=0;i<=j && i<boxX;i++,ind++) SingleSomIntegral(hypot(i*dsX,j*dsY),z,somTable+4*ind);
				else for (i=j;i<boxX;i++,ind++) SingleSomIntegral(hypot(i*dsX,j*dsY),z,somTable+4*ind);
			}
		}
	}
}

//=====================================================================================================================

static inline void CombineSomTensor(const double x, const double y,const double rho,const doublecomplex Ivals[static 4],
	doublecomplex result[static 6])
/* combine Sommerfeld integrals in the Sommerfeld tensor and add it to the result (increment)
 * x,y coordinates can be in any units (only relative values matter), rho should be equal to hypot(x,y)
 */
{
	double xr,yr; // relative (scaled by ro) transverse coordinates
	doublecomplex Irv,Izv,Irh,Iph; // values of Sommerfeld integrals

	// index tables
	Irv=Ivals[0];
	Izv=Ivals[1];
	Irh=Ivals[2];
	Iph=Ivals[3];

	if (rho==0) {
		/* particular direction doesn't matter in this case; this special case can be obtained from
		 * general case by taking, e.g. xr=1, yr=0 (since Iph+Irh=Irv=0 for rho=0).
		 */
		result[0] += Irh;
		result[3] += Irh;
		result[5] += Izv;
	}
	else {
		xr=x/rho;
		yr=y/rho;
		result[0] += xr*xr*Irh - yr*yr*Iph;
		result[1] += xr*yr*(Irh+Iph);
		result[2] += xr*Irv;
		result[3] += yr*yr*Irh - xr*xr*Iph;
		result[4] += yr*Irv;
		result[5] += Izv;
	}
}


//=====================================================================================================================

void ReflTerm_som_int(const int i,const int j,const int k,doublecomplex result[static restrict 6])
// Reflection term using the Sommerfeld integrals for integer input; arguments are described in .h file
// we do not use wrappers here, since both integer and real values are required
{
	// first, image-dipole part
	// standard variable definitions used for functions InterParams and ReflTerm_core
	double qvec[3],qmunu[6]; // distance vector (in units of d) and normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr3,kr,kr2; // |R|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3
	UnitsGridToCoordShift(i,j,k,qvec);
	InterParams(qvec,qmunu,&rr,&invr3,&kr,&kr2);
	ReflTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	// second, Sommerfeld integral part
	// compute table index
	int iT=abs(i);
	int jT=abs(j);
	size_t ind;
	// index for the table
	if (XlessY) {
		if (iT<=jT) ind=somIndex[jT]+iT;
		else ind=somIndex[iT]+jT; // effectively swap iT and jT
	}
	else {
		if (iT>=jT) ind=somIndex[jT]+iT-jT;
		else ind=somIndex[iT]+jT-iT; // effectively swap iT and jT
	}
	ind=4*(ind+k*somIndex[boxY]);
	double x=qvec[0];
	double y=qvec[1];
	double rho=hypot(x,y);
	CombineSomTensor(x,y,rho,somTable+ind,result);
	PRINT_GVAL;
}

//=====================================================================================================================

void ReflTerm_som_real(const double qvec[static restrict 3],doublecomplex result[static restrict 6])
// Reflection term using the Sommerfeld integrals for real input; arguments are described in .h file
{
	// first, image-dipole part
	// standard variable definitions used for functions InterParams and ReflTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,invr3,kr,kr2; // |R|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3
	double qv[3]; // copy of qvec
	vCopy(qvec,qv);
	InterParams(qv,qmunu,&rr,&invr3,&kr,&kr2);
	ReflTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	// second, Sommerfeld integral part
	double x=qvec[0];
	double y=qvec[1];
	double rho=hypot(x,y);
	doublecomplex Ivals[4];
	SingleSomIntegral(rho,qvec[2],Ivals);
	CombineSomTensor(x,y,rho,Ivals,result);
	PRINT_GVAL;
}

/* TO ADD NEW REFLECTION FORMULATION
 * Add here functions that actually perform the calculation of the reflection term, according to the new formulae. At
 * the end you need to have two functions according to the declarations ReflTerm_int and ReflTerm_real in
 * interaction.h. Their input and output arguments are described there as well. Precise definition of reflected Green's
 * tensor is given in the manual, in particular, it is based on CGS system of units. There are two ways to proceed:
 *
 * 1) Recommended. You create one main function with the following declaration
 * static inline void ReflTerm_<name>(double qvec[static 3],doublecomplex result[static 6])
 * which works for double input vector. After the function definition you specify "WRAPPERS_REFL(ReflTerm_<name>)",
 * which automatically produces declarations compatible to that in interaction.h. See ReflTerm_img() for example.
 *
 * 2) Create two separate functions named ReflTerm_<name>_int and ReflTerm_<name>_real with declarations described in
 * interaction.h. If the new formulation does not support arbitrary real input vector (e.g. it is based on tables), then
 * use "NO_REAL_WRAPPER(ReflTerm_<name>)" instead. See ReflTerm_som_int/real() for example.
 *
 * In any case you may benefit from existing utility functions InterParams() and ReflTerm_core(). Adhering to the
 * naming conventions is important to be able to use macros here and below in InitInteraction().
 */

//=====================================================================================================================

void InitInteraction(void)
// Initialize the interaction calculations
{
#define SET_FUNC_POINTERS(type,name) { type##_int = &type##_##name##_int; type##_real = &type##_##name##_real; }
	// set InterTerm_int (real) to point at the right functions
	switch (IntRelation) {
		case G_POINT_DIP: SET_FUNC_POINTERS(InterTerm,poi); break;
		case G_FCD: SET_FUNC_POINTERS(InterTerm,fcd); break;
		case G_FCD_ST: SET_FUNC_POINTERS(InterTerm,fcd_st); break;
		case G_IGT_SO: SET_FUNC_POINTERS(InterTerm,igt_so); break;
		case G_NLOC: SET_FUNC_POINTERS(InterTerm,nloc); break;
		case G_NLOC_AV: SET_FUNC_POINTERS(InterTerm,nloc_av); break;
#ifndef NO_FORTRAN
		case G_IGT:
			SET_FUNC_POINTERS(InterTerm,igt);
			// initialize IGT distance threshold
			if (igt_lim!=UNDEF) {
				igtLimR2=igt_lim*MAX(MAX(dsX,dsY),dsZ);
				igtLimR2*=igtLimR2;
			}
			break;
#endif
		/* TO ADD NEW INTERACTION FORMULATION
		 * Add here the assignment of function pointers for the new formulation. It is recommended to use special macro,
		 * assuming that you conformed to naming conventions for functions themselves. If new formulation requires
		 * initialization, add it here, preferably by a separate function. Initialization should honor the 'prognosis'
		 * flag. Additional memory should be counted always, but allocated only when not prognosis. If the new
		 * formulation does not support arbitrary real input vector (e.g. it is based on tables), then test for
		 * InteractionRealArgs and add an exception.
		 */
		default: LogError(ONE_POS, "Invalid interaction term calculation method: %d",(int)IntRelation);
			// no break
	}

	// Interaction through reflection from surface
	if (surface) {
#ifdef SPARSE
		local_Nz_Rm=2*boxZ-1;
#else
		local_Nz_Rm=MAX(MIN(2*local_z1,2*boxZ-1)-2*local_z0,0);
#endif
		switch (ReflRelation) {
			case GR_IMG:  SET_FUNC_POINTERS(ReflTerm,img); break;
			case GR_SOM:
				SET_FUNC_POINTERS(ReflTerm,som);
				if (!prognosis) som_init(msub*msub);
				CalcSomTable();
				break;
			/* TO ADD NEW REFLECTION FORMULATION
			 * Add here the assignment of function pointers for the new formulation. It is recommended to use special
			 * macro, assuming that you conformed to naming conventions for functions themselves. If new formulation
			 * requires initialization add it here, preferably by a separate function. Initialization should honor the
			 * 'prognosis' flag. Additional memory should be counted always, but allocated only when not prognosis. See
			 * CalcSomTable() for example. If the new formulation does not support arbitrary real input vector (e.g. it
			 * is based on tables), then test for InteractionRealArgs and add an exception.
			 */
			default: LogError(ONE_POS, "Invalid reflection term calculation method: %d",(int)ReflRelation);
				// no break
		}
		surfRCn=msubInf ? -1 : ((1-msub*msub)/(1+msub*msub));
	}

#ifdef USE_SSE3
	c1 = _mm_set_pd(1.34959795251974073996e-11,3.92582397764340914444e-14);
	c2 = _mm_set_pd(-8.86096155697856783296e-7,-3.86632385155548605680e-9);
	c3 = _mm_set_pd(1.74532925199432957214e-2,1.52308709893354299569e-4);
	zo = _mm_set_pd(0.0,1.0);
	inv_2pi = _mm_set_sd(1.0/(2*PI));
	p360 = _mm_set_sd(360.0);
	prad_to_deg = _mm_set_sd(180.0/PI);

	for (unsigned int i=0; i<=360; i++) {
		double x = (PI/180.0)*(double)i;
		exptbl[i] = _mm_set_pd(sin(x),cos(x));
	}
#endif

#undef SET_FUNC_POINTERS
}

//=====================================================================================================================

void FreeInteraction(void)
// Free buffers used for interaction calculation
{
	if (surface && ReflRelation==GR_SOM) {
		Free_general(somIndex);
		Free_cVector(somTable);
	}
	/* TO ADD NEW INTERACTION FORMULATION
	 * TO ADD NEW REFLECTION FORMULATION
	 * If you allocate any memory (for tables, etc.), free it here
	 */
}
