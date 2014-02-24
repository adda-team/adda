/* FILE : interaction.c
 * Descr: the functions used to calculate the interaction term
 *
 * Copyright (C) 2011-2014 ADDA contributors
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
#include "io.h"
#include "memory.h"
#include "vars.h"
// system headers
#include <float.h> // for DBL_EPSILON
#include <stdlib.h>

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

// tables of integrals
static double * restrict tab1,* restrict tab2,* restrict tab3,* restrict tab4,* restrict tab5,* restrict tab6,
	* restrict tab7,* restrict tab8,* restrict tab9,* restrict tab10;
/* it is preferable to declare the following as "* restrict * restrict", but it is hard to make it
* generally compatible with Free_iMatrix function syntax.
*/
static int ** restrict tab_index; // matrix for indexing of table arrays
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
// fort/propaesplibreintadda.f
void propaespacelibreintadda_(const double *Rij,const double *ka,const double *arretecube,const double *relreq,
	double *result);
#endif
// sinint.c
void cisi(double x,double *ci,double *si);
// somnec.c
void som_init(complex double epscf);
void evlua(double zphIn,double rhoIn,complex double *erv, complex double *ezv,complex double *erh, complex double *eph);

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
	vCopyIntReal(i,j,k,qvec); \
	name(qvec,result,true); }

// same as above, but calling function is different from name and accepts additional argument
# define INT_WRAPPER_INTER_3(name,func,arg) \
void name##_int(const int i,const int j,const int k,doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	vCopyIntReal(i,j,k,qvec); \
	func(qvec,result,true,arg); }

// wrapper for <name> (reflected interaction), based on integer input; arguments are described in .h file
# define INT_WRAPPER_REFL(name) \
void name##_int(const int i,const int j,const int k,doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	vCopyIntRealShift(i,j,k,qvec); \
	name(qvec,result,true); }

/* wrapper for <name>, based on real input; arguments are described in .h file
 * to keep the input argument const, it has to be duplicated since some of the formulations of InterTerms may change it
 */
# define REAL_WRAPPER_INTER(name) \
void name##_real(const double qvec_in[restrict 3],doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	vCopy(qvec_in,qvec); \
	name(qvec,result,false); }

// same as above, but calling function is different from name and accepts additional argument
# define REAL_WRAPPER_INTER_3(name,func,arg) \
void name##_real(const double qvec_in[restrict 3],doublecomplex result[static restrict 6]) { \
	double qvec[3]; \
	vCopy(qvec_in,qvec); \
	func(qvec,result,false,arg); }

// wrapper for <name>, based on real input; arguments are described in .h file
# define REAL_WRAPPER_REFL(name) \
void name##_real(const double qvec[restrict 3],doublecomplex result[static restrict 6]) \
	{ name(qvec,result,false); }

// aggregate defines
#define WRAPPERS_INTER(name) INT_WRAPPER_INTER(name) REAL_WRAPPER_INTER(name)
#define WRAPPERS_INTER_3(name,func,arg) INT_WRAPPER_INTER_3(name,func,arg) REAL_WRAPPER_INTER_3(name,func,arg)
#define WRAPPERS_REFL(name) INT_WRAPPER_REFL(name) REAL_WRAPPER_REFL(name)

/* this macro defines a void (error generating) real-input wrapper for Green's tensor formulations, which are, for
 * instance, based on tables. Doesn't generate 'int' wrapper, since the function itself is used for it.
 */
# define NO_REAL_WRAPPER(name) \
void name##_real(const double qvec_in[restrict 3] ATT_UNUSED ,doublecomplex result[static restrict 6] ATT_UNUSED ) { \
	LogError(ALL_POS,"Function "#name" to compute dipole interaction does not support real input"); }

//=====================================================================================================================
/* The following two functions are used to do common calculation parts. For simplicity it is best to call them with the
 * same set of parameters - in this respect they can be replaced by macros. But inline functions are probably easier to
 * maintain.
 */

//=====================================================================================================================

static inline void vCopyIntReal(const int i,const int j,const int k,double qvec[static 3])
// initialize real vector with integer values
{
	qvec[0]=i;
	qvec[1]=j;
	qvec[2]=k;
}

//=====================================================================================================================

static inline void InterParams(double qvec[static 3],double qmunu[static 6],double *rr,double *rn,double *invr3,
	double *kr,double *kr2,const bool unitsGrid)
/* some common variables needed by the interaction functions - needed for all except IGT
 * It will probably break down for qvec=0, so if any interaction term is required for zero argument, a new function need
 * to be created, which will avoid this singularity
 */
{
	double invrn;

	*rn=vNorm(qvec);
	invrn = 1.0/(*rn);
	vMultScal(invrn,qvec,qvec); // finalize qvec (to have unity amplitude)
	if (unitsGrid) (*rr)=(*rn)*gridspace;
	else {
		(*rr)=(*rn);
		(*rn)=(*rr)/gridspace;
	}
	*invr3=1/((*rr)*(*rr)*(*rr));
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
	 * !!! TODO: SSE3 code is a nice hack. But it should be considered carefully - is it worth it? In particular, it
	 * seems that only parts of it are really beneficial (like tabulated evaluation of imaginary exponents), and those
	 * can be incorporated into the main code (using standard C99 only).
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

static inline void InterTerm_poi(double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid)
/* Interaction term between two dipoles using the point-dipoles formulation;
 * qvec is  a distance given by either integer-valued vector (in units of d) or arbitrary-valued in real units (um),
 * controlled by unitsGrid (true of false respectively), result is for produced output
 */
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,rn,invr3,kr,kr2; // |R|, |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	InterParams(qvec,qmunu,&rr,&rn,&invr3,&kr,&kr2,unitsGrid);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);
	PRINT_GVAL;
}

WRAPPERS_INTER(InterTerm_poi)

//=====================================================================================================================

static inline void InterTerm_fcd(double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid)
/* Interaction term between two dipoles for FCD. See InterTerm_poi for more details.
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
	double rr,rn,invr3,kr,kr2; // |R|, |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double temp,kfr,ci,si,ci1,si1,ci2,si2,brd,g0,g2;
	int comp;
	doublecomplex eikfr; // exp(i*k_F*R)

	InterParams(qvec,qmunu,&rr,&rn,&invr3,&kr,&kr2,unitsGrid);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	kfr=PI*rn; // k_F*r, for FCD
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

static inline void InterTerm_fcd_st(double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid)
// Interaction term between two dipoles for static FCD (in the limit of k->inf). See InterTerm_fcd for more details.
// If needed, it can be updated to work fine for qvec==0
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,rn,invr3,kr,kr2; // |R|, |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double kfr,ci,si,brd;
	int comp;
	doublecomplex eikfr;

	InterParams(qvec,qmunu,&rr,&rn,&invr3,&kr,&kr2,unitsGrid);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	kfr=PI*rn; // k_F*r, for FCD
	eikfr=accImExp(kfr);
	// result = Gp*[3*Si(k_F*r)+k_F*r*cos(k_F*r)-4*sin(k_F*r)]*2/(3*pi)
	cisi(kfr,&ci,&si);
	brd=TWO_OVER_PI*ONE_THIRD*(3*si+kfr*creal(eikfr)-4*cimag(eikfr));
	for (comp=0;comp<NDCOMP;comp++) result[comp]*=brd;
	PRINT_GVAL;
}

WRAPPERS_INTER(InterTerm_fcd_st)

//=====================================================================================================================

static inline bool TestTableSize(const double rn)
// tests if rn fits into the table; if not, returns false and produces info message
{
	static bool warned=false;

	if (rn>TAB_RMAX) {
		if (!warned) {
			warned=true;
			LogWarning(EC_INFO,ONE_POS,"Not enough table size (available only up to R/d=%d), so O(kd^2) accuracy of "
				"Green's function is not guaranteed",TAB_RMAX);
		}
		return false;
	}
	else return true;
}

//=====================================================================================================================

void InterTerm_igt_so_int(const int i,const int j,const int k,doublecomplex result[static restrict 6])
/* Interaction term between two dipoles for approximate IGT. arguments are described in .h file
 * Can't be easily made operational for arbitrary real distance due to predefined tables.
 *
 * There is still some space for speed optimization here (e.g. move mu,nu-independent operations out of the cycles over
 * components).
 */
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,rn,invr3,kr,kr2; // |R|, |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double q2[3];
	double kd2,q4,temp,invrn2,invrn4;
	doublecomplex br,Gm0;
	int ind0,ind1,ind2,ind2m,ind3,ind4,indmunu,comp,mu,nu,mu1,nu1;
	int sigV[3],ic,sig,ivec[3],ord[3],invord[3];
	double t3q,t4q,t5tr,t6tr;
	
	vCopyIntReal(i,j,k,qvec);
	InterParams(qvec,qmunu,&rr,&rn,&invr3,&kr,&kr2,true);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	kd2=kd*kd;
	if (kr*rn < G_BOUND_CLOSE && TestTableSize(rn)) {
		//====== G close for IGT =============
		ivec[0]=i;
		ivec[1]=j;
		ivec[2]=k;
		// transformation of negative coordinates
		for (ic=0;ic<3;ic++) {
			if (ivec[ic]<0) {
				sigV[ic]=-1;
				qvec[ic]*=-1;
				ivec[ic]*=-1;
			}
			else sigV[ic]=1;
		}
		// transformation to case i>=j>=k>=0
		// building of ord; ord[x] is x-th largest coordinate (0-th - the largest)
		if (ivec[0]>=ivec[1]) { // i>=j
			if (ivec[0]>=ivec[2]) { // i>=k
				ord[0]=0;
				if (ivec[1]>=ivec[2]) { // j>=k
					ord[1]=1;
					ord[2]=2;
				}
				else {
					ord[1]=2;
					ord[2]=1;
				}
			}
			else {
				ord[0]=2;
				ord[1]=0;
				ord[2]=1;
			}
		}
		else {
			if (ivec[0]>=ivec[2]) { // i>=k
				ord[0]=1;
				ord[1]=0;
				ord[2]=2;
			}
			else {
				ord[2]=0;
				if (ivec[1]>=ivec[2]) { // j>=k
					ord[0]=1;
					ord[1]=2;
				}
				else {
					ord[0]=2;
					ord[1]=1;
				}
			}
		}
		// change parameters according to coordinate transforms
		Permutate(qvec,ord);
		Permutate_i(ivec,ord);
		// compute inverse permutation
		memcpy(invord,ord,3*sizeof(int));
		Permutate_i(invord,ord);
		if (invord[0]==0 && invord[1]==1 && invord[2]==2) memcpy(invord,ord,3*sizeof(int));
		// set some indices
		ind0=tab_index[ivec[0]][ivec[1]]+ivec[2];
		ind1=3*ind0;
		ind2m=6*ind0;
		// cycle over tensor components
		for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
			sig=sigV[mu]*sigV[nu]; // sign of some terms below
			/* indexes for tables of different dimensions based on transformed indices mu and nu '...munu' variables
			 * are invariant to permutations because both constituent vectors and indices are permutated. So this
			 * variables can be used, as precomputed above.
			 */
			mu1=invord[mu];
			nu1=invord[nu];
			/* indmunu is a number of component[mu,nu] in a symmetric matrix, but counted differently than comp.
			 * This is {{0,1,3},{1,2,4},{3,4,5}}
			 */
			indmunu=mu1+nu1;
			if (mu1==2 || nu1==2) indmunu++;
			ind2=ind2m+indmunu;
			ind3=3*ind2;
			ind4=6*ind2;
			// computing several quantities with table integrals
			t3q=DotProd(qvec,tab3+ind1);
			t4q=DotProd(qvec,tab4+ind3);
			t5tr=TrSym(tab5+ind2m);
			t6tr=TrSym(tab6+ind4);
			//====== computing Gc0 =====
			/* br = delta[mu,nu]*(-I7-I9/2-kr*(i+kr)/24+2*t3q+t5tr)
			 *    - (-3I8[mu,nu]-3I10[mu,nu]/2-qmunu*kr*(3i+kr)/24+2*t4q+t6tr)
			 */
			br = sig*(3*(tab10[ind2]/2+tab8[ind2])-2*t4q-t6tr) +(kr/24)*qmunu[comp]*(kr+I*3);
			if (dmunu[comp]) br += 2*t3q + t5tr - (kr/24)*(kr+I) - tab9[ind0]/2 - tab7[ind0];
			br*=kd2;
			// br+=I1*delta[mu,nu]*(-1+ikr+kr^2)-sig*I2[mu,nu]*(-3+3ikr+kr^2)
			br+=sig*tab2[ind2]*(3-I*3*kr-kr2);
			if (dmunu[comp]) br += tab1[ind0]*(-1+I*kr+kr2);
			// Gc0=expval*br
			result[comp]=expval*br;
		}
	}
	else {
		//====== Gfar (and part of Gmedian) for IGT =======
		// Gf0 = Gp*(1-kd^2/24)
		for (comp=0;comp<NDCOMP;comp++) result[comp]*=1-kd2/24;
		if (kr < G_BOUND_MEDIAN) {
			//===== G median for IGT ========
			vMult(qvec,qvec,q2);
			q4=DotProd(q2,q2);
			invrn2=1/(rn*rn);
			invrn4=invrn2*invrn2;
			for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
				// Gm0=expval*br*temp; temp is defined below
				temp=qmunu[comp]*(33*q4-7-12*(q2[mu]+q2[nu]));
				if (mu == nu) temp+=(1-3*q4+4*q2[mu]);
				temp*=7*invrn4/64;
				Gm0=expval*(-1+I*kr)*temp;
				// result = Gf + Gm0
				result[comp]+=Gm0;
			}
		}
	}
	PRINT_GVAL;
}

NO_REAL_WRAPPER(InterTerm_igt_so)

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

static inline void InterTerm_nloc_both(double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid,
	const bool averageH)
/* Interaction term between two dipoles using the non-local interaction;
 * qvec is  a distance given by either integer-valued vector (in units of d) or arbitrary-valued in real units (um),
 * controlled by unitsGrid (true of false respectively), result is for produced output
 * averageH specifies if the h function should be averaged over the dipole (cube) volume
 *
 * !!! Currently only static version is implemented
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
	InterParams(qvec,qmunu,&rr,&rn,&invr3,&kr,&kr2,unitsGrid);
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

void InterTerm_so_int(const int i,const int j,const int k,doublecomplex result[static restrict 6])
/* Interaction term between two dipoles with second-order corrections. arguments are described in .h file
 * Can't be easily made operational for arbitrary real distance due to predefined tables.
 *
 * There is still some space for speed optimization here (e.g. move mu,nu-independent operations out of the cycles over
 * components). But now extra time is equivalent to 2-3 main iterations. So first priority is to make something useful
 * out of SO.
 */
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,rn,invr3,kr,kr2; // |R|, |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double q2[3],qavec[3],av[3];
	double kr3,kd2,q4;
	double temp,qa,qamunu[6],invrn,invrn2,invrn3,invrn4;
	doublecomplex br,br1,m,m2,Gf1,Gm0,Gm1,Gc1,Gc2;
	int ind0,ind1,ind2,ind2m,ind3,ind4,indmunu,comp,mu,nu,mu1,nu1;
	int sigV[3],ic,sig,ivec[3],ord[3],invord[3];
	double t3q,t3a,t4q,t4a,t5tr,t5aa,t6tr,t6aa;
	const bool inter_avg=true; // temporary fixed option for SO formulation

	// next line should never happen
	if (anisotropy) LogError(ONE_POS,"Incompatibility error in InterTerm_so");

	vCopyIntReal(i,j,k,qvec);
	InterParams(qvec,qmunu,&rr,&rn,&invr3,&kr,&kr2,true);
	InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	kd2=kd*kd;
	kr3=kr2*kr;
	// only one refractive index can be used for FFT-compatible algorithm
	m=ref_index[0];
	m2=m*m;
	if (!inter_avg) {
		qa=DotProd(qvec,prop);
		for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
			// qamunu=qvec[mu]*prop[nu] + qvec[nu]*prop[mu]
			qamunu[comp]=qvec[mu]*prop[nu];
			if (dmunu[comp]) qamunu[comp]*=2;
			else qamunu[comp]+=qvec[nu]*prop[mu];
		}
	}
	if (kr*rn < G_BOUND_CLOSE && TestTableSize(rn)) {
		//====== G close =============
		// av is copy of propagation vector
		if (!inter_avg) vCopy(prop,av);
		ivec[0]=i;
		ivec[1]=j;
		ivec[2]=k;
		// transformation of negative coordinates
		for (ic=0;ic<3;ic++) {
			if (ivec[ic]<0) {
				sigV[ic]=-1;
				av[ic]*=-1;
				qvec[ic]*=-1;
				ivec[ic]*=-1;
			}
			else sigV[ic]=1;
		}
		// transformation to case i>=j>=k>=0
		// building of ord; ord[x] is x-th largest coordinate (0-th - the largest)
		if (ivec[0]>=ivec[1]) { // i>=j
			if (ivec[0]>=ivec[2]) { // i>=k
				ord[0]=0;
				if (ivec[1]>=ivec[2]) { // j>=k
					ord[1]=1;
					ord[2]=2;
				}
				else {
					ord[1]=2;
					ord[2]=1;
				}
			}
			else {
				ord[0]=2;
				ord[1]=0;
				ord[2]=1;
			}
		}
		else {
			if (ivec[0]>=ivec[2]) { // i>=k
				ord[0]=1;
				ord[1]=0;
				ord[2]=2;
			}
			else {
				ord[2]=0;
				if (ivec[1]>=ivec[2]) { // j>=k
					ord[0]=1;
					ord[1]=2;
				}
				else {
					ord[0]=2;
					ord[1]=1;
				}
			}
		}
		// change parameters according to coordinate transforms
		Permutate(qvec,ord);
		if (!inter_avg) Permutate(av,ord);
		Permutate_i(ivec,ord);
		// compute inverse permutation
		memcpy(invord,ord,3*sizeof(int));
		Permutate_i(invord,ord);
		if (invord[0]==0 && invord[1]==1 && invord[2]==2) memcpy(invord,ord,3*sizeof(int));
		// set some indices
		ind0=tab_index[ivec[0]][ivec[1]]+ivec[2];
		ind1=3*ind0;
		ind2m=6*ind0;			// cycle over tensor components
		for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
			sig=sigV[mu]*sigV[nu]; // sign of some terms below
			/* indexes for tables of different dimensions based on transformed indices mu and nu '...munu' variables
			 * are invariant to permutations because both constituent vectors and indices are permutated. So this
			 * variables can be used, as precomputed above.
			 */
			mu1=invord[mu];
			nu1=invord[nu];
			/* indmunu is a number of component[mu,nu] in a symmetric matrix, but counted differently than comp.
			 * This is {{0,1,3},{1,2,4},{3,4,5}}
			 */
			indmunu=mu1+nu1;
			if (mu1==2 || nu1==2) indmunu++;
			ind2=ind2m+indmunu;
			ind3=3*ind2;
			ind4=6*ind2;
			// computing several quantities with table integrals
			t3q=DotProd(qvec,tab3+ind1);
			t4q=DotProd(qvec,tab4+ind3);
			t5tr=TrSym(tab5+ind2m);
			t6tr=TrSym(tab6+ind4);
			if (inter_avg) {
				// <a[mu]*a[nu]>=1/3*delta[mu,nu]
				t5aa=ONE_THIRD*t5tr;
				t6aa=ONE_THIRD*t6tr;
			}
			else {
				t3a=DotProd(av,tab3+ind1);
				t4a=DotProd(av,tab4+ind3);
				t5aa=QuadForm(tab5+ind2m,av);
				t6aa=QuadForm(tab6+ind4,av);
			}
			//====== computing Gc0 =====
			/* br = delta[mu,nu]*(-I7-I9/2-kr*(i+kr)/24+2*t3q+t5tr)
			 *    - (-3I8[mu,nu]-3I10[mu,nu]/2-qmunu*kr*(3i+kr)/24+2*t4q+t6tr)
			 */
			br = sig*(3*(tab10[ind2]/2+tab8[ind2])-2*t4q-t6tr) + (kr/24)*qmunu[comp]*(kr+I*3);
			if (dmunu[comp]) br += 2*t3q + t5tr - (kr/24)*(kr+I) - tab9[ind0]/2 - tab7[ind0];
			br*=kd2;
			// br+=I1*delta[mu,nu]*(-1+ikr+kr^2)-sig*I2[mu,nu]*(-3+3ikr+kr^2)
			br+=sig*tab2[ind2]*(3-I*3*kr-kr2);
			if (dmunu[comp]) br+=tab1[ind0]*(-1+I*kr+kr2);
			// Gc0=expval*br
			result[comp]=expval*br;
			//==== computing Gc1 ======
			if (!inter_avg) {
				// br=(kd*kr/12)*(qa*(delta[mu,nu]*(-2+ikr)-qmunu*(-6+ikr))-qamunu)
				br=(6-I*kr)*qmunu[comp];
				if (dmunu[comp]) br += -2 + I*kr;
				br=(kr*kd/12)*(qa*br-qamunu[comp]);
				//  br1=(d/r)*(delta[mu,nu]*t3a*(-1+ikr)-sig*t4a*(-3+3ikr))
				br1=3*sig*t4a*(1-I*kr);
				if (dmunu[comp]) br1+=t3a*(-1+I*kr);
				br1/=rn;
				// Gc1=expval*i*m*kd*(br1+br)
				Gc1=I*expval*m*kd*(br1+br);
			}
			//==== computing Gc2 ======
			// br=delta[mu,nu]*t5aa-3*sig*t6aa-(kr/12)*(delta[mu,nu]*(i+kr)-qmunu*(3i+kr))
			br=-(kr+I*3)*qmunu[comp];
			if (dmunu[comp]) br += kr + I;
			br = -br*kr/12 - 3*sig*t6aa;
			if (dmunu[comp]) br+=t5aa;
			// Gc2=expval*(kd^2/2)*m^2*br
			Gc2=expval*(kd2/2)*m2*br;
			// result = Gc0 + [ Gc1 ] + Gc2
			if (!inter_avg) Gc2+=Gc1;
			result[comp]+=Gc2;
		}
	}
	else {
		//====== Gfar (and part of Gmedian) =======
		// br=1-(1+m^2)*kd^2/24
		br = 1 - (1+m2)*kd2/24;
		// Gf0 + Gf2 = Gp*br
		for (comp=0;comp<NDCOMP;comp++) result[comp]*=br;
		//==== compute and add Gf1 ===
		if (!inter_avg) for (comp=0;comp<NDCOMP;comp++) {
			// br = {delta[mu,nu]*(3-3ikr-2kr^2+ikr^3)-qmunu*(15-15ikr-6kr^2+ikr^3)}*qa + qamunu*(3-3ikr-kr^2)
			br = (6*kr2-15 + I*(15*kr-kr3))*qmunu[comp];
			if(dmunu[comp]) br += 3 - 2*kr2 + I*(kr3-3*kr);
			br = br*qa + (3-I*3*kr-kr2)*qamunu[comp];
			// Gf1=expval*i*m*br*kd^2/12kr
			Gf1=I*expval*m*br*kd2/(12*kr);
			// result = Gf
			result[comp]+=Gf1;
		}
		if (kr < G_BOUND_MEDIAN) {
			//===== G median ========
			vMult(qvec,qvec,q2);
			q4=DotProd(q2,q2);
			invrn=1/rn;
			invrn2=invrn*invrn;
			invrn3=invrn2*invrn;
			invrn4=invrn2*invrn2;
			for (mu=0,comp=0;mu<3;mu++) for (nu=mu;nu<3;nu++,comp++) {
				// Gm0=expval*br*temp; temp is defined below
				temp=qmunu[comp]*(33*q4-7-12*(q2[mu]+q2[nu]));
				if (mu == nu) temp+=(1-3*q4+4*q2[mu]);
				temp*=7*invrn4/64;
				Gm0=expval*(-1+I*kr)*temp;
				if (!inter_avg) {
					// Gm1=expval*i*m*temp; temp is defined below
					vMult(qvec,prop,qavec);
					temp = 3*qa*(dmunu[comp]-7*qmunu[comp]) + 6*dmunu[comp]*qvec[mu]*prop[mu]
					     - 7*(dmunu[comp]-9*qmunu[comp])*DotProd(qavec,q2)
					     + 3*(prop[mu]*qvec[nu]*(1-7*q2[mu])+prop[nu]*qvec[mu]*(1-7*q2[nu]));
					temp*=kd*invrn3/48;
					Gm1=I*expval*m*temp;
					// add Gm1 to Gm0
					Gm0+=Gm1;
				}
				// result = Gf + Gm0 + [ Gm1 ]
				result[comp]+=Gm0;
			}
		}
	}
	PRINT_GVAL;
}

NO_REAL_WRAPPER(InterTerm_so)

//=====================================================================================================================
#ifndef NO_FORTRAN

static inline void InterTerm_igt(double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid)
// Interaction term between two dipoles with integration of Green's tensor. See InterTerm_poi for more details.
{
	// standard variable definitions used for functions InterParams and InterTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rr,rn,invr3,kr,kr2; // |R|, |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3
	double tmp[12];
	int comp;

	// the following looks complicated, but should be easy to optimize by compiler
	if (igt_lim==UNDEF || DotProd(qvec,qvec)<=igt_lim*igt_lim*(unitsGrid ? 1 : (gridspace*gridspace)) ) {
		if (unitsGrid) vMultScal(gridspace,qvec,qvec);
		/* passing complex vectors from Fortran to c is not necessarily portable (at least requires extra effort in
		 * the Fortran code. So we do it through double. This is not bad for performance, since double is anyway used
		 * internally for integration in this Fortran routine.
		 */
		propaespacelibreintadda_(qvec,&WaveNum,&gridspace,&igt_eps,tmp);
		for (comp=0;comp<6;comp++) result[comp] = tmp[comp] + I*tmp[comp+6];
	}
	else {
		// The following is equivalent to InterTerm_poi, except for the 1st part of initialization performed above
		InterParams(qvec,qmunu,&rr,&rn,&invr3,&kr,&kr2,unitsGrid);
		InterTerm_core(kr,kr2,invr3,qmunu,&expval,result);
	}
	PRINT_GVAL;
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
 * static inline void InterTerm_<name>(double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid)
 * which works for double input vector and accepts additional bool argument, specifying the units.
 * After the function definition you specify "WRAPPERS_INTER(InterTerm_<name>)", which automatically produces wrappers
 * compatible with declarations in interaction.h. See InterTerm_poi() for example.
 *
 * 2) Create two separate functions named InterTerm_<name>_int and InterTerm_<name>_real. If the new formulation does
 * not support arbitrary real input vector (e.g. it is based on tables), then use
 * "NO_REAL_WRAPPER(InterTerm_<name>)" instead of the real-input declaration. See InterTerm_igt_so_int() for example.
 *
 * In any case you may benefit from existing utility functions InterParams() and InterTerm_core(). Adhering to the
 * naming conventions is important to be able to use macros here and below in InitInteraction().
 */

//=====================================================================================================================

static double * ATT_MALLOC ReadTableFile(const char * restrict sh_fname,const int size_multiplier)
{
	FILE * restrict ftab;
	double * restrict tab_n;
	int size;
	char fname[MAX_FNAME];
	int i;

	size=TAB_SIZE*size_multiplier;
	memory+=size*sizeof(double);
	if (!prognosis) {
		// allocate memory for tab_n
		MALLOC_VECTOR(tab_n,double,size,ALL);
		// open file
		SnprintfErr(ALL_POS,fname,MAX_FNAME,TAB_PATH"%s",sh_fname);
		ftab=FOpenErr(fname,"r",ALL_POS);
		// scan file
		for (i=0; i<size; i++) if (fscanf(ftab,"%lf\t",&(tab_n[i]))!=1)
			LogError(ALL_POS,"Scan error in file '%s'. Probably file is too small",fname);
		if (!feof(ftab))
			LogWarning(EC_WARN,ONE_POS,"File '%s' is longer than specified size (%d)",fname,size);
		// close file
		FCloseErr(ftab,fname,ALL_POS);
		return tab_n;
	}
	else return NULL;
}

//=====================================================================================================================

static void ReadTables(void)
{
	int i, j, ymax, Rm2, Rm2x;

	TIME_TYPE tstart=GET_TIME();
	tab1=ReadTableFile(TAB_FNAME(1),1);
	tab2=ReadTableFile(TAB_FNAME(2),6);
	tab3=ReadTableFile(TAB_FNAME(3),3);
	tab4=ReadTableFile(TAB_FNAME(4),18);
	tab5=ReadTableFile(TAB_FNAME(5),6);
	tab6=ReadTableFile(TAB_FNAME(6),36);
	tab7=ReadTableFile(TAB_FNAME(7),1);
	tab8=ReadTableFile(TAB_FNAME(8),6);
	tab9=ReadTableFile(TAB_FNAME(9),1);
	tab10=ReadTableFile(TAB_FNAME(10),6);
	Timing_FileIO += GET_TIME() - tstart;

	if (!prognosis) {
		// allocate memory for tab_index
		MALLOC_IMATRIX(tab_index,1,TAB_RMAX,0,TAB_RMAX,ALL);
		// fill tab_index
		Rm2=TAB_RMAX*TAB_RMAX;
		tab_index[1][0] = 0;
		for (i=1; i<=TAB_RMAX; i++) {
			Rm2x=Rm2-i*i;
			ymax = MIN(i,(int)floor(sqrt(Rm2x)));
			for (j=0; j<ymax; j++) tab_index[i][j+1] = tab_index[i][j] + MIN(j,(int)floor(sqrt(Rm2x-j*j)))+1;
			if (i<TAB_RMAX) tab_index[i+1][0] = tab_index[i][ymax] + MIN(ymax,(int)floor(sqrt(Rm2x-ymax*ymax)))+1;
		}
	}
	// PRINTZ("P[5,3]=%d (should be 41)\n",tab_index[5][3]);
}

//=====================================================================================================================

static void FreeTables(void)
{
	Free_iMatrix(tab_index,1,TAB_RMAX,0);
	Free_general(tab1);
	Free_general(tab2);
	Free_general(tab3);
	Free_general(tab4);
	Free_general(tab5);
	Free_general(tab6);
	Free_general(tab7);
	Free_general(tab8);
	Free_general(tab9);
	Free_general(tab10);
}

//=====================================================================================================================

static inline void vCopyIntRealShift(const int i,const int j,const int k,double qvec[static 3])
// initialize real vector with integer values
{
	qvec[0]=i;
	qvec[1]=j;
	qvec[2]=k+ZsumShift;
}

//=====================================================================================================================

static inline void ReflParams(const double qvec[static 3],double qmunu[static 6],double *invr3,double *kr,double *kr2,
	const bool unitsGrid)
// some common variables needed by the reflection functions
{
	double rr,qv[3];

	rr=vNorm(qvec);
	vMultScal(1/rr,qvec,qv); // qv is normalized qvec
	OuterSym(qv,qmunu);
	if (unitsGrid) rr*=gridspace;
	*invr3=1/(rr*rr*rr);
	*kr=WaveNum*rr;
	*kr2=(*kr)*(*kr);
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

static inline void ReflTerm_img(const double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid)
/* Reflection term using the image-dipole approximation;
 * qvec is a distance given by either integer-valued vector (in units of d) or arbitrary-valued in real units (um),
 * controlled by unitsGrid (true of false respectively). This distance is considered between probe and image of source,
 * so its z-component is the sum of heights of source and probe points above the surface.
 * result is for produced output
 */
{
	// standard variable definitions used for functions ReflParams and ReflTerm_core
	double qmunu[6]; // normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double invr3,kr,kr2; // |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	ReflParams(qvec,qmunu,&invr3,&kr,&kr2,unitsGrid);
	ReflTerm_core(kr,kr2,invr3,qmunu,&expval,result);
}

WRAPPERS_REFL(ReflTerm_img)

//=====================================================================================================================

static inline void SingleSomIntegral(double rho,const double z,doublecomplex vals[static 4])
/* computes a single Sommerfeld integral (4-element array); arguments are in real units (um)
 * currently it is a wrapper around evlua, but in the future it may be replaced by use of interpolation table
 */
{
	// TODO: these scales can be removed by changes in som_init to use proper wavenumber instead of 2*pi
	const double scale=WaveNum/TWO_PI;
	const double isc=pow(scale,3); // this is subject to under/overflow

	if (rho==0) rho=z*0.00000001; // a hack to overcome the poor precision of somnec for rho=0;
	evlua(z*scale,rho*scale,vals,vals+1,vals+2,vals+3);
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
		if (IFROOT) printf("Calculating table of Sommerfeld integrals\n");
		ind=0;
		for (k=0;k<local_Nz_Rm;k++) {
			z=(k+ZsumShift)*gridspace;
			for (j=0;j<boxY;j++) {
				if (XlessY) for (i=0;i<=j && i<boxX;i++,ind++) SingleSomIntegral(hypot(i,j)*gridspace,z,somTable+4*ind);
				else for (i=j;i<boxX;i++,ind++) SingleSomIntegral(hypot(i,j)*gridspace,z,somTable+4*ind);
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
	// standard variable definitions used for functions ReflParams and ReflTerm_core
	double qvec[3],qmunu[6]; // distance vector (in units of d) and normalized outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double invr3,kr,kr2; // |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3
	vCopyIntRealShift(i,j,k,qvec);
	ReflParams(qvec,qmunu,&invr3,&kr,&kr2,true);
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
	double invr3,kr,kr2; // |R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3
	ReflParams(qvec,qmunu,&invr3,&kr,&kr2,false);
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
 * static inline void ReflTerm_<name>(double qvec[static 3],doublecomplex result[static 6],const bool unitsGrid)
 * which works for double input vector and accepts additional bool argument, specifying the units.
 * After the function definition you specify "WRAPPERS_REFL(ReflTerm_<name>)", which automatically produces
 * declarations compatible to that in interaction.h. See ReflTerm_img() for example.
 *
 * 2) Create two separate functions named ReflTerm_<name>_int and ReflTerm_<name>_real. If the new formulation does
 * not support arbitrary real input vector (e.g. it is based on tables), then use
 * "NO_REAL_WRAPPER(ReflTerm_<name>)" instead. See ReflTerm_som_int/real() for example.
 *
 * In any case you may benefit from existing utility functions ReflParams() and ReflTerm_core(). Adhering to the
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
		case G_IGT_SO:
			if (InteractionRealArgs) PrintError("'-int igt_so' does not support calculation of interaction tensor for "
				"arbitrary real arguments");
			SET_FUNC_POINTERS(InterTerm,igt_so);
			break;
		case G_NLOC: SET_FUNC_POINTERS(InterTerm,nloc); break;
		case G_NLOC_AV: SET_FUNC_POINTERS(InterTerm,nloc_av); break;
		case G_SO:
			if (InteractionRealArgs) PrintError("'-int so' does not support calculation of interaction tensor for "
				"arbitrary real arguments");
			SET_FUNC_POINTERS(InterTerm,so);
			break;
#ifndef NO_FORTRAN
		case G_IGT: SET_FUNC_POINTERS(InterTerm,igt); break;
#endif
		/* TO ADD NEW INTERACTION FORMULATION
		 * Add here the assignment of function pointers for the new formulation. It is recommended to use special macro,
		 * assuming that you conformed to naming conventions for functions themselves. If new formulation requires
		 * initialization, add it here, preferably by a separate function. Initialization should honor the 'prognosis'
		 * flag. Additional memory should be counted always, but allocated only when not prognosis. See ReadTables()
		 * for example (called below). If the new formulation does not support arbitrary real input vector (e.g. it is
		 * based on tables), then test for InteractionRealArgs and add an exception (see G_IGT_SO for example).
		 */
		default: LogError(ONE_POS, "Invalid interaction term calculation method: %d",(int)IntRelation);
			// no break
	}
	// read tables if needed
	if (IntRelation == G_SO || IntRelation == G_IGT_SO) ReadTables();

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
			 * is based on tables), then test for InteractionRealArgs and add an exception (for example, see G_IGT_SO in
			 * IntRelation above).
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
	if (IntRelation == G_SO || IntRelation == G_IGT_SO) FreeTables();
	if (surface && ReflRelation==GR_SOM) {
		Free_general(somIndex);
		Free_cVector(somTable);
	}
	/* TO ADD NEW INTERACTION FORMULATION
	 * TO ADD NEW REFLECTION FORMULATION
	 * If you allocate any memory (for tables, etc.), free it here
	 */
}
