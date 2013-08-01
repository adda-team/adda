/* FILE : interaction.c
 * Descr: the functions used to calculate the interaction term
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
#include "interaction.h" // corresponding headers
// project headers
#include "cmplx.h"
#include "io.h"
#include "memory.h"
#include "vars.h"

// SEMI-GLOBAL VARIABLES

// defined and initialized in make_particle.c
extern double gridspace;
// defined and initialized in param.c
extern double igt_lim,igt_eps;

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

// this is used for debugging, should be empty define, when not required
#define PRINT_GVAL /*printf("%d,%d,%d: %g%+gi, %g%+gi, %g%+gi,\n%g%+gi, %g%+gi, %g%+gi\n",i,j,k,creal(result[0]),\
	cimag(result[0]),creal(result[1]),cimag(result[1]),creal(result[2]),cimag(result[2]),creal(result[3]),\
	cimag(result[3]),creal(result[4]),cimag(result[4]),creal(result[5]),cimag(result[5]));*/

//=====================================================================================================================
/* The following two functions are used to do common calculation parts. For simplicity it is best to call them with the
 * same set of parameters - in this respect they can be replaced by macros. But inline functions are probably easier to
 * maintain.
 */

static inline void CalcInterParams1(const int i,const int j,const int k,double qvec[static restrict 3],double *rn)
// some common variables needed by the interaction functions - needed always
{
	qvec[0]=i; // qvec is normalized below
	qvec[1]=j;
	qvec[2]=k;
	*rn=sqrt(DotProd(qvec,qvec)); // normalized r
}

//=====================================================================================================================

static inline void CalcInterParams2(double qvec[static restrict 3],double qmunu[static restrict 6],const double rn,
	double *invrn,double *invr3,double *kr,double *kr2)
// some common variables needed by the interaction functions - needed for all except IGT
{
	*invrn = 1.0/rn;
	vMultScal(*invrn,qvec,qvec); // finalize qvec
	double rr=rn*gridspace;
	*invr3=1/(rr*rr*rr);
	*kr=WaveNum*rr;
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

static void CalcInterTerm_core(const double kr,const double kr2,const double invr3,
	const double qmunu[static restrict 6],doublecomplex * restrict expval,doublecomplex result[static restrict 6])
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

static void CalcInterTerm_core(const double kr,const double kr2,const double invr3,
	const double qmunu[static restrict 6],doublecomplex * restrict expval,doublecomplex result[static restrict 6])
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

void CalcInterTerm_poi(const int i,const int j,const int k,doublecomplex result[static restrict 6])
/* calculates interaction term between two dipoles; given integer distance vector {i,j,k} (in units of d). All six
 * components of the symmetric matrix are computed at once. The elements in result are: [G11, G12, G13, G22, G23, G33]
 */
{
	// standard variable definitions used for functions CalcInterParams1,2 and CalcInterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rn,invrn,invr3,kr,kr2; // |R/d|, 1/|R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	CalcInterParams1(i,j,k,qvec,&rn);
	CalcInterParams2(qvec,qmunu,rn,&invrn,&invr3,&kr,&kr2);
	CalcInterTerm_core(kr,kr2,invr3,qmunu,&expval,result);
	PRINT_GVAL;
}

//=====================================================================================================================

void CalcInterTerm_fcd(const int i,const int j,const int k,doublecomplex result[static restrict 6])
/* Interaction term between two dipoles for FCD. See CalcInterTerm_poi for more details.
 *
 * FCD is based on Gay-Balmaz P., Martin O.J.F. "A library for computing the filtered and non-filtered 3D Green's tensor
 * associated with infinite homogeneous space and surfaces", Comp. Phys. Comm. 144:111-120 (2002), and
 * Piller N.B. "Increasing the performance of the coupled-dipole approximation: A spectral approach",
 * IEEE Trans.Ant.Propag. 46(8): 1126-1137. Here it differs by a factor of 4*pi*k^2.
 *
 * speed of FCD can be improved by using faster version of sici routine, using predefined tables, etc (e.g. as is
 * done in GSL library). But currently extra time for this computation is already smaller than one main iteration.
 */
{
	// standard variable definitions used for functions CalcInterParams1,2 and CalcInterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rn,invrn,invr3,kr,kr2; // |R/d|, 1/|R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double temp,kfr,ci,si,ci1,si1,ci2,si2,brd,g0,g2;
	int comp;
	doublecomplex eikfr; // exp(i*k_F*R)

	CalcInterParams1(i,j,k,qvec,&rn);
	CalcInterParams2(qvec,qmunu,rn,&invrn,&invr3,&kr,&kr2);
	CalcInterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

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

//=====================================================================================================================

void CalcInterTerm_fcd_st(const int i,const int j,const int k,doublecomplex result[static restrict 6])
// Interaction term between two dipoles for static FCD (in the limit of k->inf). See CalcInterTerm_fcd for more details.
{
	// standard variable definitions used for functions CalcInterParams1,2 and CalcInterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rn,invrn,invr3,kr,kr2; // |R/d|, 1/|R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double kfr,ci,si,brd;
	int comp;
	doublecomplex eikfr;

	CalcInterParams1(i,j,k,qvec,&rn);
	CalcInterParams2(qvec,qmunu,rn,&invrn,&invr3,&kr,&kr2);
	CalcInterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

	kfr=PI*rn; // k_F*r, for FCD
	eikfr=accImExp(kfr);
	// result = Gp*[3*Si(k_F*r)+k_F*r*cos(k_F*r)-4*sin(k_F*r)]*2/(3*pi)
	cisi(kfr,&ci,&si);
	brd=TWO_OVER_PI*ONE_THIRD*(3*si+kfr*creal(eikfr)-4*cimag(eikfr));
	for (comp=0;comp<NDCOMP;comp++) result[comp]*=brd;
	PRINT_GVAL;
}

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

void CalcInterTerm_igt_so(const int i,const int j,const int k,doublecomplex result[static restrict 6])
/* Interaction term between two dipoles for approximate IGT. See CalcInterTerm_poi for more details.
 *
 * There is still some space for speed optimization here (e.g. move mu,nu-independent operations out of the cycles over
 * components).
 */
{
	// standard variable definitions used for functions CalcInterParams1,2 and CalcInterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rn,invrn,invr3,kr,kr2; // |R/d|, 1/|R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double q2[3];
	double kd2,q4,temp,invrn2,invrn4;
	doublecomplex br,Gm0;
	int ind0,ind1,ind2,ind2m,ind3,ind4,indmunu,comp,mu,nu,mu1,nu1;
	int sigV[3],ic,sig,ivec[3],ord[3],invord[3];
	double t3q,t4q,t5tr,t6tr;
	
	CalcInterParams1(i,j,k,qvec,&rn);
	CalcInterParams2(qvec,qmunu,rn,&invrn,&invr3,&kr,&kr2);
	CalcInterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

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
			invrn2=invrn*invrn;
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

//=====================================================================================================================

void CalcInterTerm_so(const int i,const int j,const int k,doublecomplex result[static restrict 6])
/* Interaction term between two dipoles with second-order corrections. See CalcInterTerm_poi for more details.
 *
 * There is still some space for speed optimization here (e.g. move mu,nu-independent operations out of the cycles over
 * components). But now extra time is equivalent to 2-3 main iterations. So first priority is to make something useful
 * out of SO.
 */
{
	// standard variable definitions used for functions CalcInterParams1,2 and CalcInterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rn,invrn,invr3,kr,kr2; // |R/d|, 1/|R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3

	double q2[3],qavec[3],av[3];
	double kr3,kd2,q4;
	double temp,qa,qamunu[6],invrn2,invrn3,invrn4;
	doublecomplex br,br1,m,m2,Gf1,Gm0,Gm1,Gc1,Gc2;
	int ind0,ind1,ind2,ind2m,ind3,ind4,indmunu,comp,mu,nu,mu1,nu1;
	int sigV[3],ic,sig,ivec[3],ord[3],invord[3];
	double t3q,t3a,t4q,t4a,t5tr,t5aa,t6tr,t6aa;
	const bool inter_avg=true; // temporary fixed option for SO formulation

	// next line should never happen
	if (anisotropy) LogError(ONE_POS,"Incompatibility error in CalcInterTerm");

	CalcInterParams1(i,j,k,qvec,&rn);
	CalcInterParams2(qvec,qmunu,rn,&invrn,&invr3,&kr,&kr2);
	CalcInterTerm_core(kr,kr2,invr3,qmunu,&expval,result);

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
		if (!inter_avg) memcpy(av,prop,3*sizeof(double));
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

//=====================================================================================================================
#ifndef NO_FORTRAN

void CalcInterTerm_igt(const int i,const int j,const int k,doublecomplex result[static restrict 6])
/* Interaction term between two dipoles with integration of Green's tensor. See CalcInterTerm_poi for more details.
 */
{
	// standard variable definitions used for functions CalcInterParams1,2 and CalcInterTerm_core
	double qvec[3],qmunu[6]; // unit directional vector {qx,qy,qz} and its outer-product {qxx,qxy,qxz,qyy,qyz,qzz}
	double rn,invrn,invr3,kr,kr2; // |R/d|, 1/|R/d|, |R|^-3, kR, (kR)^2
	doublecomplex expval; // exp(ikR)/|R|^3
	double rtemp[3],tmp[12];
	int comp;

	CalcInterParams1(i,j,k,qvec,&rn);
	if (igt_lim==UNDEF || rn<=igt_lim) {
		vMultScal(gridspace,qvec,rtemp);
		/* passing complex vectors from Fortran to c is not necessarily portable (at least requires extra effort in
		 * the Fortran code. So we do it through double. This is not bad for performance, since double is anyway used
		 * internally for integration in this Fortran routine.
		 */
		propaespacelibreintadda_(rtemp,&WaveNum,&gridspace,&igt_eps,tmp);
		for (comp=0;comp<6;comp++) result[comp] = tmp[comp] + I*tmp[comp+6];
	}
	else {
		// The following is equivalent to CalcInterTerm_poi, except for the 1st part of initialization performed above
		CalcInterParams2(qvec,qmunu,rn,&invrn,&invr3,&kr,&kr2);
		CalcInterTerm_core(kr,kr2,invr3,qmunu,&expval,result);
	}
	PRINT_GVAL;
}

#endif

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

void InitInteraction(void)
// Initialize the interaction calculations
{
	// set CalcInterTerm to point at the right function
	switch (IntRelation) {
		case G_POINT_DIP: CalcInterTerm = &CalcInterTerm_poi; break;
		case G_FCD: CalcInterTerm = &CalcInterTerm_fcd; break;
		case G_FCD_ST: CalcInterTerm = &CalcInterTerm_fcd_st; break;
		case G_IGT_SO: CalcInterTerm = &CalcInterTerm_igt_so; break;
		case G_SO: CalcInterTerm = &CalcInterTerm_so; break;
#ifndef NO_FORTRAN
		case G_IGT: CalcInterTerm = &CalcInterTerm_igt; break;
#endif
		default: LogError(ONE_POS, "Invalid interaction term calculation method: %d",(int)IntRelation);
			// no break
	}
	// read tables if needed
	if (IntRelation == G_SO || IntRelation == G_IGT_SO) ReadTables();

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
}

//=====================================================================================================================

void FreeInteraction(void)
// Free buffers used for interaction calculation
{
	if (IntRelation == G_SO || IntRelation == G_IGT_SO) FreeTables();
}
