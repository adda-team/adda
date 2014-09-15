/* File: calculator.c
 * $Date::                            $
 * Descr: all the initialization is done here before actually calculating internal fields;
 *        includes calculation of couple constants
 *
 * Copyright (C) 2006-2010,2013-2014 ADDA contributors
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
#include "cmplx.h"
#include "comm.h"
#include "crosssec.h"
#include "debug.h"
#include "fft.h"
#include "interaction.h"
#include "io.h"
#include "memory.h"
#include "oclcore.h"
#include "Romberg.h"
#include "timing.h"
#include "vars.h"
// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>

// SEMI-GLOBAL VARIABLES

// defined and initialized in crosssec.c
extern const Parms_1D parms[2],parms_alpha;
extern const angle_set beta_int,gamma_int,theta_int,phi_int;
// defined and initialized in param.c
extern const int avg_inc_pol;
extern const double polNlocRp;
extern const char *alldir_parms,*scat_grid_parms;
// defined and initialized in timing.c
extern TIME_TYPE Timing_Init,Timing_Init_Int;
#ifdef OPENCL
extern TIME_TYPE Timing_OCL_Init;
#endif
extern size_t TotalEval;

#ifdef ACCIMEXP
extern doublecomplex * restrict imexptable;
#endif

// used in CalculateE.c
double * restrict muel_phi; // used to store values of Mueller matrix for different phi (to integrate)
double * restrict muel_phi_buf; // additional for integrating with different multipliers
	// scattered E (for scattering in the default scattering plane) for two incident polarizations
doublecomplex * restrict EplaneX, * restrict EplaneY;
doublecomplex * restrict EyzplX, * restrict EyzplY; // same for scattering in yz-plane
double dtheta_deg,dtheta_rad; // delta theta in degrees and radians
doublecomplex * restrict ampl_alphaX,* restrict ampl_alphaY; // amplitude matrix for different values of alpha
double * restrict muel_alpha; // mueller matrix for different values of alpha

// used in crosssec.c
doublecomplex * restrict E_ad; // complex field E, calculated for alldir
double * restrict E2_alldir; // square of E (scaled with msub, so ~ Poynting vector or dC/dOmega), calculated for alldir
doublecomplex cc[MAX_NMAT][3]; // couple constants
#ifndef SPARSE
doublecomplex * restrict expsX,* restrict expsY,* restrict expsZ; // arrays of exponents along 3 axes (for calc_field)
#endif
// used in iterative.c
doublecomplex *rvec;                 // current residual
doublecomplex * restrict Avecbuffer; // used to hold the result of matrix-vector products
// auxiliary vectors, used in some iterative solvers (with more meaningful names)
doublecomplex * restrict vec1,* restrict vec2,* restrict vec3,* restrict vec4;
// used in matvec.c
#ifdef SPARSE
doublecomplex * restrict arg_full; // vector to hold argvec for all dipoles
#endif

// LOCAL VARIABLES

static size_t block_theta; // size of one block of mueller matrix - 16*nTheta
static int finish_avg; // whether to stop orientation averaging; defined as int to simplify MPI casting
static double * restrict out; // used to collect both mueller matrix and integral scattering quantities when orient_avg

// EXTERNAL FUNCTIONS

// CalculateE.c
int CalculateE(enum incpol which,enum Eftype type);
bool TestExtendThetaRange(void);
void MuellerMatrix(void);
void SaveMuellerAndCS(double * restrict in);

//======================================================================================================================

static inline doublecomplex polCM(const doublecomplex m)
// computes CM polarizability from given refractive index m - (3V/4pi)*(m^2-1)/(m^2+2)
{
	doublecomplex m2=m*m;
	return (3/FOUR_PI)*dipvol*(m2-1)/(m2+2);
}

//======================================================================================================================

static inline doublecomplex polM(const doublecomplex M,const doublecomplex m)
// computes polarizability, given the M term and refractive index m
{
	doublecomplex cm=polCM(m);
	return cm/(1-(cm/dipvol)*M);
}

//======================================================================================================================

static inline doublecomplex polMplusRR(const doublecomplex M,const doublecomplex m)
// computes polarizability, given the M term (without RR part - ) and refractive index m
{
	return polM(M+I*2*kd*kd*kd/3,m);
}

//======================================================================================================================

static inline doublecomplex pol3coef(const double b1,const double b2,const double b3,const double S,
	const doublecomplex m)
/* computes polarizability, using the common formula, based on 3 coefs, S-term and refractive index m
 * M0=(b1+(b2+b3*S)*m^2)*kd^2.   Important feature is that M0 depends on m.
 */
{
	return polMplusRR((b1+(b2+b3*S)*m*m)*kd*kd,m);
}

//======================================================================================================================

static inline double ellTheta(const double a)
/* computes the following expression: a^3(theta_3(0,exp(-pi*a^2))^3-1), where
 * theta_3 is elliptic theta function of 3rd kind, in particular
 * theta_3(0,exp(-pi*a^2)) = 1 + Sum[exp(-pi*a^2*k,{k,1,inf}]
 * it obeys the following relation theta_3(0,exp(-pi*a^2)) = (1/a)*theta_3(0,exp(-pi/a^2)),
 * see e.g. J.D. Fenton and R.S. Gardiner-Garden, "Rapidly-convergent methods for evaluating elliptic integrals and
 * theta and elliptic functions," J. Austral. Math. Soc. B 24, 47-58 (1982).
 * or http://en.wikipedia.org/wiki/Theta_function#Jacobi_identities
 * so the sum need to be taken only for a>=1, then three terms are sufficient to obtain 10^-22 accuracy
 */
{
	double q,q2,q3,t,a2,a3,res;
	a2=a*a;
	a3=a*a2;
	if (a>=1) {
		q=exp(-PI*a2);
		q2=q*q;
		q3=q*q2;
		t=2*q*(1+q3*(1+q2*q3)); // t = 1 - theta_3(0,exp(-pi*a^2)) = 2*(q + q^4 + q^9)
		res=a3*t*(3+t*(3+t)); // a^3*(t^3-1)
	}
	else { // a<1, employ transformation a->1/a
		q=exp(-PI/a2);
		q2=q*q;
		q3=q*q2;
		t=1+2*q*(1+q3*(1+q2*q3)); // t = theta_3(0,exp(-pi/a^2)) = 1+ 2*(q + q^4 + q^9)
		res=t*t*t - a3; // a^3*(t^3-1)
	}
	return res;
}

//======================================================================================================================

static void CoupleConstant(doublecomplex *mrel,const enum incpol which,doublecomplex res[static 3])
/* Input is relative refractive index (mrel) - either one or three components (for anisotropic). incpol is relevant only
 * for LDR without avgpol. res is three values (diagonal of polarizability tensor).
 *
 * !!! TODO: if this function will be executed many times, it can be optimized by moving time-consuming calculation of
 * certain coefficients to one-call initialization function.
 *
 * TO ADD NEW POLARIZABILITY FORMULATION
 * This function implements calculations of polarizability. It should be updated to support new formulations.
 * The corresponding case should either be added to 'asym' list (then three components of polarizability are
 * calculated from one m) or to another one, then a scalar function is used. See comments in the code for more details.
 */
{
	double ka,kd2,S;
	int i;
	bool asym; // whether polarizability is asymmetric (for isotropic m)
	const double *incPol;
	bool pol_avg=true; // temporary fixed value for SO polarizability

	asym = (PolRelation==POL_CLDR || PolRelation==POL_SO); // whether non-scalar tensor is produced for scalar m
	// !!! this should never happen
	if (asym && anisotropy) LogError(ONE_POS,"Incompatibility error in CoupleConstant");

	kd2=kd*kd;
	if (asym) for (i=0;i<3;i++) { // loop over components of polarizability (for scalar input m)
		switch (PolRelation) {
			case POL_CLDR: res[i]=pol3coef(LDR_B1,LDR_B2,LDR_B3,prop[i]*prop[i],mrel[0]); break;
			case POL_SO: res[i]=pol3coef(SO_B1,SO_B2,SO_B3,(pol_avg ? ONE_THIRD : prop[i]*prop[i]),mrel[0]); break;
			default: LogError(ONE_POS,"Incompatibility error in CoupleConstant");
				// no break
		}
	}
	else for (i=0;i<Ncomp;i++) { // loop over components of input m
		switch (PolRelation) {
			case POL_CM: res[i]=polCM(mrel[i]); break;
			case POL_DGF: res[i]=polMplusRR(DGF_B1*kd2,mrel[i]); break;
			case POL_FCD: // M0={(4/3)kd^2+(2/3pi)log[(pi-kd)/(pi+kd)]kd^3}
				res[i]=polMplusRR(2*ONE_THIRD*kd2*(2+kd*INV_PI*log((PI-kd)/(PI+kd))),mrel[i]);
				break;
			case POL_IGT_SO: res[i]=polMplusRR(SO_B1*kd2,mrel[i]); break;
			case POL_LAK: // M=(8pi/3)[(1-ika)exp(ika)-1], a - radius of volume-equivalent (to cubical dipole) sphere
				ka=LAK_C*kd;
				res[i]=polM(2*FOUR_PI_OVER_THREE*((1-I*ka)*imExp(ka)-1),mrel[i]);
				break;
			case POL_LDR:
				if (avg_inc_pol) S=0.5*(1-DotProdSquare(prop,prop));
				else {
					if (which==INCPOL_Y) incPol=incPolY;
					else incPol=incPolX; // which==INCPOL_X
					S = DotProdSquare(prop,incPol);
				}
				res[i]=pol3coef(LDR_B1,LDR_B2,LDR_B3,S,mrel[i]);
				break;
			case POL_NLOC: // !!! additionally dynamic part should be added (if needed)
				/* Here the polarizability is derived from the condition that V_d*sum(G_h(ri))=-4pi/3, where sum is
				 * taken over the whole lattice. Then M=4pi/3+V_d*Gh(0)=V_d*sum(G_h(ri),i!=0)
				 * Moreover, the regular part (in limit Rp->0) of Green's tensor automatically sums to zero, so only the
				 * irregular part need to be considered -h(r)*4pi/3, where h(r) is a normalized Gaussian
				 */
				if (polNlocRp==0) res[i]=polCM(mrel[i]);
				else res[i]=polM(FOUR_PI_OVER_THREE*ellTheta(SQRT1_2PI*gridspace/polNlocRp),mrel[i]);
				break;
			case POL_NLOC_AV:
				if (polNlocRp==0) res[i]=polCM(mrel[i]); // polMplusRR(DGF_B1*kd2,mrel[i]); // just DGF
				else {
					double x=gridspace/(2*SQRT2*polNlocRp);
					double g0,t;
					// g0 = 1 - erf(x)^3, but careful evaluation is performed to keep precision
					if (x<1) {
						t=erf(x);
						g0=1-t*t*t;
					}
					else {
						t=erfc(x);
						g0=t*(3-3*t+t*t);
					}
					// !!! dynamic part should be added here
					res[i]=polM(FOUR_PI_OVER_THREE*g0,mrel[i]);
				}
				break;
			case POL_RRC: res[i]=polMplusRR(0,mrel[i]); break;
			default: LogError(ONE_POS,"Incompatibility error in CoupleConstant");
				// no break
		}
	}
	if (asym || anisotropy) {
		if (!orient_avg && IFROOT) PrintBoth(logfile, "CoupleConstant:"CFORM3V"\n",REIM3V(res));
	}
	else {
		res[2]=res[1]=res[0];
		if (!orient_avg && IFROOT) PrintBoth(logfile,"CoupleConstant:"CFORM"\n",REIM(res[0]));
	}
}

//======================================================================================================================

static void InitCC(const enum incpol which)
// calculate cc, cc_sqrt, and chi_inv
{
	int i,j;
	doublecomplex m;

	for(i=0;i<Nmat;i++) {
		CoupleConstant(ref_index+Ncomp*i,which,cc[i]);
		for(j=0;j<3;j++) cc_sqrt[i][j]=csqrt(cc[i][j]);
		// chi_inv=1/(V*chi)=4*PI/(V(m^2-1)); for anisotropic - by components
		for (j=0;j<Ncomp;j++) {
			m=ref_index[Ncomp*i+j];
			chi_inv[i][j]=FOUR_PI/(dipvol*(m*m-1));
		}
		// copy first component of chi_inv[i] into other two, if they are not calculated explicitly
		if (!anisotropy) chi_inv[i][2]=chi_inv[i][1]=chi_inv[i][0];
	}
#ifdef OPENCL
	/* this is done here, since InitCC can be run between different runs of the iterative solver; write is blocking to
	 * ensure completion before function end
	 */
	CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufcc_sqrt,CL_TRUE,0,sizeof(cc_sqrt),cc_sqrt,0,NULL,NULL));
#endif
}

//======================================================================================================================

static void calculate_one_orientation(double * restrict res)
// performs calculation for one orientation; may do orientation averaging and put the result in res
{
	TIME_TYPE tstart;

	if (orient_avg) {
		alph_deg=0;
		InitRotation();
		if (IFROOT) PrintBoth(logfile,"\nORIENTATION STEP beta="GFORMDEF" gamma="GFORMDEF"\n",bet_deg,gam_deg);
	}

	// calculate scattered field for y - polarized incident light
	if (IFROOT) {
		printf("\nhere we go, calc Y\n\n");
		if (!orient_avg) fprintf(logfile,"\nhere we go, calc Y\n\n");
	}
	InitCC(INCPOL_Y);
	// symR implies that prop is along z (in particle RF). Then it is fine for both definitions of scattering angles
	if (symR && !scat_grid) {
		if (CalculateE(INCPOL_Y,CE_PARPER)==CHP_EXIT) return;
	}
	else { // no rotational symmetry
		/* TODO: in case of scat_grid we run twice to get the full electric field with incoming light polarized in X and
		 * Y direction. In case of rotational symmetry this is not needed but requires lots more programming so we leave
		 * this optimization to a later time.
		 */
		if(CalculateE(INCPOL_Y,CE_NORMAL)==CHP_EXIT) return;

		if (IFROOT) {
			printf("\nhere we go, calc X\n\n");
			if (!orient_avg) fprintf(logfile,"\nhere we go, calc X\n\n");
		}
		if (PolRelation==POL_LDR && !avg_inc_pol) InitCC(INCPOL_X);
		/* TO ADD NEW POLARIZABILITY FORMULATION
		 * If new formulation depends on the incident polarization (unlikely) update the test above.
		 */

		if(CalculateE(INCPOL_X,CE_NORMAL)==CHP_EXIT) return;
	}
	D("CalculateE finished");
	MuellerMatrix();
	D("MuellerMatrix finished");
	if (IFROOT && orient_avg) {
		tstart=GET_TIME();
		if (store_mueller) printf("\nError of alpha integration (Mueller) is "GFORMDEF"\n",
			Romberg1D(parms_alpha,block_theta,muel_alpha,res+2));
		memcpy(res,muel_alpha-2,2*sizeof(double));
		D("Integration over alpha completed on root");
		Timing_Integration += GET_TIME() - tstart;
	}
	TotalEval++;
}

//======================================================================================================================

static double orient_integrand(int beta_i,int gamma_i, double * restrict res)
// function that provides interface with Romberg integration
{
	BcastOrient(&beta_i,&gamma_i,&finish_avg);
	if (finish_avg) return 0;

	bet_deg=beta_int.val[beta_i];
	gam_deg=gamma_int.val[gamma_i];
	calculate_one_orientation(res);
	return 0;
}

//======================================================================================================================

static void AllocateEverything(void)
// allocates a lot of arrays and performs memory analysis
{
	double tmp;
	size_t temp_int;
	double memmax;

	// redundant initialization to remove warnings
	temp_int=0;

	/* It may be nice to initialize all pointers to NULL here, so that any pointer, which is not initialized below, will
	 * surely stay NULL (independent of a particular compiler). But even without this forgetting to allocate a necessary
	 * vector, will surely cause segmentation fault afterwards. So we do not implement these extra tests for now.
	 */
	// allocate all the memory
	tmp=sizeof(doublecomplex)*(double)local_nRows;
	if (!prognosis) { // main 5 vectors, some of them are used in the iterative solver
		MALLOC_VECTOR(xvec,complex,local_nRows,ALL);
		MALLOC_VECTOR(rvec,complex,local_nRows,ALL);
		MALLOC_VECTOR(pvec,complex,local_nRows,ALL);
		MALLOC_VECTOR(Einc,complex,local_nRows,ALL);
		MALLOC_VECTOR(Avecbuffer,complex,local_nRows,ALL);
	}
	memory+=5*tmp;
#ifdef SPARSE
	if (!prognosis) { // overflow of 3*nvoid_Ndip is tested in MakeParticle()
		MALLOC_VECTOR(arg_full,complex,3*nvoid_Ndip,ALL);
	}
	memory+=3*nvoid_Ndip*sizeof(doublecomplex);
#endif // !SPARSE
	/* additional vectors for iterative methods. Potentially, this procedure can be fully automated for any new
	 * iterative solver, based on the information contained in structure array 'params' in file iterative.c. However,
	 * this requires different order of function calls to extract this information beforehand. So currently this part
	 * should be edited manually when needed.
	 */
	switch (IterMethod) {
		case IT_BCGS2:
			if (!prognosis) {
				MALLOC_VECTOR(vec1,complex,local_nRows,ALL);
				MALLOC_VECTOR(vec2,complex,local_nRows,ALL);
				MALLOC_VECTOR(vec3,complex,local_nRows,ALL);
				MALLOC_VECTOR(vec4,complex,local_nRows,ALL);
			}
			memory+=4*tmp;
			break;
		case IT_CGNR:
		case IT_BICG_CS:
			break;
		case IT_BICGSTAB:
		case IT_QMR_CS:
			if (!prognosis) {
				MALLOC_VECTOR(vec1,complex,local_nRows,ALL);
				MALLOC_VECTOR(vec2,complex,local_nRows,ALL);
				MALLOC_VECTOR(vec3,complex,local_nRows,ALL);
			}
			memory+=3*tmp;
			break;
		case IT_CSYM:
		case IT_QMR_CS_2:
			if (!prognosis) {
				MALLOC_VECTOR(vec1,complex,local_nRows,ALL);
				MALLOC_VECTOR(vec2,complex,local_nRows,ALL);
			}
			memory+=2*tmp;
			break;
	}
	/* TO ADD NEW ITERATIVE SOLVER
	 * Add here a case corresponding to the new iterative solver. If the new iterative solver requires any extra vectors
	 * (additionally to the default ones), i.e. number vec_N in corresponding element of structure array params in
	 * iterative.c is non-zero, then allocate memory for these vectors here. Variable memory should be incremented to
	 * reflect the total allocated memory.
	 */
#ifndef SPARSE
	MALLOC_VECTOR(expsX,complex,boxX,ALL);
	MALLOC_VECTOR(expsY,complex,boxY,ALL);
	MALLOC_VECTOR(expsZ,complex,local_Nz_unif,ALL);
#endif // !SPARSE
#ifdef ACCIMEXP
	MALLOC_VECTOR(imexptable,complex,361,ALL);
	double x;
	for(unsigned int i=0; i<=360;i++) {
		x = (PI/180.0)*(double)i;
		imexptable[i] = cos(x) + I*sin(x);
	}
#endif
	if (yzplane) {
		tmp=2*(double)nTheta;
		if (!prognosis) {
			CheckOverflow(2*tmp,ONE_POS_FUNC);
			temp_int=tmp;
			MALLOC_VECTOR(EyzplX,complex,temp_int,ALL);
			MALLOC_VECTOR(EyzplY,complex,temp_int,ALL);
		}
		memory+=2*tmp*sizeof(doublecomplex);
	}
	if (scat_plane) {
		tmp=2*(double)nTheta;
		if (!prognosis) {
			CheckOverflow(2*tmp,ONE_POS_FUNC);
			temp_int=tmp;
			MALLOC_VECTOR(EplaneX,complex,temp_int,ALL);
			MALLOC_VECTOR(EplaneY,complex,temp_int,ALL);
		}
		memory+=2*tmp*sizeof(doublecomplex);
	}
	if (all_dir) {
		ReadAlldirParms(alldir_parms);
		/* calculate size of vectors; 4 - because first it is used to store per and par components of the field, and
		 * only afterwards - squares.
		 */
		tmp=((double)theta_int.N)*phi_int.N;
		if (!prognosis) {
			CheckOverflow(2*tmp,ONE_POS_FUNC);
			temp_int=tmp;
			MALLOC_VECTOR(E_ad,complex,2*temp_int,ALL);
			MALLOC_VECTOR(E2_alldir,double,temp_int,ALL);
		}
		memory+=tmp*(sizeof(double)+2*sizeof(doublecomplex));
	}
	if (scat_grid) {
		ReadScatGridParms(scat_grid_parms);
		// calculate size of vectors - holds all per-par combinations
		tmp=2*(double)angles.N;
		if (!prognosis) {
			CheckOverflow(2*tmp,ONE_POS_FUNC);
			temp_int=tmp;
			MALLOC_VECTOR(EgridX,complex,temp_int,ALL);
			MALLOC_VECTOR(EgridY,complex,temp_int,ALL);
		}
		memory+=2*tmp*sizeof(doublecomplex);
		if (phi_integr && IFROOT) {
			tmp=16*(double)angles.phi.N;
			if (!prognosis) {
				CheckOverflow(tmp,ONE_POS_FUNC);
				temp_int=tmp;
				MALLOC_VECTOR(muel_phi,double,temp_int,ONE);
				MALLOC_VECTOR(muel_phi_buf,double,temp_int,ONE);
			}
			memory+=2*tmp*sizeof(double);
		}
	}
	if (orient_avg) {
		tmp=2*((double)nTheta)*alpha_int.N;
		if (!prognosis) {
			// this covers these 2 and next 2 malloc calls
			CheckOverflow(8*tmp+2,ONE_POS_FUNC);
			if (store_mueller) {
				temp_int=tmp;
				MALLOC_VECTOR(ampl_alphaX,complex,temp_int,ONE);
				MALLOC_VECTOR(ampl_alphaY,complex,temp_int,ONE);
			}
		}
		memory += 2*tmp*sizeof(doublecomplex);
		if (IFROOT) {
			if (!prognosis) {
				MALLOC_VECTOR(muel_alpha,double,block_theta*alpha_int.N+2,ONE);
				muel_alpha+=2;
				MALLOC_VECTOR(out,double,block_theta+2,ONE);
			}
			memory += (8*tmp*(1+1.0/alpha_int.N)+4)*sizeof(double);
		}
	}
	/* estimate of the memory (only the fastest scaling part):
	 * MatVec - (288+384nprocs/boxX [+192/nprocs])*Ndip
	 *          more exactly: gridX*gridY*gridZ*(36+48nprocs/boxX [+24/nprocs]) value in [] is only for parallel mode.
	 * For surf additionally: gridX*gridY*gridZ*(48+48nprocs/boxX)
	 * 			+ for Sommerfeld table: 128*boxZ*(boxX*boxY-(MIN(boxX,boxY))^2/2)
	 *    For OpenCL mode all MatVec part is allocated on GPU instead of main (CPU) memory (+ a few additional vectors).
	 *    However, OpenCL may additionally use up to 96*min(32,gridX)*gridY*gridZ if available.
	 * others - nvoid_Ndip*{271(CGNR,BiCG), 367(CSYM,QMR2), 415(BiCGStab,QMR), or 463(BCGS2)}
	 *          + additional 8*nvoid_Ndip for OpenCL mode and CGNR or Bi-CGSTAB
	 * PARALLEL: above is total; division over processors of MatVec is uniform, others - according to local_nvoid_Ndip
	 *
	 * Sparse mode - each processor needs (265--457, depending on iterative solver)*local_nvoid_Ndip + 60*nvoid_Ndip
	 *               and division is uniform, i.e. local_nvoid_Ndip = nvoid_Ndip/nprocs
	 *               Sommerfeld table - same as above, but it is not divided among processors.
	 *               Part of the memory is currently not distributed among processors - see issues 160,175.
	 */
	MAXIMIZE(memPeak,memory);
	double memSum=AccumulateMax(memPeak,&memmax);
	if (IFROOT) {
		PrintBoth(logfile,"Total memory usage: "FFORMM" MB\n",memSum/MBYTE);
#ifdef PARALLEL
		PrintBoth(logfile,"Maximum memory usage of single processor: "FFORMM" MB\n",memmax/MBYTE);
#endif
#ifdef OPENCL
		PrintBoth(logfile,"OpenCL memory usage: peak total - "FFORMM" MB, maximum object - "FFORMM" MB\n",
			oclMemPeak/MBYTE,oclMemMaxObj/MBYTE);
#endif
	}
}

//======================================================================================================================

void FreeEverything(void)
/* frees all allocated vectors; should not be called in prognosis mode, since arrays are not
 * actually allocated. Also called from matvec.c in PRECISE_TIMING.
 */
{
	FreeInteraction();
#ifndef SPARSE	
	Free_FFT_Dmat();
	Free_cVector(expsX);
	Free_cVector(expsY);
	Free_cVector(expsZ);
	Free_general(position); // allocated in MakeParticle();
#else	
	Free_general(position_full); // allocated in MakeParticle();
	Free_cVector(arg_full);
#endif // SPARSE
#ifdef ACCIMEXP
	Free_cVector(imexptable);
#endif
	Free_cVector(xvec);
	Free_cVector(rvec);
	Free_cVector(pvec);
	Free_cVector(Einc);
	Free_cVector(Avecbuffer);
	
	/* The following can be automated to some extent, either using the information from structure array 'params' in
	 * iterative.c or checking each vector for being NULL. However, it will anyway require manual editing if additional
	 * (e.g. fourth) vector will be added.
	 */
	switch (IterMethod) {
		case IT_BCGS2:
			Free_cVector(vec1);
			Free_cVector(vec2);
			Free_cVector(vec3);
			Free_cVector(vec4);
			break;
		case IT_CGNR:
		case IT_BICG_CS:
			break;
		case IT_BICGSTAB:
		case IT_QMR_CS:
			Free_cVector(vec1);
			Free_cVector(vec2);
			Free_cVector(vec3);
			break;
		case IT_CSYM:
		case IT_QMR_CS_2:
			Free_cVector(vec1);
			Free_cVector(vec2);
			break;
	}
	/* TO ADD NEW ITERATIVE SOLVER
	 * Add here a case corresponding to the new iterative solver. It should free the extra vectors that were allocated
	 * in AllocateEverything() above.
	 */
	if (yzplane) {
		Free_cVector(EyzplX);
		Free_cVector(EyzplY);
	}
	if (scat_plane) {
		Free_cVector(EplaneX);
		Free_cVector(EplaneY);
	}
	if (all_dir) {
		Free_general(theta_int.val);
		Free_general(phi_int.val);
		Free_cVector(E_ad);
		Free_general(E2_alldir);
	}
	if (scat_grid) {
		Free_general(angles.theta.val);
		Free_general(angles.phi.val);
		Free_cVector(EgridX);
		Free_cVector(EgridY);
		if (phi_integr && IFROOT) {
			Free_general(muel_phi);
			Free_general(muel_phi_buf);
		}
	}
	// these 2 were allocated in MakeParticle
	Free_general(DipoleCoord);
	Free_general(material);

	if (orient_avg) {
		if (IFROOT) {
			if (store_mueller) {
				Free_cVector(ampl_alphaX);
				Free_cVector(ampl_alphaY);
			}
			Free_general(muel_alpha-2);
			Free_general(out);
		}
		Free_general(alpha_int.val);
		Free_general(beta_int.val);
		Free_general(gamma_int.val);
	}
#ifdef OPENCL
	oclunload();
#endif
}

//======================================================================================================================

void Calculator (void)
{
	char fname[MAX_FNAME];

	// initialize variables
#ifdef OPENCL
	TIME_TYPE start_ocl_init=GET_TIME();
	oclinit();
	Timing_OCL_Init=GET_TIME()-start_ocl_init;
#endif

	if (nTheta!=0) {
		dtheta_deg = 180.0 / ((double)(nTheta-1));
		dtheta_rad = Deg2Rad(dtheta_deg);
		block_theta= 16*(size_t)nTheta;
		if (TestExtendThetaRange()) nTheta=2*(nTheta-1);
	}
	else dtheta_deg=dtheta_rad=block_theta=0;
	finish_avg=false;
	// Do preliminary setup for MatVec
	TIME_TYPE startInitInt=GET_TIME();
	InitInteraction();
	Timing_Init_Int=GET_TIME()-startInitInt;
#ifndef SPARSE
	// initialize D matrix (for matrix-vector multiplication)
	D("InitDmatrix started");
	InitDmatrix();
	D("InitDmatrix finished");
#endif // !SPARSE
	// allocate most (that is not already allocated; perform memory analysis
	AllocateEverything();
	// finish initialization
	if (!orient_avg) alpha_int.N=1;
	Timing_Init = GET_TIME() - tstart_main;
	// prognosis stops here
	if (prognosis) return;
	// main calculation part
	if (orient_avg) {
		if (IFROOT) {
			SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_LOG_ORAVG,directory);
			D("Romberg2D started on root");
			Romberg2D(parms,orient_integrand,block_theta+2,out,fname);
			D("Romberg2D finished on root");
			finish_avg=true;
			/* first two are dummy variables; this call corresponds to one in orient_integrand by other processors;
			 * TODO: replace by a call without unnecessary overhead
			 */
			BcastOrient(&finish_avg,&finish_avg,&finish_avg);
			SaveMuellerAndCS(out);
		}
		else while (!finish_avg) orient_integrand(0,0,NULL);
	}
	else calculate_one_orientation(NULL);
	// cleaning
	FreeEverything();
}
