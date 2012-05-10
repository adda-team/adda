/* File: calculator.c
 * $Date::                            $
 * Descr: all the initialization is done here before actually calculating internal fields; includes
 *        calculation of couple constants
 *
 * Copyright (C) 2006-2010,2012 ADDA contributors
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
#include "const.h" // keep this first
// project headers
#include "cmplx.h"
#include "comm.h"
#include "crosssec.h"
#include "debug.h"
#include "fft.h"
#include "io.h"
#include "memory.h"
#include "Romberg.h"
#include "timing.h"
#include "vars.h"
// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef OPENCL
#	include "oclcore.h"
#endif

// SEMI-GLOBAL VARIABLES

// defined and initialized in crosssec.c
extern const Parms_1D parms[2],parms_alpha;
extern const angle_set beta_int,gamma_int,theta_int,phi_int;
// defined and initialized in param.c
extern const int avg_inc_pol;
extern const char *alldir_parms,*scat_grid_parms;
// defined and initialized in timing.c
extern TIME_TYPE Timing_Init;
#ifdef OPENCL
extern TIME_TYPE Timing_OCL_Init;
#endif
extern size_t TotalEval;

// used in CalculateE.c
double * restrict muel_phi; // used to store values of Mueller matrix for different phi (to integrate)
double * restrict muel_phi_buf; // additional for integrating with different multipliers
	// scattered E (for scattering in one plane) for two incident polarizations
doublecomplex * restrict EplaneX, * restrict EplaneY;
double * restrict Eplane_buffer; // buffer to accumulate Eplane
double dtheta_deg,dtheta_rad; // delta theta in degrees and radians
	// amplitude matrix for different values of alpha
doublecomplex * restrict ampl_alphaX,* restrict ampl_alphaY;
double * restrict muel_alpha; // mueller matrix for different values of alpha
// used in fft.c
  // tables of integrals
double * restrict tab1,* restrict tab2,* restrict tab3,* restrict tab4,* restrict tab5,
	* restrict tab6,* restrict tab7,* restrict tab8,* restrict tab9,* restrict tab10;
/* it is preferable to declare the following as "* restrict * restrict", but it is hard to make it
 * generally compatible with Free_iMatrix function syntax. However, it is defined so in fft.c.
 */
int ** restrict tab_index; // matrix for indexing of table arrays
// used in crosssec.c
double * restrict E2_alldir; // square of E, calculated for alldir
double * restrict E2_alldir_buffer; // buffer to accumulate E2_alldir
doublecomplex cc[MAX_NMAT][3]; // couple constants
	// arrays of exponents along 3 axes (for calc_field)
doublecomplex * restrict expsX,* restrict expsY,* restrict expsZ;
// used in iterative.c
doublecomplex *rvec;                 // current residual
doublecomplex * restrict Avecbuffer; // used to hold the result of matrix-vector products
// auxiliary vectors, used in some iterative solvers (with more meaningful names)
doublecomplex * restrict vec1,* restrict vec2,* restrict vec3;

// LOCAL VARIABLES

static size_t block_theta; // size of one block of mueller matrix - 16*nTheta
	// whether to stop orientation averaging; defined as int to simplify MPI casting
static int finish_avg;
	// used to collect both mueller matrix and integral scattering quantities when orient_avg
static double * restrict out;

// EXTERNAL FUNCTIONS

// CalculateE.c
extern int CalculateE(enum incpol which,enum Eftype type);
extern void MuellerMatrix(void);
extern void SaveMuellerAndCS(double * restrict in);

//============================================================

static void CoupleConstant(doublecomplex *mrel,const enum incpol which,doublecomplex *res)
/* hard to maintain. It is better to separate different polarizability relations to make the
 * resulting expressions more understandable.
 */
{
	doublecomplex coup_con[3];
	doublecomplex tempa,tempb,cm,m2,t1;
	double temp,b1,b2,b3,ka,kd2;
	int i,imax,j; // counters: i is for 'asym', j is for 'anysotropy'
	double S,prop2[3];
	int asym; // whether polarizability is asymmetric (for isotropic m)
	const double *incPol;
	bool pol_avg=true; // temporary fixed value for SO polarizability

	// redundant initialization to remove warnings
	b1=b2=b3=S=0;

	asym = (PolRelation==POL_CLDR || PolRelation==POL_SO);
	// !!! this should never happen
	if (asym && anisotropy) LogError(ONE_POS,"Incompatibility error in CoupleConstant");
	if (asym) imax=3;
	else imax=1;
	if (PolRelation==POL_LDR || PolRelation==POL_CLDR) {
		b1=LDR_B1;
		b2=LDR_B2;
		b3=LDR_B3;
	}
	else if (PolRelation==POL_SO) {
		b1=SO_B1;
		b2=SO_B2;
		b3=SO_B3;
	}
	// calculate the CM couple constant CC=(3V/4pi)*(m^2-1)/(m^2+2)
	temp = 3*dipvol/FOUR_PI;
	for (j=0;j<Ncomp;j++) {
		cSquare(mrel[j],m2); // m2=m^2
		tempa[RE] = m2[RE] - 1.0;
		tempa[IM] = tempb[IM] = m2[IM];
		tempb[RE] = m2[RE] + 2.0;
		cDiv(tempa,tempb,coup_con[j]);
		coup_con[j][RE] *= temp;
		coup_con[j][IM] *= temp;

		if (PolRelation!=POL_CM) {
			if (PolRelation==POL_LDR || PolRelation==POL_CLDR || PolRelation==POL_SO) {
				// set prop_i^2
				for (i=0;i<3;i++) {
					if (pol_avg && PolRelation==POL_SO) prop2[i]=ONE_THIRD;
					else prop2[i]=prop[i]*prop[i];
				}
				// determine S coefficient for LDR
				if (PolRelation==POL_LDR) {
					if (avg_inc_pol) S=0.5*(1-DotProd(prop2,prop2));
					else {
						if (which==INCPOL_Y) incPol=incPolY;
						else incPol=incPolX; // which==INCPOL_X
						S = prop2[0]*incPol[0]*incPol[0] + prop2[1]*incPol[1]*incPol[1]
						  + prop2[2]*incPol[2]*incPol[2];
					}
				}
			}
			cEqual(coup_con[j],cm);
			for (i=0;i<imax;i++) {
				if (PolRelation==POL_LAK) { // t1=(8pi/3)[(1-ika)exp(ika)-1]
					ka=LAK_C*kd; // a - radius of volume-equivalent (to cube with size d) sphere
					imExp(ka,tempa);
					tempb[RE]=1;
					tempb[IM]=-ka;
					cMult(tempa,tempb,t1);
					t1[RE]-=1;
					cMultReal(2*FOUR_PI_OVER_THREE,t1,t1);
				}
				else { // other formulations are extensions of RR
					// RR correction
					kd2=kd*kd;
					t1[RE]=0.0;
					t1[IM]=2*kd2*kd/3; // t1=2/3*i*kd^3
					// plus more advanced corrections
					if (PolRelation==POL_DGF) t1[RE]+=DGF_B1*kd2;
					else if (PolRelation==POL_FCD) // t1+={(4/3)kd^2+(2/3pi)log[(pi-kd)/(pi+kd)]kd^3}
						t1[RE]+=2*ONE_THIRD*kd2*(2+kd*INV_PI*log((PI-kd)/(PI+kd)));
					if (PolRelation==POL_IGT_SO) t1[RE]+=SO_B1*kd2;
					else if (PolRelation==POL_LDR || PolRelation==POL_CLDR || PolRelation==POL_SO) {
						if (PolRelation!=POL_LDR) S=prop2[i];
						t1[RE]+=(b1+(b2+b3*S)*m2[RE])*kd2; // t1+=(b1+(b2+b3*S)*m^2)*kd^2
						t1[IM]+=(b2+b3*S)*m2[IM]*kd2;
					}
				}
				// CC[i]=cm/(1-(cm/V)*t1); t1 is the M-term
				cMultReal(1.0/dipvol,t1,t1);
				cMultSelf(t1,cm);
				t1[RE]=1-t1[RE];
				t1[IM]=-t1[IM];
				// 'i+j' is not robust. It assumes that only one counter is used
				cDiv(cm,t1,coup_con[i+j]);
			}
		}
	}
	if (asym || anisotropy) {
		if (!orient_avg && IFROOT) PrintBoth(logfile, "CoupleConstant:"CFORM3V"\n",coup_con[0][RE],
			coup_con[0][IM],coup_con[1][RE],coup_con[1][IM],coup_con[2][RE],coup_con[2][IM]);
	}
	else {
		cEqual(coup_con[0],coup_con[1]);
		cEqual(coup_con[0],coup_con[2]);
		if (!orient_avg && IFROOT) PrintBoth(logfile,"CoupleConstant:"CFORM"\n",
				coup_con[0][RE],coup_con[0][IM]);
	}
	memcpy(res,coup_con,3*sizeof(doublecomplex));
}

//============================================================

static void InitCC(const enum incpol which)
// calculate cc, cc_sqrt, and chi_inv
{
	int i,j;
	doublecomplex chi;

	for(i=0;i<Nmat;i++) {
		CoupleConstant(ref_index+Ncomp*i,which,cc[i]);
		for(j=0;j<3;j++) cSqrt(cc[i][j],cc_sqrt[i][j]);
		// chi_inv=1/(V*chi)=4*PI/(V(m^2-1)); for anisotropic - by components
		for (j=0;j<Ncomp;j++) {
			cSquare(ref_index[Ncomp*i+j],chi);
			chi[RE]-=1;
			cMultReal(dipvol/FOUR_PI,chi,chi);
			cInv(chi,chi_inv[i][j]);
		}
		// copy first component of chi_inv[i] into other two, if they are not calculated explicitly
		if (!anisotropy) {
			cEqual(chi_inv[i][0],chi_inv[i][1]);
			cEqual(chi_inv[i][0],chi_inv[i][2]);
		}
	}
#ifdef OPENCL
	/* this is done here, since InitCC can be run between different runs of the iterative solver
	 * write is blocking to ensure completion before function end
	 */
	CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufcc_sqrt,CL_TRUE,0,sizeof(cc_sqrt),cc_sqrt,0,
		NULL,NULL));
#endif
}

//============================================================

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

//============================================================

static void ReadTables(void)
{
	int i, j, ymax, Rm2, Rm2x;

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

	if (!prognosis) {
		// allocate memory for tab_index
		MALLOC_IMATRIX(tab_index,1,TAB_RMAX,0,TAB_RMAX,ALL);
		// fill tab_index
		Rm2=TAB_RMAX*TAB_RMAX;
		tab_index[1][0] = 0;
		for (i=1; i<=TAB_RMAX; i++) {
			Rm2x=Rm2-i*i;
			ymax = MIN(i,(int)floor(sqrt(Rm2x)));
			for (j=0; j<ymax; j++) {
				tab_index[i][j+1] = tab_index[i][j] + MIN(j,(int)floor(sqrt(Rm2x-j*j)))+1;
			}
			if (i<TAB_RMAX) tab_index[i+1][0] = tab_index[i][ymax]
			                                  + MIN(ymax,(int)floor(sqrt(Rm2x-ymax*ymax)))+1;
		}
	}
	// PRINTZ("P[5,3]=%d (should be 41)\n",tab_index[5][3]);
}

//============================================================

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

//============================================================

static void calculate_one_orientation(double * restrict res)
// performs calculation for one orientation; may do orientation averaging and put the result in res
{
	TIME_TYPE tstart;

	if (orient_avg) {
		alph_deg=0;
		InitRotation();
		if (IFROOT) PrintBoth(logfile,"\nORIENTATION STEP beta="GFORMDEF" gamma="GFORMDEF"\n",
			bet_deg,gam_deg);
	}

	// calculate scattered field for y - polarized incident light
	if (IFROOT) {
		printf("\nhere we go, calc Y\n\n");
		if (!orient_avg) fprintf(logfile,"\nhere we go, calc Y\n\n");
	}
	InitCC(INCPOL_Y);
	if (symR && !scat_grid) {
		if (CalculateE(INCPOL_Y,CE_PARPER)==CHP_EXIT) return;
	}
	else { // no rotational symmetry
		/* TODO: in case of scat_grid we run twice to get the full electric field with incoming
		 * light polarized in X and Y direction. In case of rotational symmetry this is not needed
		 * but requires lots more programming so we leave this optimization to a later time.
		 */
		if(CalculateE(INCPOL_Y,CE_NORMAL)==CHP_EXIT) return;

		if (IFROOT) {
			printf("\nhere we go, calc X\n\n");
			if (!orient_avg) fprintf(logfile,"\nhere we go, calc X\n\n");
		}
		if(PolRelation==POL_LDR && !avg_inc_pol) InitCC(INCPOL_X);

		if(CalculateE(INCPOL_X,CE_NORMAL)==CHP_EXIT) return;
	}
	D("CalculateE finished");
	MuellerMatrix();
	D("MuellerMatrix finished");
	if (IFROOT && orient_avg) {
		tstart=GET_TIME();
		/* it is more logical to use store_mueller in the following test, but for orientation
		 * averaging these two flags are identical
		 */
		if (yzplane) printf("\nError of alpha integration (Mueller) is "GFORMDEF"\n",
			Romberg1D(parms_alpha,block_theta,muel_alpha,res+2));
		memcpy(res,muel_alpha-2,2*sizeof(double));
		D("Integration over alpha completed on root");
		Timing_Integration += GET_TIME() - tstart;
	}
	TotalEval++;
}

//============================================================

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

//============================================================

static void AllocateEverything(void)
// allocates a lot of arrays and performs memory analysis
{
	double tmp;
	size_t temp_int;
	double memmax;

	// redundant initialization to remove warnings
	temp_int=0;

	/* It may be nice to initialize all pointers to NULL here, so that any pointer, which is not
	 * initialized below, will surely stay NULL (independent of a particular compiler). But even
	 * without this forgetting to allocate a necessary vector, will surely cause segmentation fault
	 * afterwards. So we do not implement these extra tests for now.
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
	/* additional vectors for iterative methods. Potentially, this procedure can be fully automated
	 * for any new iterative solver, based on the information contained in structure array 'params'
	 * in file iterative.c. However, this requires different order of function calls to extract this
	 * information beforehand. So currently this part should be edited manually when needed.
	 */
	if (IterMethod==IT_BICGSTAB || IterMethod==IT_QMR_CS) {
		if (!prognosis) {
			MALLOC_VECTOR(vec1,complex,local_nRows,ALL);
			MALLOC_VECTOR(vec2,complex,local_nRows,ALL);
			MALLOC_VECTOR(vec3,complex,local_nRows,ALL);
		}
		memory+=3*tmp;
	}
	else if (IterMethod==IT_CSYM || IterMethod==IT_QMR_CS_2) {
		if (!prognosis) {
			MALLOC_VECTOR(vec1,complex,local_nRows,ALL);
			MALLOC_VECTOR(vec2,complex,local_nRows,ALL);
		}
		memory+=2*tmp;
	}
	/* TO ADD NEW ITERATIVE SOLVER
	 * If the new iterative solver requires any extra vectors (additionally to the default ones),
	 * i.e. number vec_N in corresponding element of structure array params in iterative.c is
	 * non-zero, then change the above condition to allocate memory for these vectors. Variable
	 * memory should be incremented to reflect the total allocated memory.
	 */
	MALLOC_VECTOR(expsX,complex,boxX,ALL);
	MALLOC_VECTOR(expsY,complex,boxY,ALL);
	MALLOC_VECTOR(expsZ,complex,local_Nz_unif,ALL);
	if (yzplane) {
		tmp=2*(double)nTheta;
		if (!prognosis) {
			CheckOverflow(2*tmp,ONE_POS_FUNC);
			temp_int=tmp;
			MALLOC_VECTOR(EplaneX,complex,temp_int,ALL);
			MALLOC_VECTOR(EplaneY,complex,temp_int,ALL);
		}
		memory+=2*tmp*sizeof(doublecomplex);
#ifdef PARALLEL
		/* Buffers like Eplane_buffer are defined always (without "ifdef PARALLEL"), so that
		 * functions like Accumulate may be called even in sequential mode with buffers in their
		 * arguments. Such calls are void, but are good for generality of the code. So only
		 * allocation of these buffers is put inside "ifdef".
		 */
		if (IFROOT) { // buffer for accumulate operation
			if (!prognosis) MALLOC_VECTOR(Eplane_buffer,double,2*temp_int,ONE);
			memory+=2*tmp*sizeof(double);
		}
#endif
	}
	if (all_dir) {
		ReadAlldirParms(alldir_parms);
		/* calculate size of vectors; 4 - because first it is used to store per and par components
		 * of the field, and only afterwards - squares.
		 */
		tmp=4*((double)theta_int.N)*phi_int.N;
		if (!prognosis) {
			CheckOverflow(tmp,ONE_POS_FUNC);
			temp_int=tmp;
			MALLOC_VECTOR(E2_alldir,double,temp_int,ALL);
		}
		memory+=tmp*sizeof(double);
#ifdef PARALLEL
		if (IFROOT) { // buffer for accumulate operation
			if (!prognosis) MALLOC_VECTOR(E2_alldir_buffer,double,temp_int,ONE);
			memory+=tmp*sizeof(double);
		}
#endif
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
#ifdef PARALLEL
		if (IFROOT) { // buffer for accumulate operation
			if (!prognosis) MALLOC_VECTOR(Egrid_buffer,double,2*temp_int,ONE);
			memory+=2*tmp*sizeof(double);
		}
#endif
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
			if (yzplane) {
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
	 *          more exactly: gridX*gridY*gridZ*(36+48nprocs/boxX [+24/nprocs]) value in [] is only
	 *          for parallel mode.
	 *          For OpenCL mode all MatVec part is allocated on GPU instead of main (CPU) memory
	 * others - nvoid_Ndip*{271(CGNR,BiCG),367(CSYM,QMR2), or 415(BiCGStab,QMR)}
	 *          + additional 8*nvoid_Ndip for OpenCL mode and CGNR or Bi-CGSTAB
	 * PARALLEL: above is total; division over processors of MatVec is uniform,
	 *           others - according to local_nvoid_Ndip
	 */
	MAXIMIZE(memPeak,memory);
	double memSum=AccumulateMax(memPeak,&memmax);
	if (IFROOT) {
		PrintBoth(logfile,"Total memory usage: "FFORMM" MB\n",memSum/MBYTE);
#ifdef PARALLEL
		PrintBoth(logfile,"Maximum memory usage of single processor: "FFORMM" MB\n",memmax/MBYTE);
#endif
#ifdef OPENCL
		PrintBoth(logfile,
			"OpenCL memory usage: peak total - "FFORMM" MB, maximum object - "FFORMM" MB\n",
			oclMemPeak/MBYTE,oclMemMaxObj/MBYTE);
#endif
	}
}

//============================================================

static void FreeEverything(void)
/* frees all allocated vectors; should not be called in prognosis mode, since arrays are not
 * actually allocated.
 */
{
	if (IntRelation == G_SO || IntRelation == G_IGT_SO) FreeTables();
	Free_FFT_Dmat();
	Free_cVector(xvec);
	Free_cVector(rvec);
	Free_cVector(pvec);
	Free_cVector(Einc);
	Free_cVector(Avecbuffer);
	/* The following can be automated to some extent, either using the information from structure
	 * array 'params' in iterative.c or checking each vector for being NULL. However, it will anyway
	 * require manual editing if additional (e.g. fourth) vector will be added.
	 */
	if (IterMethod==IT_BICGSTAB || IterMethod==IT_QMR_CS) {
		Free_cVector(vec1);
		Free_cVector(vec2);
		Free_cVector(vec3);
	}
	else if (IterMethod==IT_CSYM || IterMethod==IT_QMR_CS_2) {
		Free_cVector(vec1);
		Free_cVector(vec2);
	}
	/* TO ADD NEW ITERATIVE SOLVER
	 * If the new iterative solver requires any extra vectors (so that they were allocated in
	 * function AllocateEverything() above), then this condition (immediately above) should be
	 * changed to perform freeing of these vectors.
	 */
	Free_cVector(expsX);
	Free_cVector(expsY);
	Free_cVector(expsZ);
	if (yzplane) {
		Free_cVector(EplaneX);
		Free_cVector(EplaneY);
#ifdef PARALLEL
		Free_general(Eplane_buffer);
#endif
	}
	if (all_dir) {
		Free_general(theta_int.val);
		Free_general(phi_int.val);
		Free_general(E2_alldir);
#ifdef PARALLEL
		Free_general(E2_alldir_buffer);
#endif
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
#ifdef PARALLEL
		Free_general(Egrid_buffer);
#endif
	}
	// these 3 were allocated in MakeParticle
	Free_general(DipoleCoord);
	Free_general(position);
	Free_general(material);

	if (orient_avg) {
		if (IFROOT) {
			if (yzplane) {
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

//============================================================

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
		// if not enough symmetry, calculate for +- theta (for one plane)
		if (!(symY || orient_avg)) nTheta=2*(nTheta-1);
	}
	else dtheta_deg=dtheta_rad=block_theta=0;
	finish_avg=false;
	// read tables if needed
	if (IntRelation == G_SO || IntRelation == G_IGT_SO) ReadTables();
	// initialize D matrix (for matrix-vector multiplication)
	D("InitDmatrix started");
	InitDmatrix();
	D("InitDmatrix finished");
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
			/* first two are dummy variables; this call corresponds to similar call in
			 * orient_integrand by other processors;
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
