/* File: fft.c
 * $Date::                            $
 * Descr: initialization of all FFT for matrix-vector products; and FFT procedures themselves; not used in sparse mode
 *        TODO: A lot of indirect indexing used - way to optimize.
 *
 * Copyright (C) 2006-2013 ADDA contributors
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
#include "fft.h" // corresponding header
// project headers
#include "cmplx.h"
#include "comm.h"
#include "debug.h"
#include "function.h"
#include "io.h"
#include "interaction.h"
#include "memory.h"
#include "oclcore.h"
#include "prec_time.h"
#include "vars.h"
// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef CLFFT_AMD
#	include <clAmdFft.h> //external library from AMD
#elif defined(CLFFT_APPLE)
#	include "cpp/clFFT.h" //nearly unmodified APPLE FFT header file
#	ifdef NO_CPP
#		error "Apple clFFT relies on C++ sources, hence is incompatible with NO_CPP option"
#	endif
#endif
/* standard FFT routines (FFTW3 of FFT_TEMPERTON) are required even when OpenCL is used, since
 * they are used for Fourier transform of the D-matrix
 */
#ifdef FFTW3
#	include <fftw3.h>
/* define level of planning for usual and Dmatrix (DM) FFT: FFTW_ESTIMATE (heuristics),
 * FFTW_MEASURE (default), FTW_PATIENT, or FFTW_EXHAUSTIVE
 */
#	define PLAN_FFTW FFTW_MEASURE
#	define PLAN_FFTW_DM FFTW_ESTIMATE
#	define ONLY_FOR_FFTW3 // this is used in function argument declarations
#else
#	define ONLY_FOR_FFTW3 ATT_UNUSED
#endif
// for transpose YZ
#define TR_BLOCK 64

#ifdef FFT_TEMPERTON
#	define ONLY_FOR_TEMPERTON // this is used in function argument declarations
#else
#	define ONLY_FOR_TEMPERTON ATT_UNUSED
#endif

// SEMI-GLOBAL VARIABLES

// defined and initialized in timing.c
extern TIME_TYPE Timing_FFT_Init,Timing_Dm_Init;

// used in matvec.c
doublecomplex * restrict Dmatrix; // holds FFT of the interaction matrix
#ifndef OPENCL
// holds input vector (on expanded grid) to matvec, also used as storage space in iterative.c
doublecomplex * restrict Xmatrix;
doublecomplex * restrict slices; // used in inner cycle of matvec - holds 3 components (for fixed x)
doublecomplex * restrict slices_tr; // additional storage space for slices to accelerate transpose
#endif
size_t DsizeY,DsizeZ,DsizeYZ; // size of the 'matrix' D

// used in comm.c
double * restrict BT_buffer, * restrict BT_rbuffer; // buffers for BlockTranspose

// LOCAL VARIABLES

// D2 matrix and its two slices; used only temporary for InitDmatrix
static doublecomplex * restrict slice,* restrict slice_tr,* restrict D2matrix;
static size_t D2sizeX,D2sizeY,D2sizeZ; // size of the 'matrix' D2
static size_t blockTr=TR_BLOCK;        // block size for TransposeYZ
static bool weird_nprocs;              // whether weird number of processors is used
#ifdef OPENCL
#	ifdef CLFFT_AMD
static clAmdFftPlanHandle clplanX,clplanY,clplanZ;
#	elif defined(CLFFT_APPLE)
static clFFT_Plan clplanX,clplanY,clplanZ; // clFFT plans
#	endif
#endif
#ifdef FFTW3
// FFTW3 plans: f - FFT_FORWARD; b - FFT_BACKWARD
static fftw_plan planXf_Dm,planYf_Dm,planZf_Dm;
#	ifndef OPENCL // these plans are used only if OpenCL is not used
static fftw_plan planXf,planXb,planYf,planYb,planZf,planZb;
#	endif
#elif defined(FFT_TEMPERTON)
#	ifdef NO_FORTRAN
#		error "Tempertron FFT is implemented in Fortran, hence incompatible with NO_FORTRAN option"
#	endif
#	define IFAX_SIZE 20
// arrays for Temperton FFT
static double * restrict trigsX,* restrict trigsY,* restrict trigsZ,* restrict work;
static int ifaxX[IFAX_SIZE],ifaxY[IFAX_SIZE],ifaxZ[IFAX_SIZE];
// Fortran routines from cfft99D.f
void cftfax_(const int *nn,int * restrict ifax,double * restrict trigs);
void cfft99_(double * restrict data,double * restrict _work,const double * restrict trigs,
	const int * restrict ifax,const int *inc,const int *jump,const int *nn,const int *lot,
	const int *isign);
#endif

//============================================================

INLINE size_t IndexDmatrix(const size_t x,size_t y,size_t z)
// index D matrix to store final result
{
	if (y>=DsizeY) y=gridY-y;
	if (z>=DsizeZ) z=gridZ-z;

	return(NDCOMP*(x*DsizeYZ+z*DsizeY+y));
}

//============================================================

INLINE size_t IndexGarbledD(const size_t x,int y,int z,const size_t lengthN UOIP)
// index D2 matrix after BlockTranspose
{
	if (y<0) y+=D2sizeY;
	if (z<0) z+=D2sizeZ;
#ifdef PARALLEL
	return(((z%lengthN)*D2sizeY+y)*gridX+(z/lengthN)*local_Nx+x%local_Nx);
#else
	return((z*D2sizeY+y)*gridX+x);
#endif
}

//============================================================

INLINE size_t IndexD2matrix(int x,int y,int z,const int nnn)
// index D2 matrix to store calculated elements
{
	if (x<0) x+=gridX;
	if (y<0) y+=D2sizeY;
	//  if (z<0) z+=D2sizeZ;
	return(((z-nnn*local_z0)*D2sizeY+y)*gridX+x);
}

//============================================================

INLINE size_t IndexSliceD2matrix(int y,int z)
// index slice of D2 matrix
{
	if (y<0) y+=gridY;
	if (z<0) z+=gridZ;

	return(y*gridZ+z);
}

//============================================================

INLINE size_t IndexSlice_zyD2matrix(const size_t y,const size_t z)
// index transposed slice of D2 matrix
{
	return (z*gridY+y);
}

//============================================================

void TransposeYZ(const int direction)
/* optimized routine to transpose y and z; forward: slices->slices_tr; backward: slices_tr->slices;
 * direction can be made boolean but this contradicts with existing definitions of FFT_FORWARD and
 * FFT_BACKWARD, which themselves are determined by FFT routines invocation format
 */
{
#ifdef OPENCL
	const size_t enqtglobalzy[3]={gridZ,gridY,3};
	const size_t enqtglobalyz[3]={gridY,gridZ,3};
	const size_t tblock[3]={16,16,1}; // this corresponds to BLOCK_DIM in oclkernels.cl
//TODO: test in which cases is the uncached variant faster than the cached one, to make a conditional or
// 	to remove cltransposef/b if cltransposeof/b is allways faster than cltransposef/b
	/* When calling kernels the working group size can't be smaller than the data size; hence cached
	 * kernel can be used only for large enough problems. Alternative solution is to determine the
	 * block size during ADDA runtime and pass it to kernel during its compilation. But using small
	 * block size is not efficient anyway, so falling back to noncached kernel is logical.
	 */
	bool cached=(enqtglobalzy[0]>=tblock[0] && enqtglobalzy[1]>=tblock[1]);
	cached&=(gridZ%16==0 && gridY%16==0); // this is required due to current limitation of cached kernel
	
	if (direction==FFT_FORWARD) {
		if (cached) CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,cltransposeof,3,NULL,
			enqtglobalzy,tblock,0,NULL,NULL));
		else CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,cltransposef,2,NULL,enqtglobalzy,NULL,0,
			NULL,NULL));
	}
	else {
		if (cached) CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,cltransposeob,3,NULL,
			enqtglobalyz,tblock,0,NULL,NULL));
		else CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,cltransposeb,2,NULL,enqtglobalyz,NULL,0,
			NULL,NULL));
	}
	clFinish(command_queue);
#else
	size_t y,z,Y,Z,y1,y2,z1,z2,i,j,y0,z0,Xcomp;
	doublecomplex *t0,*t1,*t2,*t3,*t4,*w0,*w1,*w2,*w3;

	if (direction==FFT_FORWARD) {
		Y=gridY;
		Z=gridZ;
		w0=slices;
		t0=slices_tr-Y;
	}
	else { // direction==FFT_BACKWARD
		Y=gridZ;
		Z=gridY;
		w0=slices_tr;
		t0=slices-Y;
	}

	y1=Y/blockTr;
	y2=Y%blockTr;
	z1=Z/blockTr;
	z2=Z%blockTr;

	for(Xcomp=0;Xcomp<3;Xcomp++) {
		w1=w0+Xcomp*gridYZ;
		t1=t0+Xcomp*gridYZ;
		for(i=0;i<=y1;i++) {
			if (i==y1) y0=y2;
			else y0=blockTr;
			w2=w1;
			t2=t1;
			for(j=0;j<=z1;j++) {
				if (j==z1) z0=z2;
				else z0=blockTr;
				w3=w2;
				t3=t2;
				for (y=0;y<y0;y++) {
					t4=t3+y;
					for (z=0;z<z0;z++) {
						cEqual(w3[z],*(t4+=Y));
					}
					w3+=Z;
				}
				w2+=blockTr;
				t2+=blockTr*Y;
			}
			w1+=blockTr*Z;
			t1+=blockTr;
		}
	}
#endif
}

//============================================================

static void transposeYZ_Dm(doublecomplex *data,doublecomplex *trans)
// optimized routine to transpose y and z for Dmatrix: data -> trans
{
	size_t y,z,Y,Z,y1,y2,z1,z2,i,j,y0,z0;
	doublecomplex *t1,*t2,*t3,*t4,*w1,*w2,*w3;

	Y=gridY;
	Z=gridZ;

	y1=Y/blockTr;
	y2=Y%blockTr;
	z1=Z/blockTr;
	z2=Z%blockTr;

	w1=data;
	t1=trans-Y;

	for(i=0;i<=y1;i++) {
		if (i==y1) y0=y2;
		else y0=blockTr;
		w2=w1;
		t2=t1;
		for(j=0;j<=z1;j++) {
			if (j==z1) z0=z2;
			else z0=blockTr;
			w3=w2;
			t3=t2;
			for (y=0;y<y0;y++) {
				t4=t3+y;
				for (z=0;z<z0;z++) {
					cEqual(w3[z],*(t4+=Y));
				}
				w3+=Z;
			}
			w2+=blockTr;
			t2+=blockTr*Y;
		}
		w1+=blockTr*Z;
		t1+=blockTr;
	}
}

//============================================================

void fftX(const int isign)
// FFT three components of (buf)Xmatrix(x) for all y,z; called from matvec
{
#ifdef OPENCL
#	ifdef CLFFT_AMD
	CL_CH_ERR(clAmdFftEnqueueTransform(clplanX,(clAmdFftDirection)isign,1,&command_queue,0,NULL,NULL,&bufXmatrix,
			NULL,NULL));
#	elif defined(CLFFT_APPLE)
	CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanX,(int)3*local_Nz*smallY,
		(clFFT_Direction)isign,bufXmatrix,bufXmatrix,0,NULL,NULL));
#	endif
	clFinish(command_queue);
#elif defined(FFTW3)
	if (isign==FFT_FORWARD) fftw_execute(planXf);
	else fftw_execute(planXb);
#elif defined(FFT_TEMPERTON)
	int nn=gridX,inc=1,jump=nn,lot=boxY;
	size_t z;

	for (z=0;z<3*local_Nz;z++)
		cfft99_((double *)(Xmatrix+z*gridX*smallY),work,trigsX,ifaxX,&inc,&jump,&nn,&lot,&isign);
#endif
}

//============================================================

void fftY(const int isign)
// FFT three components of slices_tr(y) for all z; called from matvec
{
#ifdef OPENCL
#	ifdef CLFFT_AMD
	CL_CH_ERR(clAmdFftEnqueueTransform(clplanY,(clAmdFftDirection)isign,1,&command_queue,0,NULL,NULL,&bufslices_tr,
			NULL,NULL));
#	elif defined(CLFFT_APPLE)
	CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanY,(int)3*gridZ,(clFFT_Direction)isign,
		bufslices_tr,bufslices_tr,0,NULL,NULL));
#	endif
	clFinish(command_queue);
#elif defined(FFTW3)
	if (isign==FFT_FORWARD) fftw_execute(planYf);
	else fftw_execute(planYb);
#elif defined(FFT_TEMPERTON)
	int nn=gridY,inc=1,jump=nn,lot=3*gridZ;

	cfft99_((double *)(slices_tr),work,trigsY,ifaxY,&inc,&jump,&nn,&lot,&isign);
#endif
}

//============================================================

void fftZ(const int isign)
// FFT three components of slices(z) for all y; called from matvec
{
#ifdef OPENCL
#	ifdef CLFFT_AMD
	CL_CH_ERR(clAmdFftEnqueueTransform(clplanZ,(clAmdFftDirection)isign,1,&command_queue,0,NULL,NULL,&bufslices,
			NULL,NULL));
#	elif defined(CLFFT_APPLE)
	CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanZ,(int)3*gridY,(clFFT_Direction)isign,
		bufslices,bufslices,0,NULL,NULL));
#	endif
	clFinish(command_queue);
#elif defined(FFTW3)
	if (isign==FFT_FORWARD) fftw_execute(planZf);
	else fftw_execute(planZb);
#elif defined(FFT_TEMPERTON)
	int nn=gridZ,inc=1,jump=nn,lot=boxY,Xcomp;

	for (Xcomp=0;Xcomp<3;Xcomp++)
		cfft99_((double *)(slices+gridYZ*Xcomp),work,trigsZ,ifaxZ,&inc,&jump,&nn,&lot,&isign);
#endif
}

//============================================================

static void fftX_Dm(const size_t lengthZ ONLY_FOR_TEMPERTON)
// FFT(forward) D2matrix(x) for all y,z; used for Dmatrix calculation
{
#ifdef FFTW3
	fftw_execute(planXf_Dm);
#elif defined(FFT_TEMPERTON)
	int nn=gridX,inc=1,jump=nn,lot=D2sizeY,isign=FFT_FORWARD;
	size_t z;

	for (z=0;z<lengthZ;z++)
		cfft99_((double *)(D2matrix+z*gridX*D2sizeY),work,trigsX,ifaxX,&inc,&jump,&nn,&lot,&isign);
#endif
}

//============================================================

static void fftY_Dm(void)
// FFT(forward) slice_tr(y) for all z; used for Dmatrix calculation
{
#ifdef FFTW3
	fftw_execute(planYf_Dm);
#elif defined(FFT_TEMPERTON)
	int nn=gridY,inc=1,jump=nn,lot=gridZ,isign=FFT_FORWARD;

	cfft99_((double *)slice_tr,work,trigsY,ifaxY,&inc,&jump,&nn,&lot,&isign);
#endif
}

//============================================================

static void fftZ_Dm(void)
// FFT(forward) slice(z) for all y; used for Dmatrix calculation
{
#ifdef FFTW3
	fftw_execute(planZf_Dm);
#elif defined(FFT_TEMPERTON)
	int nn=gridZ,inc=1,jump=nn,lot=gridY,isign=FFT_FORWARD;

	cfft99_((double *)slice,work,trigsZ,ifaxZ,&inc,&jump,&nn,&lot,&isign);
#endif
}

//============================================================

void CheckNprocs(void)
// checks for consistency the specified number of processors; called in the beginning from InitComm
{
	int y=nprocs;

	// initialize weird_nprocs
	weird_nprocs=false;
	// remove simple prime divisors of y
	while (y%2==0) y/=2;
#ifdef OPENCL
	/* this is redundant, since OpenCL is not currently intended to run in parallel
	 * In the future it should properly handle capabilities of AMD clFFT
	 */
	if (y!=1) PrintError("Specified number of processors (%d) is incompatible with clFFT, since "
		"the latter currently supports only FFTs with size 2^n. Please choose the number of "
		"processors to be of the same form.",nprocs);
#else
	while (y%3==0) y/=3;
	while (y%5==0) y/=5;
#	ifdef FFT_TEMPERTON
	if (y!=1) PrintError("Specified number of processors (%d) is weird (has prime divisors larger "
		"than 5). That is incompatible with Temperton FFT. Revise the number of processors "
		"(recommended) or recompile with FFTW 3 support.",nprocs);
#	elif defined(FFTW3)
	while (y%7==0) y/=7;
	// one multiplier of either 11 or 13 is allowed
	if (y%11==0) y/=11;
	else if (y%13==0) y/=13;
	if (y!=1) {
		LogWarning(EC_WARN,ONE_POS,"Specified number of processors (%d) is weird (has prime "
			"divisors larger than 13 or more than one divisor of either 11 or 13). FFTW3 will work "
			"less efficiently. It is strongly recommended to revise the number of processors.",
			nprocs);
		weird_nprocs=true;
	}
#	endif
#endif
}

//============================================================

int fftFit(int x,int divis)
/* find the first number >=x divisible by 2 only (Apple clFFT) or 2,3,5 only (Temperton FFT or clAMDFFT) or also
 * allowing 7 and one of 11 or 13 (FFTW3), and also divisible by 2 and divis.
 * If weird_nprocs is used, only the latter condition is required.
 */
{
	int y;

	if (weird_nprocs) {
		if (!IS_EVEN(divis)) divis*=2;
		return (divis*((x+divis-1)/divis));
	}
	else while (true) { // not very efficient but robust way
		if (IS_EVEN(x) && x%divis==0) {
			y=x;
			while (y%2==0) y/=2; // here Apple clFFT ends
#ifdef CLFFT_AMD // implies OPENCL
			while (y%3==0) y/=3;
			while (y%5==0) y/=5; // here AMD clFFT ends
#endif
#ifndef OPENCL
			while (y%3==0) y/=3;
			while (y%5==0) y/=5; // here Temperton FFT ends
#	ifdef FFTW3
			while (y%7==0) y/=7;
			// one multiplier of either 11 or 13 is allowed
			if (y%11==0) y/=11;
			else if (y%13==0) y/=13;
#	endif
#endif
			if (y==1) return(x);
		}
		x++;
	}
}

//=============================================================

static void fftInitBeforeD(const int lengthZ ONLY_FOR_FFTW3)
// initialize fft before initialization of Dmatrix
{
#ifdef FFTW3
	int grXint=gridX,grYint=gridY,grZint=gridZ; // this is needed to provide 'int *' to grids

	D("FFTW library version: %s\n     compiler: %s\n     codelet optimizations: %s",
		fftw_version,fftw_cc,fftw_codelet_optim);
	planYf_Dm=fftw_plan_many_dft(1,&grYint,gridZ,slice_tr,NULL,1,gridY,
		slice_tr,NULL,1,gridY,FFT_FORWARD,PLAN_FFTW_DM);
	planZf_Dm=fftw_plan_many_dft(1,&grZint,gridY,slice,NULL,1,gridZ,
		slice,NULL,1,gridZ,FFT_FORWARD,PLAN_FFTW_DM);
	planXf_Dm=fftw_plan_many_dft(1,&grXint,lengthZ*D2sizeY,D2matrix,NULL,1,gridX,
		D2matrix,NULL,1,gridX,FFT_FORWARD,PLAN_FFTW_DM);
#elif defined(FFT_TEMPERTON)
	int size,nn;

	// allocate memory
	MALLOC_VECTOR(trigsX,double,2*gridX,ALL);
	MALLOC_VECTOR(trigsY,double,2*gridY,ALL);
	MALLOC_VECTOR(trigsZ,double,2*gridZ,ALL);
	size=MAX(gridX*D2sizeY,3*gridYZ);
	MALLOC_VECTOR(work,double,2*size,ALL);
	// initialize ifax and trigs
	nn=gridX;
	cftfax_(&nn,ifaxX,trigsX);
	nn=gridY;
	cftfax_(&nn,ifaxY,trigsY);
	nn=gridZ;
	cftfax_(&nn,ifaxZ,trigsZ);
#endif
}

//============================================================

static void fftInitAfterD(void)
/* second part of fft initialization
 * completely separate code is used for OpenCL and FFTW3, because even precise-timing output is
 * significantly different. In particular, FFTW3 uses separate plans for forward and backward, while
 * clFFT (by Apple or AMD) uses one plan for both directions.
 */
{
#ifdef OPENCL
#	ifdef CLFFT_APPLE
	cl_int err; // error code
#	endif
#	ifdef PRECISE_TIMING
	SYSTEM_TIME tvp[4];
#	endif

	if (IFROOT) printf("Initializing clFFT\n");
#	ifdef PRECISE_TIMING
	GetTime(tvp);
#	endif
#	ifdef CLFFT_AMD
	CL_CH_ERR(clAmdFftSetup(NULL)); // first initialize clAmdFft
#	ifdef DEBUGFULL
	cl_uint major,minor,patch;
	CL_CH_ERR(clAmdFftGetVersion(&major,&minor,&patch));
	D("clAmdFft library version - %u.%u.%u",major,minor,patch);
#	endif
	/* Unfortunately, clAmdFft (and Apple clFFT as well) currently supports only simple regular batches of transforms
	 * (similar to fftw_plan_many_dft) but not fully flexible configurations, like offered by fftw_plan_guru_dft. So to
	 * make X transform as a single plan we have have to cycle over the whole smallY instead of (possibly smaller) boxY.
	 * This incurs a small performance hit for "non-standard" values of boxY, but should be overall faster than making
	 * an explicit loop over smaller kernels (like is now done with Temperton FFT).
	 */
	size_t xdimen[3]={gridX,1,1};
	CL_CH_ERR(clAmdFftCreateDefaultPlan(&clplanX,context,CLFFT_1D,xdimen));
	CL_CH_ERR(clAmdFftSetPlanBatchSize(clplanX,3*local_Nz*smallY));
	CL_CH_ERR(clAmdFftSetPlanPrecision(clplanX,CLFFT_DOUBLE));
	CL_CH_ERR(clAmdFftSetLayout(clplanX,CLFFT_COMPLEX_INTERLEAVED,CLFFT_COMPLEX_INTERLEAVED));
	CL_CH_ERR(clAmdFftSetPlanScale(clplanX,FFT_BACKWARD,1)); // override the default (1/N) scale for backward direction
	CL_CH_ERR(clAmdFftBakePlan(clplanX,1,&command_queue,NULL,NULL));
#	elif defined(CLFFT_APPLE)
	clFFT_Dim3 xdimen;
	xdimen.x=(unsigned int)gridX;
	xdimen.y=1;
	xdimen.z=1;
	clplanX=clFFT_CreatePlan(context,xdimen,clFFT_2D,clFFT_InterleavedComplexFormat,&err);
	CL_CH_ERR(err);
#	endif
#	ifdef PRECISE_TIMING
	GetTime(tvp+1);
#	endif
#	ifdef CLFFT_AMD
	size_t ydimen[3]={gridY,1,1};
	CL_CH_ERR(clAmdFftCreateDefaultPlan(&clplanY,context,CLFFT_1D,ydimen));
	CL_CH_ERR(clAmdFftSetPlanBatchSize(clplanY,3*gridZ));
	CL_CH_ERR(clAmdFftSetPlanPrecision(clplanY,CLFFT_DOUBLE));
	CL_CH_ERR(clAmdFftSetLayout(clplanY,CLFFT_COMPLEX_INTERLEAVED,CLFFT_COMPLEX_INTERLEAVED));
	CL_CH_ERR(clAmdFftSetPlanScale(clplanY,FFT_BACKWARD,1)); // override the default (1/N) scale for backward direction
	CL_CH_ERR(clAmdFftBakePlan(clplanY,1,&command_queue,NULL,NULL));
#	elif defined(CLFFT_APPLE)
	clFFT_Dim3 ydimen;
	ydimen.x=(unsigned int)gridY;
	ydimen.y=1;
	ydimen.z=1;
	clplanY=clFFT_CreatePlan(context,ydimen,clFFT_1D,clFFT_InterleavedComplexFormat,&err);
	CL_CH_ERR(err);
#	endif
#	ifdef PRECISE_TIMING
	GetTime(tvp+2);
#	endif
#	ifdef CLFFT_AMD
	/* Here the issue is similar to clplanX described above. However, we are using full gridY instead of boxY, which
	 * incurs at least a-factor-of-two performance hit. To solve this problem one need to execute separate plans for 3
	 * components of vectors. Unfortunately, this cannot be simply done using 3-element loop (like in Temperton FFT),
	 * because a part of cl_mem object can't be addressed independently (as can be done with simple C arrays). The only
	 * way to address this issue is to either create three separate cl_mem objects or to change the indexing of levels
	 * inside the array, so that 3 components are stored together.
	 */
	size_t zdimen[3]={gridZ,1,1};
	CL_CH_ERR(clAmdFftCreateDefaultPlan(&clplanZ,context,CLFFT_1D,zdimen));
	CL_CH_ERR(clAmdFftSetPlanBatchSize(clplanZ,3*gridY));
	CL_CH_ERR(clAmdFftSetPlanPrecision(clplanZ,CLFFT_DOUBLE));
	CL_CH_ERR(clAmdFftSetLayout(clplanZ,CLFFT_COMPLEX_INTERLEAVED,CLFFT_COMPLEX_INTERLEAVED));
	CL_CH_ERR(clAmdFftSetPlanScale(clplanZ,FFT_BACKWARD,1)); // override the default (1/N) scale for backward direction
	CL_CH_ERR(clAmdFftBakePlan(clplanZ,1,&command_queue,NULL,NULL));
#	elif defined(CLFFT_APPLE)
	clFFT_Dim3 zdimen;
	zdimen.x=(unsigned int)gridZ;
	zdimen.y=1;
	zdimen.z=1;
	clplanZ=clFFT_CreatePlan(context,zdimen,clFFT_1D,clFFT_InterleavedComplexFormat,&err);
	CL_CH_ERR(err);
#	endif
#	ifdef PRECISE_TIMING
	GetTime(tvp+3);
	// print precise timing of FFT planning
	if (IFROOT) PrintBoth(logfile,
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"         clFFT planning       \n"
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"X = "FFORMPT"  Total = "FFORMPT"\n"
		"Y = "FFORMPT"\n"
		"Z = "FFORMPT"\n\n",
		DiffSec(tvp,tvp+1),DiffSec(tvp,tvp+3),DiffSec(tvp+1,tvp+2),DiffSec(tvp+2,tvp+3));
#	endif
#elif defined(FFTW3) // this is not needed when OpenCL is used
	int lot;
	fftw_iodim dims,howmany_dims[2];
	int grYint=gridY; // this is needed to provide 'int *' to gridY
#	ifdef PRECISE_TIMING
	SYSTEM_TIME tvp[7];
#	endif
	if (IFROOT) printf("Initializing FFTW3\n");
#	ifdef PRECISE_TIMING
	GetTime(tvp);
#	endif
	lot=3*gridZ;
	planYf=fftw_plan_many_dft(1,&grYint,lot,slices_tr,NULL,1,gridY,
		slices_tr,NULL,1,gridY,FFT_FORWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GetTime(tvp+1);
#	endif
	planYb=fftw_plan_many_dft(1,&grYint,lot,slices_tr,NULL,1,gridY,
		slices_tr,NULL,1,gridY,FFT_BACKWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GetTime(tvp+2);
#	endif
	dims.n=gridZ;
	dims.is=dims.os=1;
	howmany_dims[0].n=3;
	howmany_dims[0].is=howmany_dims[0].os=gridZ*gridY;
	howmany_dims[1].n=boxY;
	howmany_dims[1].is=howmany_dims[1].os=gridZ;
	planZf=fftw_plan_guru_dft(1,&dims,2,howmany_dims,slices,slices,FFT_FORWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GetTime(tvp+3);
#	endif
	planZb=fftw_plan_guru_dft(1,&dims,2,howmany_dims,slices,slices,FFT_BACKWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GetTime(tvp+4);
#	endif
	dims.n=gridX;
	dims.is=dims.os=1;
	howmany_dims[0].n=3*local_Nz;
	howmany_dims[0].is=howmany_dims[0].os=smallY*gridX;
	howmany_dims[1].n=boxY;
	howmany_dims[1].is=howmany_dims[1].os=gridX;
	planXf=fftw_plan_guru_dft(1,&dims,2,howmany_dims,Xmatrix,Xmatrix,FFT_FORWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GetTime(tvp+5);
#	endif
	planXb=fftw_plan_guru_dft(1,&dims,2,howmany_dims,Xmatrix,Xmatrix,FFT_BACKWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GetTime(tvp+6);
	// print precise timing of FFT planning
	if (IFROOT) PrintBoth(logfile,
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"         FFTW3 planning       \n"
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"Yf = "FFORMPT"  Total = "FFORMPT"\n"
		"Yb = "FFORMPT"\n"
		"Zf = "FFORMPT"\n"
		"Zb = "FFORMPT"\n"
		"Xf = "FFORMPT"\n"
		"Xb = "FFORMPT"\n\n",
		DiffSec(tvp,tvp+1),DiffSec(tvp,tvp+6),DiffSec(tvp+1,tvp+2),DiffSec(tvp+2,tvp+3),
		DiffSec(tvp+3,tvp+4),DiffSec(tvp+4,tvp+5),DiffSec(tvp+5,tvp+6));
#	endif
	// destroy old plans
	fftw_destroy_plan(planXf_Dm);
	fftw_destroy_plan(planYf_Dm);
	fftw_destroy_plan(planZf_Dm);
#endif
}

//============================================================

void InitDmatrix(void)
/* Initializes the matrix D. D[i][j][k]=A[i1-i2][j1-j2][k1-k2]. Actually D=-FFT(G)/Ngrid.
 * Then -G.x=invFFT(D*FFT(x)) for practical implementation of FFT such that invFFT(FFT(x))=Ngrid*x.
 * G is exactly Green's tensor. The routine is called only once, so does not need to be very fast,
 * however we tried to optimize it.
 */
{
	int i,j,k,kcor,Dcomp;
	size_t x,y,z,indexfrom,indexto,ind,index,Dsize,D2sizeTot;
	double invNgrid;
	int nnn; // multiplier used for reduced_FFT or not reduced; 1 or 2
	int jstart, kstart;
	size_t lengthN;
	TIME_TYPE start,time1;
#ifdef PRECISE_TIMING
	// precise timing of the Dmatrix computation
	SYSTEM_TIME tvp[13];
	SYSTEM_TIME Timing_fftX,Timing_fftY,Timing_fftZ,Timing_Gcalc,Timing_ar1,Timing_ar2,Timing_ar3,
	Timing_BT,Timing_TYZ,Timing_beg;
	double t_fftX,t_fftY,t_fftZ,t_ar1,t_ar2,t_ar3,t_TYZ,t_beg,t_Gcalc,t_Arithm,t_FFT,t_BT,t_InitMV;

	// This should be the first occurrence of PRECISE_TIMING in the program
	SetTimerFreq();

	InitTime(&Timing_fftX);
	InitTime(&Timing_fftY);
	InitTime(&Timing_fftZ);
	InitTime(&Timing_ar1);
	InitTime(&Timing_ar2);
	InitTime(&Timing_ar3);
	InitTime(&Timing_BT);
	InitTime(&Timing_TYZ);
	GetTime(tvp);
#endif
	start=GET_TIME();

	// initialize sizes of D and D2 matrices
	D2sizeX=gridX;
	if (reduced_FFT) {
		D2sizeY=gridY/2;
		D2sizeZ=gridZ/2;
		DsizeY=gridY/2+1;
		DsizeZ=gridZ/2+1;
		nnn=1;
		jstart=0;
		kstart=0;
	}
	else {
		D2sizeY=DsizeY=gridY;
		D2sizeZ=DsizeZ=gridZ;
		nnn=2;
		jstart=1-boxY;
		kstart=1-boxZ;
	}
	// auxiliary parameters
	lengthN=nnn*local_Nz;
	DsizeYZ=DsizeY*DsizeZ;
	invNgrid=1.0/(gridX*((double)gridYZ));
	local_Nsmall=(gridX/2)*(gridYZ/(2*nprocs)); // size of X vector (for 1 component)
	// potentially this may cause unnecessary error during prognosis, but makes code cleaner
	Dsize=MultOverflow(NDCOMP*local_Nx,DsizeYZ,ONE_POS_FUNC);
	D2sizeTot=nnn*local_Nz*D2sizeY*D2sizeX; // this should be approximately equal to Dsize/NDCOMP
	if (IFROOT) fprintf(logfile,"The FFT grid is: %zux%zux%zu\n",gridX,gridY,gridZ);
#ifdef OPENCL // perform setting up of buffers and kernels
	/* create all Buffers needed on Device in MatVec
	 * When prognosis, the following code just counts required memory
	 */
	CREATE_CL_BUFFER(bufXmatrix,CL_MEM_READ_WRITE,local_Nsmall*3*sizeof(doublecomplex),NULL);
	CREATE_CL_BUFFER(bufargvec,CL_MEM_READ_WRITE,local_nRows*sizeof(doublecomplex),NULL);
	CREATE_CL_BUFFER(bufresultvec,CL_MEM_READ_WRITE,local_nRows*sizeof(doublecomplex),NULL);
	CREATE_CL_BUFFER(bufslices,CL_MEM_READ_WRITE,gridYZ*3*sizeof(doublecomplex),NULL);
	CREATE_CL_BUFFER(bufslices_tr,CL_MEM_READ_WRITE,gridYZ*3*sizeof(doublecomplex),NULL);
	/* The following are constant device buffers which are initialized with host data. But
	 * bufDmatrix is initialized in the end of this function (to be compatible with prognosis. And
	 * bufcc_sqrt is initialized in InitCC, since it may change for every run of the iterative
	 * solver.
	 */
	CREATE_CL_BUFFER(bufcc_sqrt,CL_MEM_READ_ONLY,sizeof(cc_sqrt),NULL);
	CREATE_CL_BUFFER(bufDmatrix,CL_MEM_READ_ONLY,Dsize*sizeof(*Dmatrix),NULL);
	CREATE_CL_BUFFER(bufmaterial,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		local_nvoid_Ndip*sizeof(*material),material);
	CREATE_CL_BUFFER(bufposition,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
		local_nRows*sizeof(*position),position);
	if (ipr_required) { // for inner product (only if it will be used afterwards)
		memory+=local_nvoid_Ndip*sizeof(double);
		CREATE_CL_BUFFER(bufinproduct,CL_MEM_READ_WRITE,local_nvoid_Ndip*sizeof(double),NULL);
		if (!prognosis) {
			MALLOC_VECTOR(inprodhlp,double,local_nvoid_Ndip,ALL);
			CL_CH_ERR(clSetKernelArg(clinprod,0,sizeof(cl_mem),&bufinproduct));
			CL_CH_ERR(clSetKernelArg(clinprod,1,sizeof(cl_mem),&bufresultvec));
		}
	}
	if (!prognosis) { // Setting kernel arguments which are always the same
		// for arith1
		CL_CH_ERR(clSetKernelArg(clarith1,0,sizeof(cl_mem),&bufmaterial));
		CL_CH_ERR(clSetKernelArg(clarith1,1,sizeof(cl_mem),&bufposition));
		CL_CH_ERR(clSetKernelArg(clarith1,2,sizeof(cl_mem),&bufcc_sqrt));
		CL_CH_ERR(clSetKernelArg(clarith1,3,sizeof(cl_mem),&bufargvec));
		CL_CH_ERR(clSetKernelArg(clarith1,4,sizeof(cl_mem),&bufXmatrix));
		CL_CH_ERR(clSetKernelArg(clarith1,5,sizeof(size_t),&local_Nsmall));
		CL_CH_ERR(clSetKernelArg(clarith1,6,sizeof(size_t),&smallY));
		CL_CH_ERR(clSetKernelArg(clarith1,7,sizeof(size_t),&gridX));
		// for arith2
		CL_CH_ERR(clSetKernelArg(clarith2,0,sizeof(cl_mem),&bufXmatrix));
		CL_CH_ERR(clSetKernelArg(clarith2,1,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(clarith2,2,sizeof(size_t),&gridZ));
		CL_CH_ERR(clSetKernelArg(clarith2,3,sizeof(size_t),&smallY));
		CL_CH_ERR(clSetKernelArg(clarith2,4,sizeof(size_t),&gridX));
		CL_CH_ERR(clSetKernelArg(clarith2,5,sizeof(size_t),&gridYZ));
		CL_CH_ERR(clSetKernelArg(clarith2,6,sizeof(size_t),&local_Nsmall));
		// for arith3
		CL_CH_ERR(clSetKernelArg(clarith3,0,sizeof(cl_mem),&bufslices_tr));
		CL_CH_ERR(clSetKernelArg(clarith3,1,sizeof(cl_mem),&bufDmatrix));
		CL_CH_ERR(clSetKernelArg(clarith3,2,sizeof(size_t),&local_x0));
		CL_CH_ERR(clSetKernelArg(clarith3,3,sizeof(size_t),&smallY));
		CL_CH_ERR(clSetKernelArg(clarith3,4,sizeof(size_t),&smallZ));
		CL_CH_ERR(clSetKernelArg(clarith3,5,sizeof(size_t),&gridX));
		CL_CH_ERR(clSetKernelArg(clarith3,6,sizeof(size_t),&DsizeY));
		CL_CH_ERR(clSetKernelArg(clarith3,7,sizeof(size_t),&DsizeZ));
		// for arith4
		CL_CH_ERR(clSetKernelArg(clarith4,0,sizeof(cl_mem),&bufXmatrix));
		CL_CH_ERR(clSetKernelArg(clarith4,1,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(clarith4,2,sizeof(size_t),&gridZ));
		CL_CH_ERR(clSetKernelArg(clarith4,3,sizeof(size_t),&smallY));
		CL_CH_ERR(clSetKernelArg(clarith4,4,sizeof(size_t),&gridX));
		CL_CH_ERR(clSetKernelArg(clarith4,5,sizeof(size_t),&gridYZ));
		CL_CH_ERR(clSetKernelArg(clarith4,6,sizeof(size_t),&local_Nsmall));
		// for arith5
		CL_CH_ERR(clSetKernelArg(clarith5,0,sizeof(cl_mem),&bufmaterial));
		CL_CH_ERR(clSetKernelArg(clarith5,1,sizeof(cl_mem),&bufposition));
		CL_CH_ERR(clSetKernelArg(clarith5,2,sizeof(cl_mem),&bufcc_sqrt));
		CL_CH_ERR(clSetKernelArg(clarith5,3,sizeof(cl_mem),&bufargvec));
		CL_CH_ERR(clSetKernelArg(clarith5,4,sizeof(cl_mem),&bufXmatrix));
		CL_CH_ERR(clSetKernelArg(clarith5,5,sizeof(size_t),&local_Nsmall));
		CL_CH_ERR(clSetKernelArg(clarith5,6,sizeof(size_t),&smallY));
		CL_CH_ERR(clSetKernelArg(clarith5,7,sizeof(size_t),&gridX));
		CL_CH_ERR(clSetKernelArg(clarith5,8,sizeof(cl_mem),&bufresultvec));
		// for transpose forward
		CL_CH_ERR(clSetKernelArg(cltransposef,0,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(cltransposef,1,sizeof(cl_mem),&bufslices_tr));
		CL_CH_ERR(clSetKernelArg(cltransposef,2,sizeof(size_t),&gridZ));
		CL_CH_ERR(clSetKernelArg(cltransposef,3,sizeof(size_t),&gridY));
		// for transpose backward
		CL_CH_ERR(clSetKernelArg(cltransposeb,0,sizeof(cl_mem),&bufslices_tr));
		CL_CH_ERR(clSetKernelArg(cltransposeb,1,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(cltransposeb,2,sizeof(size_t),&gridY));
		CL_CH_ERR(clSetKernelArg(cltransposeb,3,sizeof(size_t),&gridZ));
		//faster transpose kernel with cache
		//(maybe not faster for all sizes so keep the old kernel for special conditions)
		CL_CH_ERR(clSetKernelArg(cltransposeof,0,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(cltransposeof,1,sizeof(cl_mem),&bufslices_tr));
		CL_CH_ERR(clSetKernelArg(cltransposeof,2,sizeof(size_t),&gridZ));
		CL_CH_ERR(clSetKernelArg(cltransposeof,3,sizeof(size_t),&gridY));
		//setting up local cache size as 17*16*3 elements
		//note: a block is only 16*16, but 1*16 stride is needed
		//to avoid bank conflicts
		CL_CH_ERR(clSetKernelArg(cltransposeof,4,17*16*3*sizeof(doublecomplex),NULL));
		CL_CH_ERR(clSetKernelArg(cltransposeob,0,sizeof(cl_mem),&bufslices_tr));
		CL_CH_ERR(clSetKernelArg(cltransposeob,1,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(cltransposeob,2,sizeof(size_t),&gridY));
		CL_CH_ERR(clSetKernelArg(cltransposeob,3,sizeof(size_t),&gridZ));
		CL_CH_ERR(clSetKernelArg(cltransposeob,4,17*16*3*sizeof(doublecomplex),NULL));
	}
#endif
	MAXIMIZE(memPeak,memory);
	// objects which are always allocated (at least temporarily): Dmatrix,D2matrix,slice,slice_tr
	memPeak+=sizeof(doublecomplex)*((double)Dsize+D2sizeTot+2*gridYZ);
#ifndef OPENCL
	// allocated memory that is used further on, not relevant for OpenCL version
	double mem=sizeof(doublecomplex)*((double)Dsize+3*local_Nsmall+6*gridYZ);
#ifdef PARALLEL
	size_t BTsize = 6*smallY*local_Nz*local_Nx; // in doubles
	mem+=2*BTsize*sizeof(double);
#endif
	// printout some information
	if (IFROOT) {
#ifdef PARALLEL
		PrintBoth(logfile,"Memory usage for MatVec matrices (per processor): "FFORMM" MB\n",
			mem/MBYTE);
#else
		PrintBoth(logfile,"Memory usage for MatVec matrices: "FFORMM" MB\n",mem/MBYTE);
#endif
	}
	memory+=mem;
#endif
	if (prognosis) return;
	// allocate memory for Dmatrix
	MALLOC_VECTOR(Dmatrix,complex,Dsize,ALL);
	// allocate memory for D2matrix components
	MALLOC_VECTOR(D2matrix,complex,D2sizeTot,ALL);
	MALLOC_VECTOR(slice,complex,gridYZ,ALL);
	MALLOC_VECTOR(slice_tr,complex,gridYZ,ALL);
	// actually allocation of Xmatrix, slices, slices_tr is below after freeing of Dmatrix and its slice
#ifdef PARALLEL
	// allocate buffer for BlockTranspose_Dm
	size_t bufsize = 2*lengthN*D2sizeY*local_Nx;
	MALLOC_VECTOR(BT_buffer,double,bufsize,ALL);
	MALLOC_VECTOR(BT_rbuffer,double,bufsize,ALL);
#endif
	D("Initialize FFT (1st part)");
	fftInitBeforeD(lengthN);
#ifdef PRECISE_TIMING
	GetTime(tvp+1);
	Elapsed(tvp,tvp+1,&Timing_beg); // it includes a lot of OpenCL stuff
#endif
	if (IFROOT) printf("Calculating Green's function (Dmatrix)\n");
	/* Interaction matrix values are calculated all at once for performance reasons. They are stored
	 * in Dmatrix with indexing corresponding to D2matrix (to facilitate copying) but NDCOMP
	 * elements instead of one. Afterwards they are replaced by Fourier transforms (with different
	 * indexing) component-wise (in cycle over NDCOMP)
	 */

	/* fill Dmatrix with 0, this if to fill the possible gap between e.g. boxY and gridY/2; (and for R=0)
	 * probably faster than using a lot of conditionals
	 */
	for (ind=0;ind<Dsize;ind++) Dmatrix[ind][RE]=Dmatrix[ind][IM]=0;
	// fill Dmatrix with values of Green's tensor
	for(k=nnn*local_z0;k<nnn*local_z1;k++) {
		// correction of k is relevant only if reduced_FFT is not used
		if (k>(int)smallZ) kcor=k-gridZ;
		else kcor=k;
		for (j=jstart;j<boxY;j++) for (i=1-boxX;i<boxX;i++) {
			index=NDCOMP*IndexD2matrix(i,j,k,nnn);
			/* The test for zero distance is somewhat non-optimal. However, other alternatives are not perfect either:
			 * 1) complicate the loops to remove the zero element in the beginning (move tests to the upper level)
			 * 2) call the function with zero - it will produce NaN. Then set this element to zero after the loop.
			 */
			if (i!=0 || j!=0 || kcor!=0) (*CalcInterTerm)(i,j,kcor,Dmatrix+index);
		}
	} // end of i,j,k loop
	if (IFROOT) printf("Fourier transform of Dmatrix");
#ifdef PRECISE_TIMING
	GetTime(tvp+2);
	Elapsed(tvp+1,tvp+2,&Timing_Gcalc);
#endif
	for(Dcomp=0;Dcomp<NDCOMP;Dcomp++) { // main cycle over components of Dmatrix
#ifdef PRECISE_TIMING
		GetTime(tvp+2); // same as the last before cycle
#endif
		// fill D2matrix with precomputed values from Dmatrix
		for (ind=0;ind<D2sizeTot;ind++) cEqual(Dmatrix[NDCOMP*ind+Dcomp],D2matrix[ind]);
#ifdef PRECISE_TIMING
		GetTime(tvp+3);
		ElapsedInc(tvp+2,tvp+3,&Timing_ar1);
#endif
		fftX_Dm(lengthN); // fftX D2matrix
#ifdef PRECISE_TIMING
		GetTime(tvp+4);
		ElapsedInc(tvp+3,tvp+4,&Timing_fftX);
#endif
		BlockTranspose_Dm(D2matrix,D2sizeY,lengthN);
#ifdef PRECISE_TIMING
		GetTime(tvp+5);
		ElapsedInc(tvp+4,tvp+5,&Timing_BT);
#endif
		for(x=local_x0;x<local_x1;x++) {
#ifdef PRECISE_TIMING
			GetTime(tvp+6);
#endif
			for (ind=0;ind<gridYZ;ind++) slice[ind][RE]=slice[ind][IM]=0.0; // fill slice with 0.0
			for(j=jstart;j<boxY;j++) for(k=kstart;k<boxZ;k++) {
				indexfrom=IndexGarbledD(x,j,k,lengthN);
				indexto=IndexSliceD2matrix(j,k);
				cEqual(D2matrix[indexfrom],slice[indexto]);
			}
			if (reduced_FFT) {
				for(j=1;j<boxY;j++) for(k=0;k<boxZ;k++) {
					// mirror along y
					indexfrom=IndexSliceD2matrix(j,k);
					indexto=IndexSliceD2matrix(-j,k);
					if (Dcomp==1 || Dcomp==4) cInvSign2(slice[indexfrom],slice[indexto]);
					else cEqual(slice[indexfrom],slice[indexto]);
				}
				for(j=1-boxY;j<boxY;j++) for(k=1;k<boxZ;k++) {
					// mirror along z
					indexfrom=IndexSliceD2matrix(j,k);
					indexto=IndexSliceD2matrix(j,-k);
					if (Dcomp==2 || Dcomp==4) cInvSign2(slice[indexfrom],slice[indexto]);
					else cEqual(slice[indexfrom],slice[indexto]);
				}
			}
#ifdef PRECISE_TIMING
			GetTime(tvp+7);
			ElapsedInc(tvp+6,tvp+7,&Timing_ar2);
#endif
			fftZ_Dm(); // fftZ slice
#ifdef PRECISE_TIMING
			GetTime(tvp+8);
			ElapsedInc(tvp+7,tvp+8,&Timing_fftZ);
#endif
			transposeYZ_Dm(slice,slice_tr);
#ifdef PRECISE_TIMING
			GetTime(tvp+9);
			ElapsedInc(tvp+8,tvp+9,&Timing_TYZ);
#endif
			fftY_Dm(); // fftY slice_tr
#ifdef PRECISE_TIMING
			GetTime(tvp+10);
			ElapsedInc(tvp+9,tvp+10,&Timing_fftY);
#endif
			for(z=0;z<DsizeZ;z++) for(y=0;y<DsizeY;y++) {
				indexto=IndexDmatrix(x-local_x0,y,z)+Dcomp;
				indexfrom=IndexSlice_zyD2matrix(y,z);
				cMultReal(-invNgrid,slice_tr[indexfrom],Dmatrix[indexto]);
			}
#ifdef PRECISE_TIMING
			GetTime(tvp+11);
			ElapsedInc(tvp+10,tvp+11,&Timing_ar3);
#endif
		} // end slice X
		if (IFROOT) printf(".");
	} // end of Dcomp
	// free vectors used for computation of Dmatrix
	Free_cVector(D2matrix);
	Free_cVector(slice);
	Free_cVector(slice_tr);
#ifdef PARALLEL
	// deallocate buffers for BlockTranspose_Dm
	Free_general(BT_buffer);
	Free_general(BT_rbuffer);
	// allocate buffers for BlockTranspose
	MALLOC_VECTOR(BT_buffer,double,BTsize,ALL);
	MALLOC_VECTOR(BT_rbuffer,double,BTsize,ALL);
#endif
#ifdef OPENCL
	// copy Dmatrix to OpenCL buffer, blocking to ensure completion before function end
	CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufDmatrix,CL_TRUE,0,Dsize*sizeof(*Dmatrix),
		Dmatrix,0,NULL,NULL));
	Free_cVector(Dmatrix);
#else
	// allocate memory for Xmatrix, slices and slices_tr - used in matvec
	MALLOC_VECTOR(Xmatrix,complex,3*local_Nsmall,ALL);
	MALLOC_VECTOR(slices,complex,3*gridYZ,ALL);
	MALLOC_VECTOR(slices_tr,complex,3*gridYZ,ALL);
#endif
	if (IFROOT) printf("\n");
	time1=GET_TIME();
	Timing_Dm_Init=time1-start;

#ifdef PRECISE_TIMING
	GetTime(tvp+12);
	// time for extra initialization required for MatVec; it includes copying Dmatrix to GPU
	t_InitMV=DiffSec(tvp+11,tvp+12);
	// analyze and print precise timing information
	t_beg=TimerToSec(&Timing_beg);
	t_Gcalc=TimerToSec(&Timing_Gcalc);
	t_ar1=TimerToSec(&Timing_ar1);
	t_ar2=TimerToSec(&Timing_ar2);
	t_ar3=TimerToSec(&Timing_ar3);
	t_fftX=TimerToSec(&Timing_fftX);
	t_fftY=TimerToSec(&Timing_fftY);
	t_fftZ=TimerToSec(&Timing_fftZ);
	t_TYZ=TimerToSec(&Timing_TYZ);
	t_BT=TimerToSec(&Timing_BT);
	t_Arithm=t_beg+t_Gcalc+t_ar1+t_ar2+t_ar3+t_TYZ;
	t_FFT=t_fftX+t_fftY+t_fftZ;

	if (IFROOT) PrintBoth(logfile,
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"            Init Dmatrix timing            \n"
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"Begin  = "FFORMPT"    Arithmetics = "FFORMPT"\n"
		"Gcalc  = "FFORMPT"    FFT         = "FFORMPT"\n"
		"Arith1 = "FFORMPT"    Comm        = "FFORMPT"\n"
		"FFTX   = "FFORMPT"    Init MatVec = "FFORMPT"\n"
		"BT     = "FFORMPT"\n"
		"Arith2 = "FFORMPT"          Total = "FFORMPT"\n"
		"FFTZ   = "FFORMPT"\n"
		"TYZ    = "FFORMPT"\n"
		"FFTY   = "FFORMPT"\n"
		"Arith3 = "FFORMPT"\n"
		"InitMV = "FFORMPT"\n\n",
		t_beg,t_Arithm,t_Gcalc,t_FFT,t_ar1,t_BT,t_fftX,t_InitMV,t_BT,
		t_ar2,DiffSec(tvp,tvp+12),t_fftZ,t_TYZ,t_fftY,t_ar3,t_InitMV);
#endif

	fftInitAfterD();

	Timing_FFT_Init = GET_TIME()-time1;
}

//============================================================

void Free_FFT_Dmat(void)
// free all vectors that were allocated in fft.c (all used for FFT and MatVec)
{
#ifdef OPENCL
	my_clReleaseBuffer(bufXmatrix);
	my_clReleaseBuffer(bufmaterial);
	my_clReleaseBuffer(bufposition);
	my_clReleaseBuffer(bufcc_sqrt);
	my_clReleaseBuffer(bufargvec);
	my_clReleaseBuffer(bufresultvec);
	my_clReleaseBuffer(bufslices);
	my_clReleaseBuffer(bufslices_tr);
	my_clReleaseBuffer(bufDmatrix);
	if (ipr_required) {
		my_clReleaseBuffer(bufinproduct);
		Free_general(inprodhlp);
	}
#	ifdef CLFFT_AMD
	clAmdFftTeardown();
#	elif defined(CLFFT_APPLE)
	clFFT_DestroyPlan(clplanX);
	clFFT_DestroyPlan(clplanY);
	clFFT_DestroyPlan(clplanZ);
#	endif
	if (oclMem>0) LogWarning(EC_WARN,ALL_POS,
		"Possible leak of OpenCL memory (size %zu bytes) detected",oclMem);
#else
	Free_cVector(Dmatrix);
	Free_cVector(Xmatrix);
	Free_cVector(slices);
	Free_cVector(slices_tr);
#	ifdef PARALLEL
	Free_general(BT_buffer);
	Free_general(BT_rbuffer);
#	endif
#	ifdef FFTW3 // these plans are defined only when OpenCL is not used
	fftw_destroy_plan(planXf);
	fftw_destroy_plan(planXb);
	fftw_destroy_plan(planYf);
	fftw_destroy_plan(planYb);
	fftw_destroy_plan(planZf);
	fftw_destroy_plan(planZb);
#	endif
#endif
#ifdef FFT_TEMPERTON // these vectors are used even with OpenCL
	Free_general(work);
	Free_general(trigsX);
	Free_general(trigsY);
	Free_general(trigsZ);
#endif
}
