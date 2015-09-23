/* File: fft.c
 * $Date::                            $
 * Descr: initialization of all FFT for matrix-vector products; and FFT procedures themselves; not used in sparse mode
 *        TODO: A lot of indirect indexing used - way to optimize.
 *
 * Copyright (C) 2006-2014 ADDA contributors
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
	IGNORE_WARNING(-Wstrict-prototypes) // no way to change the library header
#	include <clFFT.h> //external library
	STOP_IGNORE
	// Defines precision of clFFT transforms. !!! CLFFT_DOUBLE_FAST should be tested when becomes operational
#	define PRECISION_CLFFT CLFFT_DOUBLE
#elif defined(CLFFT_APPLE)
#	include "cpp/clFFT.h" //nearly unmodified APPLE FFT header file
#	ifdef NO_CPP
#		error "Apple clFFT relies on C++ sources, hence is incompatible with NO_CPP option"
#	endif
#endif
/* standard FFT routines (FFTW3 of FFT_TEMPERTON) are required even when OpenCL is used, since they are used for Fourier
 * transform of the D-matrix
 */
#ifdef FFTW3
#	include <fftw3.h> // types.h or cmplx.h should be defined before (to match C99 complex type)
/* define level of planning for usual and Dmatrix (DM) FFT: FFTW_ESTIMATE (heuristics), FFTW_MEASURE (default),
 * FTW_PATIENT, or FFTW_EXHAUSTIVE
 */
#	define PLAN_FFTW FFTW_MEASURE
#	define PLAN_FFTW_DM FFTW_ESTIMATE
#	define ONLY_FOR_FFTW3 // this is used in function argument declarations
#else
#	define ONLY_FOR_FFTW3 ATT_UNUSED
#endif

#ifdef FFT_TEMPERTON
#	define ONLY_FOR_TEMPERTON // this is used in function argument declarations
#else
#	define ONLY_FOR_TEMPERTON ATT_UNUSED
#endif

// SEMI-GLOBAL VARIABLES

// defined and initialized in interaction.c
extern const int local_Nz_Rm;
// defined and initialized in timing.c
extern TIME_TYPE Timing_FFT_Init,Timing_Dm_Init;

// used in comm.c
double * restrict BT_buffer, * restrict BT_rbuffer; // buffers for BlockTranspose
// used in matvec.c; in OpenCL mode some of those are not used at all, others - only locally
doublecomplex * restrict Dmatrix; // holds FFT of the interaction matrix
doublecomplex * restrict Rmatrix; // holds FFT of the reflection matrix
#ifndef OPENCL
	// holds input vector (on expanded grid) to matvec, also used as storage space in iterative.c
doublecomplex * restrict Xmatrix;
doublecomplex * restrict slices; // used in inner cycle of matvec - holds 3 components (for fixed x)
doublecomplex * restrict slices_tr; // additional storage space for slices to accelerate transpose
doublecomplex * restrict slicesR,* restrict slicesR_tr; // same as above, but for reflected interaction
#endif
size_t DsizeY,DsizeZ,DsizeYZ; // size of the 'matrix' D
size_t RsizeY; // size of the 'matrix' R; in OpenCL mode it is used in oclmatvec.c
// used in oclmatvec.c
#ifdef OPENCL
/* clxslices is the number of slices required in MatVec,
 * so that more GPU memory is used and kernels can run over larger arrays
 */
size_t clxslices;
size_t local_gridX; // 'thickness' in x direction of a slice*
size_t slicesize; // total number of doublecomplex values per slice
#endif

// LOCAL VARIABLES

// D2 matrix and its two slices; used only temporary for InitDmatrix
static doublecomplex * restrict slice,* restrict slice_tr,* restrict D2matrix;
static doublecomplex * restrict R2matrix; // same for surface (slice and slice_tr are reused from Dmatrix)
static size_t D2sizeY; // size of the 'matrix' D2 (x-size is gridX), Z size is not used
static size_t R2sizeY; // size of the 'matrix' R2 (x- and z-sizes are corresponding grids)
static size_t lz_Dm,lz_Rm; // local sizes along z for D(2) and R(2) matrices
// the following two lines are defined in InitDmatrix but used in InitRmatrix, they are analogous to Dm values
static size_t Rsize,R2sizeTot; // sizes of R and R2 matrices
static int jstartR;            // starting index for y
static bool weird_nprocs;      // whether weird number of processors is used

#ifdef OPENCL
// clFFT plans
#	ifdef CLFFT_AMD
static clfftPlanHandle clplanX,clplanY,clplanZ;
static size_t clfftBufSize=0;
#	elif defined(CLFFT_APPLE)
static clFFT_Plan clplanX,clplanY,clplanZ;
#	endif
#endif
#ifdef FFTW3
// FFTW3 plans: f - FFT_FORWARD; b - FFT_BACKWARD
static fftw_plan planXf_Dm,planYf_slice,planZf_slice,planXf_Rm;
#	ifndef OPENCL // these plans are used only if OpenCL is not used
static fftw_plan planXf,planXb,planYf,planYb,planZf,planZb,planYRf,planZRf; // last two for reflected interaction
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
void cfft99_(double * restrict data,double * restrict _work,const double * restrict trigs,const int * restrict ifax,
	const int *inc,const int *jump,const int *nn,const int *lot,const int *isign);
#endif

//======================================================================================================================

static inline size_t IndexDmatrix(const size_t x,size_t y,size_t z)
// index D matrix to store final result (symmetric with respect to center for y and z)
{
	if (y>=DsizeY) y=gridY-y;
	if (z>=DsizeZ) z=gridZ-z;
	return(NDCOMP*((x*DsizeZ+z)*DsizeY+y));
}

//======================================================================================================================

static inline size_t IndexGarbledD(const size_t x,int y,int z)
// index D2 matrix after BlockTranspose (periodic over y and z)
{
	if (y<0) y+=gridY;
	if (z<0) z+=gridZ;
#ifdef PARALLEL
	return(((z%lz_Dm)*D2sizeY+y)*gridX+(z/lz_Dm)*local_Nx+x%local_Nx);
#else
	return((z*D2sizeY+y)*gridX+x);
#endif
}

//======================================================================================================================

static inline size_t IndexSliceD2matrix(int y,int z)
// index slice of D2 matrix (periodic over y and z)
{
	if (y<0) y+=gridY;
	if (z<0) z+=gridZ;
	return(y*gridZ+z);
}

//======================================================================================================================

static inline size_t Index2matrix(int x,int y,const int z,const int sizeY)
// index D2 or R2 matrix to store calculated elements (periodic over x and y), z should already be shifted
{
	if (x<0) x+=gridX;
	if (y<0) y+=gridY;
	return((z*sizeY+y)*gridX+x);
}

//======================================================================================================================

static inline size_t IndexSlice_zy(const size_t y,const size_t z)
// index transposed slice of D2 (or R2) matrix
{
	return (z*gridY+y);
}

//======================================================================================================================

static inline size_t IndexRmatrix(const size_t x,size_t y,const size_t z)
// index R matrix to store final result (symmetric with respect to center for y)
{
	if (y>=RsizeY) y=gridY-y;
	return(NDCOMP*((x*gridZ+z)*RsizeY+y));
}

//======================================================================================================================

static inline size_t IndexGarbledR(const size_t x,int y,const int z)
// index R2 matrix after BlockTranspose (periodic over y)
{
	if (y<0) y+=gridY;
#ifdef PARALLEL
	return(((z%lz_Rm)*R2sizeY+y)*gridX+(z/lz_Rm)*local_Nx+x%local_Nx);
#else
	return((z*R2sizeY+y)*gridX+x);
#endif
}

//======================================================================================================================

static inline size_t IndexSliceR2matrix(int y,const int z)
// index slice of R2 matrix (periodic over y)
{
	if (y<0) y+=gridY;
	return(y*gridZ+z);
}

//======================================================================================================================

static void transpose(const doublecomplex * restrict data,doublecomplex * restrict trans,const size_t Y,const size_t Z)
// optimized routine to transpose complex matrix with dimensions YxZ: data -> trans
{
	size_t y,z,y1,y2,z1,z2,i,j,y0,z0;
	doublecomplex *t1,*t2,*t3,*t4;
	const doublecomplex *w1,*w2,*w3;
	const size_t blockTr=64; // block size
	/* Intel compiler 11.1 seems to produce broken code for this function whenever blockTr>1, at least when the function
	 * is called in row of three from TransposeYZ(). So maybe the bug is due to incorrect inlining. This bug appears
	 * only for -O3 compilation, but not for -O2. Moreover, the bug is not present for icc 13.0.
	 * We have no idea what can be done on our side.
	 */

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
				for (z=0;z<z0;z++) *(t4+=Y)=w3[z];
				w3+=Z;
			}
			w2+=blockTr;
			t2+=blockTr*Y;
		}
		w1+=blockTr*Z;
		t1+=blockTr;
	}
}

//======================================================================================================================

void TransposeYZ(const int direction)
/* optimized routine to transpose y and z; forward: slices->slices_tr; backward: slices_tr->slices; direction can be
 * made boolean but this contradicts with existing definitions of FFT_FORWARD and FFT_BACKWARD, which themselves are
 * determined by FFT routines invocation format
 */
{
#ifdef OPENCL
	const size_t blocksize=16; //this corresponds to BLOCK_DIM in oclkernels.cl
	const size_t tblock[3]={blocksize,blocksize,1}; 
	size_t enqtglobalzy[3]={gridZ,gridY,3*local_gridX};
	size_t enqtglobalyz[3]={gridY,gridZ,3*local_gridX};

	//if the grid is not dividable by blocksize, extend it. Kernel takes care of borders
	size_t tgridZ = (gridZ%blocksize==0) ? gridZ : (gridZ/blocksize+1)*blocksize;
	size_t tgridY = (gridY%blocksize==0) ? gridY : (gridY/blocksize+1)*blocksize;
	enqtglobalzy[0]=tgridZ;
	enqtglobalzy[1]=tgridY;
	enqtglobalyz[0]=tgridY;
	enqtglobalyz[1]=tgridZ;

	if (direction==FFT_FORWARD) {
		CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,cltransposeof,3,NULL,enqtglobalzy,tblock,0,NULL,NULL));
		if (surface)
				CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,cltransposeofR,3,NULL,enqtglobalzy,tblock,0,NULL,NULL));
	}
	else CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,cltransposeob,3,NULL,enqtglobalyz,tblock,0,NULL,NULL));
#else
	size_t Xcomp,ind;

	if (direction==FFT_FORWARD) for (Xcomp=0;Xcomp<3;Xcomp++) {
		ind=Xcomp*gridYZ;
		transpose(slices+ind,slices_tr+ind,gridY,gridZ);
		if (surface) transpose(slicesR+ind,slicesR_tr+ind,gridY,gridZ);
	}
	else for (Xcomp=0;Xcomp<3;Xcomp++) { // direction==FFT_BACKWARD
		ind=Xcomp*gridYZ;
		transpose(slices_tr+ind,slices+ind,gridZ,gridY);
	}
#endif
}

//======================================================================================================================

void fftX(const int isign)
// FFT three components of (buf)Xmatrix(x) for all y,z; called from matvec
{
#ifdef OPENCL
#	ifdef CLFFT_AMD
	CL_CH_ERR(clfftEnqueueTransform(clplanX,(clfftDirection)isign,1,&command_queue,0,NULL,NULL,&bufXmatrix,NULL,
		NULL));
#	elif defined(CLFFT_APPLE)
	CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanX,(int)3*local_Nz*smallY,(clFFT_Direction)isign,bufXmatrix,
		bufXmatrix,0,NULL,NULL));
#	endif
#elif defined(FFTW3)
	if (isign==FFT_FORWARD) fftw_execute(planXf);
	else fftw_execute(planXb);
#elif defined(FFT_TEMPERTON)
	int nn=gridX,inc=1,jump=nn,lot=boxY;
	size_t z;
	/* Calls to Temperton FFT cause warnings for translation from doublecomplex to double pointers. However, such a cast
	 * is perfectly valid in C99. So we set pragmas to remove these warnings.
	 *
	 * !!! TODO: Another (ultimate) solution is to remove this routine altogether, since FFTW is perfect in all
	 * respects. This is also reasonable considering future switch to tgmath.h
	 */
	IGNORE_WARNING(-Wstrict-aliasing);
	for (z=0;z<3*local_Nz;z++) cfft99_((double *)(Xmatrix+z*gridX*smallY),work,trigsX,ifaxX,&inc,&jump,&nn,&lot,&isign);
	STOP_IGNORE;
#endif
}

//======================================================================================================================

void fftY(const int isign)
// FFT three components of slices_tr(y) for all z; called from matvec
{
#ifdef OPENCL
#	ifdef CLFFT_AMD
	CL_CH_ERR(clfftEnqueueTransform(clplanY,(clfftDirection)isign,1,&command_queue,0,NULL,NULL,&bufslices_tr,NULL,
		NULL));
	if (surface && isign==FFT_FORWARD)
		CL_CH_ERR(clfftEnqueueTransform(clplanY,(clfftDirection)isign,1,&command_queue,0,NULL,NULL,&bufslicesR_tr,
			NULL,NULL));
#	elif defined(CLFFT_APPLE)
	CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanY,(int)3*gridZ*local_gridX,(clFFT_Direction)isign,
		bufslices_tr,bufslices_tr,0,NULL,NULL));
	if (surface && isign==FFT_FORWARD)
		CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanY,(int)3*gridZ*local_gridX,(clFFT_Direction)isign,
			bufslicesR_tr,bufslicesR_tr,0,NULL,NULL));
#	endif
#elif defined(FFTW3)
	if (isign==FFT_FORWARD) {
		fftw_execute(planYf);
		if (surface) fftw_execute(planYRf);
	}
	else fftw_execute(planYb);
#elif defined(FFT_TEMPERTON)
	int nn=gridY,inc=1,jump=nn,lot=3*gridZ;

	IGNORE_WARNING(-Wstrict-aliasing);
	cfft99_((double *)(slices_tr),work,trigsY,ifaxY,&inc,&jump,&nn,&lot,&isign);
	// the same operation is applied to sliceR_tr, when required
	if (surface && isign==FFT_FORWARD) cfft99_((double *)(slicesR_tr),work,trigsY,ifaxY,&inc,&jump,&nn,&lot,&isign);
	STOP_IGNORE;
#endif
}

//======================================================================================================================

void fftZ(const int isign)
// FFT three components of slices(z) for all y; called from matvec
{
#ifdef OPENCL
#	ifdef CLFFT_AMD
	CL_CH_ERR(clfftEnqueueTransform(clplanZ,(clfftDirection)isign,1,&command_queue,0,NULL,NULL,&bufslices,NULL,
		NULL));
	if (surface && isign==FFT_FORWARD) // the same operation is applied to bufslicesR, but with inverse transform
		CL_CH_ERR(clfftEnqueueTransform(clplanZ,(clfftDirection)FFT_BACKWARD,1,&command_queue,0,NULL,NULL,
			&bufslicesR,NULL,NULL));
#	elif defined(CLFFT_APPLE)
	CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanZ,(int)3*gridY*local_gridX,(clFFT_Direction)isign,bufslices,
		bufslices,0,NULL,NULL));
	if (surface && isign==FFT_FORWARD) // the same operation is applied to bufslicesR, but with inverse transform
		CL_CH_ERR(clFFT_ExecuteInterleaved(command_queue,clplanZ,(int)3*gridY*local_gridX,(clFFT_Direction)FFT_BACKWARD,
			bufslicesR,bufslicesR,0,NULL,NULL));
#	endif
#elif defined(FFTW3)
	if (isign==FFT_FORWARD) {
		fftw_execute(planZf);
		if (surface) fftw_execute(planZRf);
	}
	else fftw_execute(planZb);
#elif defined(FFT_TEMPERTON)
	int nn=gridZ,inc=1,jump=nn,lot=boxY,Xcomp;

	IGNORE_WARNING(-Wstrict-aliasing);
	for (Xcomp=0;Xcomp<3;Xcomp++) cfft99_((double *)(slices+gridYZ*Xcomp),work,trigsZ,ifaxZ,&inc,&jump,&nn,&lot,&isign);
	if (surface && isign==FFT_FORWARD) { // the same operation is applied to slicesR, but with inverse transform
		const int invSign=FFT_BACKWARD;
		for (Xcomp=0;Xcomp<3;Xcomp++)
			cfft99_((double *)(slicesR+gridYZ*Xcomp),work,trigsZ,ifaxZ,&inc,&jump,&nn,&lot,&invSign);
	}
	STOP_IGNORE;
#endif
}

//======================================================================================================================

static void fftX_Dm(void)
// FFT(forward) D2matrix(x) for all y,z; used for Dmatrix calculation
{
#ifdef FFTW3
	fftw_execute(planXf_Dm);
#elif defined(FFT_TEMPERTON)
	int nn=gridX,inc=1,jump=nn,lot=D2sizeY,isign=FFT_FORWARD;
	size_t z;

	IGNORE_WARNING(-Wstrict-aliasing);
	for (z=0;z<lz_Dm;z++) cfft99_((double *)(D2matrix+z*gridX*D2sizeY),work,trigsX,ifaxX,&inc,&jump,&nn,&lot,&isign);
	STOP_IGNORE;
#endif
}

//======================================================================================================================

static void fftX_Rm(void)
// FFT(forward) D2matrix(x) for all y,z; used for Rmatrix calculation
{
#ifdef FFTW3
	fftw_execute(planXf_Rm);
#elif defined(FFT_TEMPERTON)
	int nn=gridX,inc=1,jump=nn,lot=R2sizeY,isign=FFT_FORWARD;
	size_t z;
	const size_t zlim=local_Nz_Rm; // can be smaller by 1 than lz_Rm

	IGNORE_WARNING(-Wstrict-aliasing);
	for (z=0;z<zlim;z++) cfft99_((double *)(R2matrix+z*gridX*R2sizeY),work,trigsX,ifaxX,&inc,&jump,&nn,&lot,&isign);
	STOP_IGNORE;
#endif
}

//======================================================================================================================

static void fftY_slice(void)
// FFT(forward) slice_tr(y) for all z; used for Dmatrix and Rmatrix calculation
{
#ifdef FFTW3
	fftw_execute(planYf_slice);
#elif defined(FFT_TEMPERTON)
	int nn=gridY,inc=1,jump=nn,lot=gridZ,isign=FFT_FORWARD;

	IGNORE_WARNING(-Wstrict-aliasing);
	cfft99_((double *)slice_tr,work,trigsY,ifaxY,&inc,&jump,&nn,&lot,&isign);
	STOP_IGNORE;
#endif
}

//======================================================================================================================

static void fftZ_slice(void)
// FFT(forward) slice(z) for all y; used for Dmatrix and Rmatrix calculation
{
#ifdef FFTW3
	fftw_execute(planZf_slice);
#elif defined(FFT_TEMPERTON)
	int nn=gridZ,inc=1,jump=nn,lot=gridY,isign=FFT_FORWARD;

	IGNORE_WARNING(-Wstrict-aliasing);
	cfft99_((double *)slice,work,trigsZ,ifaxZ,&inc,&jump,&nn,&lot,&isign);
	STOP_IGNORE;
#endif
}

//======================================================================================================================

void CheckNprocs(void)
// checks for consistency the specified number of processors; called in the beginning from InitComm
{
	int y=nprocs;

	// initialize weird_nprocs
	weird_nprocs=false;
	// remove simple prime divisors of y
	while (y%2==0) y/=2;
#ifdef OPENCL
	/* this is redundant, since OpenCL is not currently intended to run in parallel. In the future it should properly
	 * handle capabilities of AMD clFFT
	 */
	if (y!=1) PrintError("Specified number of processors (%d) is incompatible with clFFT, since the latter currently "
		"supports only FFTs with size 2^n. Please choose the number of processors to be of the same form.",nprocs);
#else
	while (y%3==0) y/=3;
	while (y%5==0) y/=5;
#	ifdef FFT_TEMPERTON
	if (y!=1) PrintError("Specified number of processors (%d) is weird (has prime divisors larger than 5). That is "
		"incompatible with Temperton FFT. Revise the number of processors (recommended) or recompile with FFTW 3 "
		"support.",nprocs);
#	elif defined(FFTW3)
	while (y%7==0) y/=7;
	// one multiplier of either 11 or 13 is allowed
	if (y%11==0) y/=11;
	else if (y%13==0) y/=13;
	if (y!=1) {
		LogWarning(EC_WARN,ONE_POS,"Specified number of processors (%d) is weird (has prime divisors larger than 13 or "
			"more than one divisor of either 11 or 13). FFTW3 will work less efficiently. It is strongly recommended "
			"to revise the number of processors.",nprocs);
		weird_nprocs=true;
	}
#	endif
#endif
}

//======================================================================================================================

int fftFit(int x,int divis)
/* find the first number >=x divisible by 2 only (Apple clFFT) or 2,3,5 only (Temperton FFT or clAMDFFT) or also
 * allowing 7 and one of 11 or 13 (FFTW3), and also divisible by 2 and divis. If weird_nprocs is used, only the latter
 * condition is required.
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

//======================================================================================================================

static void fftInitBeforeD(void)
// initialize fft before initialization of Dmatrix
{
#ifdef FFTW3
	int grXint=gridX,grYint=gridY,grZint=gridZ; // this is needed to provide 'int *' to grids

	D("FFTW library version: %s\n     compiler: %s\n     codelet optimizations: %s",fftw_version,fftw_cc,
		fftw_codelet_optim);
	planYf_slice=fftw_plan_many_dft(1,&grYint,gridZ,slice_tr,NULL,1,gridY,slice_tr,NULL,1,gridY,FFT_FORWARD,
		PLAN_FFTW_DM);
	planZf_slice=fftw_plan_many_dft(1,&grZint,gridY,slice,NULL,1,gridZ,slice,NULL,1,gridZ,FFT_FORWARD,PLAN_FFTW_DM);
	planXf_Dm=fftw_plan_many_dft(1,&grXint,lz_Dm*D2sizeY,D2matrix,NULL,1,gridX,D2matrix,NULL,1,gridX,FFT_FORWARD,
		PLAN_FFTW_DM);
	// very similar to Dm, but local_Nz_Rm can be smaller by 1 than lz_Rm
	if (surface) planXf_Rm=fftw_plan_many_dft(1,&grXint,local_Nz_Rm*R2sizeY,R2matrix,NULL,1,gridX,R2matrix,NULL,1,gridX,
		FFT_FORWARD,PLAN_FFTW_DM);
#elif defined(FFT_TEMPERTON)
	int nn;
	size_t size;

	// allocate memory
	MALLOC_VECTOR(trigsX,double,2*gridX,ALL);
	MALLOC_VECTOR(trigsY,double,2*gridY,ALL);
	MALLOC_VECTOR(trigsZ,double,2*gridZ,ALL);
	size=MAX(gridX*D2sizeY,3*gridYZ);
	if (surface) size=MAX(size,gridX*R2sizeY);
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

//======================================================================================================================

static void fftInitAfterD(void)
/* second part of fft initialization
 * completely separate code is used for OpenCL and FFTW3, because even precise-timing output is significantly different.
 * In particular, FFTW3 uses separate plans for forward and backward, while clFFT (by Apple or AMD) uses one plan for
 * both directions.
 *
 * clFft access the OpenCL buffers directly, so they are not anyhow affected by the definition of complex numbers in the
 * main part of the code (although, it is consistent with it)
 */
{
#ifdef OPENCL
#	ifdef CLFFT_APPLE
	cl_int err; // error code
#	else
	size_t bufsize;
#	endif
#	ifdef PRECISE_TIMING
	SYSTEM_TIME tvp[4];
#	endif

	if (IFROOT) printf("Initializing clFFT\n");
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp);
#	endif
#	ifdef CLFFT_AMD
	CL_CH_ERR(clfftSetup(NULL)); // first initialize clfft
#	ifdef DEBUGFULL
	cl_uint major,minor,patch;
	CL_CH_ERR(clfftGetVersion(&major,&minor,&patch));
	D("clFFT library version - %u.%u.%u",major,minor,patch);
#	endif
	/* Here and further we explicitly set all plan parameters for clFFT, even those that are equal to the default
	 * values (as recommended in clFFT manual)
	 */
	/* Unfortunately, clFFT (and Apple clFFT as well) currently supports only simple regular batches of transforms
	 * (similar to fftw_plan_many_dft) but not fully flexible configurations, like offered by fftw_plan_guru_dft.
	 * Basically the problem is due to lack of multi-dimensional (non-tightly packed) batches. So to make X transform as
	 * a single plan we have have to cycle over the whole smallY instead of (possibly smaller) boxY. This incurs a small
	 * performance hit for "non-standard" values of boxY, but should be overall faster than making an explicit loop over
	 * smaller kernels (like is now done with Temperton FFT).
	 */
	CL_CH_ERR(clfftCreateDefaultPlan(&clplanX,context,CLFFT_1D,&gridX));
	CL_CH_ERR(clfftSetPlanBatchSize(clplanX,3*local_Nz*smallY));
	CL_CH_ERR(clfftSetPlanPrecision(clplanX,PRECISION_CLFFT));
	CL_CH_ERR(clfftSetResultLocation(clplanX,CLFFT_INPLACE));
	CL_CH_ERR(clfftSetLayout(clplanX,CLFFT_COMPLEX_INTERLEAVED,CLFFT_COMPLEX_INTERLEAVED));
	CL_CH_ERR(clfftSetPlanScale(clplanX,FFT_FORWARD,1));
	CL_CH_ERR(clfftSetPlanScale(clplanX,FFT_BACKWARD,1)); // override the default (1/N) scale for backward direction
	CL_CH_ERR(clfftBakePlan(clplanX,1,&command_queue,NULL,NULL));
	CL_CH_ERR(clfftGetTmpBufSize(clplanX,&bufsize));
	clfftBufSize+=bufsize;
#	elif defined(CLFFT_APPLE)
	clFFT_Dim3 xdimen;
	xdimen.x=(unsigned int)gridX;
	xdimen.y=1;
	xdimen.z=1;
	clplanX=clFFT_CreatePlan(context,xdimen,clFFT_2D,clFFT_InterleavedComplexFormat,&err);
	CL_CH_ERR(err);
#	endif
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+1);
#	endif
#	ifdef CLFFT_AMD
	CL_CH_ERR(clfftCreateDefaultPlan(&clplanY,context,CLFFT_1D,&gridY));
	CL_CH_ERR(clfftSetPlanBatchSize(clplanY,3*gridZ*local_gridX));
	CL_CH_ERR(clfftSetPlanPrecision(clplanY,PRECISION_CLFFT));
	CL_CH_ERR(clfftSetResultLocation(clplanY,CLFFT_INPLACE));
	CL_CH_ERR(clfftSetLayout(clplanY,CLFFT_COMPLEX_INTERLEAVED,CLFFT_COMPLEX_INTERLEAVED));
	CL_CH_ERR(clfftSetPlanScale(clplanY,FFT_FORWARD,1));
	CL_CH_ERR(clfftSetPlanScale(clplanY,FFT_BACKWARD,1)); // override the default (1/N) scale for backward direction
	CL_CH_ERR(clfftBakePlan(clplanY,1,&command_queue,NULL,NULL));
	CL_CH_ERR(clfftGetTmpBufSize(clplanY,&bufsize));
	clfftBufSize+=bufsize;
#	elif defined(CLFFT_APPLE)
	clFFT_Dim3 ydimen;
	ydimen.x=(unsigned int)gridY;
	ydimen.y=1;
	ydimen.z=1;
	clplanY=clFFT_CreatePlan(context,ydimen,clFFT_1D,clFFT_InterleavedComplexFormat,&err);
	CL_CH_ERR(err);
#	endif
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+2);
#	endif
#	ifdef CLFFT_AMD
	/* Here the issue is similar to clplanX described above. However, we are using full gridY instead of boxY, which
	 * incurs at least a-factor-of-two performance hit. To solve this problem one need to execute separate plans for 3
	 * components of vectors. Unfortunately, this cannot be simply done using 3-element loop (like in Temperton FFT),
	 * because a part of cl_mem object can't be addressed independently (as can be done with simple C arrays). The only
	 * way to address this issue is to either create three separate cl_mem objects or to change the indexing of levels
	 * inside the array, so that 3 components are stored together.
	 */
	CL_CH_ERR(clfftCreateDefaultPlan(&clplanZ,context,CLFFT_1D,&gridZ));
	/* TODO: last slices can be very slightly thinner than the previous ones, but since the batchsize is part of the
	 * plan, another plan would be needed to address this. However, we ignore this currently and assume that every
	 * slice has a thickness of local_gridX.
	 * This issue also applies to clplanY.
	 */
	CL_CH_ERR(clfftSetPlanBatchSize(clplanZ,3*gridY*local_gridX));
	CL_CH_ERR(clfftSetPlanPrecision(clplanZ,PRECISION_CLFFT));
	CL_CH_ERR(clfftSetResultLocation(clplanZ,CLFFT_INPLACE));
	CL_CH_ERR(clfftSetLayout(clplanZ,CLFFT_COMPLEX_INTERLEAVED,CLFFT_COMPLEX_INTERLEAVED));
	CL_CH_ERR(clfftSetPlanScale(clplanZ,FFT_FORWARD,1));
	CL_CH_ERR(clfftSetPlanScale(clplanZ,FFT_BACKWARD,1)); // override the default (1/N) scale for backward direction
	CL_CH_ERR(clfftBakePlan(clplanZ,1,&command_queue,NULL,NULL));
	CL_CH_ERR(clfftGetTmpBufSize(clplanZ,&bufsize));
	clfftBufSize+=bufsize;
	/* In most cases clfftBufSize is zero, except some weird grid sizes like 2x2x60000. Still, we rigorously account
	 * for this memory. However, we do not update oclMemMaxObj, since even single plan is not guaranteed to allocate a
	 * single object. So we assume that clFftAmd will either handle maximum object size itself or produce a meaningful
	 * error.
	 */
	oclMem+=clfftBufSize;
	MAXIMIZE(oclMemPeak,oclMem);
#	elif defined(CLFFT_APPLE)
	clFFT_Dim3 zdimen;
	zdimen.x=(unsigned int)gridZ;
	zdimen.y=1;
	zdimen.z=1;
	clplanZ=clFFT_CreatePlan(context,zdimen,clFFT_1D,clFFT_InterleavedComplexFormat,&err);
	CL_CH_ERR(err);
#	endif
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+3);
	// print precise timing of FFT planning
	if (IFROOT) PrintBoth(logfile,
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"         clFFT planning       \n"
		"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
		"X = "FFORMPT"  Total = "FFORMPT"\n"
		"Y = "FFORMPT"\n"
		"Z = "FFORMPT"\n\n",
		DiffSystemTime(tvp,tvp+1),DiffSystemTime(tvp,tvp+3),DiffSystemTime(tvp+1,tvp+2),DiffSystemTime(tvp+2,tvp+3));
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
	GET_SYSTEM_TIME(tvp);
#	endif
	lot=3*gridZ;
	planYf=fftw_plan_many_dft(1,&grYint,lot,slices_tr,NULL,1,gridY,slices_tr,NULL,1,gridY,FFT_FORWARD,PLAN_FFTW);
	if (surface) // same operation, but applied to slicesR_tr
		planYRf=fftw_plan_many_dft(1,&grYint,lot,slicesR_tr,NULL,1,gridY,slicesR_tr,NULL,1,gridY,FFT_FORWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+1);
#	endif
	planYb=fftw_plan_many_dft(1,&grYint,lot,slices_tr,NULL,1,gridY,slices_tr,NULL,1,gridY,FFT_BACKWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+2);
#	endif
	dims.n=gridZ;
	dims.is=dims.os=1;
	howmany_dims[0].n=3;
	howmany_dims[0].is=howmany_dims[0].os=gridZ*gridY;
	howmany_dims[1].n=boxY;
	howmany_dims[1].is=howmany_dims[1].os=gridZ;
	planZf=fftw_plan_guru_dft(1,&dims,2,howmany_dims,slices,slices,FFT_FORWARD,PLAN_FFTW);
	// same operation but for slicesR and inverse transform (since correlation is computed instead of convolution)
	if (surface) planZRf=fftw_plan_guru_dft(1,&dims,2,howmany_dims,slicesR,slicesR,FFT_BACKWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+3);
#	endif
	planZb=fftw_plan_guru_dft(1,&dims,2,howmany_dims,slices,slices,FFT_BACKWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+4);
#	endif
	dims.n=gridX;
	dims.is=dims.os=1;
	howmany_dims[0].n=3*local_Nz;
	howmany_dims[0].is=howmany_dims[0].os=smallY*gridX;
	howmany_dims[1].n=boxY;
	howmany_dims[1].is=howmany_dims[1].os=gridX;
	planXf=fftw_plan_guru_dft(1,&dims,2,howmany_dims,Xmatrix,Xmatrix,FFT_FORWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+5);
#	endif
	planXb=fftw_plan_guru_dft(1,&dims,2,howmany_dims,Xmatrix,Xmatrix,FFT_BACKWARD,PLAN_FFTW);
#	ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+6);
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
		DiffSystemTime(tvp,tvp+1),DiffSystemTime(tvp,tvp+6),DiffSystemTime(tvp+1,tvp+2),DiffSystemTime(tvp+2,tvp+3),
		DiffSystemTime(tvp+3,tvp+4),DiffSystemTime(tvp+4,tvp+5),DiffSystemTime(tvp+5,tvp+6));
#	endif
#endif
#ifdef FFTW3
	// destroy old (D,R-matrix) plans; also in OpenCL mode
	fftw_destroy_plan(planXf_Dm);
	fftw_destroy_plan(planYf_slice);
	fftw_destroy_plan(planZf_slice);
	if (surface) fftw_destroy_plan(planXf_Rm);
#	ifdef OPENCL // in this case, FFTW ends here
	fftw_cleanup();
#	endif
#endif
}

//======================================================================================================================

static void InitRmatrix(const double invNgrid)
/* Initializes the matrix R. R[i][j][k]=GR[i1-i2][j1-j2][k1+k2]. Actually R=-FFT(GR)/Ngrid. Then -GR.x=invFFT(R*FFT(x))
 * for practical implementation of FFT such that invFFT(FFT(x))=Ngrid*x. GR is exactly reflected Green's tensor. The
 * routine is very similar to the corresponding part of InitDmatrix. Moreover, some initialization is delegated to the
 * latter function.
 */
{
	int i,j,k,Rcomp;
	size_t x,y,z,indexfrom,indexto,ind,index;

	// allocate memory for Rmatrix (R2matrix is allocated earlier in InitDmatrix)
	MALLOC_VECTOR(Rmatrix,complex,Rsize,ALL);
#ifdef PARALLEL
	// allocate buffer for BlockTranspose_DRm
	size_t bufsize = 2*lz_Rm*R2sizeY*local_Nx;
	MALLOC_VECTOR(BT_buffer,double,bufsize,ALL);
	MALLOC_VECTOR(BT_rbuffer,double,bufsize,ALL);
#endif
	if (IFROOT) printf("Calculating reflected Green's function (Rmatrix)\n");
	/* Interaction matrix values are calculated all at once for performance reasons. They are stored in Rmatrix with
	 * indexing corresponding to R2matrix (to facilitate copying) but NDCOMP elements instead of one. Afterwards they
	 * are replaced by Fourier transforms (with different indexing) component-wise (in cycle over NDCOMP)
	 */
	/* fill Rmatrix with 0, this if to fill the possible gap between e.g. boxY and gridY/2; (and for R=0) probably
	 * faster than using a lot of conditionals
	 */
	for (ind=0;ind<Rsize;ind++) Rmatrix[ind]=0;
	// fill Rmatrix with values of reflected Green's tensor
	for(k=0;k<local_Nz_Rm;k++) for (j=jstartR;j<boxY;j++) for (i=1-boxX;i<boxX;i++) {
			index=NDCOMP*Index2matrix(i,j,k,R2sizeY);
			(*ReflTerm_int)(i,j,k,Rmatrix+index);
	} // end of i,j,k loop
	if (IFROOT) printf("Fourier transform of Rmatrix");
	for(Rcomp=0;Rcomp<NDCOMP;Rcomp++) { // main cycle over components of Rmatrix
		// fill R2matrix with precomputed values from Rmatrix
		for (ind=0;ind<R2sizeTot;ind++) R2matrix[ind]=Rmatrix[NDCOMP*ind+Rcomp];
		fftX_Rm(); // fftX R2matrix
		BlockTranspose_DRm(R2matrix,R2sizeY,lz_Rm);
		for(x=local_x0;x<local_x1;x++) {
			for (ind=0;ind<gridYZ;ind++) slice[ind]=0.0; // fill slice with 0.0
			for(j=jstartR;j<boxY;j++) for(k=0;k<2*boxZ-1;k++) {
				indexfrom=IndexGarbledR(x,j,k);
				indexto=IndexSliceR2matrix(j,k);
				slice[indexto]=R2matrix[indexfrom];
			}
			/* here a specific symmetry is used, that elements of R depend on direction y/|rho| either as even order
			 * (0 or 2) or as odd (1) - the latter are elements 1 and 4 (see GetSomIntegral in interaction.c)
			 */
			if (reduced_FFT) for(j=1;j<boxY;j++) for(k=0;k<2*boxZ-1;k++) {
				// mirror along y
				indexfrom=IndexSliceR2matrix(j,k);
				indexto=IndexSliceR2matrix(-j,k);
				if (Rcomp==1 || Rcomp==4) slice[indexto]=-slice[indexfrom];
				else slice[indexto]=slice[indexfrom];
			}
			fftZ_slice(); // fftZ slice
			transpose(slice,slice_tr,gridY,gridZ);
			fftY_slice(); // fftY slice_tr
			for(z=0;z<gridZ;z++) for(y=0;y<RsizeY;y++) {
				indexto=IndexRmatrix(x-local_x0,y,z)+Rcomp;
				indexfrom=IndexSlice_zy(y,z);
				Rmatrix[indexto]=-invNgrid*slice_tr[indexfrom];
			}
		} // end slice X
		if (IFROOT) printf(".");
	} // end of Rcomp
	if (IFROOT) printf("\n");
#ifdef PARALLEL
	// deallocate buffers for BlockTranspose_DRm
	Free_general(BT_buffer);
	Free_general(BT_rbuffer);
#endif
#ifdef OPENCL
	// Setting kernel arguments which are always the same
	// for arith3_surface
	CL_CH_ERR(clSetKernelArg(clarith3_surface,0,sizeof(cl_mem),&bufslices_tr));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,1,sizeof(cl_mem),&bufDmatrix));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,2,sizeof(size_t),&smallY));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,3,sizeof(size_t),&smallZ));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,4,sizeof(size_t),&gridX));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,5,sizeof(size_t),&DsizeY));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,6,sizeof(size_t),&DsizeZ));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,11,sizeof(cl_mem),&bufslicesR_tr));
	CL_CH_ERR(clSetKernelArg(clarith3_surface,12,sizeof(cl_mem),&bufRmatrix));
	// for transpose forward (backward are not needed for surface)
	CL_CH_ERR(clSetKernelArg(cltransposeofR,0,sizeof(cl_mem),&bufslicesR));
	CL_CH_ERR(clSetKernelArg(cltransposeofR,1,sizeof(cl_mem),&bufslicesR_tr));
	CL_CH_ERR(clSetKernelArg(cltransposeofR,2,sizeof(size_t),&gridZ));
	CL_CH_ERR(clSetKernelArg(cltransposeofR,3,sizeof(size_t),&gridY));
	CL_CH_ERR(clSetKernelArg(cltransposeofR,4,17*16*sizeof(doublecomplex),NULL));
	// copy Rmatrix to OpenCL buffer, blocking to ensure completion before function end
	CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufRmatrix,CL_TRUE,0,Rsize*sizeof(*Rmatrix),Rmatrix,0,NULL,NULL));
	Free_cVector(Rmatrix);
#endif
}


//======================================================================================================================

void InitDmatrix(void)
/* Initializes the matrix D. D[i][j][k]=A[i1-i2][j1-j2][k1-k2]. Actually D=-FFT(G)/Ngrid. Then -G.x=invFFT(D*FFT(x)) for
 * practical implementation of FFT such that invFFT(FFT(x))=Ngrid*x. G is exactly Green's tensor. The routine is called
 * only once, so does not need to be very fast, however we tried to optimize it.
 */
{
	int i,j,k,kcor,Dcomp;
	size_t x,y,z,indexfrom,indexto,ind,index,Dsize,D2sizeTot;
	double invNgrid;
	int nnn; // multiplier used for reduced_FFT or not reduced; 1 or 2
	int jstart,kstart;
	TIME_TYPE start,time1;
#ifdef PRECISE_TIMING
	// precise timing of the Dmatrix computation
	SYSTEM_TIME tvp[15];
	SYSTEM_TIME Timing_fftX,Timing_fftY,Timing_fftZ,Timing_Gcalc,Timing_ar1,Timing_ar2,Timing_ar3,Timing_BT,Timing_TYZ,
		Timing_beg,Timing_InitMV;
	double t_fftX,t_fftY,t_fftZ,t_ar1,t_ar2,t_ar3,t_TYZ,t_beg,t_Gcalc,t_Arithm,t_FFT,t_BT,t_InitMV,t_Rm,t_Tot;

	// This should be the first occurrence of PRECISE_TIMING in the program
	SetTimerFreq();

	t_Rm=0; // redundant initialization to remove warnings
	InitTime(&Timing_fftX);
	InitTime(&Timing_fftY);
	InitTime(&Timing_fftZ);
	InitTime(&Timing_ar1);
	InitTime(&Timing_ar2);
	InitTime(&Timing_ar3);
	InitTime(&Timing_BT);
	InitTime(&Timing_TYZ);
	InitTime(&Timing_InitMV);
	GET_SYSTEM_TIME(tvp);
#endif
	start=GET_TIME();

	// initialize sizes of D and D2 matrices
	if (reduced_FFT) {
		D2sizeY=gridY/2;
		// D2sizeZ=gridZ/2;
		DsizeY=gridY/2+1;
		DsizeZ=gridZ/2+1;
		nnn=1;
		jstart=0;
		kstart=0;
	}
	else {
		D2sizeY=DsizeY=gridY;
		DsizeZ=gridZ; // also =D2sizeZ
		nnn=2;
		jstart=1-boxY;
		kstart=1-boxZ;
	}
	// auxiliary parameters
	lz_Dm=nnn*local_Nz;
	DsizeYZ=DsizeY*DsizeZ;
	invNgrid=1.0/(gridX*((double)gridYZ));
	local_Nsmall=(gridX/2)*(gridYZ/(2*nprocs)); // size of X vector (for 1 component)
	// potentially this may cause unnecessary error during prognosis, but makes code cleaner
	Dsize=MultOverflow(NDCOMP*local_Nx,DsizeYZ,ONE_POS_FUNC);
	D2sizeTot=nnn*local_Nz*D2sizeY*gridX; // this should be approximately equal to Dsize/NDCOMP
	if (IFROOT) fprintf(logfile,"The FFT grid is: %zux%zux%zu\n",gridX,gridY,gridZ);

	// part of the code for InitRmatrix is here to be compatible with prognosis and FFT init
	if (surface) {
		/* We keep option to turn off reduced_FFT for Rmatrix (at least for tests). However, its usefulness for some
		 * weird Green's tensors (non-symmetric) is limited, because certain symmetry along the z-axis is still assumed
		 * - that is R({i1,j1,k1},{i2,j2,k2})=R({i1,j1,k2},{i2,j2,k1}). Without this limitation the whole FFT part will
		 * be broken (or need significant revision).
		 *
		 * Moreover the savings by using reduced_FFT is only a factor of two for Rmatrix (in contrast to 4 for Dmatrix)
		 */
		if (reduced_FFT) {
			R2sizeY=gridY/2;
			RsizeY=gridY/2+1;
			jstartR=0;
		}
		else {
			R2sizeY=RsizeY=gridY;
			jstartR=1-boxY;
		}
		lz_Rm=2*local_Nz;
		// potentially this may cause unnecessary error during prognosis, but makes code cleaner
		Rsize=MultOverflow(NDCOMP*local_Nx,RsizeY*gridZ,ONE_POS_FUNC);
		R2sizeTot=lz_Rm*R2sizeY*gridX; // this should be approximately equal to Rsize/NDCOMP
	}
#ifdef OPENCL // perform setting up of buffers and kernels
	/* The order of allocation is such that to have all bufslices* at the end to spent whatever memory is still
	 * available on the GPU on 'thickness' of 3D slices. This enables to have a 3D FFT and longer kernel runs but less
	 * number of kernel calls.
	 */
	// create all Buffers needed on Device in MatVec; When prognosis, the following code just counts required memory
	CREATE_CL_BUFFER(bufXmatrix,CL_MEM_READ_WRITE,local_Nsmall*3*sizeof(doublecomplex),NULL);
#	ifdef OCL_BLAS
	if (IterMethod==IT_BICG_CS) { // currently, used only in one iterative solver
		// Most clBLAS functions require scratch buffer of size N, but Dznrm2 - 2N
		CREATE_CL_BUFFER(buftmp,CL_MEM_READ_WRITE,local_nRows*2*sizeof(doublecomplex),NULL);
		CREATE_CL_BUFFER(bufxvec,CL_MEM_READ_WRITE,local_nRows*sizeof(doublecomplex),NULL);
		CREATE_CL_BUFFER(bufrvec,CL_MEM_READ_WRITE,local_nRows*sizeof(doublecomplex),NULL);
	}
#	endif
	CREATE_CL_BUFFER(bufargvec,CL_MEM_READ_WRITE,local_nRows*sizeof(doublecomplex),NULL);
	CREATE_CL_BUFFER(bufresultvec,CL_MEM_READ_WRITE,local_nRows*sizeof(doublecomplex),NULL);
	if (ipr_required) { // for inner product (only if it will be used afterwards)
		memory+=local_nvoid_Ndip*sizeof(double);
		CREATE_CL_BUFFER(bufinproduct,CL_MEM_READ_WRITE,local_nvoid_Ndip*sizeof(double),NULL);
		if (!prognosis) {
			MALLOC_VECTOR(inprodhlp,double,local_nvoid_Ndip,ALL);
			CL_CH_ERR(clSetKernelArg(clinprod,0,sizeof(cl_mem),&bufinproduct));
			CL_CH_ERR(clSetKernelArg(clinprod,1,sizeof(cl_mem),&bufresultvec));
		}
	}
	/* The following are constant device buffers which are initialized with host data. They are all created here (to be
	 * compatible with prognosis), but some are initialized (filled with data) later.
	 */
	CREATE_CL_BUFFER(bufcc_sqrt,CL_MEM_READ_ONLY,sizeof(cc_sqrt),NULL);
	CREATE_CL_BUFFER(bufDmatrix,CL_MEM_READ_ONLY,Dsize*sizeof(*Dmatrix),NULL);
	if (surface) CREATE_CL_BUFFER(bufRmatrix,CL_MEM_READ_ONLY,Rsize*sizeof(*Rmatrix),NULL);
	CREATE_CL_BUFFER(bufmaterial,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,local_nvoid_Ndip*sizeof(*material),material);
	CREATE_CL_BUFFER(bufposition,CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,local_nRows*sizeof(*position),position);

	/* In the following bufslices* are allocated, based on available GPU memory.
	 * The estimation doesn't account for memory, which may be allocated further in fftInitBeforeD(), but it is zero in
	 * most cases. There are also all kind of other possible problems (e.g. with Nvidia GPUs, which gives access to the
	 * whole GPU memory), when occupying most of the GPU memory results in failures ("Out of resources") either now or
	 * during further iterations. There is little we can do with these problems, apart from using '-opt mem').
	 */
	const size_t memReserve = 100*MBYTE; // memory reserved for all other GPU needs (including desktop,etc.)
	if (save_memory) { // fall back to one-layer-at-a-time implementation
		local_gridX=1;
		clxslices=gridX;
		D("Using 1-layer x-slices (memory optimization)");
	}
	else if (prognosis) { // maximum memory, should not necessarily fit in the current GPU
		local_gridX=gridX;
		clxslices=1;
		D("Using the largest x-slices (prognosis mode)");
	}
	else if (oclMemDev<oclMem+memReserve) {
		/* insufficient GPU memory. We use minimum local_gridX, but even this is expected to cause error. But we leave
		 * the error processing to the allocation routines or OpenCL driver. It may also happen that no error will be
		 * produced, instead the execution will be slow due to memory swapping.
		 */
		local_gridX=1;
		clxslices=gridX;
		D("Using 1-layer x-slices due (insufficient memory, already occupied "FFORMM" MB)",oclMem/MBYTE);
	}
	else {
		const size_t slbufnum = surface ? 4 : 2; // number of buffers
		const size_t memAvail=oclMemDev-oclMem-memReserve; // this is guaranteed to be non-negative
		const size_t memLayer=3*gridYZ*sizeof(doublecomplex); // memory for a single layer (minimum slice)
		// local_gridX is always at least 1 (possible errors of insufficient memory are ignored here)
		local_gridX=MIN(memAvail/(slbufnum*memLayer),oclMemDevObj/memLayer);
		if (local_gridX==0) local_gridX=1;
		if (local_gridX>32) local_gridX=32;
		if (gridX%local_gridX==0) clxslices=gridX/local_gridX; // automatic uniform division
		else {
			clxslices=(gridX/local_gridX)+1;
			local_gridX=DIV_CEILING(gridX,clxslices); // adjust local_gridX to be closer to uniform division
			// if gridX<=32; the above code will set local_gridX=gridX
		}

		D("Already occupied OpenCL memory: "FFORMM" MB,\n"
			"     available for x-slices: "FFORMM" MB (excluding "FFORMM" MB reserve),\n"
			"     required for the largest x-slices: "FFORMM" MB",
			oclMem/MBYTE,memAvail/MBYTE,memReserve/MBYTE,(memLayer/MBYTE)*gridX*slbufnum);
		if (local_gridX==gridX) { // braces {} are to remove warnings
			D("Using the largest x-slices (sufficient memory)");
		}
		else {
			D("gridX (%zu) is split into %zu slices of %zu layers each (low memory)",gridX,clxslices,local_gridX);
		}
	}
	slicesize=local_gridX*gridYZ*3;
	/* TODO: Here it is possible to recover from error in memory allocation by catching this exception, decreasing the
	 * local_gridX, and trying again. However, consistency of memory sizes is not well anyway. Even if we succeed to
	 * allocate memory, some artifacts may further appear, e.g. from desktop activity.
	 */
	CREATE_CL_BUFFER(bufslices,CL_MEM_READ_WRITE,slicesize*sizeof(doublecomplex),NULL);
	CREATE_CL_BUFFER(bufslices_tr,CL_MEM_READ_WRITE,slicesize*sizeof(doublecomplex),NULL);
	if (surface) {
		CREATE_CL_BUFFER(bufslicesR,CL_MEM_READ_WRITE,slicesize*sizeof(doublecomplex),NULL);
		CREATE_CL_BUFFER(bufslicesR_tr,CL_MEM_READ_WRITE,slicesize*sizeof(doublecomplex),NULL);
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
		CL_CH_ERR(clSetKernelArg(clarith3,2,sizeof(size_t),&smallY));
		CL_CH_ERR(clSetKernelArg(clarith3,3,sizeof(size_t),&smallZ));
		CL_CH_ERR(clSetKernelArg(clarith3,4,sizeof(size_t),&gridX));
		CL_CH_ERR(clSetKernelArg(clarith3,5,sizeof(size_t),&DsizeY));
		CL_CH_ERR(clSetKernelArg(clarith3,6,sizeof(size_t),&DsizeZ));
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
		// transpose kernels, first for transpose forward
		CL_CH_ERR(clSetKernelArg(cltransposeof,0,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(cltransposeof,1,sizeof(cl_mem),&bufslices_tr));
		CL_CH_ERR(clSetKernelArg(cltransposeof,2,sizeof(size_t),&gridZ));
		CL_CH_ERR(clSetKernelArg(cltransposeof,3,sizeof(size_t),&gridY));
		/* setting up local cache size as 17*16 elements; note: a block is only 16*16, but 1*16 stride is needed to
		 * avoid bank conflicts
		 */
		CL_CH_ERR(clSetKernelArg(cltransposeof,4,17*16*sizeof(doublecomplex),NULL));
		// transpose backward
		CL_CH_ERR(clSetKernelArg(cltransposeob,0,sizeof(cl_mem),&bufslices_tr));
		CL_CH_ERR(clSetKernelArg(cltransposeob,1,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clSetKernelArg(cltransposeob,2,sizeof(size_t),&gridY));
		CL_CH_ERR(clSetKernelArg(cltransposeob,3,sizeof(size_t),&gridZ));
		CL_CH_ERR(clSetKernelArg(cltransposeob,4,17*16*sizeof(doublecomplex),NULL));
	}
#endif
	// memory estimation and exit for prognosis
	MAXIMIZE(memPeak,memory);
	/* objects which are always allocated (at least temporarily): Dmatrix,D2matrix,slice,slice_tr
	 * for surface, the peak is either by D2matrix & R2matrix, or by R2matrix & Rmatrix (the latter is mostly probable)
	 */
	memPeak+=sizeof(doublecomplex)*((double)Dsize+2*gridYZ+(surface ? (MAX(Rsize,D2sizeTot)+R2sizeTot) : D2sizeTot));
#ifndef OPENCL
	/* allocated memory that is used further on (Dmatrix,Xmatrix,slices,slices_tr), not relevant for OpenCL version;
	 * we assume that it is always larger than memPeak above (so memPeak doesn't have to be adjusted). In particular,
	 * we ignore the memory, which is temporarily allocated for BlockTranspose buffers of Dm and Rm.
	 */
	double mem=sizeof(doublecomplex)*((double)Dsize+3*local_Nsmall+6*gridYZ);
	if (surface) mem+=sizeof(doublecomplex)*((double)Rsize+6*gridYZ); // for Rmatrix, slicesR, and slicesR_tr
#ifdef PARALLEL
	const size_t BTsize = 6*smallY*local_Nz*local_Nx; // in doubles
	mem+=2*BTsize*sizeof(double);
#endif
	// printout some information
	if (IFROOT) {
#ifdef PARALLEL
		PrintBoth(logfile,"Memory usage for MatVec matrices (per processor): "FFORMM" MB\n",mem/MBYTE);
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
	/* allocate memory for R2matrix components. In principle, this can be done after D2 matrix is freed. However, this
	 * way allows us to init all FFT routines (in particular, build FFTW plans) in one go. Moreover, this should not
	 * increase the peak memory, since Rmatrix is allocated further on (see above).
	 */
	if (surface) MALLOC_VECTOR(R2matrix,complex,R2sizeTot,ALL);
	// actually allocation of Xmatrix, slices, slices_tr is below after freeing of Dmatrix and its slice
#ifdef PARALLEL
	// allocate buffer for BlockTranspose_Dm
	size_t bufsize = 2*lz_Dm*D2sizeY*local_Nx;
	MALLOC_VECTOR(BT_buffer,double,bufsize,ALL);
	MALLOC_VECTOR(BT_rbuffer,double,bufsize,ALL);
#endif
	D("Initialize FFT (1st part)");
	fftInitBeforeD();
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+1);
	Elapsed(tvp,tvp+1,&Timing_beg); // it includes a lot of OpenCL stuff
#endif
	if (IFROOT) printf("Calculating Green's function (Dmatrix)\n");
	/* Interaction matrix values are calculated all at once for performance reasons. They are stored in Dmatrix with
	 * indexing corresponding to D2matrix (to facilitate copying) but NDCOMP elements instead of one. Afterwards they
	 * are replaced by Fourier transforms (with different indexing) component-wise (in cycle over NDCOMP)
	 */
	/* fill Dmatrix with 0, this if to fill the possible gap between e.g. boxY and gridY/2; (and for R=0) probably
	 * faster than using a lot of conditionals
	 */
	for (ind=0;ind<Dsize;ind++) Dmatrix[ind]=0;
	// fill Dmatrix with values of Green's tensor
	for(k=nnn*local_z0;k<nnn*local_z1;k++) {
		// correction of k is relevant only if reduced_FFT is not used
		if (k>(int)smallZ) kcor=k-gridZ;
		else kcor=k;
		for (j=jstart;j<boxY;j++) for (i=1-boxX;i<boxX;i++) {
			index=NDCOMP*Index2matrix(i,j,k-nnn*local_z0,D2sizeY);
			/* The test for zero distance is somewhat non-optimal. However, other alternatives are not perfect either:
			 * 1) complicate the loops to remove the zero element in the beginning (move tests to the upper level)
			 * 2) call the function with zero - it will produce NaN. Then set this element to zero after the loop.
			 */
			if (i!=0 || j!=0 || kcor!=0) (*InterTerm_int)(i,j,kcor,Dmatrix+index);
		}
	} // end of i,j,k loop
	if (IFROOT) printf("Fourier transform of Dmatrix");
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+11); // same as the last time-stamp in the following loop
	Elapsed(tvp+1,tvp+11,&Timing_Gcalc);
#endif
	for(Dcomp=0;Dcomp<NDCOMP;Dcomp++) { // main cycle over components of Dmatrix
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+2);
		ElapsedInc(tvp+11,tvp+2,&Timing_InitMV);
#endif
		// fill D2matrix with precomputed values from Dmatrix
		for (ind=0;ind<D2sizeTot;ind++) D2matrix[ind]=Dmatrix[NDCOMP*ind+Dcomp];
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+3);
		ElapsedInc(tvp+2,tvp+3,&Timing_ar1);
#endif
		fftX_Dm(); // fftX D2matrix
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+4);
		ElapsedInc(tvp+3,tvp+4,&Timing_fftX);
#endif
		BlockTranspose_DRm(D2matrix,D2sizeY,lz_Dm);
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+5);
		ElapsedInc(tvp+4,tvp+5,&Timing_BT);
#endif
		for(x=local_x0;x<local_x1;x++) {
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+6);
#endif
			for (ind=0;ind<gridYZ;ind++) slice[ind]=0.0; // fill slice with 0.0
			for(j=jstart;j<boxY;j++) for(k=kstart;k<boxZ;k++) {
				indexfrom=IndexGarbledD(x,j,k);
				indexto=IndexSliceD2matrix(j,k);
				slice[indexto]=D2matrix[indexfrom];
			}
			if (reduced_FFT) { // here a specific symmetry is used, that G is a combination of tensors I and RR/|R|^2
				for(j=1;j<boxY;j++) for(k=0;k<boxZ;k++) {
					// mirror along y
					indexfrom=IndexSliceD2matrix(j,k);
					indexto=IndexSliceD2matrix(-j,k);
					if (Dcomp==1 || Dcomp==4) slice[indexto]=-slice[indexfrom];
					else slice[indexto]=slice[indexfrom];
				}
				for(j=1-boxY;j<boxY;j++) for(k=1;k<boxZ;k++) {
					// mirror along z
					indexfrom=IndexSliceD2matrix(j,k);
					indexto=IndexSliceD2matrix(j,-k);
					if (Dcomp==2 || Dcomp==4) slice[indexto]=-slice[indexfrom];
					else slice[indexto]=slice[indexfrom];
				}
			}
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+7);
			ElapsedInc(tvp+6,tvp+7,&Timing_ar2);
#endif
			fftZ_slice(); // fftZ slice
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+8);
			ElapsedInc(tvp+7,tvp+8,&Timing_fftZ);
#endif
			transpose(slice,slice_tr,gridY,gridZ);
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+9);
			ElapsedInc(tvp+8,tvp+9,&Timing_TYZ);
#endif
			fftY_slice(); // fftY slice_tr
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+10);
			ElapsedInc(tvp+9,tvp+10,&Timing_fftY);
#endif
			for(z=0;z<DsizeZ;z++) for(y=0;y<DsizeY;y++) {
				indexto=IndexDmatrix(x-local_x0,y,z)+Dcomp;
				indexfrom=IndexSlice_zy(y,z);
				Dmatrix[indexto]=-invNgrid*slice_tr[indexfrom];
			}
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+11);
			ElapsedInc(tvp+10,tvp+11,&Timing_ar3);
#endif
		} // end slice X
		if (IFROOT) printf(".");
	} // end of Dcomp
	if (IFROOT) printf("\n");
	// free vectors used for computation of Dmatrix; slice and slice_tr are freed after InitRmatrix
	Free_cVector(D2matrix);
#ifdef PARALLEL
	// deallocate buffers for BlockTranspose_DRm
	Free_general(BT_buffer);
	Free_general(BT_rbuffer);
#endif
#ifdef OPENCL
	// copy Dmatrix to OpenCL buffer, blocking to ensure completion before function end
	CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufDmatrix,CL_TRUE,0,Dsize*sizeof(*Dmatrix),Dmatrix,0,NULL,NULL));
	Free_cVector(Dmatrix);
#endif
	if (surface) { // only the total execution time of InitRmatrix is timed
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+12);
#endif
			InitRmatrix(invNgrid);
			Free_cVector(R2matrix); // free it here since it was allocated above
#ifdef PRECISE_TIMING
			GET_SYSTEM_TIME(tvp+13);
			t_Rm=DiffSystemTime(tvp+12,tvp+13);
#endif
	}
	Free_cVector(slice);
	Free_cVector(slice_tr);
#ifdef PARALLEL
	// allocate buffers for BlockTranspose
	MALLOC_VECTOR(BT_buffer,double,BTsize,ALL);
	MALLOC_VECTOR(BT_rbuffer,double,BTsize,ALL);
#endif
#ifndef OPENCL
	// allocate memory for Xmatrix, slices and slices_tr - used in matvec
	MALLOC_VECTOR(Xmatrix,complex,3*local_Nsmall,ALL);
	MALLOC_VECTOR(slices,complex,3*gridYZ,ALL);
	MALLOC_VECTOR(slices_tr,complex,3*gridYZ,ALL);
	if (surface) { // additional slices for reflection interaction
		MALLOC_VECTOR(slicesR,complex,3*gridYZ,ALL);
		MALLOC_VECTOR(slicesR_tr,complex,3*gridYZ,ALL);
	}
#endif
	time1=GET_TIME();
	Timing_Dm_Init=time1-start;

#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+14);
	/* time for extra initialization required for MatVec; it includes copying Dmatrix to GPU. Earlier we included in it
	 * a few printfs in the loop
	 */
	ElapsedInc(tvp+11,tvp+14,&Timing_InitMV);
	t_InitMV=TimerToSec(&Timing_InitMV);
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
	t_Tot=DiffSystemTime(tvp,tvp+14);

	if (surface) { // correct InitMV and Total time by that of InitRmatrix
		t_InitMV-=t_Rm;
		t_Tot-=t_Rm;
	}

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
		t_beg,t_Arithm,t_Gcalc,t_FFT,t_ar1,t_BT,t_fftX,t_InitMV,t_BT,t_ar2,t_Tot,t_fftZ,t_TYZ,t_fftY,
		t_ar3,t_InitMV);
	if (surface && IFROOT) PrintBoth(logfile,"Additionally time for initialization of Rmatrix = "FFORMPT"\n\n",t_Rm);
#endif

	fftInitAfterD();

	Timing_FFT_Init = GET_TIME()-time1;
}

//======================================================================================================================

void Free_FFT_Dmat(void)
// free all vectors that were allocated in fft.c (all used for FFT and MatVec)
{
#ifdef OPENCL
#	ifdef OCL_BLAS
	if (IterMethod==IT_BICG_CS) { // currently, used only in one iterative solver
		my_clReleaseBuffer(buftmp);
		my_clReleaseBuffer(bufxvec);
		my_clReleaseBuffer(bufrvec);
	}
#	endif
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
	if (surface) {
		my_clReleaseBuffer(bufRmatrix);
		my_clReleaseBuffer(bufslicesR);
		my_clReleaseBuffer(bufslicesR_tr);
	}
#	ifdef CLFFT_AMD
	CL_CH_ERR(clfftTeardown());
	oclMem-=clfftBufSize;
#	elif defined(CLFFT_APPLE)
	// the following do not return error status
	clFFT_DestroyPlan(clplanX);
	clFFT_DestroyPlan(clplanY);
	clFFT_DestroyPlan(clplanZ);
#	endif
	if (oclMem>0) LogWarning(EC_WARN,ALL_POS,"Possible leak of OpenCL memory (size %zu bytes) detected",oclMem);
#else
	Free_cVector(Dmatrix);
	Free_cVector(Xmatrix);
	Free_cVector(slices);
	Free_cVector(slices_tr);
	if (surface) {
		Free_cVector(Rmatrix);
		Free_cVector(slicesR);
		Free_cVector(slicesR_tr);
	}
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
	if (surface) {
		fftw_destroy_plan(planYRf);
		fftw_destroy_plan(planZRf);
	}
	fftw_cleanup();
#	endif
#endif
#ifdef FFT_TEMPERTON // these vectors are used even with OpenCL
	Free_general(work);
	Free_general(trigsX);
	Free_general(trigsY);
	Free_general(trigsZ);
#endif
}
