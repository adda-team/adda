/* File: matvec.c
 * $Date::                            $
 * Descr: calculate local matrix vector product of decomposed interaction matrix with r_k or p_k, using a FFT-based
 *        convolution algorithm. Also contains code for SPARSE (non-FFT) mode.
 *
 * Copyright (C) 2006-2013 ADDA contributors
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
#include "fft.h"
#include "io.h"
#include "interaction.h"
#include "linalg.h"
#include "prec_time.h"
#include "sparse_ops.h"
#include "vars.h"

// SEMI-GLOBAL VARIABLES

#ifdef SPARSE
// defined and initialized in calculator.c
extern doublecomplex * restrict arg_full;
#else
// defined and initialized in fft.c
extern const doublecomplex * restrict Dmatrix,* restrict Rmatrix;
extern doublecomplex * restrict Xmatrix,* restrict slices,* restrict slices_tr,* restrict slicesR,* restrict slicesR_tr;
extern const size_t DsizeY,DsizeZ;
#endif // !SPARSE
extern const size_t RsizeY;
// defined and initialized in timing.c
extern size_t TotalMatVec;

// EXTERNAL FUNCTIONS

#ifdef PRECISE_TIMING
// calculator.c
void FreeEverything(void); // for proper finalization
#endif

#ifndef SPARSE
//======================================================================================================================

static inline size_t IndexSliceZY(const size_t y,const size_t z)
{
	return z*gridY+y;
}

//======================================================================================================================

static inline size_t IndexSliceYZ(const size_t y,const size_t z)
{
	return y*gridZ+z;
}

//======================================================================================================================

static inline size_t IndexGarbledX(const size_t x,const size_t y,const size_t z)
{
#ifdef PARALLEL
	return ((z%local_Nz)*smallY+y)*gridX+(z/local_Nz)*local_Nx+x%local_Nx;
#else
	return (z*smallY+y)*gridX+x;
#endif
}

//======================================================================================================================

static inline size_t IndexXmatrix(const size_t x,const size_t y,const size_t z)
{
	return (z*smallY+y)*gridX+x;
}

//======================================================================================================================

static inline size_t IndexDmatrix_mv(size_t x,size_t y,size_t z,const bool transposed)
{
	if (transposed) { // used only for G_SO
		/* reflection along the x-axis can't work in parallel mode, since the corresponding values are generally stored
		 * on a different processor. A rearrangement of memory distribution is required to remove this limitation.
		 */
		if (x>0) x=gridX-x;
		if (y>0) y=gridY-y;
		if (z>0) z=gridZ-z;
	}
	else {
		if (y>=DsizeY) y=gridY-y;
		if (z>=DsizeZ) z=gridZ-z;
	}

	return NDCOMP*((x*DsizeZ+z)*DsizeY+y);
}

//======================================================================================================================

static inline size_t IndexRmatrix_mv(size_t x,size_t y,size_t z,const bool transposed)
{
	if (transposed) { // used only for G_SO !!!
		/* reflection along the x-axis can't work in parallel mode, since the corresponding values are generally stored
		 * on a different processor. A rearrangement of memory distribution is required to remove this limitation.
		 */
		if (x>0) x=gridX-x;
		if (y>0) y=gridY-y;
	}
	else {
		if (y>=RsizeY) y=gridY-y;
	}

	return NDCOMP*((x*gridZ+z)*RsizeY+y);
}

#endif // !SPARSE

//======================================================================================================================

#ifndef SPARSE
void MatVec (doublecomplex * restrict argvec,    // the argument vector
             doublecomplex * restrict resultvec, // the result vector
             double *inprod,         // the resulting inner product
             const bool her,         // whether Hermitian transpose of the matrix is used
             TIME_TYPE *timing,      // this variable is incremented by total time
             TIME_TYPE *comm_timing) // this variable is incremented by communication time
/* This function implements matrix-vector product. If we want to calculate the inner product as well, we pass 'inprod'
 * as a non-NULL pointer. if 'inprod' is NULL, we don't calculate it. 'argvec' always remains unchanged afterwards,
 * however it is not strictly const - some manipulations may occur during the execution. comm_timing can be NULL, then
 * it is ignored.
 */
{
	size_t j,x;
	bool ipr,transposed;
	size_t boxY_st=boxY,boxZ_st=boxZ; // copies with different type
	size_t i;
	doublecomplex fmat[6],xv[3],yv[3],xvR[3],yvR[3];
	size_t index,y,z,Xcomp;
	unsigned char mat;
#ifdef PRECISE_TIMING
	SYSTEM_TIME tvp[18];
	SYSTEM_TIME Timing_FFTXf,Timing_FFTYf,Timing_FFTZf,Timing_FFTXb,Timing_FFTYb,Timing_FFTZb,Timing_Mult1,Timing_Mult2,
		Timing_Mult3,Timing_Mult4,Timing_Mult5,Timing_BTf,Timing_BTb,Timing_TYZf,Timing_TYZb,Timing_ipr;
	double t_FFTXf,t_FFTYf,t_FFTZf,t_FFTXb,t_FFTYb,t_FFTZb,t_Mult1,t_Mult2,t_Mult3,t_Mult4,t_Mult5,t_ipr,t_BTf,t_BTb,
		t_TYZf,t_TYZb,t_Arithm,t_FFT,t_Comm;

#endif

	/* A = I + S.D.S
	 * S = sqrt(C)
	 * A.x = x + S.D.(S.x)
	 * A(H).x = x + (S(T).D(T).S(T).x(*))(*)
	 * C,S - diagonal => symmetric
	 * (!! will change if tensor (non-diagonal) polarizability is used !!)
	 * D - symmetric (except for G_SO)
	 *
	 * D.x=F(-1)(F(D).F(X))
	 * F(D) is just a vector
	 *
	 * G_SO: F(D(T)) (k) =  F(D) (-k)
	 *       k - vector index
	 *
	 * For reflected matrix the situation is similar.
	 * R.x=F(-1)(F(R).H(X)), where R is a vector, similar with G, where R[i,j,k=0] is for interaction of two bottom
	 * dipoles. H(X) is FxFy(Fz^(-1)(X)), where Fx,... are Fourier transforms along corresponding coordinates. It can be
	 * computed along with F(X).
	 * Matrix R is symmetric (as a whole), but not in small parts, so R(i,j)=R(j,i)(T). Hence, in contrast to D, for
	 * 'transpose' actual transpose (changing sign of a few elements) of 3x3 submatrix is required along with addressing
	 * different elements of F(R).
	 *
	 * For (her) three additional operations of nConj are used. Should not be a problem, but can be avoided by a more
	 * complex code.
	 */
	TIME_TYPE tstart=GET_TIME();
	transposed=(!reduced_FFT) && her;
	ipr=(inprod!=NULL);
	if (ipr && !ipr_required) LogError(ONE_POS,"Incompatibility error in MatVec");
#ifdef PRECISE_TIMING
	InitTime(&Timing_FFTYf);
	InitTime(&Timing_FFTZf);
	InitTime(&Timing_FFTYb);
	InitTime(&Timing_FFTZb);
	InitTime(&Timing_Mult2);
	InitTime(&Timing_Mult3);
	InitTime(&Timing_Mult4);
	InitTime(&Timing_TYZf);
	InitTime(&Timing_TYZb);
	GET_SYSTEM_TIME(tvp);
#endif
	// FFT_matvec code
	if (ipr) *inprod = 0.0;
	// fill Xmatrix with 0.0
	for (i=0;i<3*local_Nsmall;i++) Xmatrix[i]=0.0;

	// transform from coordinates to grid and multiply with coupling constant
	if (her) nConj(argvec); // conjugated back afterwards

	for (i=0;i<local_nvoid_Ndip;i++) {
		// fill grid with argvec*sqrt_cc
		j=3*i;
		mat=material[i];
		index=IndexXmatrix(position[j],position[j+1],position[j+2]);
		// Xmat=cc_sqrt*argvec
		for (Xcomp=0;Xcomp<3;Xcomp++) Xmatrix[index+Xcomp*local_Nsmall]=cc_sqrt[mat][Xcomp]*argvec[j+Xcomp];
	}
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+1);
	Elapsed(tvp,tvp+1,&Timing_Mult1);
#endif
	// FFT X
	fftX(FFT_FORWARD); // fftX (buf)Xmatrix
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+2);
	Elapsed(tvp+1,tvp+2,&Timing_FFTXf);
#endif
#ifdef PARALLEL
	BlockTranspose(Xmatrix,comm_timing);
#endif
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+3);
	Elapsed(tvp+2,tvp+3,&Timing_BTf);
#endif
	// following is done by slices
	for(x=local_x0;x<local_x1;x++) {
		/* TODO: if z and y FFTs are interchanged, then computing reflected interaction can be optimized even further.
		 * Moreover, the typical situation of particles near surfaces, like large particulate slabs, correspond to the
		 * smallest dimension along z, which will also benefit from such interchange (then gridZ do not have to divide
		 * nprocs) - issue 177
		 */
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+4);
#endif
		// clear slice
		for(i=0;i<3*gridYZ;i++) slices[i]=0.0;
		// fill slices with values from Xmatrix
		for(y=0;y<boxY_st;y++) for(z=0;z<boxZ_st;z++) {
			i=IndexSliceYZ(y,z);
			j=IndexGarbledX(x,y,z);
			for (Xcomp=0;Xcomp<3;Xcomp++) slices[i+Xcomp*gridYZ]=Xmatrix[j+Xcomp*local_Nsmall];
		}
		// create a copy of slice, which is further transformed differently
		if (surface) memcpy(slicesR,slices,3*gridYZ*sizeof(doublecomplex));
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+5);
		ElapsedInc(tvp+4,tvp+5,&Timing_Mult2);
#endif
		// FFT z&y
		fftZ(FFT_FORWARD); // fftZ (buf)slices (and reflected terms)
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+6);
		ElapsedInc(tvp+5,tvp+6,&Timing_FFTZf);
#endif
		TransposeYZ(FFT_FORWARD); // including reflecting terms
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+7);
		ElapsedInc(tvp+6,tvp+7,&Timing_TYZf);
#endif
		fftY(FFT_FORWARD); // fftY (buf)slices_tr (and reflected terms)
#ifdef PRECISE_TIMING//
		GET_SYSTEM_TIME(tvp+8);
		ElapsedInc(tvp+7,tvp+8,&Timing_FFTYf);
#endif//
		// do the product D~*X~  and R~*X'~
		for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
			i=IndexSliceZY(y,z);
			for (Xcomp=0;Xcomp<3;Xcomp++) xv[Xcomp]=slices_tr[i+Xcomp*gridYZ];
			j=IndexDmatrix_mv(x-local_x0,y,z,transposed);
			memcpy(fmat,Dmatrix+j,6*sizeof(doublecomplex));
			if (reduced_FFT) { // symmetry with respect to reflection (x_i -> x_2N-i) is the same as in r-space
				if (y>=DsizeY) { // we assume that compiler will optimize x*=-1 into negation of sign
					fmat[1]*=-1;
					if (z>=DsizeZ) fmat[2]*=-1;
					else fmat[4]*=-1;
				}
				else if (z>=DsizeZ) {
					fmat[2]*=-1;
					fmat[4]*=-1;
				}
			}
			cSymMatrVec(fmat,xv,yv); // yv=fmat.xv
			if (surface) {
				for (Xcomp=0;Xcomp<3;Xcomp++) xvR[Xcomp]=slicesR_tr[i+Xcomp*gridYZ];
				j=IndexRmatrix_mv(x-local_x0,y,z,transposed);
				memcpy(fmat,Rmatrix+j,6*sizeof(doublecomplex));
				if (reduced_FFT && y>=RsizeY) {
					fmat[1]*=-1;
					fmat[4]*=-1;
				}
				if (transposed) { // corresponds to transpose of 3x3 matrix
					fmat[2]*=-1;
					fmat[4]*=-1;
				}
				// yv+=fmat.xvR
				cReflMatrVec(fmat,xvR,yvR);
				cvAdd(yvR,yv,yv);
			}
			for (Xcomp=0;Xcomp<3;Xcomp++) slices_tr[i+Xcomp*gridYZ]=yv[Xcomp];
		}
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+9);
		ElapsedInc(tvp+8,tvp+9,&Timing_Mult3);
#endif
		// inverse FFT y&z
		fftY(FFT_BACKWARD); // fftY (buf)slices_tr
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+10);
		ElapsedInc(tvp+9,tvp+10,&Timing_FFTYb);
#endif
		TransposeYZ(FFT_BACKWARD);
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+11);
		ElapsedInc(tvp+10,tvp+11,&Timing_TYZb);
#endif
		fftZ(FFT_BACKWARD); // fftZ (buf)slices
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+12);
		ElapsedInc(tvp+11,tvp+12,&Timing_FFTZb);
#endif
		//arith4 on host
		// copy slice back to Xmatrix
		for(y=0;y<boxY_st;y++) for(z=0;z<boxZ_st;z++) {
			i=IndexSliceYZ(y,z);
			j=IndexGarbledX(x,y,z);
			for (Xcomp=0;Xcomp<3;Xcomp++) Xmatrix[j+Xcomp*local_Nsmall]=slices[i+Xcomp*gridYZ];
		}
#ifdef PRECISE_TIMING
		GET_SYSTEM_TIME(tvp+13);
		ElapsedInc(tvp+12,tvp+13,&Timing_Mult4);
#endif
	} // end of loop over slices
	// FFT-X back the result
#ifdef PARALLEL
	BlockTranspose(Xmatrix,comm_timing);
#endif
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+14);
	Elapsed(tvp+13,tvp+14,&Timing_BTb);
#endif
	fftX(FFT_BACKWARD); // fftX (buf)Xmatrix
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+15);
	Elapsed(tvp+14,tvp+15,&Timing_FFTXb);
#endif
	// fill resultvec
	for (i=0;i<local_nvoid_Ndip;i++) {
		j=3*i;
		mat=material[i];
		index=IndexXmatrix(position[j],position[j+1],position[j+2]);
		for (Xcomp=0;Xcomp<3;Xcomp++) // result=argvec+cc_sqrt*Xmat
			resultvec[j+Xcomp]=argvec[j+Xcomp]+cc_sqrt[mat][Xcomp]*Xmatrix[index+Xcomp*local_Nsmall];
		// norm is unaffected by conjugation, hence can be computed here
		if (ipr) *inprod+=cvNorm2(resultvec+j);
	}
	if (her) {
		nConj(resultvec);
		nConj(argvec); // conjugate back argvec, so it remains unchanged after MatVec
	}
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+16);
	Elapsed(tvp+15,tvp+16,&Timing_Mult5);
#endif
	if (ipr) MyInnerProduct(inprod,double_type,1,comm_timing);
#ifdef PRECISE_TIMING
	GET_SYSTEM_TIME(tvp+17);
	Elapsed(tvp+16,tvp+17,&Timing_ipr);

	t_Mult1=TimerToSec(&Timing_Mult1);
	t_Mult2=TimerToSec(&Timing_Mult2);
	t_Mult3=TimerToSec(&Timing_Mult3);
	t_Mult4=TimerToSec(&Timing_Mult4);
	t_Mult5=TimerToSec(&Timing_Mult5);
	t_TYZf=TimerToSec(&Timing_TYZf);
	t_TYZb=TimerToSec(&Timing_TYZb);
	t_BTf=TimerToSec(&Timing_BTf);
	t_BTb=TimerToSec(&Timing_BTb);
	t_FFTXf=TimerToSec(&Timing_FFTXf);
	t_FFTXb=TimerToSec(&Timing_FFTXb);
	t_FFTYf=TimerToSec(&Timing_FFTYf);
	t_FFTYb=TimerToSec(&Timing_FFTYb);
	t_FFTZf=TimerToSec(&Timing_FFTZf);
	t_FFTZb=TimerToSec(&Timing_FFTZb);
	t_ipr=TimerToSec(&Timing_ipr);

	t_Arithm=t_Mult1+t_Mult2+t_Mult3+t_Mult4+t_Mult5+t_TYZf+t_TYZb;
	t_FFT=t_FFTXf+t_FFTYf+t_FFTZf+t_FFTXb+t_FFTYb+t_FFTZb;
	t_Comm=t_BTf+t_BTb+t_ipr;

	if (IFROOT) {
		PrintBoth(logfile,
			"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
			"                MatVec timing              \n"
			"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
			"Arith1 = "FFORMPT"    Arithmetics = "FFORMPT"\n"
			"FFTXf  = "FFORMPT"    FFT         = "FFORMPT"\n"
			"BTf    = "FFORMPT"    Comm        = "FFORMPT"\n"
			"Arith2 = "FFORMPT"\n"
			"FFTZf  = "FFORMPT"          Total = "FFORMPT"\n"
			"TYZf   = "FFORMPT"\n"
			"FFTYf  = "FFORMPT"\n"
			"Arith3 = "FFORMPT"\n"
			"FFTYb  = "FFORMPT"\n"
			"TYZb   = "FFORMPT"\n"
			"FFTZb  = "FFORMPT"\n"
			"Arith4 = "FFORMPT"\n"
			"BTb    = "FFORMPT"\n"
			"FFTXb  = "FFORMPT"\n"
			"Arith5 = "FFORMPT"\n"
			"InProd = "FFORMPT"\n\n",
			t_Mult1,t_Arithm,t_FFTXf,t_FFT,t_BTf,t_Comm,t_Mult2,t_FFTZf,DiffSystemTime(tvp,tvp+17),t_TYZf,t_FFTYf,
			t_Mult3,t_FFTYb,t_TYZb,t_FFTZb,t_Mult4,t_BTb,t_FFTXb,t_Mult5,t_ipr);
		printf("\nPrecise timing is complete. Finishing execution.\n");
	}
	FreeEverything();
	Stop(EXIT_SUCCESS);
#endif
	(*timing) += GET_TIME() - tstart;
	TotalMatVec++;
}

#else // SPARSE is defined

//======================================================================================================================

/* The sparse MatVec is implemented completely separately from the non-sparse version. Although there is some code
 * duplication, this probably makes the both versions easier to maintain.
*/
void MatVec (doublecomplex * restrict argvec,    // the argument vector
             doublecomplex * restrict resultvec, // the result vector
             double *inprod,         // the resulting inner product
             const bool her,         // whether Hermitian transpose of the matrix is used
             TIME_TYPE *timing,      // this variable is incremented by total time
             TIME_TYPE *comm_timing) // this variable is incremented by communication time
{
	const bool ipr = (inprod != NULL);
	size_t i,j,i3;

	TIME_TYPE tstart=GET_TIME();
	if (her) nConj(argvec);
	// TODO: can be replaced by nMult_mat
	for (j=0; j<local_nvoid_Ndip; j++) CcMul(argvec,arg_full+3*local_nvoid_d0,j);
#	ifdef PARALLEL
	AllGather(NULL,arg_full,cmplx3_type,comm_timing);
#	endif
	for (i=0; i<local_nvoid_Ndip; i++) {
		i3 = 3*i;
		cvInit(resultvec+i3);
		for (j=0; j<nvoid_Ndip; j++) AijProd(arg_full,resultvec,i,j);
	}
	// TODO: can be replaced by a specially designed function from linalg.c
	for (i=0; i<local_nvoid_Ndip; i++) DiagProd(argvec,resultvec,i);
	if (her) {
		nConj(resultvec);
		nConj(argvec);
	}
	if (ipr) (*inprod)=nNorm2(resultvec,comm_timing);
	(*timing) += GET_TIME() - tstart;
	TotalMatVec++;
}

#endif // SPARSE
