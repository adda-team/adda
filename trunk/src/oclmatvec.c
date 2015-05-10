/* File: oclmatvec.c
 * $Date::                            $
 * Descr: calculate local matrix vector product of decomposed interaction matrix with r_k or p_k, using a FFT-based
 *        convolution algorithm. This code is special for OpenCL mode.
 *
 * Copyright (C) 2014 ADDA contributors
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
#include "comm.h"
#include "fft.h"
#include "io.h"
#include "oclcore.h"
#include "vars.h"

// SEMI-GLOBAL VARIABLES

// defined and initialized in fft.c
extern const size_t RsizeY,clxslices,local_gridX,slicesize;
// defined and initialized in timing.c
extern size_t TotalMatVec;

//======================================================================================================================

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
	size_t j;
	bool ipr,transposed;
	size_t boxY_st=boxY,boxZ_st=boxZ; // copies with different type

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
	// FFT_matvec code
	if (ipr) *inprod = 0.0;
	const cl_char ndcomp=NDCOMP;
	// little workaround for kernel cannot take bool arguments
	const cl_char transp=(cl_char)transposed;
	const cl_char redfft=(cl_char)reduced_FFT;

	/* following two calls to clSetKernelArg can be moved to fft.c, since the arguments are constant. However, this
	 * requires setting auxiliary variables redfft and ndcomp as globals, since the kernel is called below.
	 */
	CL_CH_ERR(clSetKernelArg(clarith3,7,sizeof(cl_char),&ndcomp));
	CL_CH_ERR(clSetKernelArg(clarith3,8,sizeof(cl_char),&redfft));
	CL_CH_ERR(clSetKernelArg(clarith3,9,sizeof(cl_char),&transp));
	if (surface) { // arguments for surface-version of the arith3 kernel
		CL_CH_ERR(clSetKernelArg(clarith3_surface,7,sizeof(cl_char),&ndcomp));
		CL_CH_ERR(clSetKernelArg(clarith3_surface,8,sizeof(cl_char),&redfft));
		CL_CH_ERR(clSetKernelArg(clarith3_surface,9,sizeof(cl_char),&transp));
		CL_CH_ERR(clSetKernelArg(clarith3_surface,10,sizeof(size_t),&RsizeY));
	}
	// write into buffers e.g. upload to device; non-blocking
	if (bufupload) CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufargvec,CL_FALSE,0,local_nRows*sizeof(doublecomplex),
		argvec,0,NULL,NULL));

	size_t xmsize=local_Nsmall*3;
	if (her) {
		CL_CH_ERR(clSetKernelArg(clnConj,0,sizeof(cl_mem),&bufargvec));
		CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clnConj,1,NULL,&local_nRows,NULL,0,NULL,NULL));
	}
	// setting (buf)Xmatrix with zeros (on device)
	CL_CH_ERR(clSetKernelArg(clzero,0,sizeof(cl_mem),&bufXmatrix));
	CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clzero,1,NULL,&xmsize,NULL,0,NULL,NULL));
	CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clarith1,1,NULL,&local_nvoid_Ndip,NULL,0,NULL,NULL));
	// FFT X
	fftX(FFT_FORWARD); // fftX (buf)Xmatrix

	/* In OpenCL mode the free memory on the GPU was determined during fft.c and if enough memory is available slices
	 * contain the full fft grid. If not, FFT grid is split into "clxslices" parts with "local_gridX" length and kernels
	 * run with offsets inside a loop. Setting arith2, arith4 and arith3 run slices global work sizes and offsets.
	 */
	size_t gwsarith24[3]={local_gridX,boxZ_st,boxY_st};
	size_t gwsclarith3[3]={gridY,gridZ,local_gridX};
	size_t gwo24[3]={0,0,0};
	size_t gwo3[3]={0,0,0};

	for (size_t xsect=0; xsect<clxslices; xsect++) {
		// global work size offset in x axis, when not in first slice
		gwo24[0]=xsect*local_gridX;
		gwo3[2]=xsect*local_gridX;
		// last bunch of slices reached; ensure that arith2,3,4 do not run over border of gridX
		if (xsect+1==clxslices && gridX%local_gridX!=0) {
			gwsarith24[0]=gridX%local_gridX; 
			gwsclarith3[2]=gridX%local_gridX; 
		}
		// every index and argument prepared, starting arith2
		CL_CH_ERR(clSetKernelArg(clzero,0,sizeof(cl_mem),&bufslices));
		CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clzero,1,NULL,&slicesize,NULL,0,NULL,NULL));
		CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clarith2,3,gwo24,gwsarith24,NULL,0,NULL,NULL));
		if (surface) CL_CH_ERR(clEnqueueCopyBuffer(command_queue,bufslices,bufslicesR,0,0,
			slicesize*sizeof(doublecomplex),0,NULL,NULL));

		fftZ(FFT_FORWARD); // fftZ (buf)slices (and reflected terms)
		TransposeYZ(FFT_FORWARD); // including reflecting terms
		fftY(FFT_FORWARD); // fftY (buf)slices_tr (and reflected terms)
		// arith3 on Device
		if (surface) 
			CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clarith3_surface,3,gwo3,gwsclarith3,NULL,0,NULL,NULL));
		else 
			CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clarith3,3,gwo3,gwsclarith3,NULL,0,NULL,NULL));
		// inverse FFT y&z
		fftY(FFT_BACKWARD); // fftY (buf)slices_tr
		TransposeYZ(FFT_BACKWARD);
		fftZ(FFT_BACKWARD); // fftZ (buf)slices

		CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clarith4,3,gwo24,gwsarith24,NULL,0,NULL,NULL));
	}

	// FFT-X back the result
	fftX(FFT_BACKWARD); // fftX (buf)Xmatrix
	CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clarith5,1,NULL,&local_nvoid_Ndip,NULL,0,NULL,NULL));
	if (ipr) {
		/* calculating inner product in OpenCL is more complicated than usually. The norm for each element is calculated
		 * inside GPU, but the sum is taken by CPU afterwards. Hence, additional large buffers are required.
		 * Potentially, this can be optimized.
		 */
		CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clinprod,1,NULL,&local_nvoid_Ndip,NULL,0,NULL,NULL));
		CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufinproduct,CL_TRUE,0,local_nvoid_Ndip*sizeof(double),inprodhlp,0,
			NULL,NULL));
		// sum up on the CPU after calculating the norm on GPU; hence the read above is blocking
		for (j=0;j<local_nvoid_Ndip;j++) *inprod+=inprodhlp[j];
	}
	if (her) {
		CL_CH_ERR(clSetKernelArg(clnConj,0,sizeof(cl_mem),&bufresultvec));
		CL_CH_ERR(clEnqueueNDRangeKernel(command_queue,clnConj,1,NULL,&local_nRows,NULL,0,NULL,NULL));
	}
	// blocking read to finalize queue
	if (bufupload) CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufresultvec,CL_TRUE,0,local_nRows*sizeof(doublecomplex),
		resultvec,0,NULL,NULL));
	if (ipr) MyInnerProduct(inprod,double_type,1,comm_timing);
	(*timing) += GET_TIME() - tstart;
	TotalMatVec++;
}
