/* File: oclcore.h
 * $Date::                            $
 * Descr: all common OpenCL variables and functions; void in non-OpenCL mode
 *
 * Copyright (C) 2010-2014 ADDA contributors
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
#ifdef OPENCL

#ifndef __oclcore_h
#define __oclcore_h

// project headers
#include "function.h"
#include "io.h"
#include <stdbool.h> // for bool
// system headers
#ifdef __APPLE__
#	include <OpenCL/cl.h>
#else
#	include <CL/cl.h>
#endif

// This is redundant test for now, but may be easily updated in the future
#ifndef CL_VERSION_1_0
#	error "OpenCL version at least 1.0 is required"
#endif

// global OpenCL variables; names should not interfere with other parts of the code
extern cl_context context;
extern bool bufupload;
extern cl_command_queue command_queue;
extern cl_kernel clzero,clarith1,clarith2,clarith3,clarith3_surface,clarith4,clarith5,clnConj,clinprod,cltransposeof,
	cltransposeob,cltransposeofR;
extern cl_mem bufXmatrix,bufmaterial,bufposition,bufcc_sqrt,bufargvec,bufresultvec,bufslices,bufslices_tr,bufDmatrix,
	bufinproduct;
#ifdef OCL_BLAS
extern cl_mem buftmp,bufrvec,bufxvec;
#endif
extern cl_mem bufRmatrix,bufslicesR,bufslicesR_tr;
extern double *inprodhlp;
extern size_t oclMem,oclMemPeak,oclMemMaxObj;
extern cl_ulong oclMemDev,oclMemDevObj;
extern int gpuInd;

/* checks error status of CL functions; can either be used as a wrapper that returns error status or applied to return
 * error value. It is defined as a macro to incorporate exact position in a source file where it was called.
 */
#define CL_CH_ERR(a) CheckCLErr(a,ALL_POS,NULL)
#define CREATE_CL_BUFFER(buf,flags,size,ptr) buf=my_clCreateBuffer(flags,size,ptr,ALL_POS,#buf)
void oclinit(void);
void oclunload(void);
void PrintCLErr(cl_int err,ERR_LOC_DECL,const char * restrict msg) ATT_NORETURN;
cl_mem my_clCreateBuffer(cl_mem_flags mem_flags,size_t size,void *host_ptr,ERR_LOC_DECL,const char *name);
void my_clReleaseBuffer(cl_mem buffer);

//======================================================================================================================

static inline void CheckCLErr(const cl_int err,ERR_LOC_DECL,const char * restrict msg)
/* Checks error code and prints error if necessary. It is an inline wrapper, so it can be called after each CL function
 * without worrying about performance. Optional argument msg is added to the error message, if not NULL. Incorporating
 * it into the macro CL_CH_ERR is not trivial, since the argument can be a call to the function.
 */
{
	if (err != CL_SUCCESS) PrintCLErr(err,ERR_LOC_CALL,msg);
}

#endif // __oclcore_h

#endif // OPENCL
