/* File: oclcore.h
 * $Date::                            $
 * Descr: all common OpenCL variables and functions
 *
 * Copyright (C) 2010-2012 ADDA contributors
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

#ifndef __oclvars_h
#define __oclvars_h

// system headers
#ifdef __APPLE__
#	include <OpenCL/cl.h>
#else
#	include <CL/cl.h>
#endif

// global OCL variables; names should not interfere with other parts of the code
extern cl_context context;
extern cl_command_queue command_queue;
extern cl_kernel clzero,clarith1,clarith2,clarith3,clarith4,clarith5,clnConj,clinprod,cltransposef,
	cltransposeb,cltransposeof,cltransposeob;
extern cl_mem bufXmatrix,bufmaterial,bufposition,bufcc_sqrt,bufargvec,bufresultvec,bufslices,
	bufslices_tr,bufDmatrix,bufinproduct;
extern double *inprodhlp;

void checkErr(cl_int err,const char * name);

#endif // __oclvars_h
