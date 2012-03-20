/* File: oclcore.c
 * $Date::                            $
 * Descr: core parts of OpenCL code
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
#include "const.h" // keep this first
#include "oclcore.h" // corresponding header
// project headers
#include "memory.h"
// system headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// GLOBAL VARIABLES
/* used in fft.c and matvec.c. Some of them are used only in one file, but we do not differentiate
 * them for simplicity (to keep all extern declarations in one place)
 */
cl_context context;
cl_command_queue command_queue;
cl_kernel clarith1,clarith2,clarith3,clarith4,clarith5,clzero,clinprod,clnConj,cltransposef,
	cltransposeb,cltransposeof,cltransposeob;
cl_mem bufXmatrix,bufmaterial,bufposition,bufcc_sqrt,bufargvec,bufresultvec,bufslices,bufslices_tr,
	bufDmatrix,bufinproduct;
double *inprodhlp; // extra buffer (on CPU) for calculating inner product in MatVec

// SEMI-GLOBAL VARIABLES

// used in cpp/fft_setup.cpp
const char *coptions;

// LOCAL VARIABLES

/* platform name. Defines whether the device is by NVIDIA, AMD, or something else ...
 * It is used for compiler options and (potentially) for conditional kernels
 */
enum platname {
	PN_NVIDIA, // Nvidia GPU
	PN_AMD,    // AMD GPU
	PN_UNDEF   // Unknown OpenCL platform
};
enum platname platformname;
// ids of program, platform and device
cl_program program;
cl_platform_id used_platform_id;
cl_device_id device_id;

// The kernel source is either loaded from oclkernels.cl at runtime or included at compile time
//#define READ_CL_SOURCE_AT_RUNTIME

// For some reason Eclipse points out a syntax error in the following block. Just ignore it.
#ifndef READ_CL_SOURCE_AT_RUNTIME
	const char stringifiedSourceCL[]=
	// the following is a pure string generated automatically from oclkernels.cl at compile time
#	include "ocl/oclkernels.clstr"
	;
#endif

//========================================================================

static const char *print_cl_errstring(cl_int err)
// produces meaningful error message from the error code
{
	switch (err) {
		case CL_SUCCESS:                         return "Success!";
		case CL_DEVICE_NOT_FOUND:                return "Device not found.";
		case CL_DEVICE_NOT_AVAILABLE:            return "Device not available";
		case CL_COMPILER_NOT_AVAILABLE:          return "Compiler not available";
		case CL_MEM_OBJECT_ALLOCATION_FAILURE:   return "Memory object allocation failure";
		case CL_OUT_OF_RESOURCES:                return "Out of resources";
		case CL_OUT_OF_HOST_MEMORY:              return "Out of host memory";
		case CL_PROFILING_INFO_NOT_AVAILABLE:    return "Profiling information not available";
		case CL_MEM_COPY_OVERLAP:                return "Memory copy overlap";
		case CL_IMAGE_FORMAT_MISMATCH:           return "Image format mismatch";
		case CL_IMAGE_FORMAT_NOT_SUPPORTED:      return "Image format not supported";
		case CL_BUILD_PROGRAM_FAILURE:           return "Program build failure";
		case CL_MAP_FAILURE:                     return "Map failure";
		case CL_INVALID_VALUE:                   return "Invalid value";
		case CL_INVALID_DEVICE_TYPE:             return "Invalid device type";
		case CL_INVALID_PLATFORM:                return "Invalid platform";
		case CL_INVALID_DEVICE:                  return "Invalid device";
		case CL_INVALID_CONTEXT:                 return "Invalid context";
		case CL_INVALID_QUEUE_PROPERTIES:        return "Invalid queue properties";
		case CL_INVALID_COMMAND_QUEUE:           return "Invalid command queue";
		case CL_INVALID_HOST_PTR:                return "Invalid host pointer";
		case CL_INVALID_MEM_OBJECT:              return "Invalid memory object";
		case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: return "Invalid image format descriptor";
		case CL_INVALID_IMAGE_SIZE:              return "Invalid image size";
		case CL_INVALID_SAMPLER:                 return "Invalid sampler";
		case CL_INVALID_BINARY:                  return "Invalid binary";
		case CL_INVALID_BUILD_OPTIONS:           return "Invalid build options";
		case CL_INVALID_PROGRAM:                 return "Invalid program";
		case CL_INVALID_PROGRAM_EXECUTABLE:      return "Invalid program executable";
		case CL_INVALID_KERNEL_NAME:             return "Invalid kernel name";
		case CL_INVALID_KERNEL_DEFINITION:       return "Invalid kernel definition";
		case CL_INVALID_KERNEL:                  return "Invalid kernel";
		case CL_INVALID_ARG_INDEX:               return "Invalid argument index";
		case CL_INVALID_ARG_VALUE:               return "Invalid argument value";
		case CL_INVALID_ARG_SIZE:                return "Invalid argument size";
		case CL_INVALID_KERNEL_ARGS:             return "Invalid kernel arguments";
		case CL_INVALID_WORK_DIMENSION:          return "Invalid work dimension";
		case CL_INVALID_WORK_GROUP_SIZE:         return "Invalid work group size";
		case CL_INVALID_WORK_ITEM_SIZE:          return "Invalid work item size";
		case CL_INVALID_GLOBAL_OFFSET:           return "Invalid global offset";
		case CL_INVALID_EVENT_WAIT_LIST:         return "Invalid event wait list";
		case CL_INVALID_EVENT:                   return "Invalid event";
		case CL_INVALID_OPERATION:               return "Invalid operation";
		case CL_INVALID_GL_OBJECT:               return "Invalid OpenGL object";
		case CL_INVALID_BUFFER_SIZE:             return "Invalid buffer size";
		case CL_INVALID_MIP_LEVEL:               return "Invalid mip-map level";
		default:                                 return "Unknown";
	}
}

//========================================================================

void checkErr(cl_int err,const char *name)
// checks error code and stops if necessary
{
	if (err != CL_SUCCESS) {
		printf("ERROR: %s / Error code (%i : %s)\n",name,err,print_cl_errstring(err));
		exit(EXIT_FAILURE);
	}
}

//========================================================================

static void GetDevice(void)
// set OpenCL device number and related variable
{
	cl_int err; // error code
	cl_uint num_of_platforms,num_of_devices;
	cl_platform_id *platform_id;
	char pname[MAX_LINE],devicename[MAX_LINE];
	cl_int devtype=CL_DEVICE_TYPE_GPU; // set preferred device type

	// little trick to get just the number of the Platforms
	err=clGetPlatformIDs(0,NULL,&num_of_platforms);
	checkErr(err,"Get number of platforms");
	// dynamic array of platformids creation at runtime, stored in heap
	platform_id=(cl_platform_id *)
		voidVector(num_of_platforms*sizeof(platform_id),ALL_POS,"platform_id");
	err=clGetPlatformIDs(num_of_platforms,platform_id,NULL);
	checkErr(err,"Get platforms IDs");

	for (unsigned int i=0;i<num_of_platforms;i++) {
		/* Some errors are allowed here, since some platforms/devices may be not supported.
		 * The execution continues normally if at least one supported device is found.
		 * So errors are handled directly here, not through
		 */
		clGetPlatformInfo(platform_id[i],CL_PLATFORM_NAME,sizeof(pname),pname,NULL);
		err=clGetDeviceIDs(platform_id[i],devtype,1,&device_id,&num_of_devices);
		if (err==CL_SUCCESS) {
			if (strcmp(pname,"NVIDIA CUDA")==0) platformname=PN_NVIDIA;
			else if (strcmp(pname,"ATI Stream")==0) platformname=PN_AMD;
			/* the program can potentially work with unknown compatible device, but performance is,
			 * to some extent, unpredictable.
			 */
			else platformname=PN_UNDEF;
			err=clGetDeviceInfo(device_id,CL_DEVICE_NAME,sizeof(devicename),devicename,NULL);
			checkErr(err,"Get device name");
			printf("Found OpenCL device %s, based on %s.\n",devicename,pname);
			used_platform_id=platform_id[i];
			break;
		}
	}
	// if all platforms of the above cycle fail, then error is produced
	checkErr(err,"No valid preferred OpenCL device found");
	Free_general(platform_id);
}

//========================================================================

void oclinit(void)
// initialize OpenCL environment
{
	cl_int err; // error code
	const char *cssPtr[1]; // pointer to (array of) cSourceString, required to avoid warnings

	// getting the first OpenCL device which is a GPU
	GetDevice();
	/* cl_context_properties is a strange list of item:
	 * first comes the name of the property as next element the corresponding value
	 */
	cl_context_properties properties[3];
	properties[0]=CL_CONTEXT_PLATFORM;
	properties[1]=(cl_context_properties)used_platform_id;
	properties[2]=0; // last one must be zero
	context=clCreateContext(
		properties,
		1, //number of devices to use
		&device_id,
		NULL, // error info as string
		NULL, // user data
		&err
	);
	checkErr(err,"Create Context");

	command_queue=clCreateCommandQueue(
		context,
		device_id,
		CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, // command queue properties
		&err
	);
	checkErr(err,"Command_queue");

#ifdef READ_CL_SOURCE_AT_RUNTIME
	size_t sourceStrSize;
	char *cSourceString;
	FILE *file=FOpenErr("oclkernels.cl","rb",ALL_POS);
	fseek(file,0,SEEK_END);
	sourceStrSize=ftell(file);
	fseek(file,0,SEEK_SET);
	MALLOC_VECTOR(cSourceString,char,sourceStrSize+1,ALL);
	fread(cSourceString,sourceStrSize,1,file);
	fclose(file);
	cssPtr[0]=(const char *)cSourceString;
#else
	cssPtr[0]=stringifiedSourceCL;
#endif
	program=clCreateProgramWithSource(
		context, // valid Context
		1,       // number of strings in next parameter
		cssPtr,  // array of source strings
		NULL,    // length of each string
		&err     // error code
	);
	checkErr(err,"Creating program");

	if (platformname==PN_NVIDIA) coptions="-cl-mad-enable"; //for NVIDIA GPUs
	else if (platformname==PN_AMD) coptions="-DAMD -cl-mad-enable"; //for AMD GPUs
	else if (platformname==PN_UNDEF) coptions=""; // no options for unknown OpenCL platform

	// special error handling to enable debugging of OpenCL kernels
	if (clBuildProgram(program,0,NULL,coptions,NULL,NULL) != CL_SUCCESS) {
		printf("Error building OpenCL program\n");
		char buffer[8192];
		clGetProgramBuildInfo(
			program,              // valid program object
			device_id,            // valid device_id that executable was built
			CL_PROGRAM_BUILD_LOG, // indicate to retrieve build log
			sizeof(buffer),       // size of the buffer to write log to
			buffer,               // the actual buffer to write log to
			NULL);                // the actual size in bytes of data copied to buffer
		printf("%s\n",buffer);
		exit(EXIT_FAILURE);
	}

	clzero=clCreateKernel(program,"clzero",&err);
	checkErr(err,"creating kernel clzero");
	clarith1=clCreateKernel(program,"arith1",&err);
	checkErr(err,"creating kernel clarith1");
	clarith2=clCreateKernel(program,"arith2",&err);
	checkErr(err,"creating kernel clartih2");
	clarith3=clCreateKernel(program,"arith3",&err);
	checkErr(err,"creating kernel clartih3");
	clarith4=clCreateKernel(program,"arith4",&err);
	checkErr(err,"creating kernel clarith4");
	clarith5=clCreateKernel(program,"arith5",&err);
	checkErr(err,"creating kernel clarith5");
	clnConj=clCreateKernel(program,"nConj",&err);
	checkErr(err,"creating kernel clnConj");
	clinprod=clCreateKernel(program,"inpr",&err);
	checkErr(err,"creating kernel clinprod");
	cltransposef=clCreateKernel(program,"transpose",&err);
	checkErr(err,"creating kernel cltransposef");
	cltransposeb=clCreateKernel(program,"transpose",&err);
	checkErr(err,"creating kernel cltransposeb");
	cltransposeof=clCreateKernel(program,"transposeo",&err);
	checkErr(err,"creating kernel cltransposeof");
	cltransposeob=clCreateKernel(program,"transposeo",&err);
	checkErr(err,"creating kernel cltransposeob");
#ifdef READ_CL_SOURCE_AT_RUNTIME
	Free_general(cSourceString);
#endif

}

//========================================================================

void oclunload(void)
// unload all OpenCL kernels and similar stuff; buffers are released in Free_FFT_Dmat()
{
	clReleaseProgram(program);
	clReleaseKernel(clzero);
	clReleaseKernel(clarith1);
	clReleaseKernel(clarith2);
	clReleaseKernel(clarith3);
	clReleaseKernel(clarith4);
	clReleaseKernel(clarith5);
	clReleaseKernel(clnConj);
	clReleaseKernel(clinprod);
	clReleaseKernel(cltransposef);
	clReleaseKernel(cltransposeb);
	clReleaseKernel(cltransposeof);
	clReleaseKernel(cltransposeob);
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);
}
