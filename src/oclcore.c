/* File: oclcore.c
 * $Date::                            $
 * Descr: core parts of OpenCL code
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

/* This file should be compiled only in OpenCL mode, hence the following declaration is redundant. However, it helps
 * proper syntax checking in IDE, such as Eclipse.
 */
#ifndef OPENCL
#  define OPENCL
#endif

#include "const.h" // keep this first
#include "oclcore.h" // corresponding header
// project headers
#include "debug.h"
#include "memory.h"
#include "vars.h"
// system headers
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// GLOBAL VARIABLES
/* used in fft.c and matvec.c. Some of them are used only in one file, but we do not differentiate them for simplicity
 * (to keep all extern declarations in one place)
 */
cl_context context;
cl_command_queue command_queue;
cl_kernel clarith1,clarith2,clarith3,clarith3_surface,clarith4,clarith5,clzero,clinprod,clnConj,cltransposeof,
	cltransposeob,cltransposeofR;
cl_mem bufXmatrix,bufmaterial,bufposition,bufcc_sqrt,bufargvec,bufresultvec,bufslices,bufslices_tr,bufDmatrix,
	bufinproduct;

/* defines if bufargvec and bufresultvec are to be uploaded in the beginning of MatVec
 * it is false only in case when the iterative solver handles the buffer upload
 */
bool bufupload=true;

#ifdef OCL_BLAS
cl_mem buftmp;  // temporary buffer for dot products and Norm (as required by clBLAS library)
cl_mem bufrvec; // buffer used in iterative solver
cl_mem bufxvec; // buffer used in iterative solver
#endif

cl_mem bufRmatrix,bufslicesR,bufslicesR_tr; //for surface
double *inprodhlp; // extra buffer (on CPU) for calculating inner product in MatVec
// OpenCL memory counts (current, peak, and maximum for a single object)
size_t oclMem,oclMemPeak,oclMemMaxObj;
// OpenCL memory available at device (total and for a single object); cl_ulong should be compatible with size_t
cl_ulong oclMemDev,oclMemDevObj;
int gpuInd; // index of GPU to use (starting from 0)

// SEMI-GLOBAL VARIABLES

// used in cpp/fft_setup.cpp
const char *coptions;

// LOCAL VARIABLES

/* platform name. Defines whether the device is by NVIDIA, AMD, or something else ... It is used for compiler options
 * and (potentially) for conditional kernels
 */
enum platname {
	PN_NVIDIA, // Nvidia GPU
	PN_AMD,    // AMD GPU
	PN_UNDEF   // Unknown OpenCL platform
};
enum platname platformname; // also used in cpp/fft_setup.cpp

// ids of program, platform and device
static cl_program program;
static cl_platform_id used_platform_id;
static cl_device_id device_id;

struct string {
	char *text; // text of string
	size_t len; // length
	size_t mem; // amount of allocated memory
};

// The kernel source is either loaded from oclkernels.cl at runtime or included at compile time
//#define OCL_READ_SOURCE_RUNTIME

#ifndef OCL_READ_SOURCE_RUNTIME
	const char stringifiedSourceCL[]=""
	// the following is a pure string generated automatically from oclkernels.cl at compile time
#	include "ocl/oclkernels.clstr"
	;
#endif

//======================================================================================================================

#define STR_SIZE_STEP 100
static struct string StrInit(void)
// initializes string structure
{
	struct string str;
	str.len=0;
	str.mem=STR_SIZE_STEP;
	MALLOC_VECTOR(str.text,char,str.mem,ALL);
	str.text[0]=0;
	return str;
}

//======================================================================================================================

static void StrCatSpace(struct string *dest,const char *src)
/* concatenates dest string (described by structure) and src string into the dest. Additionally, space is added between
 * the strings, if dest is originally not empty. If needed additional memory is allocated automatically.
 */
{
	if (dest->len>0) {
		dest->text[dest->len]=' ';
		dest->len++;
	}
	size_t slen=strlen(src);
	if (dest->len+slen>=dest->mem) {
		do dest->mem+=STR_SIZE_STEP; while (dest->len+slen>=dest->mem);
		REALLOC_VECTOR(dest->text,char,dest->mem,ALL);
	}
	strcpy(dest->text+dest->len,src);
	dest->len+=slen;
}
#undef STR_SIZE_STEP

//======================================================================================================================

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

//======================================================================================================================

void PrintCLErr(cl_int err,ERR_LOC_DECL,const char * restrict msg)
// Prints explicit information about CL error. Optional argument msg is added to the error message, if not NULL.
{
	if (msg==NULL) LogError(ERR_LOC_CALL,"CL error code %d: %s\n",err,print_cl_errstring(err));
	else LogError(ERR_LOC_CALL,"%s (CL error code %d: %s)\n",msg,err,print_cl_errstring(err));
}

//======================================================================================================================

static char* ATT_MALLOC dyn_clGetPlatformInfo(cl_platform_id plat_id,cl_platform_info param_name)
/* wrapper for clGetPlatformInfo with string return value, it automatically allocates the string to hold the result and
 * returns it to the caller. All error checks are performed inside.
 */
{
	size_t size;
	char *res;
	CL_CH_ERR(clGetPlatformInfo(plat_id,param_name,0,NULL,&size));
	MALLOC_VECTOR(res,char,size,ALL);
	CL_CH_ERR(clGetPlatformInfo(plat_id,param_name,size,res,NULL));
	return res;
}

//======================================================================================================================

static char* ATT_MALLOC dyn_clGetDeviceInfo(cl_device_id dev_id,cl_device_info param_name)
/* wrapper for clGetDeviceInfo with string return value, it automatically allocates the string to hold the result and
 * returns it to the caller. All error checks are performed inside.
 */
{
	size_t size;
	char *res;
	CL_CH_ERR(clGetDeviceInfo(dev_id,param_name,0,NULL,&size));
	MALLOC_VECTOR(res,char,size,ALL);
	CL_CH_ERR(clGetDeviceInfo(dev_id,param_name,size,res,NULL));
	return res;
}

//======================================================================================================================

static char* ATT_MALLOC dyn_clGetProgramBuildInfo(cl_program prog,cl_device_id dev_id,cl_program_build_info param_name)
/* wrapper for clGetProgramBuildInfo with string return value, it automatically allocates the string to hold the result
 * and returns it to the caller. All error checks are performed inside.
 */
{
	size_t size;
	char *res;
	CL_CH_ERR(clGetProgramBuildInfo(prog,dev_id,param_name,0,NULL,&size));
	MALLOC_VECTOR(res,char,size,ALL);
	CL_CH_ERR(clGetProgramBuildInfo(prog,dev_id,param_name,size,res,NULL));
	return res;
}
//======================================================================================================================

static void GetDevice(struct string *copt_ptr)
// set OpenCL device number and related variable
{
	cl_int err; // error code
	cl_uint num_of_platforms,num_of_devices;
	cl_platform_id *platform_id;
	cl_device_id *devices;
	const cl_int devtype=CL_DEVICE_TYPE_GPU; // set preferred device type
	int gpuN=0;

	printf("Searching for OpenCL devices\n");
	// little trick to get just the number of the Platforms
	CL_CH_ERR(clGetPlatformIDs(0,NULL,&num_of_platforms));
	/* OpenCL standard is somewhat unclear whether clGetPlatformIDs can return zero num_of_platforms with successful
	 * return status. So we additionally test it.
	 */
	if (num_of_platforms==0) LogError(ALL_POS,"No OpenCL platform was found");
	// dynamic array of platformids creation at runtime, stored in heap
	platform_id=(cl_platform_id *)voidVector(num_of_platforms*sizeof(platform_id),ALL_POS,"platform_id");
	CL_CH_ERR(clGetPlatformIDs(num_of_platforms,platform_id,NULL));

	for (unsigned int i=0;i<num_of_platforms;i++) {
		/* Some errors are allowed here, since some platforms/devices may be not supported. The execution continues
		 * normally if at least one (or gpuInd+1) supported device is found. So errors are handled directly here, not
		 * through standard error handler
		 */
		err=clGetDeviceIDs(platform_id[i],devtype,0,NULL,&num_of_devices);
		if (err==CL_SUCCESS && (gpuN+=num_of_devices)>gpuInd) {
			// choose specific device based on gpuInd
			devices=(cl_device_id *)voidVector(num_of_devices*sizeof(device_id),ALL_POS,"devices");
			CL_CH_ERR(clGetDeviceIDs(platform_id[i],devtype,num_of_devices,devices,NULL));
			device_id=devices[gpuInd-gpuN+num_of_devices];
			Free_general(devices);
			// get platform and device name
			used_platform_id=platform_id[i];
			char *pname=dyn_clGetPlatformInfo(platform_id[i],CL_PLATFORM_NAME);
			if (strcmp(pname,"NVIDIA CUDA")==0) platformname=PN_NVIDIA;
			else if (strcmp(pname,"ATI Stream")==0) platformname=PN_AMD;
			// the program can potentially work with unknown compatible device, but performance is unpredictable.
			else platformname=PN_UNDEF;
			char *devicename=dyn_clGetDeviceInfo(device_id,CL_DEVICE_NAME);
			PrintBoth(logfile,"Using OpenCL device %s, based on %s.\n",devicename,pname);
			Free_general(pname);
			Free_general(devicename);
			// get further device info
#ifdef DEBUGFULL
			char *dev_vers=dyn_clGetDeviceInfo(device_id,CL_DEVICE_VERSION);
			char *dr_vers=dyn_clGetDeviceInfo(device_id,CL_DRIVER_VERSION);
			D("Device version: %s. Driver version: %s.",dev_vers,dr_vers);
			Free_general(dev_vers);
			Free_general(dr_vers);
#endif
			CL_CH_ERR(clGetDeviceInfo(device_id,CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(oclMemDev),&oclMemDev,NULL));
			CL_CH_ERR(clGetDeviceInfo(device_id,CL_DEVICE_MAX_MEM_ALLOC_SIZE,sizeof(oclMemDevObj),&oclMemDevObj,NULL));
			// round numbers are expected so .0f is used, float is just for convenience
			PrintBoth(logfile,"Device memory: total - %.0f MB, maximum object - %.0f MB\n",oclMemDev/MBYTE,
				oclMemDevObj/MBYTE);
			char *dev_ext=dyn_clGetDeviceInfo(device_id,CL_DEVICE_EXTENSIONS);
			StrCatSpace(copt_ptr,"-DUSE_DOUBLE");
			if (strstr(dev_ext,"cl_khr_fp64")==NULL) {
				// fallback for old AMD GPUs
				if (strstr(dev_ext,"cl_amd_fp64")!=NULL) StrCatSpace(copt_ptr,"-DDOUBLE_AMD");
				else LogError(ALL_POS,"Double precision is not supported by OpenCL device");
			}
			Free_general(dev_ext);
			break;
		}
	}
	// if all platforms of the above cycle fail, then error is produced with the last error code
	if (gpuN==0) LogError(ALL_POS,"No OpenCL-compatible GPU found");
	else if (gpuN<=gpuInd) LogError(ALL_POS,"The specified GPU index (%d) must be less than the total number of "
		"available GPUs (%d)",gpuInd,gpuN);
	Free_general(platform_id);
}

//======================================================================================================================

void oclinit(void)
// initialize OpenCL environment
/* The whole OpenCL part relies on correspondence between complex and two-part, which is mandated by C99. The data
 * transfers are done through void pointers, so compiler doesn't generate warnings for casts. So that should be both
 * convenient and portable.
 */
{
	cl_int err; // error code
	const char *cssPtr[1]; // pointer to (array of) cSourceString, required to avoid warnings

	D("Starting OpenCL init");
	struct string cl_opt=StrInit();
	// getting the first OpenCL device which is a GPU
	GetDevice(&cl_opt);
	/* cl_context_properties is a strange list of item: first comes the name of the property and next element is the
	 * corresponding value
	 */
	cl_context_properties properties[3];
	properties[0]=CL_CONTEXT_PLATFORM;
	properties[1]=(cl_context_properties)used_platform_id;
	properties[2]=0; // last one must be zero
	context=clCreateContext(properties,1,&device_id,NULL,NULL,&err);
	CL_CH_ERR(err);
	// since we use only one context the following variables are implicitly associated with it
	oclMem=oclMemPeak=oclMemMaxObj=0;

	// for now we use in-order execution only, since it is much safer
	command_queue=clCreateCommandQueue(context,device_id,0,&err);
	CL_CH_ERR(err);

#ifdef OCL_READ_SOURCE_RUNTIME
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
	D("Creating CL program");
	program=clCreateProgramWithSource(context,1,cssPtr,NULL,&err);
	CL_CH_ERR(err);
	// finalize build options and point coptions (used for sharing with clFFT) to it
	StrCatSpace(&cl_opt,"-cl-mad-enable -cl-fast-relaxed-math");
	/* OpenCL standard does define size_t-analogue type to avoid confusion. But we want to use home-made definition in
	 * OpenCL source, since size_t is used a lot in ADDA to indicate the largest problem solvable on current hardware.
	 * Since the proper OpenCL type is chosen during compilation of the main ADDA source, there should be no portability
	 * problems, even though the OpenCL sources may be compiled (at ADDA runtime) on a different hardware than ADDA
	 * itself.
	 * In the following it is more logical to compare against CL_UINT_MAX and CL_ULONG_MAX, but at least the latter
	 * contains typecast "(cl_ulong)" which is badly interpreted by the preprocessor. The values used should always be
	 * the same according to OpenCL standard.
	 */
#if (SIZE_MAX==UINT32_MAX)
	StrCatSpace(&cl_opt,"-DSIZET_UINT");
#elif (SIZE_MAX==UINT64_MAX)
	StrCatSpace(&cl_opt,"-DSIZET_ULONG");
#else
#	error "No OpenCL alternative for size_t. Create an issue at http://code.google.com/p/a-dda/issues/"
#endif
	coptions=cl_opt.text;
	D("Building CL program with options: '%s'",cl_opt.text);
	// special error handling to enable debugging of OpenCL kernels
	err=clBuildProgram(program,0,NULL,cl_opt.text,NULL,NULL);
	if (err!=CL_SUCCESS) {
		char *buffer=dyn_clGetProgramBuildInfo(program,device_id,CL_PROGRAM_BUILD_LOG);
		printf("Following errors occurred during building of OpenCL program:\n%s\n",buffer);
		Free_general(buffer);
		CL_CH_ERR(err);
	}
	clzero=clCreateKernel(program,"clzero",&err);
	CL_CH_ERR(err);
	clarith1=clCreateKernel(program,"arith1",&err);
	CL_CH_ERR(err);
	clarith2=clCreateKernel(program,"arith2",&err);
	CL_CH_ERR(err);
	clarith3=clCreateKernel(program,"arith3",&err);
	CL_CH_ERR(err);
	clarith4=clCreateKernel(program,"arith4",&err);
	CL_CH_ERR(err);
	clarith5=clCreateKernel(program,"arith5",&err);
	CL_CH_ERR(err);
	clnConj=clCreateKernel(program,"nConj",&err);
	CL_CH_ERR(err);
	clinprod=clCreateKernel(program,"inpr",&err);
	CL_CH_ERR(err);
	/* In principle a single kernel can be used for all transpose operations, including surface ones. However, this will
	 * require setting the kernel arguments just before the execution. Thus, using multiple kernels seem to be a bit
	 * faster.
	 */
	cltransposeof=clCreateKernel(program,"transposeo",&err);
	CL_CH_ERR(err);
	cltransposeob=clCreateKernel(program,"transposeo",&err);
	CL_CH_ERR(err);
	if (surface) { // kernels for surface
		clarith3_surface=clCreateKernel(program,"arith3_surface",&err);
		CL_CH_ERR(err);
		cltransposeofR=clCreateKernel(program,"transposeo",&err);
		CL_CH_ERR(err);
	}
#ifdef OCL_READ_SOURCE_RUNTIME
	Free_general(cSourceString);
#endif
	CL_CH_ERR(clUnloadCompiler());
	D("OpenCL init complete");
}

//======================================================================================================================

cl_mem my_clCreateBuffer(cl_mem_flags mem_flags,size_t size,void *host_ptr,ERR_LOC_DECL,const char *name)
// wrapper to create buffer, which also adjusts memory counts and takes care of errors
{
	cl_mem buf=NULL; // default value to return during prognosis
	oclMem+=size;
	MAXIMIZE(oclMemPeak,oclMem);
	MAXIMIZE(oclMemMaxObj,size);
	if (!prognosis) {
		cl_int err;
		buf=clCreateBuffer(context,mem_flags,size,host_ptr,&err);
		if (err!=CL_SUCCESS)
			PrintCLErr(err,ERR_LOC_CALL,dyn_sprintf("Failed to allocate OpenCL object '%s'",name));
	}
	return buf;
}

//======================================================================================================================

void my_clReleaseBuffer(cl_mem buffer)
// wrapper to release buffer and decrease memory count
{
	cl_uint count;
	CL_CH_ERR(clGetMemObjectInfo(buffer,CL_MEM_REFERENCE_COUNT,sizeof(count),&count,NULL));
	if (count==1) {
		size_t size;
		CL_CH_ERR(clGetMemObjectInfo(buffer,CL_MEM_SIZE,sizeof(size),&size,NULL));
		if (oclMem<size) LogWarning(EC_WARN,ALL_POS,"Inconsistency detected in handling OpenCL memory: remaining "
			"memory (%zu) is smaller than size of the object to be freed (%zu)",oclMem,size);
		else oclMem-=size;
	}
	CL_CH_ERR(clReleaseMemObject(buffer));
}

//======================================================================================================================

void oclunload(void)
// unload all OpenCL kernels and similar stuff; buffers are released in Free_FFT_Dmat()
{
	CL_CH_ERR(clReleaseProgram(program));
	CL_CH_ERR(clReleaseKernel(clzero));
	CL_CH_ERR(clReleaseKernel(clarith1));
	CL_CH_ERR(clReleaseKernel(clarith2));
	CL_CH_ERR(clReleaseKernel(clarith3));
	CL_CH_ERR(clReleaseKernel(clarith4));
	CL_CH_ERR(clReleaseKernel(clarith5));
	CL_CH_ERR(clReleaseKernel(clnConj));
	CL_CH_ERR(clReleaseKernel(clinprod));
	CL_CH_ERR(clReleaseKernel(cltransposeof));
	CL_CH_ERR(clReleaseKernel(cltransposeob));
	if (surface) { // kernels for surface
		CL_CH_ERR(clReleaseKernel(clarith3_surface));
		CL_CH_ERR(clReleaseKernel(cltransposeofR));
	}
	CL_CH_ERR(clReleaseCommandQueue(command_queue));
	CL_CH_ERR(clReleaseContext(context));
}
