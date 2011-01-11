/* File: oclvars.h
 * $Author: yurkin $
 * $Date:: 2010-09-30 19:52:58 +0200 #$
 * Descr: initialization and unload of all OpenCL variables
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * Copyright (C) 2009,2010 Institute of Chemical Kinetics and Combustion & University of Amsterdam
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

#include<stdio.h>
#include<CL/cl.h>
#include "cpp/clFFT.h"

cl_int err; //initialize error code variable for OpenCL
cl_device_id device_id; //create device and platform variables 
cl_platform_id used_platform_id;

cl_int devtype=CL_DEVICE_TYPE_GPU; //set preffered device type
//char * coptions="-DAMD"; //for AMD GPUs
char * coptions="";

cl_context context;
cl_command_queue command_queue;
cl_kernel clzero, clarith1,clarith2,clarith3, clarith4, clarith5, clnConj, clinprod, cltransposef, cltransposeb;
cl_program program;
cl_mem bufXmatrix, bufmaterial, bufposition, bufcc_sqrt, bufargvec, bufresultvec, bufslices, bufslices_tr, bufDmatrix;




void oclinit(void)
{
    err = GetDevice(
            &used_platform_id,
            &device_id,
            devtype
            );
    checkErr(err, "No Valid prefered OpenCL device found");
    //cl_context_properties is a strange list of item:
    //first comes the name of the property as next element the corrosponding value
    cl_context_properties properties[3];
    properties[0]=CL_CONTEXT_PLATFORM;
    properties[1]=(cl_context_properties)used_platform_id;
    properties[2]=0; //last one must be zero
    context = clCreateContext(
            properties,
            1, //number of devices to use
            &device_id,
            NULL, //errorinfo as string
            NULL, //userdata
            &err
            );

    checkErr(err, "Create Context");
    command_queue = clCreateCommandQueue(
            context,
            device_id,
            CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE, //command queue properties
            &err
            ); 
    checkErr(err, "Command_queue");

    FILE* file= fopen("oclkernels.cl","r");
    size_t sourceStrSize;
    fseek(file, 0, SEEK_END);
    sourceStrSize = ftell(file);
    fseek(file, 0, SEEK_SET);
    char* cSourceString = (char *)malloc(sourceStrSize+1);
    fread(cSourceString, sourceStrSize,1, file);
    program = clCreateProgramWithSource(
            context, //valid Context
            1, //number of strings in next parameter
            (const char **) &cSourceString, //array of source strings
            NULL, //length of each string
            &err
            );

    checkErr(err, "building program");
    if (clBuildProgram(program, 0, NULL, coptions, NULL, NULL) != CL_SUCCESS)
    {
        printf("Error building program\n");
        char buffer[8192];
        size_t length=8192;
        clGetProgramBuildInfo(
                program,              // valid program object
                device_id,            // valid device_id that executable was built
                CL_PROGRAM_BUILD_LOG, // indicate to retrieve build log
                sizeof(buffer),       // size of the buffer to write log to
                buffer,               // the actual buffer to write log to
                &length);             // the actual size in bytes of data copied to buffer
        printf("%s\n",buffer);
        exit(1);
    }

    clzero = clCreateKernel( program, "clzero", &err);
    checkErr(err, "creating kernel clzero");

    clarith1 = clCreateKernel( program, "arith1", &err);
    checkErr(err, "creating kernel clarith1");
    
    clarith2 = clCreateKernel( program, "arith2", &err);
    checkErr(err, "creating kernel clartih2");

    clarith3 = clCreateKernel( program, "arith3", &err);
    checkErr(err, "creating kernel clartih3");
    
    clarith4 = clCreateKernel( program, "arith4", &err);
    checkErr(err, "creating kernel clarith4");
    
    clarith5 = clCreateKernel( program, "arith5", &err);
    checkErr(err, "creating kernel clarith5");


    clnConj = clCreateKernel( program, "nConj", &err);
    checkErr(err, "creating kernel clnConj");
    
    clinprod = clCreateKernel( program, "inpr", &err);
    checkErr(err, "creating kernel clinprod");
    

    cltransposef = clCreateKernel( program, "transpose", &err);
    checkErr(err, "creating kernel cltransposef");
    
    cltransposeb = clCreateKernel( program, "transpose", &err);
    checkErr(err, "creating kernel cltransposeb");
    
    free(cSourceString);
}

void oclunload(void)
{
    clReleaseMemObject(bufXmatrix);
    clReleaseMemObject(bufDmatrix);
    clReleaseMemObject(bufmaterial);
    clReleaseMemObject(bufposition);
    clReleaseMemObject(bufcc_sqrt);
    clReleaseMemObject(bufargvec);
    clReleaseMemObject(bufresultvec);
    clReleaseMemObject(bufslices);
    clReleaseMemObject(bufslices_tr);

    clReleaseProgram(program);
    clReleaseKernel(cltransposef);
    clReleaseKernel(cltransposeb);
    clReleaseKernel(clzero);
    clReleaseKernel(clarith1);
    clReleaseKernel(clarith2);
    clReleaseKernel(clarith3);
    clReleaseKernel(clarith4);
    clReleaseKernel(clarith5);
    clReleaseCommandQueue(command_queue);
    clReleaseContext(context);
}
