/* File: oclgetdev.c
 * $Date::                            $
 * Descr: function to get an OpenCL device
 *
 * Copyright (C) 2010-2011 ADDA contributors
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

#include <stdio.h>
#include <string.h>
#include <CL/cl.h>

cl_int GetDevice(cl_platform_id* used_platform_id, cl_device_id* device_id, cl_int* platformname,cl_int devtype)
{
    cl_int err;
    cl_uint num_of_platforms;
    //little trick to get just the number of the Platforms
    err = clGetPlatformIDs(
                           0,
                           NULL,
                           &num_of_platforms
                           );
    if (err!=CL_SUCCESS) return err; 
                            
    // dynamic array of platformids creation at runtime, stored in heap   
    cl_platform_id* platform_id = (cl_platform_id *)malloc(num_of_platforms);
    char pname[1024];
    err = clGetPlatformIDs(
                           num_of_platforms,
                           platform_id,
                           NULL
                           );
    
    if (err!=CL_SUCCESS) return err; 
    cl_uint num_of_devices; 
    int i;
    for (i=0; i<num_of_platforms; i++)
    {
        clGetPlatformInfo(
                          platform_id[i],
                          CL_PLATFORM_NAME,
                          sizeof(pname),
                          pname,
                          NULL
                          );

        err = clGetDeviceIDs(
                             platform_id[i],
                             devtype,
                             1,
                             device_id,
                             &num_of_devices
                             );
        if (err==CL_SUCCESS)
        {
			if (strcmp(pname,"ATI Stream")==0)
				*platformname=1;
			if (strcmp(pname,"NVIDIA CUDA")==0)
				*platformname=0;
			*used_platform_id = platform_id[i];
			break; 
        }
    }
    free(platform_id);
    return err;
}
