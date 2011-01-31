/* File: oclvars.h
 * $Date::                            $
 * Descr: initialization and unload of all OpenCL variables
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

#ifndef __oclvars_h
#define __oclvars_h
extern cl_int err;
extern cl_device_id device_id;
extern cl_platform_id used_platform_id;
extern cl_context context;
extern cl_command_queue command_queue;
extern cl_kernel clzero, clarith1,clarith2,clarith3, clarith4, clarith5, clnConj, clinprod, cltransposef, cltransposeb;
extern cl_mem bufXmatrix, bufmaterial, bufposition, bufcc_sqrt, bufargvec, bufresultvec, bufslices, bufslices_tr, bufDmatrix;
extern cl_program program;

extern void oclinit(void);
extern void oclunload(void);
#endif
