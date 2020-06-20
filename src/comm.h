/* File: comm.h
 * $Date::                            $
 * Descr: definitions of communication global variables and routines
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
#ifndef __comm_h
#define __comm_h

// project headers
#include "types.h"    // needed for doublecomplex
#include "function.h" // for function attributes
#include "timing.h"   // for TIME_TYPE

// UOIP - Used Only In Parallel; to remove spurious 'unused' warnings in sequential mode
#ifdef PARALLEL
#	define UOIP
#else
#	define UOIP ATT_UNUSED
#endif

typedef enum {uchar_type,int_type,int3_type,sizet_type,double_type,double3_type,cmplx_type,cmplx3_type} var_type;

void Stop(int) ATT_NORETURN;
void Synchronize(void);
double AccumulateMax(double data,double *max);
void Accumulate(void * restrict data UOIP,const var_type type UOIP,size_t n UOIP,TIME_TYPE *timing UOIP);
void MyInnerProduct(void * restrict data,const var_type type,size_t n,TIME_TYPE *timing);
void InitComm(int *argc_p,char ***argv_p);
void ParSetup(void);
void SetupLocalD(void);
void MyBcast(void * restrict data,const var_type type,const size_t n_elem,TIME_TYPE *timing);
void BcastOrient(int *i,int *j,int *k);
void ReadField(const char * restrict fname,doublecomplex *restrict field);

#ifndef SPARSE
void BlockTranspose(doublecomplex * restrict X,TIME_TYPE *timing);
void BlockTranspose_DRm(doublecomplex * restrict X,size_t lengthY,size_t lengthZ);
// used by granule generator
void SetGranulComm(double z0,double z1,double gdZ,int gZ,size_t gXY,size_t buf_size,int *lz0,int *lz1,int sm_gr);
void CollectDomainGranul(unsigned char * restrict dom,size_t gXY,int lz0,int locgZ,TIME_TYPE *timing);
void FreeGranulComm(int sm_gr);
void ExchangeFits(bool * restrict data,const size_t n,TIME_TYPE *timing);
#endif // !SPARSE

#ifdef PARALLEL
// these functions are defined only in parallel mode
void CatNFiles(const char * restrict dir,const char * restrict tmpl,const char * restrict dest);
bool ExchangePhaseShifts(doublecomplex * restrict bottom, doublecomplex * restrict top,TIME_TYPE *timing);
void AllGather(void * restrict x_from,void * restrict x_to,var_type type,TIME_TYPE *timing);

/* The advantage of using this define is that compiler may remove an unnecessary test in sequential mode. The define do
 * not include common 'if', etc. to make the structure of the code (in the main text) immediately visible.
 */
#	define IFROOT (ringid==ADDA_ROOT)
#else
#	define IFROOT (true)
#endif

#endif // __comm_h
