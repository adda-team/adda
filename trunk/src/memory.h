/* File: memory.h
 * $Date::                            $
 * Descr: definitions of functions for memory allocation and freeing; also includes overflows checks
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
#ifndef __memory_h
#define __memory_h

// project headers
#include "const.h"    // for enum types
#include "function.h" // for function attributes
#include "io.h"
#include "types.h"    // for doublecomplex
// system headers
#include <stddef.h>   // for size_t

#define MBYTE 1048576.0
#define FFORMM "%.1f" // format for memory footprint (estimates) in MB
// for conciseness
#define OTHER_ARGUMENTS ERR_LOC_DECL,const char *name

void CheckOverflow(double size,OTHER_ARGUMENTS);
size_t MultOverflow(size_t a,size_t b,OTHER_ARGUMENTS);
// allocate
doublecomplex *complexVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
double **doubleMatrix(size_t rows,size_t cols,OTHER_ARGUMENTS) ATT_MALLOC;
double *doubleVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
double *doubleVector2(size_t nl,size_t nh,OTHER_ARGUMENTS) ATT_MALLOC;
int **intMatrix (size_t nrl,size_t nrh,size_t ncl,size_t nch,OTHER_ARGUMENTS) ATT_MALLOC;
int *intVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
unsigned short *ushortVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
char *charVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
unsigned char *ucharVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
bool *boolVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
size_t *sizetVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
void *voidVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
// reallocate; only two for now, more can be easily added
double *doubleRealloc(double *ptr,const size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
char *charRealloc(char *ptr,const size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
// free
void Free_cVector(doublecomplex * restrict v);
void Free_dMatrix(double ** restrict m,size_t rows);
void Free_dVector2(double * restrict v,size_t nl);
void Free_iMatrix(int ** restrict m,size_t nrl,size_t nrh,size_t ncl);
void Free_general(void * restrict v);

// macros to use for allocation and reallocation
#define MALLOC_VECTOR(vec,type,size,who) vec=type##Vector(size,who,POSIT,#vec)
#define MALLOC_DVECTOR2(vec,nl,nh,who) vec=doubleVector2(nl,nh,who,POSIT,#vec)
#define MALLOC_DMATRIX(vec,rows,cols,who) vec=doubleMatrix(rows,cols,who,POSIT,#vec)
#define MALLOC_IMATRIX(vec,nrl,nrh,ncl,nch,who) vec=intMatrix(nrl,nrh,ncl,nch,who,POSIT,#vec)
#define REALLOC_VECTOR(vec,type,size,who) vec=type##Realloc(vec,size,who,POSIT,#vec)

#endif //__memory_h
