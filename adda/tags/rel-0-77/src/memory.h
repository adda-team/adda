/* FILE: memory.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        memory allocation and freeing
 *        also includes overflows checks
 *
 * Copyright (C) 2006-2007 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __memory_h
#define __memory_h

#include <stddef.h>    /* for size_t */
#include "function.h"  /* for function attributes */

#define MBYTE 1048576.0
/* for conciseness */
#define OTHER_ARGUMENTS const int who,const char *fname,const int line,const char *name

void CheckOverflow(double size,OTHER_ARGUMENTS);
size_t MultOverflow(size_t a,size_t b,OTHER_ARGUMENTS);
/* allocate */
doublecomplex *complexVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
double **doubleMatrix(size_t rows,size_t cols,OTHER_ARGUMENTS) ATT_MALLOC;
double *doubleVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
double *doubleVector2(size_t nl,size_t nh,OTHER_ARGUMENTS) ATT_MALLOC;
int **intMatrix (size_t nrl,size_t nrh,size_t ncl,size_t nch,OTHER_ARGUMENTS) ATT_MALLOC;
int *intVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
unsigned short *ushortVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
char *charVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
unsigned char *ucharVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
void *voidVector(size_t size,OTHER_ARGUMENTS) ATT_MALLOC;
/* free */
void Free_cVector(doublecomplex *v);
void Free_dMatrix(double **m,size_t rows);
void Free_dVector2(double *v,size_t nl);
void Free_iMatrix(int **m,size_t nrl,size_t nrh,size_t ncl);
void Free_general(void *v);

/* macros to use for allocation */
#define MALLOC_VECTOR(vec,type,size,who) vec=type##Vector(size,who,POSIT,#vec)
#define MALLOC_DVECTOR2(vec,nl,nh,who) vec=doubleVector2(nl,nh,who,POSIT,#vec)
#define MALLOC_DMATRIX(vec,rows,cols,who) vec=doubleMatrix(rows,cols,who,POSIT,#vec)
#define MALLOC_IMATRIX(vec,nrl,nrh,ncl,nch,who) vec=intMatrix(nrl,nrh,ncl,nch,who,POSIT,#vec)

#endif /*__memory_h*/
