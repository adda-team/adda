/* FILE: memory.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        memory allocation and freeing
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __memory_h
#define __memory_h

#include "function.h"  /* for function attributes */

#define MBYTE 1048576.0
/* for conciseness */
#define OTHER_ARGUMENTS const int who,const char *fname,const int line,const char *name

/* allocate */
doublecomplex *complexVector(int size,OTHER_ARGUMENTS) ATT_MALLOC;
double **doubleMatrix(int rows,int cols,OTHER_ARGUMENTS) ATT_MALLOC;
double *doubleVector(int size,OTHER_ARGUMENTS) ATT_MALLOC;
double *doubleVector2(int nl,int nh,OTHER_ARGUMENTS) ATT_MALLOC;
int **intMatrix (int nrl,int nrh,int ncl,int nch,OTHER_ARGUMENTS) ATT_MALLOC;
int *intVector(int size,OTHER_ARGUMENTS) ATT_MALLOC;
unsigned short *ushortVector(int size,OTHER_ARGUMENTS) ATT_MALLOC;
char *charVector(int size,OTHER_ARGUMENTS) ATT_MALLOC;
unsigned char *ucharVector(int size,OTHER_ARGUMENTS) ATT_MALLOC;
void *voidVector(int size,OTHER_ARGUMENTS) ATT_MALLOC;
/* free */
void Free_cVector(doublecomplex *v);
void Free_dMatrix(double **m,int rows);
void Free_dVector2(double *v,int nl);
void Free_iMatrix(int **m,int nrl,int nrh,int ncl);
void Free_general(void *v);

/* macros to use for allocation */
#define MALLOC_VECTOR(vec,type,size,who) vec=type##Vector(size,who,POSIT,#vec)
#define MALLOC_DVECTOR2(vec,nl,nh,who) vec=doubleVector2(nl,nh,who,POSIT,#vec)
#define MALLOC_DMATRIX(vec,rows,cols,who) vec=doubleMatrix(rows,cols,who,POSIT,#vec)
#define MALLOC_IMATRIX(vec,nrl,nrh,ncl,nch,who) vec=intMatrix(nrl,nrh,ncl,nch,who,POSIT,#vec)

#endif /*__memory_h*/
