/* FILE: comm.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of communication global variables
 *        and routines
 */
#ifndef __comm_h
#define __comm_h

#include "cmplx.h"  /* needed for doublecomplex */

/* Using ROOT!=0 should work, however it was not thoroughly tested.
   Hence do not change without necessity */
#define ROOT 0  /* ringid of root processor */

#if defined(MPI)
#define PARALLEL
#endif

typedef enum {char_type,int_type,double_type,cmplx_type} var_type;

void Stop(int);
void Synchronize(void);
void BlockTranspose(doublecomplex *X);
void BlockTranspose_Dm(doublecomplex *X,int lengthY,int lengthZ);
void AccumulateMax(double *data,double *max);
void Accumulate(double *data,int size,double *buffer);
void MyInnerProduct(void *a,var_type type,int n_elem);
void AllGather(void *x_from,void *x_to,var_type type,int n_elem);
void InitComm(int argc,char **argv);
void ParSetup(void);
void BcastOrient(int *i,int *j,int *k);

#ifdef PARALLEL
#define printz if (ringid==ROOT) printf
#define fprintz if (ringid==ROOT) fprintf
#define sprintz if (ringid==ROOT) sprintf
#define systemz if (ringid==ROOT) system
#define fclosez if (ringid==ROOT) fclose
#define fflushz if (ringid==ROOT) fflush
#else
#define fclosez fclose
#define systemz system
#define printz printf
#define fprintz fprintf
#define sprintz sprintf
#define fflushz fflush
#endif


#endif /*__comm_h*/
