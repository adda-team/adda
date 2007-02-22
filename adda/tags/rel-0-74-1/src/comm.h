/* FILE: comm.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of communication global variables
 *        and routines
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
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
void InitComm(int *argc_p,char ***argv_p);
void ParSetup(void);
void MyBcast(void *data,var_type type,int n_elem);
void BcastOrient(int *i,int *j,int *k);

#ifdef PARALLEL
/* this functions are defined only in parallel mode */
void CatNFiles(const char *dir,const char *tmpl,const char *dest);

/* analogs of frequently used functions that should be executed only by the ROOT processor */
/* not safe if used in constructions like { if PRINTZ(...); else } */
# define PRINTZ if (ringid==ROOT) printf
# define FPRINTZ if (ringid==ROOT) fprintf
# define SPRINTZ if (ringid==ROOT) sprintf
# define SYSTEMZ if (ringid==ROOT) system
# define FCLOSEZ if (ringid==ROOT) fclose
# define FFLUSHZ if (ringid==ROOT) fflush
# define PRINTBOTHZ if (ringid==ROOT) PrintBoth
#else
# define FCLOSEZ fclose
# define SYSTEMZ system
# define PRINTZ printf
# define FPRINTZ fprintf
# define SPRINTZ sprintf
# define FFLUSHZ fflush
# define PRINTBOTHZ PrintBoth
#endif


#endif /*__comm_h*/
