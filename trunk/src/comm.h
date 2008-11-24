/* FILE: comm.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of communication global variables
 *        and routines
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __comm_h
#define __comm_h

#include "types.h"    // needed for doublecomplex
#include "function.h" // for function attributes
#include "timing.h"   // for TIME_TYPE

typedef enum {char_type,int_type,double_type,cmplx_type} var_type;

void Stop(int) ATT_NORETURN;
void Synchronize(void);
void BlockTranspose(doublecomplex *X);
void BlockTranspose_Dm(doublecomplex *X,size_t lengthY,size_t lengthZ);
void AccumulateMax(double *data,double *max);
void Accumulate(double *data,size_t size,double *buffer,TIME_TYPE *timing);
void MyInnerProduct(void *a,var_type type,size_t n_elem,TIME_TYPE *timing);
void AllGather(void *x_from,void *x_to,var_type type,size_t n_elem);
void InitComm(int *argc_p,char ***argv_p);
void ParSetup(void);
void MyBcast(void *data,var_type type,size_t n_elem,TIME_TYPE *timing);
void BcastOrient(int *i,int *j,int *k);
// used by granule generator
void SetGranulComm(double z0,double z1,double gdZ,int gZ,size_t gXY,size_t buf_size,int *lz0,
                   int *lz1,int sm_gr);
void CollectDomainGranul(unsigned char *dom,size_t gXY,int lz0,int locgZ,TIME_TYPE *timing);
void FreeGranulComm(int sm_gr);
void ExchangeFits(char *data,const size_t n,TIME_TYPE *timing);

#ifdef PARALLEL
// this functions are defined only in parallel mode
void CatNFiles(const char *dir,const char *tmpl,const char *dest);

/* analogs of frequently used functions that should be executed only by the ROOT processor
 * !!! not safe if used in constructions like { if (...) PRINTZ(...); else }
 */
#	define PRINTZ if (ringid==ROOT) printf
#	define FPRINTZ if (ringid==ROOT) fprintf
#	define SPRINTZ if (ringid==ROOT) sprintf
#	define STRCPYZ if (ringid==ROOT) strcpy
#	define FCLOSEZ if (ringid==ROOT) fclose
#	define FFLUSHZ if (ringid==ROOT) fflush
#	define PRINTBOTHZ if (ringid==ROOT) PrintBoth
#else
#	define PRINTZ printf
#	define FPRINTZ fprintf
#	define SPRINTZ sprintf
#	define STRCPYZ strcpy
#	define FCLOSEZ fclose
#	define FFLUSHZ fflush
#	define PRINTBOTHZ PrintBoth
#endif

#endif // __comm_h
