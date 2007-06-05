/* FILE: comm.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of communication global variables
 *        and routines
 */
#ifndef __comm_h
#define __comm_h

/* ringid of root processor */
#define ROOT 0

#if defined(MPI)
#define PARALLEL
#endif

#include "types.h"   /* needed to define Parms_1D */

typedef enum {char_type,int_type,double_type,cmplx_type} var_type;

/* initialized in init_comm (comm.c) */
extern int nprocs,ringid;
/* initialized in par_setup (comm.c) */
extern int local_d0,local_d1,local_Ndip,Ndip,local_z0,local_z1,local_Nz,
           local_x0,local_x1,local_Nx,local_z1_f,local_Nz_f;
/* initialized in make_particle */
extern int local_nvoid_Ndip,nvoid_Ndip;


void LogError(int ErrCode, int who, char *FName, int line, char *Format,...);

void init_comm(int *argc,char ***argv);
void stop(int);
void par_setup(void);
void synchronize(void);

void Bcast_parms(Parms_1D *parms);
void Bcast_orient(int *i,int *j,int *k);

void block_transpose(doublecomplex *X);
void block_transpose_Dm(doublecomplex *X,int lengthY,int lengthZ);
void AccumulateMax(double *data,double *max);
void accumulate(double *data,int size,double *buffer);
void my_inner_product(void *a,var_type type,int n_elem);

void all_gather(void *x_from,void *x_to,var_type type,int n_elem);

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
