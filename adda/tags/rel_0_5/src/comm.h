#ifndef __comm_h
#define __comm_h

#define LOCAL 0
#define GLOBAL 1
#define GLOB_LOC 2
 
#if defined(MPI)
#define PARALLEL
#endif

extern int local_d0,local_d1,local_Ndip,Ndip;
extern int local_z0,local_z1,local_Nz,local_x0,local_x1,local_Nx,local_z1_f,local_Nz_f;

extern int nprocs;

extern int ringid;

void LogError(int ErrCode, int who, char *FName, int line, char *Format,...);

void init_comm(int *argc,char ***argv);
void stop(int);
void par_setup(void);
void synchronize(void);

void Bcast_parms(Parms_1D *parms);
void Bcast_orient(int *i,int *j,int *k);

void block_transpose(double *X, int nnn);
void accumulate(double *,int);
void my_inner_product(double *a);

void all_gather(void *x_from,void *x_to,char data_type[],int n_elem);
void all_gather_dcomplex(void *x_from,void *x_to,int n_elem);
void all_gather_int(void *x_from,void *x_to,int n_elem);
void all_gather_REAL(void *x_from,void *x_to,int n_elem);

#ifdef PARALLEL
#define printz if (ringid==0) printf
#define fprintz if (ringid==0) fprintf
#define sprintz if (ringid==0) sprintf
#define systemz if (ringid==0) system
#define fclosez if (ringid==0) fclose
#else
#define fclosez fclose
#define systemz system
#define printz printf
#define fprintz fprintf
#define sprintz sprintf
#endif


#endif /*__comm_h*/