/* FILE : comm.c
 * AUTH : Martijn Frijlink 
 * DESCR: The main intention of this library is to incorporate all
 *        parallelisation related code, so most of it is directly
 *        involved in or closely related to inter-proces communication,
 *        hence its name. The only acception is the LogError-function,
 *        originally in its own file, which is merely defined here
 *        because I could not think of a better place (and wanted to
 *        reduce the number of source-files.
 *
 *        Currently is developed by Maxim Yurkin
 */
#include <stdio.h>
#include "stdarg.h"
#include <time.h>
#include "debug.h"
#include "cmplx.h"
#include "types.h"
#include "comm.h"
#include "const.h"

#if defined(MPI)
 #include <mpi.h>
 MPI_Datatype mpi_dcomplex;
#endif

int nprocs;

int local_z0,local_z1,local_Nz,local_x0,local_x1,local_Nx,local_z1_f,local_Nz_f;
int local_d0,local_d1,local_Ndip,Ndip;
int ringid;

char buffer1[40960];

/*===========================================*/

void init_comm(int *argc,char ***argv)
{
#if defined(MPI)
  int 
    dcmplx_blocklength[2]={1,1},dcmplx_displs[2] = {0,1};
  MPI_Datatype
    dcmplx_type[2];

  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&ringid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);

 
  /* Create MPI-type for sending dcomplex-variables */
  dcmplx_type[0] = MPI_DOUBLE; dcmplx_type[1] = MPI_DOUBLE;
  MPI_Type_struct(2,dcmplx_blocklength,dcmplx_displs,dcmplx_type,
                  &mpi_dcomplex);
  MPI_Type_commit(&mpi_dcomplex);
#elif !defined(PARALLEL)
  nprocs=1; ringid=0;
#endif
}

/*===========================================*/

void stop(int code)
{
  /* exit the program */
  printf("end %i\n",ringid);
#if defined(MPI)
  if (code)
    {
      printf("Aborting proces %d\n",ringid);
      MPI_Abort(MPI_COMM_WORLD,code);
    }
  else
    {
      MPI_Finalize();
    }
#endif
  exit(code);
}

/*===========================================*/

void synchronize(void)
{
#if defined(MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/*===========================================*/

void Bcast_orient(int *i, int *j, int *k)
{
#ifdef MPI
  int buf[3];

  if (ringid==0) {  
    buf[0]=*i; 
    buf[1]=*j;
    buf[2]=*k;
  }
  synchronize();
  MPI_Bcast(buf,3,MPI_INT,0,MPI_COMM_WORLD);
  
  if (ringid!=0) {
    *i=buf[0]; 
    *j=buf[1];
    *k=buf[2];
  }
#endif
}

/*===========================================*/

void Bcast_parms(Parms_1D *parms)
{
#ifdef MPI
  int root=0;  
  int
    len[6] = {1,1,1,1,1,1};
  MPI_Datatype
    mpi_parms,
    type[6];
  MPI_Aint
    disp[6];

  disp[0] = 0;                     type[0] = MPI_DOUBLE;
  disp[1] = sizeof(double);          type[1] = MPI_INT;
  disp[2] = disp[1] + sizeof(int); type[2] = MPI_INT;
  disp[3] = disp[2] + sizeof(int); type[3] = MPI_DOUBLE;
  disp[4] = disp[3] + sizeof(double);type[4] = MPI_DOUBLE;
  disp[5] = disp[4] + sizeof(double);type[5] = MPI_INT;

  MPI_Type_struct(6,len,disp,type,&mpi_parms);
  MPI_Type_commit(&mpi_parms);

  MPI_Bcast(parms,2,mpi_parms,root,MPI_COMM_WORLD);
#endif
}

/*===========================================*/
 
void accumulate(/* gather and add scattered fields on proces root */
                double *vector,int n)
{
#if defined(MPI)
  int root=0;
  int i;
  double *temp[GLOB_LOC];
 
  temp[LOCAL] = (double *) malloc(n*sizeof(double));
  temp[GLOBAL] = (double *) calloc(n,sizeof(double));
 
  memcpy(temp[LOCAL],vector,n*sizeof(double));
  MPI_Reduce(temp[LOCAL],temp[GLOBAL],n,MPI_DOUBLE,
             MPI_SUM,root,MPI_COMM_WORLD);
  memcpy(vector,temp[GLOBAL],n*sizeof(double));
 
  free(temp[LOCAL]);
  free(temp[GLOBAL]);
#endif
}

/*===========================================*/

void accumulate_new(/* gather and add scattered fields on proces root */
                double *vector,int n)
     /* we noticed some problems witht the MPI_Reduce call running in MPICH on the */
     /* UvA Blue Beowulf cluster. So here we replace MPI_Reduce with plain point- */
     /* to-point communication. */
{
  int root=0; /* processor 0 will accumulate everything */

  double *eper_buffer;
  int i,j;

#if defined(MPI)
  MPI_Status status;

  MPI_Barrier(MPI_COMM_WORLD);
  if (ringid!=0) {
    MPI_Send(vector,n,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
  }
  else {
    if ((eper_buffer = (double *)malloc(n*sizeof(double))) == NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc eper_buffer");
    for(i=1; i<nprocs; i++){
      MPI_Recv(eper_buffer,n,MPI_DOUBLE,i,0,MPI_COMM_WORLD,&status);
      for(j=0;j<n;j++) vector[j]+=eper_buffer[j];
    }
    free(eper_buffer);
  }     
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/*===========================================*/

void my_inner_product(/* gather local inproducts of a distributed vector or
			 sum, stored in *a, add them and return them in *a */
                      double *a)
{
#if defined(MPI)
  int root=0;
  double temp[GLOB_LOC];
 
  temp[LOCAL] = *a;
  MPI_Allreduce(&temp[LOCAL],&temp[GLOBAL],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  *a = temp[GLOBAL];
#endif
}

/*===========================================*/

void block_transpose(/* Do the data-transposition, i.e. exchange,
                        between fftX and fftY&fftZ */
                      double *X, int nnn)
{
  int tstart;  
  extern clock_t Timing_OneIterComm;  
#if defined(MPI)  
  double
    *buffer, *rbuffer;
  int
    bufsize,
    x,y,z,
    transmission,part;
  extern int smallY,smallZ;
  int position;
  MPI_Status status;
  bufsize = 2*nnn*local_Nz_f*smallY*local_Nx; 
  buffer = (double *) malloc(bufsize*sizeof(double));
  rbuffer = (double *) malloc(bufsize*sizeof(double));
#endif
  tstart=clock();
#if defined(MPI)
  for(transmission=1;transmission<nprocs;transmission++)
    {
      part=ringid ^ transmission;
      position = 0;
      for(z=nnn*local_z0;z<nnn*local_z1_f;z++) for(y=0;y<smallY;y++) {
        memcpy(&buffer[position],
               &X[2*index_Xmatrix(local_Nx*part,y,z,nnn)],
               2*local_Nx*sizeof(double));
        if ( (position += 2*local_Nx) > bufsize)
          LogError(EC_ERROR,ALL,POSIT,"(transm=%d) buffer-overflow: bufsize=%d position=%d",transmission,bufsize,position);
        }
      /* it looks like that the original code by Martijn, see below, may block in some situations, e.g. -grid 160 on */
      /* 32 processors. So, lets see what happens if we change all this with MPI_Sendrecv routines */
     /*  MPI_Send(buffer,bufsize,MPI_DOUBLE,part,0,MPI_COMM_WORLD); */
     /* MPI_Barrier(MPI_COMM_WORLD); */
 
      /* MPI_Recv(buffer,bufsize,MPI_DOUBLE,part,0,MPI_COMM_WORLD,&status); */

      MPI_Sendrecv(buffer, bufsize, MPI_DOUBLE, part, 0,
		  rbuffer, bufsize, MPI_DOUBLE, part, 0,
		  MPI_COMM_WORLD,&status);
  
      position = 0;
      for(z=nnn*local_z0;z<nnn*local_z1_f;z++) for(y=0;y<smallY;y++) {
        memcpy(&X[2*index_Xmatrix(local_Nx*part,y,z,nnn)],
               &rbuffer[position],
               2*local_Nx*sizeof(double));
        if ( (position += 2*local_Nx) > bufsize)
          LogError(EC_ERROR,ALL,POSIT,"(transm=%d) buffer-overflow: bufsize=%d position=%d",transmission,bufsize,position);
      }
      MPI_Barrier(MPI_COMM_WORLD);

      /* free(buffer);*/
    }
  free(buffer);
  free(rbuffer);
#endif
  Timing_OneIterComm += clock() - tstart;
}

/*===========================================*/

void par_setup(void)
{
  extern int boxX,boxY,boxZ,smallZ,smallX;
  int unitZ,unitX;

#ifdef PARALLEL
  unitZ=(float) smallZ/nprocs+.99999;
  local_z0=ringid*unitZ;
  local_z1_f=(ringid+1)*unitZ;
  if (local_z1_f > boxZ) local_z1=boxZ;
  else local_z1=local_z1_f;
  unitX=(float) smallX/nprocs+.99999;
  local_x0=ringid*unitX;
  local_x1=(ringid+1)*unitX;
#else
  local_z0=0;
  local_z1_f=smallZ;
  local_z1=boxZ;
  local_x0=0;
  local_x1=smallX;
#endif
  if (local_z1<=local_z0) {
    LogError(EC_WARN,ALL,POSIT,"No real dipoles are assigned");
    local_z1=local_z0;
  }
  local_Nz_f=local_z1_f-local_z0;
  local_Nz=local_z1-local_z0;
  local_Nx=local_x1-local_x0;
  local_d0=boxX*boxY*local_z0;
  local_d1=boxX*boxY*local_z1;
  local_Ndip=local_d1-local_d0;
  Ndip=boxX*boxY*boxZ;
  printf("%i :  %i %i %i %i %i %i\n",ringid,local_z0,local_z1,local_z1_f,local_Ndip,Ndip,local_Nx);
}

/*===========================================*/

void all_gather(void *x_from,void *x_to,
		char data_type[],
		int n_elem)       
     /* Gather distributed arrays */
{
#if defined(MPI)
  MPI_Datatype mes_type;

  if (strcmp(data_type,"int")==0)
    mes_type = MPI_INT;
  if (strcmp(data_type,"double")==0)
    mes_type = MPI_DOUBLE;
  if (strcmp(data_type,"doublecomplex")==0)
    {
      mes_type = MPI_DOUBLE;
      n_elem *= 2;
    }

  MPI_Allgather(x_from,n_elem,mes_type,
		x_to,n_elem,mes_type,MPI_COMM_WORLD);
#endif
}

/*============================================*/

/* LogError
 * we use sprintf a couple of times, because we want each node to
 * generate an atomic message, not a couple of messages after
 * each other, since other nodes may then interfere with our
 * output!
 */
void LogError(int code, int who, char *fname, int lineN, char *fmt, ... )
   /* performs output of error specified by code at fname:lineN                  */
   /* fmt + args (...) specify error message                                     */
   /* who specifies whether 1 (ringid=0) or all processors should produce output */
   /* if code is EC_ERROR program aborts after output                            */
{
  va_list args;
  char line[255];
  char id_str[30];
 
  if (who==ALL || ringid==0) {
    va_start(args, fmt);
    
    strcpy(id_str,"");
    #ifdef PARALLEL 
    if (who==ALL) sprintf(id_str," - ringID=%d",ringid);
    #endif

    if (code==EC_ERROR) strcpy(line,"ERROR:");
    else if (code==EC_WARN) strcpy(line,"WARNING:");
    else sprintf(line,"Error code=%d:",code);
    sprintf(line+strlen(line)," (%s:%d%s) ", fname, lineN, id_str);
    vsprintf(line+strlen(line), fmt, args);
    strcat(line,"\n");
    fprintf(stderr, line);
    fflush(stderr);
    va_end(args);
  }
  if (code==EC_ERROR) stop(1);
}