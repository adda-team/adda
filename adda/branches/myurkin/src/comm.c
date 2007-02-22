/* File : comm.c
   Author : Martijn Frijlink 
   Descr : The main intention of this library is to incorporate all
           parallelisation related code, so most of it is directly 
           involved in or closely related to inter-proces communication,
           hence its name. The only acception is the LogError-function,
	   originally in its own file, which is merely defined here 
	   because I could not think of a better place (and wanted to 
	   reduce the number of source-files.
*/
#include <stdio.h>
#include "stdarg.h"
#include <sys/time.h>
#include "debug.h"
#include "cmplx.h"
#include "types.h"
#include "comm.h"

#if defined(PVM)
#include <pvm3.h>
#elif defined(MPI)
#include <mpi.h>
#endif

#if defined(MPI)
MPI_Datatype mpi_dcomplex;
MPI_Datatype mpi_REAL;
#endif

int nprocs=1;
int mytid=0;
int tids[256];            /* array of task id */
int me;

int *all_nodes;
int __fake_nnodes;
int total_nodes;
static int myid;

int local_z0,local_z1,local_Nz;
int local_d0,local_d1,local_Ndip,Ndip;
int ringid;
int local_xs_unit;

char buffer1[40960];

extern unsigned long extime(void);

#ifndef PowerPVM
int TimeNow(void) {return(clock());}
#endif

void init_comm(int *argc,char ***argv)
{
#if defined(PVM)
  int i,info,narch;
  struct pvmhostinfo *hostp;
 
  mytid = pvm_mytid();
  info = pvm_config( &nprocs, &narch, &hostp );

  me = pvm_joingroup( "foo" );
  if (pvm_parent() == PvmNoParent){
     me = 0;
   }
  else me=1;
 
  printf("Nprocs= %i me = %d mytid = %d\n",nprocs,me,mytid);
 
  if( me == 0 ) {
    pvm_spawn((*argv)[0], &(*argv)[1], 0, "", nprocs-1, &tids[1]);
    tids[0]=mytid;
  }
 
  if (me==0) {
    pvm_initsend(PvmDataDefault);
    pvm_pkint(&nprocs, 1, 1);
    pvm_pkint(tids, nprocs, 1);
    pvm_mcast(&tids[1], nprocs-1, 0);
  }
  else {
    pvm_recv( -1, 0);
    pvm_upkint(&nprocs, 1, 1);
    pvm_upkint(tids, nprocs, 1);
  }
  for(i=0;i<nprocs;i++) if (tids[i]==mytid) ringid=i;
#elif defined(MPI)
  int 
    dcmplx_blocklength[2]={1,1},dcmplx_displs[2] = {0,1};
  MPI_Datatype
    dcmplx_type[2];

  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&ringid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  me = (ringid==0) ? 0 : 1; 
  mytid = 0;

  /* Create MPI-type for sending REAL's */ 
#if SINGLE
  mpi_REAL = MPI_FLOAT;
#elif DOUBLE
  mpi_REAL = MPI_DOUBLE;
#endif
 
  /* Create MPI-type for sending dcomplex-variables */
  dcmplx_type[0] = mpi_REAL; dcmplx_type[1] = mpi_REAL;
  MPI_Type_struct(2,dcmplx_blocklength,dcmplx_displs,dcmplx_type,
                  &mpi_dcomplex);
  MPI_Type_commit(&mpi_dcomplex);
#elif !defined(PARALLEL)
  nprocs=1; ringid=0;
  me=0; mytid=0;
#endif
}

void stop(int code)
{
  /* exit the program */
  printf("end %i\n",ringid);
#if defined(PVM)
  pvm_lvgroup( "foo" );
  pvm_exit();
#elif defined(MPI)
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

void synchronize(void)
{
#if defined(PVM)
  pvm_barrier("foo",nprocs);
#elif defined(MPI)
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void Bcast_parms(Parms_1D *parms)
{
  int root=0;
#ifdef MPI
  int
    len[6] = {1,1,1,1,1,1};
  MPI_Datatype
    mpi_parms,
    type[6];
  MPI_Aint
    disp[6];

  disp[0] = 0;                     type[0] = mpi_REAL;
  disp[1] = sizeof(REAL);          type[1] = MPI_INT;
  disp[2] = disp[1] + sizeof(int); type[2] = MPI_INT;
  disp[3] = disp[2] + sizeof(int); type[3] = mpi_REAL;
  disp[4] = disp[3] + sizeof(REAL);type[4] = mpi_REAL;
  disp[5] = disp[4] + sizeof(REAL);type[5] = MPI_INT;

  MPI_Type_struct(6,len,disp,type,&mpi_parms);
  MPI_Type_commit(&mpi_parms);

  MPI_Bcast(parms,2,mpi_parms,root,MPI_COMM_WORLD);
#endif
}
 
void accumulate(/* gather and add scattered fields on proces root */
                REAL *vector,int n)
{
  int root=0;
#if defined(PVM)
  REAL *eper_buffer;
  int i,j;
 
  if ((eper_buffer = (REAL *)malloc(n*sizeof(REAL))) == NULL) {
    LogError (EC_ERROR,"comm.c",
              "processor %d, ringID %d, could not malloc Eperb",
              MyProcId, RingId);
    stop(1);
  }
 
  /* send data to proc root */
  pvm_barrier("foo", nprocs );
  if (ringid!=root) {
    pvm_initsend(PvmDataDefault);
    pvm_pkfloat(vector,n,1);
    pvm_send(tids[root],1);
  }
  else { /* wait for data */
    for( i=0 ; i<nprocs-1 ; i++ ) {
      pvm_recv( -1, 1);
      pvm_upkfloat(eper_buffer,n,1);
      for(j=0;j<n;j++) vector[j]+=eper_buffer[j];
    }
  }
 
  free(eper_buffer);
#elif defined(MPI)
  int i;
  REAL *temp[GLOB_LOC];
 
  temp[LOCAL] = (REAL *) malloc(n*sizeof(REAL));
  temp[GLOBAL] = (REAL *) calloc(n,sizeof(REAL));
 
  memcpy(temp[LOCAL],vector,n*sizeof(REAL));
  MPI_Reduce(temp[LOCAL],temp[GLOBAL],n,mpi_REAL,
             MPI_SUM,root,MPI_COMM_WORLD);
  memcpy(vector,temp[GLOBAL],n*sizeof(REAL));
 
  free(temp[LOCAL]);
  free(temp[GLOBAL]);
#endif
}

void accumulate_new(/* gather and add scattered fields on proces root */
                REAL *vector,int n)
     /* we noticed some problems witht the MPI_Reduce call running in MPICH on the */
     /* UvA Blue Beowulf cluster. So here we replace MPI_Reduce with plain point- */
     /* to-point communication. */
{
  int root=0; /* processor 0 will accumulate everything */

  REAL *eper_buffer;
  int i,j;

 #if defined(MPI)
  MPI_Status status;
#endif
 
 #if defined(PVM)
  /* send data to proc root */
  pvm_barrier("foo", nprocs );
  if (ringid!=root) {
    pvm_initsend(PvmDataDefault);
    pvm_pkfloat(vector,n,1);
    pvm_send(tids[root],1);
  }
  else { /* wait for data */
  if ((eper_buffer = (REAL *)malloc(n*sizeof(REAL))) == NULL) {
    LogError (EC_ERROR,"comm.c",
              "processor %d, ringID %d, could not malloc Eperb",
              MyProcId, RingId);
    stop(1);
  }
    for( i=0 ; i<nprocs-1 ; i++ ) {
      pvm_recv( -1, 1);
      pvm_upkfloat(eper_buffer,n,1);
      for(j=0;j<n;j++) vector[j]+=eper_buffer[j];
    }
 free(eper_buffer);
  }
 
#elif defined(MPI)
  MPI_Barrier(MPI_COMM_WORLD);
  if (ringid!=0) {
    MPI_Send(vector,n,mpi_REAL,0,0,MPI_COMM_WORLD);
  }
  else {
  if ((eper_buffer = (REAL *)malloc(n*sizeof(REAL))) == NULL) {
    LogError (EC_ERROR,"comm.c",
              "processor %d, ringID %d, could not malloc Eperb",
              MyProcId, RingId);
    stop(1);
  }
    for(i=1; i<nprocs; i++){
      MPI_Recv(eper_buffer,n,mpi_REAL,i,0,MPI_COMM_WORLD,&status);
      for(j=0;j<n;j++) vector[j]+=eper_buffer[j];
    }
 free(eper_buffer);
  }     
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void my_inner_product(/* gather local inproducts of a distributed vector or
			 sum, stored in *a, add them and return them in *a */
                      double *a)
{
  int root=0;
#if defined(PVM)
  double temp;
  int i;
 
  pvm_barrier("foo",nprocs);
 
  if (ringid!=root) {
    pvm_initsend(PvmDataDefault);
    pvm_pkdouble(a,1,1);
    pvm_send(tids[root],24);
  }
  else { /* wait for data */
    for( i=0 ; i<nprocs-1 ; i++ ) {
      pvm_recv( -1, 24);
      pvm_upkdouble(&temp,1,1);
      *a+=temp;
    }
  }
  /* processor root holds a; broadcast to all others */
  if (ringid==root) {
    pvm_initsend(PvmDataDefault);
    pvm_pkdouble(a,1,1);
    pvm_mcast(tids,nprocs,25);
  }
  else {
    pvm_recv(tids[root],25);
    pvm_upkdouble(a,1,1);
  }
#elif defined(MPI)
  double temp[GLOB_LOC];
 
  temp[LOCAL] = *a;
  MPI_Allreduce(&temp[LOCAL],&temp[GLOBAL],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  *a = temp[GLOBAL];
#endif
}

void block_transpose(/* Do the data-transposition, i.e. exchange,
                        between fftX and fftY&fftZ */
                      REAL *X)
{
  REAL
    *buffer, *rbuffer;
  int
    bufsize,
    x,y,z,
    transmission,part;
  extern int smallY,smallZ;
  extern unsigned long Timing_OneIterComm;
  int tstart,tstop;
#if defined(MPI)
  int position;
  MPI_Status status;
  bufsize = 2*local_Nz*smallY*local_xs_unit; 
  buffer = (REAL *) malloc(bufsize*sizeof(REAL));
  rbuffer = (REAL *) malloc(bufsize*sizeof(REAL));
#endif
 
  tstart=extime();

  for(transmission=1;transmission<nprocs;transmission++)
    {
      part=ringid ^ transmission;
#if defined(PVM)
      pvm_initsend(PvmDataDefault);
 
      for(z=local_z0;z<local_z1;z++) for(y=0;y<smallY;y++) {
        buffer = &X[2*index_Xmatrix(local_xs_unit*part,y,z)];
        pvm_pkfloat(buffer,2*local_xs_unit,1);
        }
      pvm_send(tids[part],200+transmission);
 
      pvm_barrier("foo",nprocs);
 
      pvm_recv(tids[part],200+transmission);
      for(z=local_z0;z<local_z1;z++) for(y=0;y<smallY;y++) {
        buffer = &X[2*index_Xmatrix(local_xs_unit*part,y,z)];
        pvm_upkfloat(buffer,2*local_xs_unit,1);
      }
 
      pvm_barrier("foo",nprocs);
#elif defined(MPI)
 
      position = 0;
      for(z=local_z0;z<local_z1;z++) for(y=0;y<smallY;y++) {
        memcpy(&buffer[position],
               &X[2*index_Xmatrix(local_xs_unit*part,y,z)],
               2*local_xs_unit*sizeof(REAL));
        if ( (position += 2*local_xs_unit) > bufsize)
          printf("block_transpose, ringid=%d transm=%d\n buffer-overflow: bufsize=%d position=%d",ringid,transmission,bufsize,position);
        }
      /* it looks like that the original code by Martijn, see below, may block in some situations, e.g. -grid 160 on */
      /* 32 processors. So, lets see what happens if we change all this with MPI_Sendrecv routines */
     /*  MPI_Send(buffer,bufsize,mpi_REAL,part,0,MPI_COMM_WORLD); */
     /* MPI_Barrier(MPI_COMM_WORLD); */
 
      /* MPI_Recv(buffer,bufsize,mpi_REAL,part,0,MPI_COMM_WORLD,&status); */

      MPI_Sendrecv(buffer, bufsize, mpi_REAL, part, 0,
		  rbuffer, bufsize, mpi_REAL, part, 0,
		  MPI_COMM_WORLD,&status);

      position = 0;
      for(z=local_z0;z<local_z1;z++) for(y=0;y<smallY;y++) {
        memcpy(&X[2*index_Xmatrix(local_xs_unit*part,y,z)],
               &rbuffer[position],
               2*local_xs_unit*sizeof(REAL));
        if ( (position += 2*local_xs_unit) > bufsize)
          printf("block_transpose, ringid=%d transm=%d\n buffer-overflow: bufsize=%d position=%d",ringid,transmission,bufsize,position);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      /* free(buffer);*/
#endif
    }
  free(buffer);
  free(rbuffer);
  tstop=extime();
  Timing_OneIterComm += tstop - tstart;
}

void par_setup(void)
{
  extern int boxX,boxY,boxZ,smallZ,smallX;
  int unit;

#ifdef PARALLEL
  unit=(float) smallZ/nprocs+.99999;
  local_z0=ringid*unit;
  local_z1=(ringid+1)*unit;
  if (local_z1 > boxZ) local_z1=boxZ;
  local_Nz=local_z1-local_z0;
  local_d0=boxX*boxY*local_z0;
  local_d1=boxX*boxY*local_z1;
  local_Ndip=local_d1-local_d0;
  local_xs_unit=smallX/nprocs;
#else
  local_d0=0;
  local_d1=boxX*boxY*boxZ;
  local_z0=0;
  local_z1=smallZ;
  local_Nz=local_z1-local_z0;
  local_Ndip=local_d1-local_d0;
  local_xs_unit=smallX/nprocs;
#endif
  Ndip=boxX*boxY*boxZ;
  printf("%i :  %i %i %i %i %i %i\n",ringid,local_d0,local_d1,local_z0,local_z1,local_Ndip,Ndip);
}

void all_gather(void *x_from,void *x_to,
		char data_type[],
		int n_elem)       
     /* Gather distributed arrays */
{
#if defined(MPI)
  MPI_Datatype mes_type;

  if (strcmp(data_type,"int")==0)
    mes_type = MPI_INT;
  if (strcmp(data_type,"REAL")==0)
    mes_type = mpi_REAL;
  if (strcmp(data_type,"dcomplex")==0)
    {
      mes_type = mpi_REAL;
      n_elem *= 2;
    }

  MPI_Allgather(x_from,n_elem,mes_type,
		x_to,n_elem,mes_type,MPI_COMM_WORLD);
#endif
}

unsigned long extime(void)
{
  return (unsigned long) TimeNow();
}






/* LogError
 * we use sprintf a couple of times, because we want each node to
 * generate an atomic message, not a couple of messages after
 * each other, since other nodes may then interfere with our
 * output!
 */
PUBLIC void
LogError( int code, char *fname, char *fmt, ... )
{
  va_list args;
  char    line[255];
 
  va_start(args, fmt);
  sprintf(line, "LogError: error %d, function %s: ", code, fname);
  vsprintf(line+strlen(line), fmt, args);
  sprintf(line+strlen(line), "\n");
  fprintf(stderr, line);
  fflush(stderr);
  va_end(args);
}
 
PUBLIC void
AbortServer(int status)
{
  stop(status);
}
