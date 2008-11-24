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
#include <stdlib.h>
#include "stdarg.h"
#include <time.h>
#include "debug.h"
#include "cmplx.h"
#include "comm.h"
#include "const.h"
/* for getch in stop */
#ifdef _WIN32
# include "conio.h"
#endif

#ifdef MPI
 #include <mpi.h>
 MPI_Datatype mpi_dcomplex;
#endif

int nprocs;   /* total number of processes */
int ringid;   /* id of current process */
int Ntrans;   /* number of transmissions; used in calc_partner */

int local_z0,local_z1,local_Nz,local_x0,       /* common variables; initialized in par_setup() */
    local_x1,local_Nx,local_z1_f,local_Nz_f,
    local_d0,local_d1,local_Ndip,Ndip;



/*===========================================*/

void init_comm(int *argc,char ***argv)
  /* initialize communications in the beginning of the program */
{
#ifdef MPI
  int dcmplx_blocklength[2]={1,1};
  MPI_Aint dcmplx_displs[2] = {0,1};
  MPI_Datatype dcmplx_type[2];

  MPI_Init(argc,argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&ringid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  /* initialize Ntrans */
  if ((nprocs%2)==0) Ntrans=nprocs-1;
  else Ntrans=nprocs;

  /* Create MPI-type for sending dcomplex-variables */
  dcmplx_type[0] = MPI_DOUBLE; dcmplx_type[1] = MPI_DOUBLE;
  MPI_Type_struct(2,dcmplx_blocklength,dcmplx_displs,dcmplx_type,&mpi_dcomplex);
  MPI_Type_commit(&mpi_dcomplex);
#elif !defined(PARALLEL)
  nprocs=1; ringid=ROOT;
#endif
}

/*===========================================*/

void stop(int code)
  /* stops the programm with code */
{
  /* exit the program */
  printf("end %i\n",ringid);
#ifdef MPI
  if (code) {
    printf("Aborting process %d\n",ringid);
    MPI_Abort(MPI_COMM_WORLD,code);
  }
  else MPI_Finalize();
#endif
#ifdef _WIN32
  printz("\nProgram has finished execution.\nPress any key to close window...");
  /* waits for pressed key */
  getch();
#endif
  exit(code);
}

/*===========================================*/

void synchronize(void)
  /* synchronizes all processes */
{
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/*===========================================*/

void Bcast_orient(int *i, int *j, int *k)
    /* cast current orientation angle (in orientation averaging)
       to all processes from root */
{
#ifdef MPI
  int buf[3];

  if (ringid==ROOT) {
    buf[0]=*i;
    buf[1]=*j;
    buf[2]=*k;
  }

  MPI_Bcast(buf,3,MPI_INT,ROOT,MPI_COMM_WORLD);

  if (ringid!=ROOT) {
    *i=buf[0];
    *j=buf[1];
    *k=buf[2];
  }
#endif
}

/*===========================================*/

void Bcast_parms(Parms_1D *parms)
   /* cast integration parameters to all processes
      !!! each processor should read it from file - will be faster !!! */
{
#ifdef MPI
  int len[6] = {1,1,1,1,1,1};
  MPI_Datatype mpi_parms,type[6];
  MPI_Aint disp[6];

  disp[0] = 0;                     type[0] = MPI_DOUBLE;
  disp[1] = sizeof(double);          type[1] = MPI_INT;
  disp[2] = disp[1] + sizeof(int); type[2] = MPI_INT;
  disp[3] = disp[2] + sizeof(int); type[3] = MPI_DOUBLE;
  disp[4] = disp[3] + sizeof(double);type[4] = MPI_DOUBLE;
  disp[5] = disp[4] + sizeof(double);type[5] = MPI_INT;

  MPI_Type_struct(6,len,disp,type,&mpi_parms);
  MPI_Type_commit(&mpi_parms);

  MPI_Bcast(parms,2,mpi_parms,ROOT,MPI_COMM_WORLD);
#endif
}

/*===========================================*/

void AccumulateMax(double *data,double *max)
    /* given element, accumulates their sum and maximum on root processor */
{
#ifdef MPI
  double buf;
                   /* potentially can be optimized by combining into one operation */
  MPI_Reduce(data,&buf,1,MPI_DOUBLE,MPI_SUM,ROOT,MPI_COMM_WORLD);
  MPI_Reduce(data,max,1,MPI_DOUBLE,MPI_MAX,ROOT,MPI_COMM_WORLD);
  *data=buf;
#endif
}

/*===========================================*/

void accumulate(double *vector,int n,double *buf)
    /* gather and add scattered fields on proces root */
{
#ifdef MPI
  MPI_Reduce(vector,buf,n,MPI_DOUBLE,MPI_SUM,ROOT,MPI_COMM_WORLD);
  memcpy(vector,buf,n*sizeof(double));
#endif
}

/*===========================================*/

void my_inner_product(void *data,var_type type,int n_elem)
   /* gather values stored in *a, add them and return them in *a
      works for all types; not optimized for long data (allocates memory at every call) */
{
#ifdef MPI
  int size;
  MPI_Datatype mes_type;
  void *temp;

  if (type==int_type) {
    mes_type=MPI_INT;
    size=n_elem*sizeof(int);
  }
  else if (type==double_type) {
    mes_type=MPI_DOUBLE;
    size=n_elem*sizeof(double);
  }
  else if (type==cmplx_type) {
    mes_type=MPI_DOUBLE;
    n_elem*=2;
    size=n_elem*sizeof(double);
  }
  else LogError(EC_ERROR,ONE,POSIT,"my_inner_product: variable type %d is not supported",type);

  if ((temp=malloc(size))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"could not malloc temp");
  MPI_Allreduce(data,temp,n_elem,mes_type,MPI_SUM,MPI_COMM_WORLD);
  memcpy(data,temp,size);
  free(temp);
#endif
}

/*============================================================*/
/* following two functions are used only for block_transpose */
#ifdef PARALLEL

INLINE int index_block(int x,int y,int z,int lengthY)
   /* index block; used in block_transpose */
{
  extern int gridX;

  return((z*lengthY+y)*gridX+x);
}

/*============================================================*/

INLINE int calc_partner(int tran)
   /* calculate ringid of partner processor for current transmission;
      used in block_transpose
      many different implementations are possible; the only requirements are
      1) f(tran,f(tran,ringid))=ringid
      2) f({1,2,Ntrans},ringid)={0,1,Ntrans}\ringid
      where f=nprocs <=> skip this transmission (for odd nprocs) */
{
  int part;

  if (ringid==0) part=tran;
  else if (ringid==tran) part=0;
  else {
    part=2*tran-ringid;
    if (part<=0) part+=Ntrans;
    else if (part>Ntrans) part-=Ntrans;
  }
  return part;
}

#endif
/*============================================================*/

void block_transpose(doublecomplex *X)
  /* do the data-transposition, i.e. exchange, between fftX and fftY&fftZ
     specializes at Xmatrix; do 3 components in one message */
{
#ifdef MPI
  clock_t tstart;
  extern clock_t Timing_OneIterComm;
  extern double *BT_buffer, *BT_rbuffer;
  extern int smallY,local_Nsmall;
  int bufsize,y,z,transmission,part;
  int posit,step,Xpos,msize,Xcomp;
  MPI_Status status;

  tstart=clock();

  step=2*local_Nx;
  msize=local_Nx*sizeof(doublecomplex);
  bufsize=6*local_Nz_f*smallY*local_Nx;

  for(transmission=1;transmission<=Ntrans;transmission++) {
    if ((part=calc_partner(transmission))!=nprocs) {   /* if part==nprocs then scip this transmission */
      posit=0;
      Xpos=local_Nx*part;
      for(Xcomp=0;Xcomp<3;Xcomp++) for(z=0;z<local_Nz_f;z++) for(y=0;y<smallY;y++) {
        memcpy(BT_buffer+posit,X+Xcomp*local_Nsmall+index_block(Xpos,y,z,smallY),msize);
        posit+=step;
      }

      MPI_Sendrecv(BT_buffer, bufsize, MPI_DOUBLE, part, 0,
        	  BT_rbuffer, bufsize, MPI_DOUBLE, part, 0,
        	  MPI_COMM_WORLD,&status);

      posit=0;
      Xpos=local_Nx*part;
      for(Xcomp=0;Xcomp<3;Xcomp++) for(z=0;z<local_Nz_f;z++) for(y=0;y<smallY;y++) {
        memcpy(X+Xcomp*local_Nsmall+index_block(Xpos,y,z,smallY),BT_rbuffer+posit,msize);
        posit+=step;
      }
    }
  }
  Timing_OneIterComm += clock() - tstart;
#endif
}

/*===========================================*/

void block_transpose_Dm(doublecomplex *X, int lengthY, int lengthZ)
  /* do the data-transposition, i.e. exchange, between fftX and fftY&fftZ
     specialized for D matrix */
{
#ifdef MPI
  extern double *BT_buffer, *BT_rbuffer;
  int bufsize,y,z,transmission,part;
  int posit,step,Xpos,msize;
  MPI_Status status;

  step=2*local_Nx;
  msize=local_Nx*sizeof(doublecomplex);
  bufsize = 2*lengthZ*lengthY*local_Nx;

  for(transmission=1;transmission<=Ntrans;transmission++) {
    if ((part=calc_partner(transmission))!=nprocs) {
      posit=0;
      Xpos=local_Nx*part;
      for(z=0;z<lengthZ;z++) for(y=0;y<lengthY;y++) {
        memcpy(BT_buffer+posit,X+index_block(Xpos,y,z,lengthY),msize);
        posit+=step;
      }

      MPI_Sendrecv(BT_buffer,bufsize,MPI_DOUBLE,part,0,
        	   BT_rbuffer,bufsize,MPI_DOUBLE,part,0,
        	   MPI_COMM_WORLD,&status);

      posit=0;
      Xpos=local_Nx*part;
      for(z=0;z<lengthZ;z++) for(y=0;y<lengthY;y++) {
        memcpy(X+index_block(Xpos,y,z,lengthY),BT_rbuffer+posit,msize);
        posit+=step;
      }
    }
  }
#endif
}

/*===========================================*/

void par_setup(void)
   /* initialize common parameters; need to do in the beginning to enable call to make_particle */
{
  extern int boxX,boxY,boxZ,smallY,smallZ;
  extern int gridX,gridY,gridZ,gridYZ,gridXY,gridXZ;

  extern int fft_fit(int size, int div);
  /* calculate size of 3d grid */
  gridX=fft_fit(2*boxX,nprocs);
  gridY=fft_fit(2*boxY,1);
  gridZ=fft_fit(2*boxZ,2*nprocs);
  /* initialise some variables */
  smallY=gridY/2;
  smallZ=gridZ/2;
  gridXY=gridX*gridY;
  gridYZ=gridY*gridZ;
  gridXZ=gridX*gridZ;
#ifdef PARALLEL
  int unitZ,unitX;

  unitZ=(float) smallZ/nprocs+.99999;
  local_z0=ringid*unitZ;
  local_z1_f=(ringid+1)*unitZ;
  if (local_z1_f > boxZ) local_z1=boxZ;
  else local_z1=local_z1_f;
  unitX=(float) gridX/nprocs+.99999;
  local_x0=ringid*unitX;
  local_x1=(ringid+1)*unitX;
#else
  local_z0=0;
  local_z1_f=smallZ;
  local_z1=boxZ;
  local_x0=0;
  local_x1=gridX;
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

void all_gather(void *x_from,void *x_to,var_type type,int n_elem)
     /* Gather distributed arrays; works for all types */
{
#ifdef MPI
  /* need to be rewritten when n_elem are unequal on each processor */
  MPI_Datatype mes_type;

  if (type==int_type) mes_type = MPI_INT;
  else if (type==double_type) mes_type = MPI_DOUBLE;
  else if (type==cmplx_type) {
    mes_type = MPI_DOUBLE;
    n_elem *= 2;
  }
  else LogError(EC_ERROR,ONE,POSIT,"all_gather: variable type %d is not supported",type);

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
   /* who specifies whether 1 (ringid=ROOT) or all processors should produce output */
   /* if code is EC_ERROR program aborts after output                            */
{
  va_list args;
  char line[255];
  char id_str[30];
  extern char logname[];
  extern FILE *logfile;

  if (who==ALL || ringid==ROOT) {
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
    fprintf(stderr,line);
    fflush(stderr);
    va_end(args);
  }
  if (code==EC_ERROR) {
    /* duplicate error message in logfile */
    if (logname[0]!=0) {  /* otherwise can't do anything */
      if (ringid==ROOT) {
        /* logfile==NULL with logname!=0 may only occur when error is
           in opening of the logfile itself */
        if (logfile!=NULL) {
          fprintf(logfile,line);
          fclose(logfile);
        }
      }
      else if (who==ALL) {
        logfile=fopen(logname,"a");
        fprintf(logfile,line);
        fclose(logfile);
      }
    }
    stop(1);
  }
}



