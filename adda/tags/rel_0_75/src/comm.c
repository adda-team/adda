/* FILE : comm.c
 * AUTH : Maxim Yurkin 
 * DESCR: The main intention of this library is to incorporate all
 *        parallelisation related code, so most of it is directly
 *        involved in or closely related to inter-proces communication,
 *        hence its name.
 *
 *        Previous versions were by Martijn Frijlink
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include "vars.h"
#include "debug.h"
#include "cmplx.h"
#include "comm.h"
#include "const.h"
#include "io.h"
#include "fft.h"
/* for getch in Stop */
#ifdef RUN_BCB
# include <conio.h>
#endif

#ifdef MPI
 #include <mpi.h>
 MPI_Datatype mpi_dcomplex;
#endif

/* SEMI-GLOBAL VARIABLES */

/* defined and allocated in fft.c */
extern double *BT_buffer, *BT_rbuffer;

/* LOCAL VARIABLES */

#ifdef PARALLEL
static int Ntrans;   /* number of transmissions; used in CalcPartner */

/* First funtions that are defined only in parallel mode */
/*===========================================*/

static void RecoverCommandLine(int *argc_p,char ***argv_p)
  /* eliminate all NULL pointers from argv, shift the rest, and
     adjust argc accordingly. Used in InitComm */
{
   int i,j;

   for (i=0,j=0;i<(*argc_p);i++) {
     if ((*argv_p)[i]==NULL) j++;
     else if (j!=0) (*argv_p)[i-j]=(*argv_p)[i];
   }
   (*argc_p)-=j;
}
/*============================================================*/

INLINE int IndexBlock(const int x,const int y,const int z,const int lengthY)
   /* index block; used in BlockTranspose */
{
  return((z*lengthY+y)*gridX+x);
}

/*============================================================*/

INLINE int CalcPartner(const int tran)
   /* calculate ringid of partner processor for current transmission;
      used in BlockTranspose
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

/*===========================================*/

void CatNFiles(const char *dir,const char *tmpl,const char *dest)
   /* cat several temporary files (one for each processor, names defines by the template 'temp'
      that should contain %d to be replaced by ringid). Files are located in directory 'dir'.
      Combined into 'dest' in the same directory. Afterwards temporary files are removed. */
{
  int i,c;
  FILE *in,*out;
  char fname_out[MAX_TMP_FNAME],fname_in[MAX_TMP_FNAME];

  /* produce full path of destination file and open it */
  sprintf(fname_out,"%s/%s",directory,dest);
  out=FOpenErr(fname_out,"w",ONE_POS);
  for (i=0;i<nprocs;i++) {
    /* produce full path of tmp file and open it */
    sprintf(fname_in,"%s/",directory);
    sprintf(fname_in+strlen(fname_in),tmpl,i);
    in=FOpenErr(fname_in,"r",ONE_POS);
    /* copy file in to out */
    while((c=getc(in))!=EOF) putc(c,out);
    /* close and remove tmp file */
    FCloseErr(in,fname_in,ONE_POS);
    RemoveErr(fname_in,ONE_POS);
  }
  /* close destination file */
  FCloseErr(out,fname_out,ONE_POS);
}
#endif

/* Further routines are defined always, although in sequential mode they are blank */
/*===========================================*/

void InitComm(int *argc_p,char ***argv_p)
  /* initialize communications in the beginning of the program */
{
#ifdef MPI
  int dcmplx_blocklength[2]={1,1};
  MPI_Aint dcmplx_displs[2] = {0,1};
  MPI_Datatype dcmplx_type[2];

  /* MPI_Init may alter argc and argv and interfere with normal parsing of command
     line parameters. The way of altering is implementation depending. MPI searches for
     mpi parameters in command line and removes them (we assume some kind of removing does
     take place - otherwise ADDA will produce error 'unknown parameter'). The best would be
     to change argc and argv so that they look like no special command line arguments are
     present. However, MPICH 1.2.5, for example, just replaces corresponding parameters by
     NULLs. To incorporate it we introduce special function to restore the command line */
  MPI_Init(argc_p,argv_p);
  RecoverCommandLine(argc_p,argv_p);
  /* initialize ringid and nprocs */
  MPI_Comm_rank(MPI_COMM_WORLD,&ringid);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  /* initialize Ntrans */
  if ((nprocs%2)==0) Ntrans=nprocs-1;
  else Ntrans=nprocs;
  /* Create MPI-type for sending dcomplex-variables */
  dcmplx_type[0] = MPI_DOUBLE; dcmplx_type[1] = MPI_DOUBLE;
  MPI_Type_struct(2,dcmplx_blocklength,dcmplx_displs,dcmplx_type,&mpi_dcomplex);
  MPI_Type_commit(&mpi_dcomplex);
  /* if MPI crashes, it happens here */
  Synchronize();
#elif !defined(PARALLEL)
  nprocs=1;
  ringid=ROOT;
#endif
  /* check if wierd number of processors is specified
     called even in sequenatial mode to initialize weird_nprocs */
  CheckNprocs();
}

/*===========================================*/

void Stop(const int code)
  /* stops the programm with code */
{
#ifdef MPI
  if (code) {  /* error occured */
    fflush(stdout);
    fprintf(stderr,"Aborting process %d\n",ringid);
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD,code);
  }
  else {      /* regular termination */
    /* wait for all processors */
    fflush(stdout);
    Synchronize();
    /* finalize MPI communications */
    MPI_Finalize();
  }
#endif
/* if run under Borland C++ Builder, don't close the window automatically */
#ifdef RUN_BCB
  PRINTZ("\nProgram has finished execution.\nPress any key to close window...");
  /* waits for pressed key */
  getch();
#endif
  /* in the code '0' corresponds to success, and '1' to failure */
  if (code==0) exit(EXIT_SUCCESS);
  else if (code==1) exit(EXIT_FAILURE);
  else exit(code);
}

/*===========================================*/

void Synchronize(void)
  /* synchronizes all processes */
{
#ifdef MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
}

/*===========================================*/

void MyBcast(void *data,const var_type type,int n_elem)
   /* casts values stored in *data from ROOT processor to all other;
      works for all types; not optimized for long data (allocates memory at every call) */
{
#ifdef MPI
  int size;
  MPI_Datatype mes_type;
  void *temp;

  if (type==char_type) {
    mes_type=MPI_CHAR;
    size=n_elem*sizeof(char);
  }
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
  else LogError(EC_ERROR,ONE_POS,"MyBcast: variable type %u is not supported",type);
  /* allocate buffer */
  if ((temp=malloc(size))==NULL)
    LogError(EC_ERROR,ALL_POS,"could not malloc temp");
  /* fill buffer on ROOT processor */
  if (ringid==ROOT) memcpy(temp,data,size);
  /* cast */
  MPI_Bcast(temp,n_elem,mes_type,ROOT,MPI_COMM_WORLD);
  /* copy buffer on other processors */
  if (ringid!=ROOT) memcpy(data,temp,size);

  free(temp);
#endif
}

/*===========================================*/

void BcastOrient(int *i, int *j, int *k)
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

void Accumulate(double *vector,const int n,double *buf)
    /* gather and add scattered fields on proces root */
{
#ifdef MPI
  MPI_Reduce(vector,buf,n,MPI_DOUBLE,MPI_SUM,ROOT,MPI_COMM_WORLD);
  memcpy(vector,buf,n*sizeof(double));
#endif
}

/*===========================================*/

void MyInnerProduct(void *data,const var_type type,int n_elem)
   /* gather values stored in *data, add them and return them in *data
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
  else LogError(EC_ERROR,ONE_POS,"MyInnerProduct: variable type %u is not supported",type);

  if ((temp=malloc(size))==NULL)
    LogError(EC_ERROR,ALL_POS,"could not malloc temp");
  MPI_Allreduce(data,temp,n_elem,mes_type,MPI_SUM,MPI_COMM_WORLD);
  memcpy(data,temp,size);
  free(temp);
#endif
}

/*============================================================*/

void BlockTranspose(doublecomplex *X)
  /* do the data-transposition, i.e. exchange, between fftX and fftY&fftZ
     specializes at Xmatrix; do 3 components in one message */
{
#ifdef MPI
  clock_t tstart;
  int bufsize,y,z,transmission,part;
  int posit,step,Xpos,msize,Xcomp;
  MPI_Status status;

  tstart=clock();

  step=2*local_Nx;
  msize=local_Nx*sizeof(doublecomplex);
  bufsize=6*local_Nz*smallY*local_Nx;

  for(transmission=1;transmission<=Ntrans;transmission++) {
    /* if part==nprocs then scip this transmission */
    if ((part=CalcPartner(transmission))!=nprocs) {
      posit=0;
      Xpos=local_Nx*part;
      for(Xcomp=0;Xcomp<3;Xcomp++) for(z=0;z<local_Nz;z++) for(y=0;y<smallY;y++) {
        memcpy(BT_buffer+posit,X+Xcomp*local_Nsmall+IndexBlock(Xpos,y,z,smallY),msize);
        posit+=step;
      }

      MPI_Sendrecv(BT_buffer, bufsize, MPI_DOUBLE, part, 0,
        	  BT_rbuffer, bufsize, MPI_DOUBLE, part, 0,
        	  MPI_COMM_WORLD,&status);

      posit=0;
      Xpos=local_Nx*part;
      for(Xcomp=0;Xcomp<3;Xcomp++) for(z=0;z<local_Nz;z++) for(y=0;y<smallY;y++) {
        memcpy(X+Xcomp*local_Nsmall+IndexBlock(Xpos,y,z,smallY),BT_rbuffer+posit,msize);
        posit+=step;
      }
    }
  }
  Timing_OneIterComm += clock() - tstart;
#endif
}

/*===========================================*/

void BlockTranspose_Dm(doublecomplex *X,const int lengthY,const int lengthZ)
  /* do the data-transposition, i.e. exchange, between fftX and fftY&fftZ
     specialized for D matrix */
{
#ifdef MPI
  int bufsize,y,z,transmission,part;
  int posit,step,Xpos,msize;
  MPI_Status status;

  step=2*local_Nx;
  msize=local_Nx*sizeof(doublecomplex);
  bufsize = 2*lengthZ*lengthY*local_Nx;

  for(transmission=1;transmission<=Ntrans;transmission++) {
    if ((part=CalcPartner(transmission))!=nprocs) {
      posit=0;
      Xpos=local_Nx*part;
      for(z=0;z<lengthZ;z++) for(y=0;y<lengthY;y++) {
        memcpy(BT_buffer+posit,X+IndexBlock(Xpos,y,z,lengthY),msize);
        posit+=step;
      }

      MPI_Sendrecv(BT_buffer,bufsize,MPI_DOUBLE,part,0,
        	   BT_rbuffer,bufsize,MPI_DOUBLE,part,0,
        	   MPI_COMM_WORLD,&status);

      posit=0;
      Xpos=local_Nx*part;
      for(z=0;z<lengthZ;z++) for(y=0;y<lengthY;y++) {
        memcpy(X+IndexBlock(Xpos,y,z,lengthY),BT_rbuffer+posit,msize);
        posit+=step;
      }
    }
  }
#endif
}

/*===========================================*/

void ParSetup(void)
   /* initialize common parameters; need to do in the beginning to enable call to MakeParticle */
{
#ifdef PARALLEL
  int unitZ,unitX;
#endif
  /* calculate size of 3d grid */
  gridX=fftFit(2*boxX,nprocs);
  gridY=fftFit(2*boxY,1);
  gridZ=fftFit(2*boxZ,2*nprocs);
  /* initialise some variables */
  smallY=gridY/2;
  smallZ=gridZ/2;
  gridYZ=gridY*gridZ;
#ifdef PARALLEL
  unitZ=(int)ceil(((double)smallZ)/((double)nprocs));
  local_z0=ringid*unitZ;
  local_z1=(ringid+1)*unitZ;
  if (local_z1 > boxZ) local_z1_coer=boxZ;
  else local_z1_coer=local_z1;
  unitX=(int)ceil(((double)gridX)/((double)nprocs));
  local_x0=ringid*unitX;
  local_x1=(ringid+1)*unitX;
#else
  local_z0=0;
  local_z1=smallZ;
  local_z1_coer=boxZ;
  local_x0=0;
  local_x1=gridX;
#endif
  if (local_z1_coer<=local_z0) {
    LogError(EC_INFO,ONE_POS,"No real dipoles are assigned");
    local_z1_coer=local_z0;
  }
  local_Nz=local_z1-local_z0;
  local_Nx=local_x1-local_x0;
  local_d0=boxX*boxY*local_z0;
  local_d1=boxX*boxY*local_z1_coer;
  local_Ndip=local_d1-local_d0;
  printf("%i :  %i %i %i %i %i \n",
         ringid,local_z0,local_z1_coer,local_z1,local_Ndip,local_Nx);
}

/*===========================================*/

void AllGather(void *x_from,void *x_to,const var_type type,int n_elem)
     /* Gather distributed arrays; works for all types */
{
#ifdef MPI
  /* need to be rewritten when n_elem are unequal on each processor */
  MPI_Datatype mes_type;

  if (type==char_type) mes_type = MPI_CHAR;
  else if (type==int_type) mes_type = MPI_INT;
  else if (type==double_type) mes_type = MPI_DOUBLE;
  else if (type==cmplx_type) {
    mes_type = MPI_DOUBLE;
    n_elem *= 2;
  }
  else LogError(EC_ERROR,ONE_POS,"AllGather: variable type %u is not supported",type);

  MPI_Allgather(x_from,n_elem,mes_type,
		x_to,n_elem,mes_type,MPI_COMM_WORLD);
#endif
}




