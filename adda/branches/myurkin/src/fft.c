/* FILE: fft_init.c
 * AUTH: Michel Grimminck, using old code of Alfons Hoekstra
 * DATE: januari 1995
 */


#include <stdio.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "types.h"
#include "comm.h"

#ifdef PowerPVM
#define cfft99_ cfft99
#define cftfax_ cftfax
#endif

#define FFT NEW

extern dcomplex *dCvector(int nl,int nh);

extern unsigned long Timing_FFT_Init;

extern int RingId;           /* ID of this processor in the Ring */
extern int MyProcId;
extern int symD;
extern double gridspaceX,gridspaceY,gridspaceZ;
extern int fsc,prognose;

dcomplex *Dmatrix;           /* size: Dmatrix[gridX][gridY][gridZ][6] */
dcomplex *D2matrix;
int DsizeX,DsizeY,DsizeZ;    /* the size of the 'matrix' D */ 
int gridX,gridY,gridZ,Ngrid; /* the size of the 'matrix' X */
int gridXY,gridYZ,gridXZ;
int D2sizeX,D2sizeY,D2sizeZ;
int DsizeXY,DsizeYZ,DsizeXZ;
int smallX,smallY,smallZ,Nsmall;    /* the size of the reduced matrix X */
int local_Nsmall;
int NDcomp;                  /* number of components of D */
int Dmatrix_nel;
int minY,maxY,minZ,maxZ;

#define signed 

signed char *useZ,*useX,*useY;

/**************************************************************************/
int index_Dmatrix(int x,int y,int z)
{
  extern int abs(int);
  
  if (z<0 && fsc==true) z=DsizeZ-z;
  if (x>=DsizeX) x=gridX-x;
  if (y>=DsizeY) y=gridY-y;
  if (z>=DsizeZ && fsc==false) z=gridZ-z;
  
  return(NDcomp*(z*DsizeX*DsizeY+y*DsizeX+x));
}
/**************************************************************************/
int index_tDmatrix(int x,int y,int z)
{
  extern int abs(int);
  
  if (z<0 && fsc==true) z=DsizeZ-z;
  if (x>=DsizeX) x=gridX-x;
  if (y>=DsizeY) y=gridY-y;
  if (z>=DsizeZ && fsc==false) z=gridZ-z;
  
  if (x>=DsizeX || y>=DsizeY || z>=DsizeZ) printf("ERROR:index_tDmatrix:%i %i %i\n",x,y,z);
  
  return(NDcomp*(x*DsizeYZ+z*DsizeY+y));
}
/**************************************************************************/
int index_D2matrix(int x,int y,int z)
{
  extern int abs(int);
  
  if (x<0) x=D2sizeX+x;
  if (y<0) y=D2sizeY+y;
  if (z<0) z=D2sizeZ+z;
  
  if (x>=gridX || y>=gridY || z>=gridZ || x<0 || y<0 || z<0) printf("ERROR:index_D2matrix:%i %i %i\n",x,y,z);
  
  return((z*D2sizeX*D2sizeY+y*D2sizeX+x));
}
/**************************************************************************/
int index_sliceD2matrix(int y,int z)
{
  if (y<0) y=gridY+y;
  if (z<0) z=gridZ+z;
  if (y*gridZ+z<0 || y*gridZ+z>gridY*gridZ) printf("boudn error\n %i %i",y,z);
  return((y*gridZ+z));
}
/**************************************************************************/
int index_slice_zyD2matrix(int y,int z)
{
  if (y<0) y=gridY+y;
  if (z<0) z=gridZ+z;
  if (y*gridZ+z<0 || y*gridZ+z>gridY*gridZ) printf("boudn error\n %i %i",y,z);
  
  return((z*gridY+y));
}
/**************************************************************************/
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(data,nn,isign)
     REAL data[];
     int nn,isign;
{
  int n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  
  data-=1;
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2) {
    if (j > i) {
      SWAP(data[j],data[i]);
      SWAP(data[j+1],data[i+1]);
    }
    m=n >> 1;
    while (m >= 2 && j > m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }

  mmax=2;
  while (n > mmax) {
    istep=2*mmax;
    theta=6.28318530717959/(isign*mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
	j=i+mmax;
	tempr=wr*data[j]-wi*data[j+1];
	tempi=wr*data[j+1]+wi*data[j];
	data[j]=data[i]-tempr;
	data[j+1]=data[i+1]-tempi;
	data[i] += tempr;
	data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}
/**************************************************************************/
void transposeYZ(REAL *data,int X,int Y,int Z,int start,int end)
     /* optimised routine to transpose y and z.
      * We know smallX and smallY are even, so we can do
      * two points at the same time
      */
{
  int x,y,z,a,b;
  REAL *temp,*work,*nowtemp,*nowwork;
  
  temp=(REAL *) malloc(2*sizeof(REAL)*Z*Y);
  if (temp==0) {
    printf("malloc temp failed\n");
    stop(1);
  }
  a=2*Y;b=a+1;
  for(x=start;x<end;x++) {
    work=data+2*x*Y*Z;
    for(y=0;y<2*Y;y+=2) {
      nowwork=work+y*Z;
      nowtemp=temp+y;
      for(z=0;z<2*Z;z+=4) {
	nowtemp[0]=nowwork[z];
	nowtemp[1]=nowwork[z+1];
	nowtemp[a]=nowwork[z+2];
	nowtemp[b]=nowwork[z+3];
	nowtemp+=4*Y;
      }
    }
    for(z=0;z<Z;z++) {
      nowwork=work+2*z*Y;
      nowtemp=temp+2*z*Y;
      for(y=0;y<2*Y;y+=4) {
	nowwork[y]=nowtemp[y];
	nowwork[y+1]=nowtemp[y+1];
	nowwork[y+2]=nowtemp[y+2];
	nowwork[y+3]=nowtemp[y+3];
      }
    }
  }
  free(temp);
}
/**************************************************************************/
void fftX(REAL *data,int isign,dcomplex *dummy,int z0,int z1)
{
  int nn=smallX,dim=1,iform=1,is,x,y,z,inc=1,jump=nn,lot=1;
  int ifax[20];
  REAL *trigs,*work;
  
  if ((trigs = (REAL *) malloc(2*nn*sizeof(REAL))) == NULL) {
    printf("malloc trigs\n"); AbortServer (1);
  }
  if ((work = (REAL *) malloc(2*nn*lot*sizeof(REAL))) == NULL) {
    printf("malloc work\n"); AbortServer (1);
  }
  
  is=-isign;
#ifdef FFT
  
  cftfax_ (&nn,ifax,trigs);
  for (z=z0;z<z1;z++) for(y=0;y<smallY;y++) {
    
    cfft99_ (data+2*(z*smallX*smallY+y*smallX),work,trigs,ifax,&inc,&jump,&nn,
	     &lot,&is);
  }
  
#else
  for (z=z0;z<z1;z++) for(y=0;y<smallY;y++) {
    four1(data+2*(z*smallX*smallY+y*smallX),nn,is);
  }
#endif
  free(work);
  free(trigs);
}
/**************************************************************************/
void fftY(REAL *data,int isign,int start,int end,dcomplex *dummy)
{
  int nn=gridY,dim=1,iform=1,is,x,y,z,inc=1,jump=nn,lot=gridZ/2;
  int ifax[20];
  REAL *trigs,*work;
  
  if ((trigs = (REAL *) malloc(2*nn*sizeof(REAL))) == NULL) {
    printf("malloc trigs\n"); AbortServer (1);
  }
  if ((work = (REAL *) malloc(2*nn*lot*sizeof(REAL))) == NULL) {
    printf("malloc work\n"); AbortServer (1);
  }
  
  is=-isign;
#ifdef FFT
  cftfax_ (&nn,ifax,trigs);
  for(x=start;x<end;x++) for (z=0;z<gridZ;z+=gridZ/2) {
    cfft99_ (data+2*(x*gridZ*gridY+z*gridY),work,trigs,ifax,&inc,&jump,&nn,
	     &lot,&is);
  }
#else
  for(x=start;x<end;x++) for (z=0;z<gridZ;z++) {
    four1(data+2*(x*gridZ*gridY+z*gridY),nn,is);
  }
#endif
  free(work);
  free(trigs);
}
/**************************************************************************/
void fftZ(REAL *data,int isign,int start,int end,dcomplex *dummy)
{
  int nn=gridZ,dim=1,iform=1,is,x,y,z,inc=1,jump=nn,lot=1;
  int ifax[20];
  REAL *trigs,*work;
  
  if ((trigs = (REAL *) malloc(2*nn*sizeof(REAL))) == NULL) {
    printf("malloc trigs\n"); AbortServer (1);
  }
  if ((work = (REAL *) malloc(2*nn*lot*sizeof(REAL))) == NULL) {
    printf("malloc work\n"); AbortServer (1);
  }
    
  is=-isign;
#ifdef FFT
  cftfax_ (&nn,ifax,trigs);
  for(y=0;y<gridY;y++) for(x=start;x<end;x++) if (useY[y]==true) {
    cfft99_ (data+2*(x*gridZ*gridY+y*gridZ),work,trigs,ifax,&inc,&jump,&nn,
	     &lot,&is);
  }
#else
  for(y=0;y<gridY;y++) for(x=start;x<end;x++) if (useY[y]==true) {
    four1(data+2*(x*gridZ*gridY+y*gridZ),nn,is);
  }
#endif
  free(work);
  free(trigs);
}
#undef SWAP
/**************************************************************************/
double puls(double x)
{
  if (x==0) return(1);
  return(sin(x)/x);
}
/**************************************************************************/
int fft_fit(int x)
     /* find the first number >=x divisible by 2,3 and 5 only
      * the resulting number is divisible by 2
      */
{
  int y;
  
  do {
    y=x;
    while (y%2==0) y/=2;
    while (y%3==0) y/=3;
    while (y%5==0) y/=5;
    if (y==1 && x%2==0) return(x);
    x++;
  } while(true);
}
/**************************************************************************/
void fill_use()
{
  int i;

  /* allocate memory for use arrays */
  if ((useZ = (signed char *) malloc(sizeof(signed char)*gridZ)) == NULL) {
    LogError (EC_ERROR, "fft",
              "processor %d, ringID %d, could not malloc useZ",
              MyProcId, RingId);
    AbortServer (1);
  }
  /* allocate memory for use arrays */
  if ((useY = (signed char *) malloc(sizeof(signed char)*gridY)) == NULL) {
    LogError (EC_ERROR, "fft",
              "processor %d, ringID %d, could not malloc useY",
              MyProcId, RingId);
    AbortServer (1);
  }
  /* allocate memory for use arrays */
  if ((useX = (signed char *) malloc(sizeof(signed char)*gridX)) == NULL) {
    LogError (EC_ERROR, "fft",
              "processor %d, ringID %d, could not malloc useX",
              MyProcId, RingId);
    AbortServer (1);
  }
  for(i=0;i<gridZ;i++) useZ[i]=true;
  for(i=0;i<gridY;i++) useY[i]=true;
  for(i=0;i<gridX;i++) useZ[i]=true;
}
/**************************************************************************/

void setup_use(REAL **dip_coord,int ndip,double invgridspaceX,double invgridspaceY,double invgridspaceZ)
{
  int i,j,x,y,z;
  double rtemp0,rtemp1,rtemp2;
  int minX,maxX;
  extern short int *position;
  extern int boxX,boxY,boxZ;

  minY=100000,minZ=100000,maxY=-100000,maxZ=-100000;
  
  for(i=0;i<gridZ;i++) useZ[i]=false;
  for(i=0;i<gridY;i++) useY[i]=false;
  for(i=0;i<gridX;i++) useZ[i]=false;
  
  minX=-boxX/2; maxX=minX+boxX-1;
  minY=-boxY/2; maxY=minY+boxY-1;
  minZ=0; maxZ=boxZ-1; /* ??? */
  for(i=minX;i<=maxX;i++) {
    j=i; if (j<0) j+=gridX;
    useX[j]=true;
  }
  for(i=minY;i<=maxY;i++) {
    j=i; if (j<0) j+=gridY;
    useY[j]=true;
  }
  for(i=minZ;i<=maxZ;i++) {
    j=i; if (j<0) j+=gridZ;
    useZ[j]=true;
  }
  printf("setup_use: %i, %i %i %i %i\n",ringid,minY,maxY,minZ,maxZ);
}

/**************************************************************************/
void init_Dmatrix(int X,int Y,int Z,
		  double lambda,
		  double dips_lambda,
		  double wavenum)
     /* initialises the matrix D. D[i][j][k]=A[i1-i2][j1-j2][k1-k2]
	The routine is called only once, so needs not to be fast.
	*/
{
  int x,y,z,i,j,k,index,mem,Dcomp,indexfrom,indexto;
  char *Fname = __FILE__;
  dcomplex temp_minus[6],temp_plus[6];
  double nx,ny,nz,correction;
  double invNgrid;
  REAL work[4000];
  dcomplex *slice;
  
  extern FILE *logfile;
  extern int memory;
  extern double dplX,dplY,dplZ;
  
  /* these variables from function couple_matrix */
  double f00r, f01r, f02r, f10r, f11r, f12r, f20r, f21r, f22r;
  double f00i, f01i, f02i, f10i, f11i, f12i, f20i, f21i, f22i;
  double rr, rtemp0, rtemp1, rtemp2;
  double rt00, rt01, rt02, rt11, rt12, rt22;
  double rrsq, rrcb, kr;
  double temp;
  int start;

  /* calculate size of 3d grid */ 
  gridX=fft_fit(2*X);
  gridY=fft_fit(2*Y);
  gridZ=fft_fit(2*Z);
  
  /* initialise small X,Y,Z as its big counterparts.
   * After the FFT of D set them to their correct values
   */
  
  smallX=gridX; smallZ=gridZ/2;par_setup();
  smallX=gridX; smallY=gridY; smallZ=gridZ;
  gridXY=gridX*gridY; gridYZ=gridY*gridZ; gridXZ=gridX*gridZ;
#ifdef PARALLEL
  DsizeX=gridX;
#else
  DsizeX=gridX/2+1;
#endif
  
  DsizeY=gridY/2+1;
  D2sizeX=gridX; D2sizeY=gridY/2; D2sizeZ=gridZ/2;
  
  if (fsc==true) DsizeZ=gridZ;
  else DsizeZ=gridZ/2+1;
  DsizeXY=DsizeX*DsizeY; DsizeYZ=DsizeY*DsizeZ; DsizeXZ=DsizeX*DsizeZ;
  
  if (DsizeX==DsizeY && DsizeX==DsizeZ && gridspaceX==gridspaceY && gridspaceX==gridspaceZ && fsc==false) symD=true; else symD=false;
#ifdef PARALLEL
  symD=false;
#endif
  Ngrid=gridX*gridY*gridZ; 
  invNgrid=1.0/Ngrid;
  Nsmall=smallX*smallY*smallZ/4;
  local_Nsmall=Nsmall/nprocs; /* ??*/
  if (symD==true) NDcomp=2; else NDcomp=6;
  
  fprintz(logfile,"The 3d grid is:%ix%ix%i\n",gridX,gridY,gridZ);
  mem=sizeof(dcomplex)*(3*local_Nsmall/4+NDcomp*(2*local_Nz)*DsizeY*DsizeZ);
  memory+=mem;
  printz("Memory usage for X and D matrix:%.1f Mb\n",mem/1048576.0);
  fprintz(logfile,"Memory usage for X and D matrix:%.1f Mb\n",mem/1048576.0);
  if (ringid==0) fflush(logfile);
  
  if (prognose==true) return;
  /* allocate memory for Dmatrix */
  if ((Dmatrix = dCvector(0, Dmatrix_nel=NDcomp*(2*local_Nz)*DsizeY*DsizeZ)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc Dmatrix",
              MyProcId, RingId);
    AbortServer (1);
  }
  
  /* allocate memory for D2matrix components */
  if ((D2matrix = dCvector(0, (2*local_Nz)*D2sizeY*D2sizeZ)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc D2matrix",
              MyProcId, RingId);
    AbortServer (1);
  }
  if ((slice = dCvector(0, gridY*gridZ)) == NULL) {
    LogError (EC_ERROR, Fname,
              "processor %d, ringID %d, could not malloc slice",
              0, 0);
    AbortServer (1);
  }
  
  fill_use();
  smallX=D2sizeX; smallY=D2sizeY; smallZ=D2sizeZ;
  start=clock(); 
  for(Dcomp=0;Dcomp<NDcomp;Dcomp++) {
    
    /* fill Dmatrix with 0.0 */
    
    for(z=local_z0;z<local_z1;z++) for(y=0;y<D2sizeY;y++) for(x=0;x<D2sizeX;x++) {
      i=index_Xmatrix(x,y,z);
      D2matrix[i].r=D2matrix[i].i=0.0;
    }
    for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
      i=index_sliceD2matrix(y,z);
      slice[i].r=slice[i].i=0.0;
    }
    /* fill D (F'i-j) */
    for(k=local_z0;k<local_z1;k++)  for(j=0;j<Y;j++) for(i=1-X;i<X;i++) {
      double ar, ai, br, bi, expr, expi;
      
      if (i==0 && j==0 && k==0) { /* self interaction */  
	index=index_Xmatrix(i,j,k);
	D2matrix[index].r=0.0;
	D2matrix[index].i=0.0;
	continue;
      }
      
      rtemp0=i*gridspaceX;
      rtemp1=j*gridspaceY;
      rtemp2=k*gridspaceZ;
      
      /* start of old code to calculate Fij */
      rt00 = rtemp0 * rtemp0;
      rt11 = rtemp1 * rtemp1;
      rt22 = rtemp2 * rtemp2;
      
      rr = sqrt(rt00 + rt11 + rt22);
      nx=rtemp0/rr;
      ny=rtemp1/rr;
      nz=rtemp2/rr;
      
      if(rr == 0.0) {
	LogError(EC_ERROR, Fname,
		 "dipole and observer at the same location\n");
	AbortServer(1);
      }
      
      rrsq = rr * rr;
      rrcb = rrsq * rr;
      
      kr = wavenum * rr;
      rtemp0 /= rr;
      rtemp1 /= rr;
      rtemp2 /= rr;
      
      rt00 = rtemp0 * rtemp0;
      rt01 = rtemp0 * rtemp1;
      rt02 = rtemp0 * rtemp2;
      rt11 = rtemp1 * rtemp1;
      rt12 = rtemp1 * rtemp2;
      rt22 = rtemp2 * rtemp2;
      
      br = wavenum * wavenum / rr - 1.0 / rrcb;
      bi = wavenum / rrsq;
      
      temp = ar = 2.0 / rrcb - br;
      ai = -3.0 * bi;
      expr = cos (kr);
      
      while (kr > 2 * PI)
	kr -= 2 * PI;
      
      expi = sqrt(1.0 - expr * expr);
      if (kr > PI)
	expi *= -1.0;
      
      ar = temp * expr - ai * expi;
      ai = temp * expi + ai * expr;
      temp = br;
      br = temp * expr - bi * expi;
      bi = temp * expi + bi * expr;

      if (fsc==true) {
	correction=puls(PI/dplZ*(1+nz))*puls(PI/dplX*nx)*puls(PI/dplY*ny);
	ar*=correction; ai*=correction;
	br*=correction; bi*=correction;
      }

      f00r = rt00 * ar + br;
      f00i = rt00 * ai + bi;
      f01r = f10r = rt01 * ar;
      f01i = f10i = rt01 * ai;
      f02r = f20r = rt02 * ar;
      f02i = f20i = rt02 * ai;
      f11r = rt11 * ar + br;
      f11i = rt11 * ai + bi;
      f12r = f21r = rt12 * ar;
      f12i = f21i = rt12 * ai;
      f22r = rt22 * ar + br;
      f22i = rt22 * ai + bi;
      /* end of old code */
      
      /* store results */
      index=index_Xmatrix(i,j,k);
      switch((char) Dcomp) {
      case 0: {
	D2matrix[index].r=f00r;
	D2matrix[index].i=f00i;
	break;
      }
      case 1: {
	D2matrix[index].r=f01r;
	D2matrix[index].i=f01i;
	break;
      }
      case 2: {
	D2matrix[index].r=f02r;
	D2matrix[index].i=f02i;
	break;
      }
      case 3: {
	D2matrix[index].r=f11r;
	D2matrix[index].i=f11i;
	break;
      }
      case 4: {
	D2matrix[index].r=f12r;
	D2matrix[index].i=f12i;
	break;
      }
      case 5: {
	D2matrix[index].r=f22r;
	D2matrix[index].i=f22i;
	break;
      }
      } /* end of switch */
    } /* end of i,j,k loop */
    fftX((REAL *) D2matrix,1,work,0,local_Nz);
    block_transpose((REAL*) D2matrix);
    
    for(i=2*local_z0;i<2*local_z1;i++) {
      x=i;
      for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
	index=index_sliceD2matrix(y,z);
	slice[index].r=slice[index].i=0.0;
      }
      
      for(k=1-Z;k<Z;k++) for(j=1-Y;j<Y;j++) {
	indexfrom=index_garbledX(i,j,k);
	indexto=index_sliceD2matrix(j,k);
	
	slice[indexto].r=D2matrix[indexfrom].r;
	slice[indexto].i=D2matrix[indexfrom].i;
      }
      
      for(k=0;k<Z;k++) for(j=1;j<Y;j++) {
	indexfrom=index_sliceD2matrix(j,k);
	
	/* mirror along y */
	indexto=index_sliceD2matrix(-j,k);
	if (Dcomp==1 || Dcomp==4) {
	  slice[indexto].r=-slice[indexfrom].r;
	  slice[indexto].i=-slice[indexfrom].i;
	}
	else {
	  slice[indexto].r=slice[indexfrom].r;
	  slice[indexto].i=slice[indexfrom].i;
	}
      }
      for(k=1;k<Z;k++) for(j=1-Y;j<Y;j++) {
	indexfrom=index_sliceD2matrix(j,k);
	
	/* mirror along z */
	indexto=index_sliceD2matrix(j,-k);
	if (Dcomp==2 || Dcomp==4) {
	  slice[indexto].r=-slice[indexfrom].r;
	  slice[indexto].i=-slice[indexfrom].i;
	}
	else {
	  slice[indexto].r=slice[indexfrom].r;
	  slice[indexto].i=slice[indexfrom].i;
	}
      }  
      
      fftZ(slice,1,0,1,work);
      transposeYZ(slice,0,gridY,gridZ,0,1);
      fftY(slice,1,0,1,work);
      
      for(z=0;z<DsizeZ;z++) for(y=0;y<DsizeY;y++) {
	indexto=index_tDmatrix(x-2*local_z0,y,z);
	indexfrom=index_slice_zyD2matrix(y,z);
	Dmatrix[indexto+Dcomp].r=-slice[indexfrom].r*invNgrid;
	Dmatrix[indexto+Dcomp].i=-slice[indexfrom].i*invNgrid;
	
      }
      
      
    } /* end slice X */
    /* FFT D */
    printz("."); fflush(stdout);
    /*test_interrupt();*/
  } /* end of Dcomp */
  
  free(D2matrix);
  
  Timing_FFT_Init = clock()-start;
  /* restore some variables to their correct value:
   */
  
  smallY=gridY/2; smallZ=gridZ/2;
  Nsmall=smallX*smallY*smallZ;
  free(slice);
}     
 
void free_FFT_Dmat(void)
{
  free(Dmatrix);
}
