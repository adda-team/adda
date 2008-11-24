/* FILE: matvec.c
 * AUTH: Michel Grimminck
 * DATE: Februari 1995
 */

/* calculate local matrix vector product of decomposed interaction
 * matrix with rk or pk, using a FFT based convolution algorithm
 */

#include <stdio.h>
#include <sys/time.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "types.h"
#include "comm.h"

extern dcomplex *dCvector(int nl,int nh);

extern dcomplex *Dmatrix;

extern int profile[10];
static char *Fname = "MatVec_nim";

void test_interrupt(void)
{
  struct tm *nowtime;
  FILE *test;
  char filename[40];
  time_t t;
#ifdef PARALLEL
  return;
#else
  time(&t);
  nowtime=localtime(&t);
  sprintf(filename,"/tmp/scatter.lock.%i",nowtime->tm_mday);
  do {
    test=fopen(filename,"r");
    if (test!=NULL) {
      fclose(test);
      printf("Scatter is locked!\n");
      system("ls -l /tmp/scatter.lock.*");
      system("sleep 1200");
    }
  } while(test!=NULL);
#endif
}

int round(double x)
     /* round a double to nearest integer */
{
  if (x<0) x-=1.0;
  return((int) x);
}

int index_slice_yz(int y,int z)
{
  extern int gridY,gridZ;
  
  if (y<0) y+=gridY;
  if (z<0) z+=gridZ;
  
  return(z*gridY+y);
}

int index_slice_zy(int y,int z)
{
  extern int gridY,gridZ;
  
  if (y<0) y+=gridY;
  if (z<0) z+=gridZ;
  
  return(y*gridZ+z);
}

int index_garbledX(int x,int y,int z)
{
  extern smallX,smallY,smallZ;
  int zbox;
  
  if (x<0) x+=smallX;
  if (y<0) y+=smallY;
  if (z<0) z+=smallZ;
  
  zbox=z/local_Nz;
  
  x=x%local_xs_unit;
  if (abs(x)>=smallX || abs(y)>=smallY || abs(z)>=smallZ) printf("ERROR:index_Xmatrix: %i %i %i\n",x,y,z);
  return((z%local_Nz)*smallX*smallY+zbox*local_xs_unit+y*smallX+x);
}

int index_Xmatrix(int x,int y,int z)
{
  extern int abs(int);
  extern smallX,smallY,smallZ;
  
  if (x<0) x+=smallX;
  if (y<0) y+=smallY;
  if (z<0) z+=smallZ;
  z-=local_z0;
  if (abs(x)>=smallX || abs(y)>=smallY || abs(z)>=smallZ) printf("ERROR:index_Xmatrix: %i %i %i\n",x,y,z);
  return((z*smallX*smallY+y*smallX+x));
}

/************* MatVec{Inp}_nim **************/

/* This function implements both MatVec_nim and MatVecAndInp_nim.
 * the difference is that when we want to calculate the inproduct
 * as well, we pass 'inprod' as a non-null pointer. if 'inprod' is
 * a NULL, we don't calculate it.
 */
void
both_MatVec  (dcomplex *argvec,     /* the argument vector */
	      dcomplex *resultvec,  /* the result vector */
	      double *inprod,       /* the result inner product */
	      int nldip,            /* number of local dipoles */
	      int ndip,             /* total number of dipoles */
	      int procid,           /* processor ring number */
	      int nproc,            /* total number of processors */
	      double wavenum,       /* wavenumber */
	      int her)              /* 0 for non-hermitic, 1 for hermetic */
{
  extern int smallX,smallY,smallZ;
  int i, j;
  int base;                     /* base vector */
  int reminder = ndip % nproc;  /* number of proc containing one dip more */
  
  double f00r, f01r, f02r, f11r, f12r, f22r;
  double f00i, f01i, f02i, f11i, f12i, f22i;
  double temp,in;
  double rtemp0,rtemp1,rtemp2;
  dcomplex *slice0,*slice1,*slice2;
  /* FFT variables: */
  double invgridspaceX,invgridspaceY,invgridspaceZ;
  dcomplex *Xmatrix0,*Xmatrix1,*Xmatrix2,*swap,x0,x1,x2;
  int index,x,y,z,Xcomp,sliceX,endofX,mat;
  REAL work[100];
  extern int gridX,gridY,gridZ,Ngrid,Nsmall;
  extern int local_Nsmall;
  extern int symD,symX;
  extern double gridspaceX,gridspaceY,gridspaceZ;
  extern int fsc;
  extern short int *position;
  extern int Dmatrix_nel;
  extern int minY,maxY,minZ,maxZ;
  extern doublecomplex cc[10];
  extern int *material;
  extern int Nmat;
  extern unsigned long Timing_OneIterComm;
  int tstart,tstop;
  static mvcounter=0;
  int start;
  
  
  test_interrupt();
  profile[0]-=extime();
  /* FFT_matvec code */
  if (inprod)
    *inprod = 0.0;
  invgridspaceX=1.0/gridspaceX;
  invgridspaceY=1.0/gridspaceY;
  invgridspaceZ=1.0/gridspaceZ;
  
  /* allocate memory for Xmatrix */
  if ((Xmatrix0 = dCvector(0, local_Nsmall)) == NULL) {
    LogError (EC_ERROR, Fname,
	      "processor %d, ringID %d, could not malloc Xmatrix",
	      0, 0);
    AbortServer (1);
  }
  /* allocate memory for Xmatrix */
  if ((Xmatrix1 = dCvector(0, local_Nsmall)) == NULL) {
    LogError (EC_ERROR, Fname,
	      "processor %d, ringID %d, could not malloc Xmatrix",
	      0, 0);
    AbortServer (1);
  }
  if ((Xmatrix2 = dCvector(0, local_Nsmall)) == NULL) {
    LogError (EC_ERROR, Fname,
	      "processor %d, ringID %d, could not malloc Xmatrix",
	      0, 0);
    AbortServer (1);
  }
  if ((slice0 = dCvector(0, gridY*gridZ)) == NULL) {
    LogError (EC_ERROR, Fname,
	      "processor %d, ringID %d, could not malloc slice",
	      0, 0);
    AbortServer (1);
  }
  if ((slice1 = dCvector(0, gridY*gridZ)) == NULL) {
    LogError (EC_ERROR, Fname,
	      "processor %d, ringID %d, could not malloc slice",
	      0, 0);
    AbortServer (1);
  }
  if ((slice2 = dCvector(0, gridY*gridZ)) == NULL) {
    LogError (EC_ERROR, Fname,
	      "processor %d, ringID %d, could not malloc slice",
	      0, 0);
    AbortServer (1);
  }
  
  for(Xcomp=0;Xcomp<3;Xcomp++) { 
    /* for each x,y,z component do:
     * X=0;
     * fill X with cc*argvec;
     * FFT x
     */
    
    /* fill Xmatrix with 0.0 */
    if (Xcomp==0) for(i=0;i<local_Nsmall;i++) Xmatrix0[i].r=Xmatrix0[i].i=0.0;
    if (Xcomp==1) for(i=0;i<local_Nsmall;i++) Xmatrix1[i].r=Xmatrix1[i].i=0.0;
    if (Xcomp==2) for(i=0;i<local_Nsmall;i++) Xmatrix2[i].r=Xmatrix2[i].i=0.0;
    
    /* transform from coordinates to grid and multiply with coupling constant */
    
    for (i = local_d0; i < local_d1; i++) {
      /* fill grid with E*coupleconstant */
      j=3*(i-local_d0); index=index_Xmatrix(position[j],
					    position[j+1],
					    position[j+2]);
      mat=material[i-local_d0];
      switch(Xcomp) {
      case 0: { /* x field */
	if (!her)
	  {
	    Xmatrix0[index].r=cc[mat].r*argvec[j].r - cc[mat].i*argvec[j].i;
	    Xmatrix0[index].i=cc[mat].r*argvec[j].i + cc[mat].i*argvec[j].r;
	  }
	else
	  {
	    Xmatrix0[index].r=argvec[j].r;
	    Xmatrix0[index].i=-argvec[j].i;
	  }
	break;
      }
      case 1: { /* y field */
	if (!her)
	  {
	    Xmatrix1[index].r=cc[mat].r*argvec[j+1].r - cc[mat].i*argvec[j+1].i;
	    Xmatrix1[index].i=cc[mat].r*argvec[j+1].i + cc[mat].i*argvec[j+1].r;
	  }
	else
	  {
	    Xmatrix1[index].r=argvec[j+1].r;
	    Xmatrix1[index].i=-argvec[j+1].i;
	  }
	break;
      }
      case 2: { /* z field */
	if (!her)
	  {
	    Xmatrix2[index].r=cc[mat].r*argvec[j+2].r - cc[mat].i*argvec[j+2].i;
	    Xmatrix2[index].i=cc[mat].r*argvec[j+2].i + cc[mat].i*argvec[j+2].r;
	  }
	else
	  {
	    Xmatrix2[index].r=argvec[j+2].r;
	    Xmatrix2[index].i=-argvec[j+2].i;
	  }
	break;
      }
      }
    }
    
    /* FFT X */
    if (Xcomp==0) {
      fftX((REAL *) Xmatrix0,1,work,0,local_Nz);
      block_transpose((REAL*) Xmatrix0);
    }
    else if (Xcomp==1) {
      fftX((REAL *) Xmatrix1,1,work,0,local_Nz);
      block_transpose((REAL*)Xmatrix1);
    }
    else if (Xcomp==2) {
      fftX((REAL *) Xmatrix2,1,work,0,local_Nz);
      block_transpose((REAL*)Xmatrix2);
    }
  }
  
  /* do the product D~*X~ */
  for(sliceX=2*local_z0;sliceX<2*local_z1;sliceX++) {
    profile[1]-=extime();
    x=sliceX;
    /* clear slice */
    for(i=0;i<gridZ*gridY;i++) {
      slice0[i].r=slice0[i].i=0.0;
      slice1[i].r=slice1[i].i=0.0;
      slice2[i].r=slice2[i].i=0.0;
    }
#ifndef PARALLEL
    for(y=minY;y<=maxY;y++) for(z=minZ;z<=maxZ;z++) {
      i=index_slice_zy(y,z);
      j=index_Xmatrix(sliceX,y,z);
      slice0[i].r=Xmatrix0[j].r;  slice0[i].i=Xmatrix0[j].i;
      slice1[i].r=Xmatrix1[j].r;  slice1[i].i=Xmatrix1[j].i;
      slice2[i].r=Xmatrix2[j].r;  slice2[i].i=Xmatrix2[j].i;
    }
#else
    for(y=minY;y<=maxY;y++) for(z=minZ;z<=maxZ;z++) {
      i=index_slice_zy(y,z);
      j=index_garbledX(sliceX,y,z);
      slice0[i].r=Xmatrix0[j].r;  slice0[i].i=Xmatrix0[j].i;
      slice1[i].r=Xmatrix1[j].r;  slice1[i].i=Xmatrix1[j].i;
      slice2[i].r=Xmatrix2[j].r;  slice2[i].i=Xmatrix2[j].i;
    }
#endif
    profile[1]+=extime();
    
    fftZ(slice0,1,0,1,work);
    transposeYZ(slice0,0,gridY,gridZ,0,1);
    fftY(slice0,1,0,1,work);
    
    fftZ(slice1,1,0,1,work);
    transposeYZ(slice1,0,gridY,gridZ,0,1);
    fftY(slice1,1,0,1,work);
    
    fftZ(slice2,1,0,1,work);
    transposeYZ(slice2,0,gridY,gridZ,0,1);
    fftY(slice2,1,0,1,work);
    profile[2]-=extime();
    
    ram_memory(Dmatrix,Dmatrix_nel*sizeof(dcomplex));
    for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
      /* store old X */
      i=index_slice_yz(y,z);
      x0.r=slice0[i].r; x0.i=slice0[i].i;
      x1.r=slice1[i].r; x1.i=slice1[i].i;
      x2.r=slice2[i].r; x2.i=slice2[i].i;
      
      if (symD==true) { /* symmetric D, only 2 components in memory */
	int j1,j2,j3,j4,j5;
	j1=index_tDmatrix(x,y,z); j2=index_tDmatrix(y,x,z);
	j3=index_tDmatrix(x,z,y); j4=index_tDmatrix(z,y,x);
	j5=index_tDmatrix(y,z,x);
	/* get D~ in f00r...f22i */
	f00r=Dmatrix[j1].r;   f00i=Dmatrix[j1].i;
	f01r=Dmatrix[j1+1].r; f01i=Dmatrix[j1+1].i;   
	f02r=Dmatrix[j3+1].r; f02i=Dmatrix[j3+1].i;   
	f11r=Dmatrix[j2].r; f11i=Dmatrix[j2].i;
	f12r=Dmatrix[j5+1].r; f12i=Dmatrix[j5+1].i;   
	f22r=Dmatrix[j4].r; f22i=Dmatrix[j4].i;
      }
      else { /* assym. D, 6 components in memory */
	j=index_tDmatrix(x-2*local_z0,y,z);
	f00r=Dmatrix[j].r;        f00i=Dmatrix[j].i;
	f01r=Dmatrix[j+1].r;      f01i=Dmatrix[j+1].i;
	f02r=Dmatrix[j+2].r;      f02i=Dmatrix[j+2].i;
	f11r=Dmatrix[j+3].r;      f11i=Dmatrix[j+3].i;
	f12r=Dmatrix[j+4].r;      f12i=Dmatrix[j+4].i;
	f22r=Dmatrix[j+5].r;      f22i=Dmatrix[j+5].i;
      }
#ifndef PARALLEL
      if (x>gridX/2) {
	f01r=-f01r; f01i=-f01i;
	f02r=-f02r; f02i=-f02i;
      }
#endif
      if (y>gridY/2) {
	f01r=-f01r; f01i=-f01i;
	f12r=-f12r; f12i=-f12i;
      }
      if (z>gridZ/2 && fsc==false) {
	f02r=-f02r; f02i=-f02i;
	f12r=-f12r; f12i=-f12i;
      }
      
      
      slice0[i].r  = f00r * x0.r - f00i * x0.i + f01r * x1.r - f01i * x1.i +
	f02r * x2.r - f02i * x2.i;
      slice0[i].i  = f00r * x0.i + f00i * x0.r + f01r * x1.i + f01i * x1.r +
	f02r * x2.i + f02i * x2.r;
      
      slice1[i].r  = f01r * x0.r - f01i * x0.i + f11r * x1.r - f11i * x1.i +
	f12r * x2.r - f12i * x2.i;
      slice1[i].i  = f01r * x0.i + f01i * x0.r + f11r * x1.i + f11i * x1.r +
	f12r * x2.i + f12i * x2.r;
      
      slice2[i].r  = f02r * x0.r - f02i * x0.i + f12r * x1.r - f12i * x1.i +
	f22r * x2.r - f22i * x2.i;
      slice2[i].i  = f02r * x0.i + f02i * x0.r + f12r * x1.i + f12i * x1.r +
	f22r * x2.i + f22i * x2.r;
    }
    profile[2]+=extime();
    fftY(slice0,-1,0,1,work);
    transposeYZ(slice0,0,gridY,gridZ,0,1);
    fftZ(slice0,-1,0,1,work);
    
    fftY(slice1,-1,0,1,work);
    transposeYZ(slice1,0,gridY,gridZ,0,1);
    fftZ(slice1,-1,0,1,work);
    
    fftY(slice2,-1,0,1,work);
    transposeYZ(slice2,0,gridY,gridZ,0,1);
    fftZ(slice2,-1,0,1,work);
    
    /* copy slice back to Xmatrix */
#ifndef PARALLEL
    for(y=minY;y<=maxY;y++) for(z=minZ;z<=maxZ;z++) {
      i=index_slice_zy(y,z);
      j=index_Xmatrix(sliceX,y,z);
      Xmatrix0[j].r=slice0[i].r;  Xmatrix0[j].i=slice0[i].i;
      Xmatrix1[j].r=slice1[i].r;  Xmatrix1[j].i=slice1[i].i;
      Xmatrix2[j].r=slice2[i].r;  Xmatrix2[j].i=slice2[i].i;
    }
#else*
    for(y=minY;y<=maxY;y++) for(z=minZ;z<=maxZ;z++) {
      i=index_slice_zy(y,z);
      j=index_garbledX(sliceX,y,z);
      Xmatrix0[j].r=slice0[i].r;  Xmatrix0[j].i=slice0[i].i;
      Xmatrix1[j].r=slice1[i].r;  Xmatrix1[j].i=slice1[i].i;
      Xmatrix2[j].r=slice2[i].r;  Xmatrix2[j].i=slice2[i].i;
    }
#endif
  } /* end of sliceX */
  
  profile[3]-=extime();
  for(Xcomp=0;Xcomp<3;Xcomp++) { 
    /* for each x,y,z component do:
     * FFT x[]
     * restore dipole array
     */
    
    /* FFT back the result and store results in resultvec */
    if (Xcomp==0) block_transpose((REAL*)Xmatrix0);
    if (Xcomp==1) block_transpose((REAL*)Xmatrix1);
    if (Xcomp==2) block_transpose((REAL*)Xmatrix2);
    
    if (Xcomp==0) {ram_memory(Xmatrix0,local_Nsmall*sizeof(dcomplex)); fftX(Xmatrix0,-1,work,0,local_Nz);}
    if (Xcomp==1) {ram_memory(Xmatrix1,local_Nsmall*sizeof(dcomplex)); fftX(Xmatrix1,-1,work,0,local_Nz);}
    if (Xcomp==2) {ram_memory(Xmatrix2,local_Nsmall*sizeof(dcomplex)); fftX(Xmatrix2,-1,work,0,local_Nz);}
    
    for (i = local_d0; i < local_d1; i++) if (material[i-local_d0]<Nmat-1) {
      /* restore resultvec from grid */
      
      j=3*(i-local_d0); index=index_Xmatrix(position[j],
					    position[j+1],
					    position[j+2]);
      
      mat=material[i-local_d0];
      if (Xcomp==0) {
	if (!her)
	  {
	    resultvec[j].r=argvec[j].r + Xmatrix0[index].r;
	    resultvec[j].i=argvec[j].i + Xmatrix0[index].i;
	  }
	else
	  {
	    resultvec[j].r=argvec[j].r + cc[mat].r*Xmatrix0[index].r - cc[mat].i*Xmatrix0[index].i;
	    resultvec[j].i=argvec[j].i - cc[mat].r*Xmatrix0[index].i - cc[mat].i*Xmatrix0[index].r;
	  }
      }
      else if (Xcomp==1) {
	if (!her)
	  {
	    resultvec[j+1].r=argvec[j+1].r + Xmatrix1[index].r;
	    resultvec[j+1].i=argvec[j+1].i + Xmatrix1[index].i;
	  }
	else
	  {
	    resultvec[j+1].r=argvec[j+1].r + cc[mat].r*Xmatrix1[index].r - cc[mat].i*Xmatrix1[index].i;
	    resultvec[j+1].i=argvec[j+1].i - cc[mat].r*Xmatrix1[index].i - cc[mat].i*Xmatrix1[index].r;
	  }
      }
      else if (Xcomp==2) {
	if (!her)
	  {
	    resultvec[j+2].r=argvec[j+2].r+Xmatrix2[index].r;
	    resultvec[j+2].i=argvec[j+2].i+Xmatrix2[index].i;
	  }
	else
	  {
	    resultvec[j+2].r=argvec[j+2].r + cc[mat].r*Xmatrix2[index].r - cc[mat].i*Xmatrix2[index].i;
	    resultvec[j+2].i=argvec[j+2].i - cc[mat].r*Xmatrix2[index].i - cc[mat].i*Xmatrix2[index].r;
	  }
      }
      if (inprod) 
	*inprod+=resultvec[j+Xcomp].r*resultvec[j+Xcomp].r + resultvec[j+Xcomp].i*resultvec[j+Xcomp].i;
    } 
    
  }profile[3]+=extime();

  free(slice2); free(slice1); free(slice0);
  free(Xmatrix2); free(Xmatrix1); free(Xmatrix0);
  
  profile[0]+=extime();
  if (inprod)
    {
      tstart=extime();
      my_inner_product(inprod);
      tstop=extime();
      Timing_OneIterComm+=tstop-tstart;
    }
  mvcounter++;
}







