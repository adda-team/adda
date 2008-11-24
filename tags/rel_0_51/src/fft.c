/* FILE: fft.c
 * AUTH: Michel Grimminck (using old code by Alfons Hoekstra)
 * DESCR: Initialization of all FFT for matrix-vector products
 *        and FFT itself     
 *        A lot of indirect indexing used - way to optimize
 *
 *        Currently is developed by Maxim Yurkin
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "types.h"
#include "comm.h"

#define FFT

extern doublecomplex *dCvector(int nl,int nh);

extern clock_t Timing_FFT_Init;

extern double gridspace;
extern int prognose;

doublecomplex *Dmatrix;      /* size: Dmatrix[gridX][gridY][gridZ][6] */
doublecomplex *D2matrix;
int DsizeX,DsizeY,DsizeZ;    /* the size of the 'matrix' D */ 
int gridX,gridY,gridZ,Ngrid; /* the size of the 'matrix' X */
int gridXY,gridYZ,gridXZ;
int D2sizeX,D2sizeY,D2sizeZ;
int DsizeXY,DsizeYZ,DsizeXZ;
int smallX,smallY,smallZ;    /* the size of the reduced matrix X */
int local_Nsmall;
int NDcomp;                  /* number of components of D */
int Dmatrix_nel;
int minY,maxY,minZ,maxZ;

char *useY;

/*============================================================*/

int index_tDmatrix(int x,int y,int z)
{
  extern int abs(int);
  
  if (x>=DsizeX) x=gridX-x;
  if (y>=DsizeY) y=gridY-y;
  if (z>=DsizeZ) z=gridZ-z;
  
  if (x>=DsizeX || y>=DsizeY || z>=DsizeZ) 
    LogError(EC_ERROR,ALL,POSIT,"index_tDmatrix:%i %i %i",x,y,z);
 
  return(NDcomp*(x*DsizeYZ+z*DsizeY+y));
}

/*============================================================*/

int index_D2matrix(int x,int y,int z)
{
  extern int abs(int);
  
  if (x<0) x=D2sizeX+x;
  if (y<0) y=D2sizeY+y;
  if (z<0) z=D2sizeZ+z;
  
  if (x>=gridX || y>=gridY || z>=gridZ)
    LogError(EC_ERROR,ALL,POSIT,"index_D2matrix:%i %i %i",x,y,z);
  
  return(z*D2sizeX*D2sizeY+y*D2sizeX+x);
}

/*============================================================*/

int index_sliceD2matrix(int y,int z)
{
  if (y<0) y=gridY+y;
  if (z<0) z=gridZ+z;
  if (y>=gridY || z>=gridZ)
    LogError(EC_ERROR,ALL,POSIT,"index_sliceD2matrix:%i %i",y,z);
  
  return(y*gridZ+z);
}

/*============================================================*/

int index_slice_zyD2matrix(int y,int z)
{
  if (y<0) y=gridY+y;
  if (z<0) z=gridZ+z;
  if (y>=gridY || z>=gridZ) 
    LogError(EC_ERROR,ALL,POSIT,"index_slicezyD2matrix:%i %i",y,z);
  
  return(z*gridY+y);
}

/*============================================================*/

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void four1(data,nn,isign)
     double data[];
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
    theta=TWO_PI/(isign*mmax);
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

/*============================================================*/

void transposeYZ(double *data,int X,int Y,int Z,int start,int end)
     /* optimised routine to transpose y and z.
      * We know smallX and smallY are even, so we can do
      * two points at the same time
      */
{
  int x,y,z,a,b;
  double *temp,*work,*nowtemp,*nowwork;
  
  if ((temp=(double *) malloc(2*sizeof(double)*Z*Y))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"transposeYZ: could not malloc temp");
  
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
    memcpy(work,temp,2*sizeof(double)*Z*Y);
  }
  free(temp);
}

/*============================================================*/

void fftX(double *data,int isign,doublecomplex *dummy,int z0,int z1)
{
  int nn=smallX,dim=1,iform=1,is,x,y,z,inc=1,jump=nn,lot=1;
  int ifax[20];
  double *trigs,*work;
  
  if ((trigs = (double *) malloc(2*nn*sizeof(double))) == NULL) 
    LogError(EC_ERROR,ALL,POSIT,"fftX: could not malloc trigs");
  if ((work = (double *) malloc(2*nn*lot*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"fftX: could not malloc work");
  
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

/*============================================================*/

void fftY(double *data,int isign,int start,int end,doublecomplex *dummy)
{
  int nn=gridY,dim=1,iform=1,is,x,y,z,inc=1,jump=nn,lot=gridZ/2;
  int ifax[20];
  double *trigs,*work;
  
  if ((trigs = (double *) malloc(2*nn*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"fftY: could not malloc trigs");
  if ((work = (double *) malloc(2*nn*lot*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"fftY: could not malloc work");
  
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

/*============================================================*/

void fftZ(double *data,int isign,int start,int end,doublecomplex *dummy)
{
  int nn=gridZ,dim=1,iform=1,is,x,y,z,inc=1,jump=nn,lot=1;
  int ifax[20];
  double *trigs,*work;
  
  if ((trigs = (double *) malloc(2*nn*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"fftZ: could not malloc trigs");
  if ((work = (double *) malloc(2*nn*lot*sizeof(double))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"fftZ: could not malloc work");
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

/*============================================================*/

void test_fft () 
{ /* perform FFT many times and measures time */
  
  double *data, t;
  int i, NN=200, ofs;
  clock_t tstart;
  
  data=dvector(0,2*gridY*gridZ+NN);
  for (i=0;i<=2*gridY*gridZ+NN;i++) {
    data[i]=(double)rand();
  }
  
  tstart=clock();
  for (i=0;i<NN;i++) {
    ofs=floor(NN*((double)random())/RAND_MAX);
    fftZ(&data[ofs],1,0,1,NULL);
  }
  t=TO_SEC(clock()-tstart);
  
  printz("%d complex FFTs of size %d\n"\
         "total time = %g\n"\
         "per 1 FFT  = %g\n",gridY*NN,gridZ,t,t/(gridY*NN));
}
                                                               
/*=============================================================*/

int fft_fit(int x, int div)
     /* find the first number >=x divisible by 2,3 and 5 only
      * the resulting number is divisible by 2 and div
      */
{
  int y;
  
  do {
    y=x;
    while (y%2==0) y/=2;
    while (y%3==0) y/=3;
    while (y%5==0) y/=5;
    if (y==1 && x%2==0 && x%div==0) return(x);
    x++;
  } while(true);
}

/*============================================================*/

void fill_use()
{
  int i;

  /* allocate memory for use arrays */
  if ((useY = (char *) malloc(sizeof(char)*gridY)) == NULL) 
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc useY");
  for(i=0;i<gridY;i++) useY[i]=true;
}

/*============================================================*/

void setup_use(double **dip_coord,int ndip)
{
  int i,j,x,y,z;
  double rtemp0,rtemp1,rtemp2;
  int minX,maxX;
  extern int boxX,boxY,boxZ;

  minY=100000,minZ=100000,maxY=-100000,maxZ=-100000;
  
  for(i=0;i<gridY;i++) useY[i]=false;
  
  minX=-boxX/2; maxX=minX+boxX-1;
  minY=-boxY/2; maxY=minY+boxY-1;
  minZ=0; maxZ=boxZ-1; /* ??? */
  for(i=minY;i<=maxY;i++) {
    j=i; if (j<0) j+=gridY;
    useY[j]=true;
  }
  printf("setup_use: %i, %i %i %i %i\n",ringid,minY,maxY,minZ,maxZ);
}

/**************************************************************************/

void CalcInterTerm(int i, int j, int k, int mu, int nu, doublecomplex *result, double wavenum)
{
  extern double *tab1,*tab2,*tab3,*tab4,*tab5,*tab6,*tab7,*tab8,*tab9,*tab10; /* tables of integrals */
  extern double dpl;
  extern double prop[3];
  extern int IntRelation;
  extern double kd;
  extern doublecomplex ref_index[10];
  extern int **tab_index;

  double rr, rtemp[3], qvec[3], q2[3], invr, invr3, qavec[3], av[3];
  double rr2, kr, kr2, kr3, kd2, q4, rn, rn2, rn4;
  double temp, qmunu, qa, qamunu, invrn, invrn2, invrn3, invrn4, dmunu;
  doublecomplex expval, br, br1, m, m2, Gf1, Gm0, Gm1, Gc1, Gc2;
  int ind0, ind1, ind2, ind2m, ind3, ind4, ind5, indmunu;
  int sigV[3], ic, sig, ivec[3], ord[3], invord[3];
  double t3q, t3a, t4q, t4a, t5tr, t5aa, t6tr, t6aa;
 
  if (i==0 && j==0 && k==0) { /* self interaction */  
    result->r=0.0;
    result->i=0.0;
    return;
  }
  
/*  int pr;
  pr=(i==5 && j==4 && k==3);           
  if (pr) printz("%d,%d: ",mu,nu);  */   /* for debugging */
  
  rtemp[0]=i*gridspace;
  rtemp[1]=j*gridspace;
  rtemp[2]=k*gridspace;
  
  rr2 = DotProd(rtemp,rtemp);
  rr = sqrt(rr2);
  invr = 1/rr;
  invr3 = invr*invr*invr;
  MultScal(invr,rtemp,qvec);
  kr = wavenum * rr;
  kr2 = kr*kr;
  qmunu=qvec[mu]*qvec[nu];
  
  br.r=(3-kr2)*qmunu;
  br.i=-3*kr*qmunu;
  if(mu==nu) {                       /* br=delta[mu,nu]*(-1+ikr+kr^2)  */
    br.r+=kr2-1;                     /*    -qmunu*(-3+3ikr+kr^2)       */
    br.i+=kr;
  }
  
  cExp(kr,&expval);
  cMultReal(invr3,expval,&expval);     /* expval=Exp(ikr)/rr^3 */
  cMult(br,expval,result);             /* result=Gp            */
  
  if (IntRelation == SOrd) {
    kd2=kd*kd;
    kr3=kr2*kr;
    rn=rr/gridspace;              /* normalized r */
    m=ref_index[0];       /* only one refractive index can be used for FFT-compatible algorithm !!! */  
    cSquare(m,&m2);
    qa=DotProd(qvec,prop);
    
    qamunu=qvec[mu]*prop[nu];          /* qamunu=qvec[mu]*prop[nu] + qvec[nu]*prop[mu] */
    if (mu==nu) qamunu*=2;
    else qamunu+=qvec[nu]*prop[mu];
    
    if (kr*rn < 1) {       /* G close */
      memcpy(av,prop,3*sizeof(double));   /* av is copy of propagation vector */
      ivec[0]=i;
      ivec[1]=j;
      ivec[2]=k;
      /* transformation of negative coordinates */
      for (ic=0;ic<3;ic++) {
        if (ivec[ic]<0) {
          sigV[ic]=-1;
          av[ic]*=-1;
          qvec[ic]*=-1;
          ivec[ic]*=-1;
        } 
        else sigV[ic]=1;
      }
      i=ivec[0];
      j=ivec[1];
      k=ivec[2];
      sig=sigV[mu]*sigV[nu];           /* sign of some terms below */
      /* transformation to case i>=j>=k>=0 */
        /* building of ord */
      if (i>=j) {
        if (i>=k) {
          ord[0]=0;
          if (j>=k) { 
            ord[1]=1; 
            ord[2]=2;
          }
          else { 
            ord[1]=2; 
            ord[2]=1; 
          }
        }
        else {
          ord[0]=2;
          ord[1]=0; 
          ord[2]=1;          /* ord[x] is x-th largest coordinate (0-th - the largest) */
        }
      }
      else {
        if (i>=k) {
          ord[0]=1;
          ord[1]=0; 
          ord[2]=2;
        }
        else {
          ord[2]=0;
          if (j>=k) { 
            ord[0]=1; 
            ord[1]=2;
          }
          else { 
            ord[0]=2; 
            ord[1]=1; 
          }
        }
      }
      /* change parameters according to coordinate transforms */
      permutate(qvec,ord);
      permutate(av,ord);
      permutate_i(ivec,ord);
      i=ivec[0];
      j=ivec[1];
      k=ivec[2];
      /* compute inverse permutation */
      memcpy(invord,ord,3*sizeof(int));
      permutate_i(invord,ord);
      if (invord[0]==0 && invord[1]==1 && invord[2]==2) memcpy(invord,ord,3*sizeof(int));
      /* compute transformed indices mu and nu */
      mu=invord[mu];
      nu=invord[nu];
      /* indexes for tables of different dimensions */
      indmunu=mu+nu;
      if (mu==2 || nu==2) indmunu++;  /* indmunu is a number of component[mu,nu] in symmetric matrix */
            
      ind0=tab_index[i][j]+k;
      ind1=3*ind0;
      ind2m=6*ind0;
      ind2=ind2m+indmunu;
      ind3=3*ind2;
      ind4=6*ind2;
      /* computing several quantities with table integrals */
      t3q=DotProd(qvec,&(tab3[ind1]));
      t3a=DotProd(av,&(tab3[ind1]));
      t4q=DotProd(qvec,&(tab4[ind3]));
      t4a=DotProd(av,&(tab4[ind3]));
      t5tr=TrSym(&(tab5[ind2m]));
      t5aa=QuadForm(&(tab5[ind2m]),av);
      t6tr=TrSym(&(tab6[ind4]));
      t6aa=QuadForm(&(tab6[ind4]),av);
      /* computing Gc0 */
      temp=kr/24;                      /* temp = kr/12 */
      
      br.r=sig*(3*(tab10[ind2]/2+tab8[ind2])-2*t4q-t6tr)+temp*qmunu*kr;
      br.i=3*temp*qmunu;                                   /* br=delta[mu,nu]*(-I7-I9/2-kr*(i+kr)/24+2*t3q+t5tr) */
      if (mu==nu) {                                        /*   -(-3I8[mu,nu]-3I10[mu,nu]/2-qmunu*kr*(i+kr)/24   */
        br.r+=2*t3q+t5tr-temp*kr-tab9[ind0]/2-tab7[ind0];  /*     +2*t4q+t6tr)                                   */
        br.i-=temp;
      }
      
      cMultReal(kd2,br,&br);          /* br*=kd^2 */
      
      br.r+=sig*tab2[ind2]*(3-kr2);   /* br+=I1*delta[mu,nu]*(-1+ikr+kr^2) */
      br.i-=sig*tab2[ind2]*3*kr;      /*     -sig*I2[mu,nu]*(-3+3ikr+kr^2) */
      if (mu==nu) {                   
        br.r+=tab1[ind0]*(kr2-1);
        br.i+=tab1[ind0]*kr;
      }
      
      cMult(expval,br,result);        /* Gc0=expval*br */
      /* computing Gc1 */
      br.r=6*qmunu;        
      br.i=-kr*qmunu;      /* br=(kd*kr/24)*(qa*(delta[mu,nu]*(-2+ikr)-qmunu*(-6+ikr)) */
      if (mu==nu) {        /*                -qamunu)                                  */
        br.r-=2;
        br.i+=kr;
      }
      cMultReal(qa,br,&br);
      br.r-=qamunu;
      cMultReal(2*temp*kd,br,&br);

      br1.r=3*sig*t4a;
      br1.i=-kr*br1.r;
      if (mu==nu) {        /*  br1=(d/r)*(delta[mu,nu]*t3h*(-1+ikr)-sig*t4h*(-3+3ikr)) */
        br1.r-=t3a;
        br1.i+=t3a*kr;
      }
      cMultReal(1/rn,br1,&br1);
     
      cAdd(br,br1,&Gc1);
      cMult(m,Gc1,&Gc1);
      cMultReal(kd,Gc1,&Gc1);
      cMult(expval,Gc1,&Gc1);
      cMult_i(&Gc1);                   /* Gc1=expval*i*m*kd*(br1+br) */
      /* computing Gc2 */
      br.r=-kr*qmunu;        
      br.i=-3*qmunu;      /* br=delta[mu,nu]*t5aa-3*sig*t6aa                 */
      if (mu==nu) {       /*    -(kr/12)*(delta[mu,nu]*(i+kr)-qmunu*(3i+kr)) */
        br.r+=kr;
        br.i+=1;
      }
      cMultReal(-2*temp,br,&br);
      br.r-=3*sig*t6aa;
      if (mu==nu) br.r+=t5aa;

      cMult(m2,br,&Gc2);
      cMultReal(kd2/2,Gc2,&Gc2);     /* Gc2=expval*(kd^2/2)*m^2*br */
      cMult(expval,Gc2,&Gc2);
      
      cAdd(Gc1,Gc2,&Gc1);
      cAdd(Gc1,*result,result);        /* result = Gc0 + Gc1 + Gc2 */
    }
    else {              /* G median and G far */
      temp=kd2/24;             /* temp=kd^2/24;        */
      br.r=1-(1+m2.r)*temp;    /* br=1-(1+m^2)*kd^2/24 */
      br.i=-m2.i*temp;
      cMult(*result,br,result);   /* result = Gp*br */

      br.r=(6*kr2-15)*qmunu;
      br.i=(15*kr-kr3)*qmunu;
      if(mu==nu) {                       /* br={delta[mu,nu]*(3-3ikr-2kr^2+ikr^3)  */
        br.r+=3-2*kr2;                   /*    -qmunu*(15-15ikr-6kr^2+ikr^3)}*qa   */
        br.i+=kr3-3*kr;                  /*    +qamunu*(3-3ikr-kr^2)               */
      }
      cMultReal(qa,br,&br);
      br.r+=(3-kr2)*qamunu;
      br.i-=3*kr*qamunu;

      temp*=2/kr;                /* temp = kd^2/(12*kr) */
      cMult(m,br,&Gf1);
      cMultReal(temp,Gf1,&Gf1);
      cMult(expval,Gf1,&Gf1);
      cMult_i(&Gf1);                     /* Gf1=expval*i*m*temp*br */
      cAdd(Gf1,*result,result);          /* result = Gf  */
      
      if (kr < 1) {                      /* G median */
        vMult(qvec,qvec,q2);
        q4=DotProd(q2,q2);
        invrn=1/rn;
        invrn2=invrn*invrn;
        invrn3=invrn2*invrn;
        invrn4=invrn2*invrn2;

        temp=qmunu*(33*q4-7-12*(q2[mu]+q2[nu]));     
        if (mu == nu) temp+=(1-3*q4+4*q2[mu]);
        temp*=7*invrn4/64;
        br.r=-1;
        br.i=kr;
        cMultReal(temp,br,&Gm0);
        cMult(expval,Gm0,&Gm0);                     /* Gm0=expval*br*temp */
        
        vMult(qvec,prop,qavec);
        if (mu == nu) dmunu=1;
        else dmunu=0;
        temp=3*qa*(dmunu-7*qmunu)+6*dmunu*qvec[mu]*prop[mu]-7*(dmunu-9*qmunu)*DotProd(qavec,q2)+
          3*(prop[mu]*qvec[nu]*(1-7*q2[mu])+prop[nu]*qvec[mu]*(1-7*q2[nu]));
        temp*=kd*invrn3/48;
        cMultReal(temp,m,&Gm1);
        cMult_i(&Gm1);
        cMult(expval,Gm1,&Gm1);    /* Gm1=expval*i*m*temp */
        
        cAdd(Gm0,Gm1,&Gm0);
        cAdd(Gm0,*result,result);   /* result=Gf+Gm0+Gm1 */
      }

    }
  }
  /* if (pr) printz("%d,%d: %f+%fi\n",mu,nu,result->r,result->i); */
}

/*============================================================*/

void init_Dmatrix(int X,int Y,int Z,
		  double lambda,
		  double dips_lambda,
		  double wavenum)
     /* initial./ises the matrix D. D[i][j][k]=A[i1-i2][j1-j2][k1-k2]
	The routine is called only once, so needs not to be fast.
	*/
{
  int x,y,z,i,j,k,kcor,index,mem,Dcomp,indexfrom,indexto;
  doublecomplex temp_minus[6],temp_plus[6];
  double invNgrid;
  double work[4000];
  doublecomplex *slice;
  
  extern FILE *logfile;
  extern int memory;
  extern double dpl;
  extern int reduced_FFT;

  int nnn; /* multiplier used for reduced_FFT or not reduced */
  int jstart, kstart;

  /* these variables from function couple_matrix */
  int start, mu, nu;

  start=clock(); 

  /* calculate size of 3d grid */ 
  gridX=fft_fit(2*X,nprocs);
  gridY=fft_fit(2*Y,1);
  gridZ=fft_fit(2*Z,2*nprocs);

  /* initialise small X,Y,Z as its big counterparts.
   * After the FFT of D set them to their correct values
   */
  
  smallX=gridX; smallZ=gridZ/2; par_setup();
  smallX=gridX; smallY=gridY; smallZ=gridZ; 
  gridXY=gridX*gridY; gridYZ=gridY*gridZ; gridXZ=gridX*gridZ;
/*#ifdef PARALLEL*/
  DsizeX=gridX;
/*#else
  DsizeX=gridX/2+1;
#endif*/
  fill_use();

  /* temporary - test FFT */
/*  test_fft();
  stop(1);  */
  
  D2sizeX=gridX; 
  
  if (reduced_FFT) {
    D2sizeY=gridY/2; 
    D2sizeZ=gridZ/2;
    DsizeY=gridY/2+1;
    DsizeZ=gridZ/2+1;
    nnn=1;
    jstart=0;
    kstart=0;
  }
  else {
    D2sizeY=DsizeY=gridY;
    D2sizeZ=DsizeZ=gridZ;
    nnn=2;
    jstart=1-Y;
    kstart=1-Z;
  }

  DsizeXY=DsizeX*DsizeY; DsizeYZ=DsizeY*DsizeZ; DsizeXZ=DsizeX*DsizeZ;
  
  Ngrid=gridX*gridY*gridZ; 
  invNgrid=1.0/Ngrid;
  local_Nsmall=smallX*smallY*smallZ/(4*nprocs); /* size of X vector */
  NDcomp=6;
  
  fprintz(logfile,"The 3d grid is:%ix%ix%i\n",gridX,gridY,gridZ);
  mem=sizeof(doublecomplex)*(3*local_Nsmall+NDcomp*local_Nx*DsizeY*DsizeZ);
  memory+=mem;
  printz("Memory usage for X and D matrix:%.1f Mb\n",mem/1048576.0);
  fprintz(logfile,"Memory usage for X and D matrix:%.1f Mb\n",mem/1048576.0);
  if (ringid==0) fflush(logfile);
  
  if (prognose==true) return;
  /* allocate memory for Dmatrix */
  if ((Dmatrix = dCvector(0, Dmatrix_nel=NDcomp*local_Nx*DsizeY*DsizeZ)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc Dmatrix");
  
  /* allocate memory for D2matrix components */
  if ((D2matrix = dCvector(0, (nnn*local_Nz_f)*D2sizeY*D2sizeX)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc D2matrix");
  if ((slice = dCvector(0, gridY*gridZ)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc slice");
  
  smallX=D2sizeX; smallY=D2sizeY; smallZ=D2sizeZ;

  for(Dcomp=0;Dcomp<NDcomp;Dcomp++) {

    switch((char) Dcomp) {    /* determine mu,nu */
      case 0: {
        mu=0;
        nu=0;
	break;
      }
      case 1: {
        mu=0;
        nu=1;
	break;
      }
      case 2: {
        mu=0;
        nu=2;
	break;
      }
      case 3: {
        mu=1;
        nu=1;
	break;
      }
      case 4: {
        mu=1;
        nu=2;
	break;
      }
      case 5: {
        mu=2;
        nu=2;
	break;
      }
    } /* end of switch */
    
    /* fill Dmatrix with 0.0 */ /* actually only needed for few indices */
    for(z=nnn*local_z0;z<nnn*local_z1_f;z++) for(y=0;y<D2sizeY;y++) for(x=0;x<D2sizeX;x++) {
      i=index_Xmatrix(x,y,z,nnn);
      D2matrix[i].r=D2matrix[i].i=0.0;
    }      
    for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
      i=index_sliceD2matrix(y,z);
      slice[i].r=slice[i].i=0.0;
    }
    
    /* fill D (F'i-j) */
    for(k=nnn*local_z0;k<nnn*local_z1_f;k++) {
      if(k>(gridZ/2)) kcor=k-gridZ;
      else kcor=k;
      for(j=jstart;j<Y;j++) for(i=1-X;i<X;i++) {
        index=index_Xmatrix(i,j,k,nnn);                              
        CalcInterTerm(i,j,kcor,mu,nu,&(D2matrix[index]),wavenum);   /* calculate F[mu][nu] */
      }
    } /* end of i,j,k loop */
    synchronize();
    fftX((double *) D2matrix,1,work,0,nnn*local_Nz_f);
    block_transpose((double*) D2matrix,nnn); /* should be changed for different boxes and MPI */

    for(i=local_x0;i<local_x1;i++) {
      x=i;
      for(z=0;z<gridZ;z++) for(y=0;y<gridY;y++) {
	index=index_sliceD2matrix(y,z);        /* actually only needed for few indices */
	slice[index].r=slice[index].i=0.0;
      }  
      
      for(k=kstart;k<Z;k++) for(j=jstart;j<Y;j++) {
	indexfrom=index_garbledX(i,j,k,nnn);
	indexto=index_sliceD2matrix(j,k);
	
	slice[indexto].r=D2matrix[indexfrom].r;
	slice[indexto].i=D2matrix[indexfrom].i;
      }
      
      if (reduced_FFT) {
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
      }
      fftZ(slice,1,0,1,work);
      transposeYZ(slice,0,gridY,gridZ,0,1);
      fftY(slice,1,0,1,work);
      
      for(z=0;z<DsizeZ;z++) for(y=0;y<DsizeY;y++) {
	indexto=index_tDmatrix(x-local_x0,y,z);
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
  free(slice);
}     

/*============================================================*/

void free_FFT_Dmat(void)
{
  free(Dmatrix);
}
