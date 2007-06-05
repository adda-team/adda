/* FILE: make_particlce.c
 * AUTH: Alfons Hoekstra    
 * DESCR: This module calculates the coordinates of dipoles
 *        for a sphere, for the 'fractional packing.
 *
 *        rewriten,
 *        Michel Grimminck 1995
 *        -----------------------------------------------------------
 *        included ellipsoidal particle for work with Victor Babenko
 *        september 2002
 *        --------------------------------------------------------
 *        included many more new particles:
 *        leucocyte, stick, rotatable oblate spheroid, lymphocyte,
 *        rotatable RBC, etc etc etc
 *        December 2003 - February 2004, by Konstantin Semianov
 * 
 *        Currently is developed by Maxim Yurkin
 */

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <string.h>
#include "const.h"
#include "cmplx.h"
#include "types.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"

char geom_format[] = "%d %d %d\n";  /* format of the geom file */

double volume_ratio;
double dradius, radius, dthick;

int local_nvoid_Ndip,nvoid_Ndip;   /* number of local and total non-void dipoles -
                                     used for effetive radius calculation */

extern int boxX,boxY,boxZ;
extern int nDip; 		/* defined in calculator */
extern short int *position;
extern int symX,symY,symZ,symR;
extern int nlocalRows;
extern int memory;
extern double gridspace;
extern double lambda;
extern double sizeX;

extern FILE *logfile;
extern double dpl;
extern char aggregate_file[];
extern int NoSymmetry;
extern int symmetry_enforced;
extern int volcor;

extern int Nmat;
extern int mat_count[MAXNMAT+1];

extern char *material;
extern double coat_ratio,coat_x,coat_y,coat_z;
extern double append_ratio, ratio_x, ratio_y, ratio_z, cratio_x, cratio_y, cratio_z;
extern double xc0,yc0,zc0,xc1,yc1,zc1,xc2,yc2,zc2,xc3,yc3,zc3,xc4,yc4,zc4,
ratio0, ratio1, ratio2, ratio3, ratio4, ratio5, delta, inc_ratio;
extern double diskratio;
extern double ellipsX,ellipsY,ellipsZ;
extern double aspect_r, betaY, betaZ;
extern double aspect_ratio, alphaY, alphaX;

extern clock_t Timing_FileIO, Timing_Particle;
extern char save_geom_fname[];
extern int jagged;
extern int shape;
extern double *DipoleCoord;
extern int save_geom;
extern char directory[];

/* temporary arrays */
char *material_tmp;
double *DipoleCoord_tmp;
short int *position_tmp;

/*============================================================*/

void SaveGeometry(void)
  /* saves dipole configuration to .geom file */
{
  char fname[200], *buf1, *buf2;
  FILE *geom;
  int i,j,buf_size;

#ifdef PARALLEL
  sprintf(fname,"%s/%d.tmp",directory,ringid);
#else
  strcpy(fname,save_geom_fname);
#endif  
  if ((geom=fopen(fname,"w"))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Failed to open file '%s'",fname);
  for(i=0;i<local_nvoid_Ndip;i++) {
    j=3*i;
    fprintf(geom,geom_format,position[j],position[j+1],position[j+2]);
  }
  fclose(geom);
#ifdef PARALLEL
  /* wait for all processes to save their part of geometry */
  synchronize();
  if (ringid==ROOT) {
    /* combine all files into one */
    buf_size=(10+strlen(directory))*nprocs+30;
    if ((buf1=(char *)malloc(buf_size*sizeof(char))) == NULL)
      LogError(EC_ERROR,ONE,POSIT,"Could not malloc buf1");
    if ((buf2=(char *)malloc(buf_size*sizeof(char))) == NULL)
      LogError(EC_ERROR,ONE,POSIT,"Could not malloc buf2");
    strcpy(buf1,"");
    for (i=0;i<nprocs;i++) sprintf(buf1+strlen(buf1)," %s/%d.tmp",directory,i);
    sprintf(buf2,"cat%s >> %s",buf1,save_geom_fname);
    system(buf2);
    free(buf1);
    free(buf2);
  }
  remove(fname);
#endif
}

/*===========================================================*/

void InitDipFile(char *filename,int *maxX, int *maxY, int *maxZ)
   /* read dipole file first to determine box sizes */
{
  FILE *input;
  int x, y, z;

  if ((input=fopen(filename,"r"))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not open dipole file - %s",filename);

  *maxX=*maxY=*maxZ=0;
  while(!feof(input)) {
    if (fscanf(input,geom_format,&x,&y,&z)!=3)
      LogError(EC_ERROR,ONE,POSIT,"Could not scan from dipole file - %s - offset=%d",filename,ftell(input));
    
    if (x>*maxX) *maxX=x;
    if (y>*maxY) *maxY=y;
    if (z>*maxZ) *maxZ=z;
  }
  *maxX=jagged*(*maxX+1);
  *maxY=jagged*(*maxY+1);
  *maxZ=jagged*(*maxZ+1);
  fclose(input);
}

/*===========================================================*/

void ReadDipFile(char *filename)
   /* read dipole file */
{
  FILE *input;
  int x, y, z, x0, y0, z0;
  int index;
  
  if ((input=fopen(filename,"r"))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not open dipole file - %s",filename);
  
  while(!feof(input)) {
    if (fscanf(input,geom_format,&x0,&y0,&z0)!=3)
      LogError(EC_ERROR,ONE,POSIT,"Could not scan from dipole file - %s - offset=%d",filename,ftell(input));
    
    for (z=jagged*z0;z<jagged*(z0+1);z++) if (z>=local_z0 && z<local_z1) 
      for (x=jagged*x0;x<jagged*(x0+1);x++) for (y=jagged*y0;y<jagged*(y0+1);y++) {
        index=(z-local_z0)*boxX*boxY+y*boxX+x;
        material_tmp[index]=0;
    }
  }
  fclose(input);
}

/*==========================================================*/

int fitbox(int box)
   /* finds the smallest value for which program would work (should be even and divide jagged) */
{
  if (jagged%2==0) return (jagged*((box+jagged-1)/jagged));
  else return (2*jagged*((box+2*jagged-1)/(2*jagged)));;
}

/*==========================================================*/

void init_shape(void)
   /* perform of initialization of symmetries and boxY, boxZ
    * Estimate the volume of the particle, when not discretisized.
    * Check whether enough refractive indices is specified
    */
{
  int n_boxX, n_boxYi, n_boxZi, temp;             /* new values for dimensions */
  double n_boxY, n_boxZ;
  clock_t tstart;
  int Nmat_need;

  tstart=clock();

  volume_ratio=UNDEF; /* for some shapes it is initialized below; if not, volume correction is not used */
  if (boxX==UNDEF) {
    if (shape!=READ) boxX=16*jagged; /* default value for boxX */
  }
  else {
    temp=boxX;
    if ((boxX=fitbox(boxX))!=temp) LogError(EC_WARN,ONE,POSIT,"boxX has been adjusted from %i to %i",temp,boxX);
  }
  n_boxX=boxX;

  if (shape==SPHERE) {
    symX=symY=symZ=symR=true;
    volume_ratio=PI/6;
    n_boxY=n_boxZ=boxX;
    Nmat_need=1;
  }
  else if (shape==ELLIPSOIDAL) {
    symX=symY=symZ=true;
    if (1==ellipsY) symR=true; else symR=false;
    volume_ratio=PI/6*ellipsY*ellipsZ;
    n_boxY=ellipsY*boxX;
    n_boxZ=ellipsZ*boxX;
    Nmat_need=1;
  }
  else if(shape==DISK) {
    symX=symY=symZ=symR=true;
    volume_ratio=PI/4*diskratio;
    n_boxY=boxX;
    n_boxZ=diskratio*boxX;
    Nmat_need=1;
  }
/*  else if(shape==RBC) {
    symX=symY=symZ=false;
    symR=false;
    volume_ratio=UNDEF;
  }
  else if(shape==RBC_ROT) {
    symX=symY=symZ=false;
    symR=false;
    volume_ratio=boxX*boxX*boxX;
  }
  else if(shape==SDISK_ROT) {
    symX=symY=symZ=false;
    symR=false;
    volume_ratio=boxX*boxX*boxX;
  }
  else if(shape==SPHEROID_ROT) {
    symX=symY=symZ=false;
    symR=false;
    volume_ratio=boxX*boxX*boxX;
  }
  else if(shape==CYLINDER) {
    symX=symY=symZ=false;
    symR=false;
    volume_ratio=boxX*boxX*boxX;
  }
  else if(shape==STICK) {
    symX=symY=symZ=true;
    symR=false;
    if ((alphaX==0) && (alphaY==0)) {symR=true;}
    volume_ratio=boxX*boxX*boxX;
  }
  else if (shape==LYMPHOCYTE1) {
    if (Nmat<2) printz("WARNING: only one refraction index is given\n");
    symX=symY=symZ=false;
    symR=false;
    volume_ratio=boxX*boxX*boxX;
  }
  else if (shape==LYMPHOCYTE2) {
    if (Nmat<2) printz("WARNING: only one refraction index is given\n");
    symX=symY=symZ=true;
    symR=true;
    if (coat_x!=0 || cratio_x!=1 || ratio_x!=1) {symX=false; symR=false;}
    if (coat_y!=0 || cratio_y!=1 || ratio_y!=1) {symY=false; symR=false;}
    if (coat_z!=0 || cratio_z!=1 || ratio_z!=1) symZ=false;
    volume_ratio=boxX*boxX*boxX;
  }
  else if (shape==LEUCOCYTE2) {
    if (Nmat<3) printz("WARNING: only one refraction index is given\n");
    symX=symY=symZ=false;
    symR=false;
    volume_ratio=boxX*boxX*boxX;
  }  */

  else if (shape==COATED) {
    symX=symY=symZ=symR=true;
    volume_ratio=PI/6;
    if (coat_x!=0) {symX=false; symR=false;}
    if (coat_y!=0) {symY=false; symR=false;}
    if (coat_z!=0) symZ=false;
    n_boxY=n_boxZ=boxX;
    Nmat_need=2;
  }
  else if (shape==SPHEREBOX) {
    symX=symY=symZ=true;
    if (boxX==boxY) symR=true; else symR=false;
    if (coat_x!=0) {symX=false; symR=false;}
    if (coat_y!=0) {symY=false; symR=false;}
    if (coat_z!=0) symZ=false;
    n_boxY=n_boxZ=boxX;
    Nmat_need=2;
  }
  else if (shape==BOX) {
    symX=symY=symZ=true;
    if (boxY==UNDEF || boxX==boxY) symR=true;
    else symR=false;
    n_boxY=n_boxZ=boxX;
    Nmat_need=1;
  }
/*  else if (shape==PRISMA) {
    symX=true;
    symY=symZ=false;
    symR=false;
    volume_ratio=.5*boxX*boxY*boxZ;
  }  */
  else if (shape==LINE) {
    symX=true;
    symY=symZ=true;
    symR=false;
    n_boxY=n_boxZ=2;
    Nmat_need=1;
  }
  else if (shape==READ) {
    symX=symY=symZ=symR=false;
    InitDipFile(aggregate_file,&n_boxX,&n_boxYi,&n_boxZi);
    n_boxY=n_boxYi;
    n_boxZ=n_boxZi;
    Nmat_need=1;
  }
  /* check if enough refr. indices */
  if (Nmat<Nmat_need)
    LogError(EC_ERROR,ONE,POSIT,"Only %d refraction indices is given. %d is required",Nmat,Nmat_need);

  if (symmetry_enforced==true) symX=symY=symZ=symR=true;
  else if (NoSymmetry==true) symX=symY=symZ=symR=false;

  if (boxX==UNDEF) boxX=fitbox(n_boxX);
  else if (n_boxX>boxX)
    LogError(EC_ERROR,ONE,POSIT,"Particle (boxX=%d) does not fit into specified boxX=%d", n_boxX, boxX);

  n_boxY=ceil(n_boxY);
  n_boxZ=ceil(n_boxZ);
  if (boxY==UNDEF) {    /* assumed that boxY and boxZ are either both defined or both not defined */
    boxY=fitbox(n_boxY);
    boxZ=fitbox(n_boxZ);
  }
  else {
    temp=boxY;
    if ((boxY=fitbox(boxY))!=temp) LogError(EC_WARN,ONE,POSIT,"boxY has been adjusted from %i to %i",temp,boxY);
    temp=boxZ;
    if ((boxZ=fitbox(boxZ))!=temp) LogError(EC_WARN,ONE,POSIT,"boxZ has been adjusted from %i to %i",temp,boxZ);

    if (n_boxY>boxY || n_boxZ>boxZ)
      LogError(EC_ERROR,ONE,POSIT,"Particle (boxY,Z={%d,%d}) does not fit into specified boxY,Z={%d,%d}",
               (int) n_boxY, (int) n_boxZ, boxY, boxZ);
  }
  Timing_Particle = clock() - tstart;
}

/*==========================================================*/

void make_particle(void)
  /* creates a particle; initializes all dipoles counts, dpl, gridspace */
{
  int i, j, k, dipcount,index,InPart;
  double x, y, z;
  double radiussq;
  double xr,yr,zr,xcoat,ycoat,zcoat;
  double r_eps, delta1, thet_, rr, phi_, deltaphi;
  int n, ni, N, ii,ij,ik;
  double xr_,yr_,zr_,CX,CY,SX,SY,CZ,SZ, bradius, cradius, clength, cthick, x_,y_,z_;
  double rrc0, rrc1, rrc2, rrc3, rrc4, rrc5, radius0, radius1, radius2, radius3, radius4, radius5,
         x_incl, y_incl, z_incl;
  double centreX,centreY,centreZ;
  int xj,yj,zj;
  int mat;
  clock_t tstart;

  tstart=clock();

  dipcount=index=0;
  radiussq=pow((boxX/2)*gridspace,2);
  /* assumed that box's are even */
  centreX=centreY=centreZ=jagged/2.0;
  /* allocate temporary memory; even if prognose, since they are needed for exact estimation
     they will be reallocated afterwards (when nlocalRows is known)*/
  if ((material_tmp = (char *) malloc(local_Ndip*sizeof(char))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc material_tmp");
  if ((DipoleCoord_tmp = dvector(0,3*local_Ndip-1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc DipoleCoord_tmp");
  if ((position_tmp = (short int *) malloc(3*sizeof(short int)*local_Ndip)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc position_tmp");

  /* way to optimize !!! */
  for(k=local_z0-boxZ/2;k<local_z1-boxZ/2;k++)
    for(j=-boxY/2;j<boxY/2;j++)
      for(i=-boxX/2;i<boxX/2;i++) {
        x=(i+centreX);
        y=(j+centreY);
        z=(k+centreZ);

        xj=jagged*((i+boxX/2)/jagged)-boxX/2;
        yj=jagged*((j+boxY/2)/jagged)-boxY/2;
        zj=jagged*((k+boxZ/2)/jagged)-boxZ/2;

        xr=(xj+centreX)/(boxX);
        yr=(yj+centreY)/(boxX);
        zr=(zj+centreZ)/(boxX);

        mat=Nmat;  /* corresponds to void */

        if (shape==BOX) mat=0;
        else if (shape==DISK) {
          dradius = boxX/2.0;
          dthick = diskratio*dradius;
          if(abs(x) <= dthick/2) if(y*y + z*z <= dradius*dradius) mat = 0;
        }
     /* else if (shape==RBC) { /* here we assume units of micrometers */
      /*  dradius = (boxX*Wgridspace)/2.0;  /* assuming X, Y, Z is equal */
      /*  radius = y*y + z*z;
          if(radius <= dradius*dradius) {
            radius /= (dradius*dradius);
            dthick = (dradius/3.91)*sqrt(1.0-radius) *
            (0.81 + 7.83*radius - 4.39*radius*radius);
            if(x<dthick/2.0 && x>-dthick/2.0) mat = 0;
          }
        } */
     /* else if (shape==RBC_ROT) { /* rottatable erithrocyte around two axes Y and Z */
      /*  dradius = (boxX*gridspace)/2.0;  /* assuming X, Y, Z is equal */
      /*  CY=cos(PI/180*betaY);  SY=sin(PI/180*betaY);
          CZ=cos(PI/180*betaZ);  SZ=sin(PI/180*betaZ);
          x_=x*CY*CZ-y*SZ-z*SY*CZ;
          y_=x*CY*SZ+y*CZ-z*SY*SZ;
          z_=x*SY+z*CY;
          radius=y_*y_+z_*z_;
          if (radius <dradius*dradius) {
            radius /= (dradius*dradius);
            dthick=(dradius*aspect_r)* sqrt(1.0-radius) *
            (0.81 + 7.83*radius - 4.39*radius*radius); /* for previous aspect_r=1/3.91 */
      /*    if (x_<dthick/2.0 && x_>-dthick/2.0) mat = 0;
          }
        }
        else if  (shape==SDISK_ROT) {
          xr= (i+centreX)/(boxX);
          yr= (j+centreY)/(boxY);
          zr= (k+centreZ)/(boxZ);
          CY=cos(PI/180*betaY);  SY=sin(PI/180*betaY);
          CZ=cos(PI/180*betaZ);  SZ=sin(PI/180*betaZ);
          xr_=xr*CY*CZ-yr*SZ-zr*SY*CZ;
          yr_=xr*CY*SZ+yr*CZ-zr*SY*SZ;
          zr_=xr*SY+zr*CY;
          cthick=(aspect_r)/2.0;
          radius=(xr_*xr_+yr_*yr_+zr_*zr_)/0.25;
          if (radius <1.0000001 && xr_<cthick && xr_>-cthick) mat=0;
        }
        else if  (shape==SPHEROID_ROT) {
          xr= (i+centreX)/(boxX);
          yr= (j+centreY)/(boxY);
          zr= (k+centreZ)/(boxZ);
          CY=cos(PI/180*betaY);  SY=sin(PI/180*betaY);
          CZ=cos(PI/180*betaZ);  SZ=sin(PI/180*betaZ);
          xr_=xr*CY*CZ-yr*SZ-zr*SY*CZ;
          yr_=xr*CY*SZ+yr*CZ-zr*SY*SZ;
          zr_=xr*SY+zr*CY;
          xr_/=aspect_r;
          radius=(xr_*xr_+yr_*yr_+zr_*zr_)/0.25;
          if (radius <1.0000001) mat=0;
        }
        else if (shape==STICK) {
          xr= (i+centreX)/(boxX);
          yr= (j+centreY)/(boxY);
          zr= (k+centreZ)/(boxZ);
          CY=cos(PI/180*alphaY);  SY=sin(PI/180*alphaY);
          CX=cos(PI/180*alphaX);  SX=sin(PI/180*alphaX);
          xr_=xr*CY - zr*SY;
          yr_=-xr*SY*SX+yr*CX-zr*CY*SX;
          zr_=xr*SY*CX+yr*SX+zr*CY*CX;
          bradius=append_ratio*append_ratio*aspect_ratio*aspect_ratio/3.99999999999;
          clength=append_ratio*(1-aspect_ratio)/2.0;
          radius=(yr_*yr_+xr_*xr_);
          if ((radius <bradius && zr_<clength && zr_>-clength) ||
   	        (radius+(zr_-clength)*(zr_-clength))< bradius ||
   	        (radius+(zr_+clength)*(zr_+clength))<bradius) mat = 0;
        }  */
    /*  else if (shape==CYLINDER) {
           xr= (i+centreX)/(boxX);
           yr= (j+centreY)/(boxY);
           zr= (k+centreZ)/(boxZ);
           bradius=aspect_ratio*aspect_ratio*append_ratio*append_ratio/3.99999999999;
           clength=append_ratio/2.0;
           radius=(yr*yr+xr*xr);
           if (radius <bradius && zr<clength && zr>-clength) mat=0;
        }
        else if (shape==LYMPHOCYTE1) {
          xr= (i+centreX)/(boxX);
          yr= (j+centreY)/(boxY);
          zr= (k+centreZ)/(boxZ);
          xcoat=xr-coat_x/2;
          ycoat=yr-coat_y/2;
          zcoat=zr-coat_z/2;
          r_eps=0.25*0.25*(1-append_ratio)*(1-append_ratio);
          if (xr*xr+yr*yr+zr*zr < .250001*append_ratio*append_ratio) mat=0;
          if (xcoat*xcoat+ycoat*ycoat+zcoat*zcoat <0.250001*coat_ratio*coat_ratio*append_ratio*append_ratio) mat=1;
          if ((xr*xr+yr*yr+zr*zr >0.25000*append_ratio*append_ratio) && (xr*xr+yr*yr+zr*zr < 0.250001)) {
            delta1=PI/12.0;
            for (ni=0; ni<=12; ni++) {
              thet_=ni*delta1; rr=(append_ratio+1)/4.0*sin(thet_); z_=(append_ratio+1)/4.0*cos(thet_);
    	        deltaphi=2.0*PI/(23.0*(1-abs(1.0-2.0*thet_/PI))+1);
    	        phi_=0;
    	        while (phi_<2.0*PI) {
    	          x_=rr*cos(phi_);  y_=rr*sin(phi_);
    	          bradius=(xr-x_)*(xr-x_)+(yr-y_)*(yr-y_)+(zr-z_)*(zr-z_);
    	          if (bradius<r_eps) mat=0;
    	          phi_=phi_+deltaphi;
    	        }
            }
          }
        }
        else if (shape==LYMPHOCYTE2) {
          xr= (i+centreX)/(boxX);
          yr= (j+centreY)/(boxY);
          zr= (k+centreZ)/(boxZ);
          xcoat=xr-coat_x/2;
          ycoat=yr-coat_y/2;
          zcoat=zr-coat_z/2;
          radius=(xr*xr/ratio_x/ratio_x+yr*yr/ratio_y/ratio_y+zr*zr/ratio_z/ratio_z)/append_ratio/append_ratio;
          cradius=(xcoat*xcoat/cratio_x/cratio_x+ycoat*ycoat/cratio_y/cratio_y+zcoat*zcoat/cratio_z/cratio_z)/append_ratio/append_ratio;
          if (radius < 0.250001) mat=0;
          if (cradius<0.250001*coat_ratio*coat_ratio) mat=1;
        }  */
   /*   else if (shape==LEUCOCYTE2) {
          xr= (i+centreX)/(boxX);
          yr= (j+centreY)/(boxY);
          zr= (k+centreZ)/(boxZ);
          N=(1/delta);
          rr=(1-inc_ratio/2.0)*(1-inc_ratio/2.0)*0.25;
          r_eps=0.25*inc_ratio*inc_ratio;
          radius=xr*xr+yr*yr+zr*zr;
          if (xc0>=1) xc0=-(xc0-1.0);
          if (xc1>=1) xc1=-(xc1-1.0);
          if (xc2>=1) xc2=-(xc2-1.0);
          if (xc3>=1) xc3=-(xc3-1.0);
          if (xc4>=1) xc4=-(xc4-1.0);
          if (yc0>=1) yc0=-(yc0-1.0);
          if (yc1>=1) yc1=-(yc1-1.0);
          if (yc2>=1) yc2=-(yc2-1.0);
          if (yc3>=1) yc3=-(yc3-1.0);
          if (yc4>=1) yc4=-(yc4-1.0);
          if (zc0>=1) zc0=-(zc0-1.0);
          if (zc1>=1) zc1=-(zc1-1.0);
          if (zc2>=1) zc2=-(xc2-1.0);
          if (zc3>=1) zc3=-(zc3-1.0);
          if (zc4>=1) zc4=-(zc4-1.0);
          rrc0=0.25*(ratio0+inc_ratio/2.0)*(ratio0+inc_ratio/2.0);
          radius0= (xc0-xr)*(xc0-xr)+(yc0-yr)*(yc0-yr)+(zc0-zr)*(zc0-zr);
          rrc1=0.25*(ratio1+inc_ratio/2.0)*(ratio1+inc_ratio/2.0);
          radius1= (xc1-xr)*(xc1-xr)+(yc1-yr)*(yc1-yr)+(zc1-zr)*(zc1-zr);
          rrc2=0.25*(ratio2+inc_ratio/2.0)*(ratio2+inc_ratio/2.0);
          radius2= (xc2-xr)*(xc2-xr)+(yc2-yr)*(yc2-yr)+(zc2-zr)*(zc2-zr);
          rrc3=0.25*(ratio3+inc_ratio/2.0)*(ratio3+inc_ratio/2.0);
          radius3= (xc3-xr)*(xc3-xr)+(yc3-yr)*(yc3-yr)+(zc3-zr)*(zc3-zr);
          rrc4=0.25*(ratio4+inc_ratio/2.0)*(ratio4+inc_ratio/2.0);
          radius4= (xc4-xr)*(xc4-xr)+(yc4-yr)*(yc4-yr)+(zc4-zr)*(zc4-zr);
          if (radius< 0.250001) {
            mat=0;
            if ( radius0<0.250001*ratio0*ratio0)  mat=1;
            else if (radius1<0.250001*ratio1*ratio1)  mat=1;
            else if (radius2<0.250001*ratio2*ratio2)  mat=1;
            else if (radius3<0.250001*ratio3*ratio3)  mat=1;
            else if (radius4<0.250001*ratio4*ratio4)  mat=1;
            else if (radius <rr) {
              for (ik=1; ik<N; ik++) {
                z_incl=(zr+0.5-delta*ik)*(zr+0.5-delta*ik);
                for (ij=1; ij<N; ij++) {
                  y_incl=(yr+0.5-delta*ij)*(yr+0.5-delta*ij);
                  for (ii=1; ii<N; ii++) {
         	       x_incl=(xr+0.5-delta*ii)*(xr+0.5-delta*ii);
         	       if (x_incl+y_incl+z_incl<r_eps) mat=1;
         	     }
                }
              }
            }
          }
        }
        else if (shape==PRISMA && y*boxZ>=-z*boxY) mat=0;*/

        else if (shape==SPHERE) {
          if (xr*xr+yr*yr+zr*zr < .250001) mat=0;
        }
        else if (shape==ELLIPSOIDAL) {
          if (x*x/(ellipsX*ellipsX) +
    	  y*y/(ellipsY*ellipsY) +
    	  z*z/(ellipsZ*ellipsZ) < 0.25 * boxX*boxX) mat = 0;
        }
        else if (shape==COATED) {
          xcoat=xr-coat_x/2;
          ycoat=yr-coat_y/2;
          zcoat=zr-coat_z/2;
          if (xr*xr+yr*yr+zr*zr < .250001) mat=0;
          if (xcoat*xcoat+ycoat*ycoat+zcoat*zcoat < .250001*coat_ratio*coat_ratio) mat=1;
        }
        else if (shape==SPHEREBOX) {
          if (xr*xr+yr*yr+zr*zr < .250001*coat_ratio*coat_ratio) mat=1;
          else mat=0;
        }
        else if (shape==LINE) {
          if (j==0 && k==0) mat=0;
        }
        position_tmp[3*index]=(short int)(i+boxX/2);
        position_tmp[3*index+1]=(short int)(j+boxY/2);
        position_tmp[3*index+2]=(short int)(k+boxZ/2);
        /* afterwards multiplied by gridspace */
        DipoleCoord_tmp[3*index] = x;
        DipoleCoord_tmp[3*index+1] = y;
        DipoleCoord_tmp[3*index+2] = z;
        material_tmp[index]=(char)mat;
        index++;
  } /* End box loop */
  if (shape==READ) ReadDipFile(aggregate_file);

  /* initialization of mat_count and dipoles counts */
  for(i=0;i<=Nmat;i++) mat_count[i]=0;
  for(i=0;i<local_Ndip;i++) mat_count[material_tmp[i]]++;
  local_nvoid_Ndip=local_Ndip-mat_count[Nmat];
  my_inner_product(mat_count,int_type,Nmat+1);
  nvoid_Ndip=Ndip-mat_count[Nmat];
  nlocalRows=3*local_nvoid_Ndip;
  /* initialize dpl and gridspace */
  if (dpl==UNDEF) {
    /* dpl is corrected to give correct volume */
    if (volcor && volume_ratio!=UNDEF)
      dpl=lambda*pow(nvoid_Ndip/volume_ratio,.333333333333)/sizeX;
    else dpl=lambda*boxX/sizeX;
  }
  gridspace=lambda/dpl;
  /* allocate main particle arrays, using precise nlocalRows
     even when prognose to enable save_geom afterwards */
  if ((material = (char *) malloc(local_nvoid_Ndip*sizeof(char))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc material");
  if ((DipoleCoord = dvector(0,nlocalRows-1)) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc DipoleCoord");
  if ((position = (short int *) malloc(nlocalRows*sizeof(short int))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc position");
  memory+=(3*(sizeof(short int)+sizeof(double))+sizeof(char))*local_nvoid_Ndip;
  /* copy nontrivial part of arrays */
  index=0;
  for (i=0;i<local_Ndip;i++) if (material_tmp[i]<Nmat) {
    material[index]=material_tmp[i];
    MultScal(gridspace,DipoleCoord_tmp+3*i,DipoleCoord+3*index);  /* DipoleCoord=gridspace*DipoleCoord_tmp */
    memcpy(position+3*index,position_tmp+3*i,3*sizeof(short int));
    index++;
  }
  /* free temporary memory */
  free(material_tmp);
  free_dvector(DipoleCoord_tmp,0);
  free(position_tmp);

  /* save geometry */
  if (save_geom) SaveGeometry();

  Timing_Particle += clock() - tstart;
}

