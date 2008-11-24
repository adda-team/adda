/* FILE: CalculateE.c
 * AUTH: Alfons Hoekstra
 * DATE: -
 *
 * January 2004 : include module to compute full Mueller Matrix over
 * full space angle, not very efficient, must be improved (A. Hoekstra)
 */

/* the module that will calculate the E field, it is very messy with
   global variables etc. That must be straightened out later, (no way)
first we
   concentrate on the harness of the implementation...
   */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#include "cmplx.h"
#include "types.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "crosssec.h"

extern int nDip; 		/* defined in calculator */		
extern int nTheta;		/* defined in calculator */
extern int nRing;		/* defined in calculator */
extern int RingId;		/* defined in calculator */
extern int nlocalDip;		/* defined in calculator */
extern int nlocalRows;		/* defined in calculator */
extern REAL **DipoleCoord;	/* defined in calculator */
extern dcomplex *Avecbuffer;	/* defined in calculator */
extern double WaveNum;		/* defined in calculator */
extern double LAMBDA;           /* defined in calculator */
extern double eps;		/* defined in calculator */

extern int yzplane;		/* defined in main */
extern int xzplane;		/* defined in main */
extern int store_all_dir;
extern int all_dir;		/* defined in main */
extern int calc_Cext;		/* defined in main */ 
extern int calc_Cabs; 		/* defined in main */
extern int calc_Csca;        	/* defined in main */
extern int calc_Csca_diff;      /* defined in main */
extern int calc_vec;            /* defined in main */
extern int calc_asym; 		/* defined in main */
extern int calc_mat_force;      /* defined in main */
extern int store_force;

extern dcomplex *x;		/* defined in main */
extern dcomplex *r;		/* defined in main */
extern dcomplex *p;		/* defined in main */
extern dcomplex *buffer;	/* defined in main */
extern dcomplex *Eper, *Epar;   /* defined in main */
extern dcomplex *Einc;          /* defined in main */
extern char directory[200];     /* defined in main */
extern FILE *logfile;           /* defined in main */
extern int store_field;
extern int store_int_field;
extern double gridspaceX,gridspaceY,gridspaceZ;

extern int *material;
extern int Nmat;
extern int mat_count[10];
extern doublecomplex ref_index[10];

double epsB;          /* the stopping criterium */
double inprodB;       /* innerproduct of b */
double inprodR;       /* innerproduct of rk */

extern dcomplex *dCvector ( );
extern dcomplex **dCmatrix ( );
extern double *dvector ( );
extern void GenerateB ( );
extern void MatVec_nim ( );
extern void MatVecHer_nim ( );
extern void couple_matrix ( );
extern int CGNR ( );
extern int symR;

void save_in_mie(void)
{
  extern double beam_x0,beam_y0,beam_z0,coat_ratio;
  extern double LAMBDA;
  extern doublecomplex ref_index[10];
  extern int mat_count[10];
  extern int shape;
  FILE *out;
  double r,rin;
  extern double dplX,dplY,dplZ;
  char buffer[1024];
  int i,Nmie,inner,outer;
  extern int beamtype;
  extern double beam_w0,beam_x0,beam_y0,beam_z0;

  inner=0; outer=0;
  if (shape==COATED) {
    inner=1;
    outer=0;
  }
  if (beamtype==PLANE) Nmie=1; else Nmie=2;
  r=2.0*PI*pow((0.75/PI)*nDip/dplX/dplY/dplZ,.333333333333);
  
  inner=0; outer=0; rin=r;
  if (shape==COATED) {
    inner=1;
    outer=0;
    rin=coat_ratio*r;
    rin=2.0*PI*pow((0.75/PI)*mat_count[1]/dplX/dplY/dplZ,.333333333333);
  }

  for(i=0;i<Nmie;i++) {
    sprintz(buffer,"%s/in_mie%i",directory,i);
    out = fopen (buffer, "wA");
    if (out==NULL) {
      printz("File 'in_mie' write error\n");
      stop(1);
    }
    if (beamtype==PLANE || (beam_x0==0.0 && beam_y0==0.0 && beam_z0==0.0)) {
      fprintz(out,"%i\n",nTheta+1);
      fprintz(out,"%g\n%g\n",rin,r);
      fprintz(out,"%f\n%f\n%f\n%f\n",
	      ref_index[inner].r,ref_index[inner].i,
	      ref_index[outer].r,ref_index[outer].i);
      fprintz(out,"1.0\n%.5E\n",LAMBDA);
      if (beamtype==PLANE) {
	fprintz(out,"%.5E\n",1e6*LAMBDA);
	fprintz(out,"0\n");
      }
      else {
	fprintz(out,"%.5E\n",beam_w0);
	fprintz(out,"99\n");
      }
    } 
    else { /* focussed beam */
      fprintz(out,"%.7E\n",LAMBDA);
      fprintz(out,"%.7E\n",beam_w0);
      fprintz(out,"%.7E\n",r/PI*LAMBDA);
      fprintz(out,"%.7E\n",ref_index[0].r);
      fprintz(out,"%.7E\n1.0\n",ref_index[0].i);
      if (i==1) fprintz(out,"%.7E\n%.7E\n%.7E\n",beam_x0,beam_y0,beam_z0);
      else fprintz(out,"%.7E\n%.7E\n%.7E\n",beam_y0,-beam_x0,beam_z0);
      fprintz(out,"%.1f\n0.0\n",90.0*i);
      fprintz(out,"%.15E\n",180.0/nTheta);
      fprintz(out,"%.7E\n",(double) nTheta);
      fprintz(out,"12.0\n0\n");
    }
    fclose(out);
  }
}

double cmulth(doublecomplex a,doublecomplex b,int part)
     /* complex multiplycation; returns real(a*b_tranposed) or imag(a*b_tranposed) */
{
  doublecomplex temp;
  
  temp.r=a.r*b.r + a.i*b.i;
  temp.i=a.i*b.r - a.r*b.i;
  
  if (part==real) return(temp.r);
  else return(temp.i);
}

void compute_mueller_matrix(double matrix[4][4], doublecomplex s1, doublecomplex s2, doublecomplex s3, doublecomplex s4)
/* computer mueller matrix from scattering matrix elements s1, s2, s3, s4, accoording to formula 3.16 */
/* from Bohren and Huffman, with our own corrections for S41, S42, S43 */
{
	matrix[0][0] = 0.5*(cmulth(s1,s1,real)+cmulth(s2,s2,real)+
	                    cmulth(s3,s3,real)+cmulth(s4,s4,real));
	matrix[0][1] = 0.5*(cmulth(s2,s2,real)-cmulth(s1,s1,real)+
	                    cmulth(s4,s4,real)-cmulth(s3,s3,real));
	matrix[0][2] = cmulth(s2,s3,real)+cmulth(s1,s4,real);
	matrix[0][3] = cmulth(s2,s3,imag)-cmulth(s1,s4,imag);
	
	matrix[1][0] = 0.5*(cmulth(s2,s2,real)-cmulth(s1,s1,real)+
	                    cmulth(s3,s3,real)-cmulth(s4,s4,real));
	matrix[1][1] = 0.5*(cmulth(s2,s2,real)+cmulth(s1,s1,real)+
	                    -cmulth(s3,s3,real)-cmulth(s4,s4,real));
	matrix[1][2] = cmulth(s2,s3,real)-cmulth(s1,s4,real);
	matrix[1][3] = cmulth(s2,s3,imag)+cmulth(s1,s4,imag);

	matrix[2][0] = cmulth(s2,s4,real)+cmulth(s1,s3,real);
	matrix[2][1] = cmulth(s2,s4,real)-cmulth(s1,s3,real);
	matrix[2][2] = cmulth(s1,s2,real)+cmulth(s3,s4,real);
	matrix[2][3] = cmulth(s2,s1,imag)+cmulth(s4,s3,imag);

	matrix[3][0] = cmulth(s4,s2,imag)+cmulth(s1,s3,imag);
	matrix[3][1] = cmulth(s4,s2,imag)-cmulth(s1,s3,imag);
	matrix[3][2] = cmulth(s1,s2,imag)-cmulth(s3,s4,imag);
	matrix[3][3] = cmulth(s1,s2,real)-cmulth(s3,s4,real);
}
void compute_mueller_matrix_norm(double matrix[4][4], doublecomplex s1, doublecomplex s2, doublecomplex s3, doublecomplex s4)
/* computer mueller matrix from scattering matrix elements s1, s2, s3, s4, accoording to formula 3.16 */
/* from Bohren and Huffman, with our own corrections for S41, S42, S43 */
{
	matrix[0][0] = 0.5*(cmulth(s1,s1,real)+cmulth(s2,s2,real)+
	                    cmulth(s3,s3,real)+cmulth(s4,s4,real));
	matrix[0][1] = 0.5*(cmulth(s2,s2,real)-cmulth(s1,s1,real)+
	                    cmulth(s4,s4,real)-cmulth(s3,s3,real))/matrix[0][0];
	matrix[0][2] = (cmulth(s2,s3,real)+cmulth(s1,s4,real))/matrix[0][0];
	matrix[0][3] = (cmulth(s2,s3,imag)-cmulth(s1,s4,imag))/matrix[0][0];
	
	matrix[1][0] = 0.5*(cmulth(s2,s2,real)-cmulth(s1,s1,real)+
	                    cmulth(s3,s3,real)-cmulth(s4,s4,real))/matrix[0][0];
	matrix[1][1] = 0.5*(cmulth(s2,s2,real)+cmulth(s1,s1,real)+
	                    -cmulth(s3,s3,real)-cmulth(s4,s4,real))/matrix[0][0];
	matrix[1][2] = (cmulth(s2,s3,real)-cmulth(s1,s4,real))/matrix[0][0];
	matrix[1][3] = (cmulth(s2,s3,imag)+cmulth(s1,s4,imag))/matrix[0][0];

	matrix[2][0] = (cmulth(s2,s4,real)+cmulth(s1,s3,real))/matrix[0][0];
	matrix[2][1] = (cmulth(s2,s4,real)-cmulth(s1,s3,real))/matrix[0][0];
	matrix[2][2] = (cmulth(s1,s2,real)+cmulth(s3,s4,real))/matrix[0][0];
	matrix[2][3] = (cmulth(s2,s1,imag)+cmulth(s4,s3,imag))/matrix[0][0];

	matrix[3][0] = (cmulth(s2,s4,imag)+cmulth(s1,s3,imag))/matrix[0][0];
	matrix[3][1] = (cmulth(s2,s4,imag)-cmulth(s1,s3,imag))/matrix[0][0];
	matrix[3][2] = (cmulth(s2,s1,imag)-cmulth(s3,s4,imag))/matrix[0][0];
	matrix[3][3] = (cmulth(s1,s2,real)-cmulth(s3,s4,real))/matrix[0][0];
}

void mueller_matrix(void)
{
  FILE *mueller,*parpar,*parper,*perper,*perpar, *E_Y, *E_X;
  double matrix[4][4];
  double theta, phi, tempr, tempi;
  double exr, exi, eyr, eyi, ezr, ezi;
  double sphi, cphi, stheta, ctheta;
  double exper_r, exper_i, expar_r, expar_i;
  double eyper_r, eyper_i, eypar_r, eypar_i;
  doublecomplex s1,s2,s3,s4;
  char stringbuffer[500];
  
  if (yzplane)
    {
      if (me!=0) return;
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-par");
      if(xzplane==true) strcat(stringbuffer,"-yz");
      perpar = fopen (stringbuffer, "r");
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-per");
      if(xzplane==true) strcat(stringbuffer,"-yz");
      perper = fopen (stringbuffer, "r");
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-par");
      if(xzplane==true) strcat(stringbuffer,"-yz");
      parpar = fopen (stringbuffer, "r");
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-per");
      if(xzplane==true) strcat(stringbuffer,"-yz");
      parper = fopen (stringbuffer, "r");
      
      if (parpar==NULL || parper==NULL || perper==NULL || perpar==NULL) {
	printz("File read error\n");
	stop(1);
      }
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/mueller");
      if(xzplane==true) strcat(stringbuffer,"-yz");
      mueller = fopen (stringbuffer, "wA");
      if (mueller==NULL) {
	printz("File 'mueller' write error\n");
	stop(1);
      }
      fprintz(mueller,"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");
      
      while(!feof(parpar)) { /* not debugged */
	fscanf(parpar,"%lf %lf %lf\n",&theta,&s2.r,&s2.i);
	fscanf(parper,"%lf %lf %lf\n",&theta,&s4.r,&s4.i);
	fscanf(perper,"%lf %lf %lf\n",&theta,&s1.r,&s1.i);
	fscanf(perpar,"%lf %lf %lf\n",&theta,&s3.r,&s3.i);
	
	compute_mueller_matrix(matrix,s1,s2,s3,s4);

	fprintz(mueller,
		"%.2f %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E\n",
		theta,
		matrix[0][0],
		matrix[0][1],
		matrix[0][2],
		matrix[0][3],
		matrix[1][0],
		matrix[1][1],
		matrix[1][2],
		matrix[1][3],
		matrix[2][0],
		matrix[2][1],
		matrix[2][2],
		matrix[2][3],
		matrix[3][0],
		matrix[3][1],
		matrix[3][2],
		matrix[3][3]);
      }
      fclose(parpar); fclose(parper);
      fclose(perper); fclose(perpar);
      fclose(mueller);
      
      /* compress data */
      sprintz(stringbuffer,"gzip -9 %s/par-par",directory);
      if(xzplane==true) strcat(stringbuffer,"-yz"); systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/par-per",directory); 
      if(xzplane==true) strcat(stringbuffer,"-yz"); systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/per-par",directory);
      if(xzplane==true) strcat(stringbuffer,"-yz"); systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/per-per",directory);
      if(xzplane==true) strcat(stringbuffer,"-yz"); systemz(stringbuffer);
    }

  if(xzplane==true) 
    {
      /* also calculate mueller in xz plane */
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-par-xz");
      perpar = fopen (stringbuffer, "r");
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-per-xz");
      perper = fopen (stringbuffer, "r");
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-par-xz");
      parpar = fopen (stringbuffer, "r");
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-per-xz");
      parper = fopen (stringbuffer, "r");
      
      if (parpar==NULL || parper==NULL || perper==NULL || perpar==NULL) {
	printz("File read error\n");
	stop(1);
      }
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/mueller-xz");
      mueller = fopen (stringbuffer, "wA");
      if (mueller==NULL) {
	printz("File 'mueller' write error\n");
	stop(1);
      }

      fprintz(mueller,"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");
      
      while(!feof(parpar)) { /* not debugged */
	fscanf(parpar,"%lf %lf %lf\n",&theta,&s2.r,&s2.i);
	fscanf(parper,"%lf %lf %lf\n",&theta,&s4.r,&s4.i);
	fscanf(perper,"%lf %lf %lf\n",&theta,&s1.r,&s1.i);
	fscanf(perpar,"%lf %lf %lf\n",&theta,&s3.r,&s3.i);
	
	compute_mueller_matrix(matrix,s1,s2,s3,s4);

	fprintz(mueller,
		"%.2f %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E\n",
		theta,
		matrix[0][0],
		matrix[0][1],
		matrix[0][2],
		matrix[0][3],
		matrix[1][0],
		matrix[1][1],
		matrix[1][2],
		matrix[1][3],
		matrix[2][0],
		matrix[2][1],
		matrix[2][2],
		matrix[2][3],
		matrix[3][0],
		matrix[3][1],
		matrix[3][2],
		matrix[3][3]);

      }
      fclose(parpar); fclose(parper);
      fclose(perper); fclose(perpar);
      fclose(mueller);
      
      /* compress data */
      sprintz(stringbuffer,"gzip -9 %s/par-par-xz",directory);
      systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/par-per-xz",directory); 
      systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/per-par-xz",directory);
      systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/per-per-xz",directory);
      systemz(stringbuffer);
    }

    if(all_dir && store_all_dir){
      /* compute Mueller Matrix in full space angle. E-fields are stored in file EgridX and EgridY */
      /* for incoming X and Y polarized light. From this compute the Mueller matrix. */
      /* this is done by first computing the scattered fields in the par-per frame, using */
      /* Espar = cos(theta)cos(phi)Ex + cos(theta)sin(phi)Ey - sin(theta)Ez */
      /* Esper = sin(phi)Ex - cos(phi)Ey */
      /* Next, this is converted to the scattering matrix elements (see e.g Bohren and Huffman) : */
      /* s2 = cos(phi)E'X'par + sin(phi)E'Y'par */
      /* s3 = sin(phi)E'X'par - cos(phi)E'Y'par */
      /* s4 = cos(phi)E'X'per + sin(phi)E'Y'per */
      /* s1 = sin(phi)E'X'per - cos(phi)E'Y'per */
      /* from this the mueller matrix elements are computed */
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/EgridX");
      E_X = fopen (stringbuffer, "r");
      strcpy(stringbuffer,directory); strcat(stringbuffer,"/EgridY");
      E_Y = fopen (stringbuffer, "r");

      if (E_X==NULL || E_Y==NULL) {
	printz("File read error E_X or E_Y\n");
	stop(1);
      }

      strcpy(stringbuffer,directory); strcat(stringbuffer,"/mueller-AllDir");
      mueller = fopen (stringbuffer, "wA");
      if (mueller==NULL) {
	printz("File 'mueller' write error\n");
	stop(1);
      }

      while(!feof(E_Y)) { /* not debugged */
	fscanf(E_X,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&theta,&phi,&exr,&exi,&eyr,&eyi,&ezr,&ezi);
	sphi = sin(phi*PI/180.0); cphi = cos(phi*PI/180.0);
	stheta = sin(theta*PI/180.0); ctheta = cos(theta*PI/180.0);
	expar_r = ctheta*cphi*exr + ctheta*sphi*eyr - stheta*ezr;
	expar_i = ctheta*cphi*exi + ctheta*sphi*eyi - stheta*ezi;
	exper_r = sphi*exr - cphi*eyr;
	exper_i = sphi*exi - cphi*eyi;
	fscanf(E_Y,"%lf %lf %lf %lf %lf %lf %lf %lf\n",&theta,&phi,&exr,&exi,&eyr,&eyi,&ezr,&ezi);
	eypar_r = ctheta*cphi*exr + ctheta*sphi*eyr - stheta*ezr;
	eypar_i = ctheta*cphi*exi + ctheta*sphi*eyi - stheta*ezi;
	eyper_r = sphi*exr - cphi*eyr;
	eyper_i = sphi*exi - cphi*eyi;

	s2.r = cphi*expar_r + sphi*eypar_r;
	s2.i = cphi*expar_i + sphi*eypar_i;
	s3.r = sphi*expar_r - cphi*eypar_r;
	s3.i = sphi*expar_i - cphi*eypar_i;
	s4.r = cphi*exper_r + sphi*eyper_r;
	s4.i = cphi*exper_i + sphi*eyper_i;
	s1.r = sphi*exper_r - cphi*eyper_r;
	s1.i = sphi*exper_i - cphi*eyper_i;
	
	compute_mueller_matrix(matrix,s1,s2,s3,s4);

	fprintz(mueller,
		"%.2f %.2f %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E %.7E\n",
		theta,
		phi,
		matrix[0][0],
		matrix[0][1],
		matrix[0][2],
		matrix[0][3],
		matrix[1][0],
		matrix[1][1],
		matrix[1][2],
		matrix[1][3],
		matrix[2][0],
		matrix[2][1],
		matrix[2][2],
		matrix[2][3],
		matrix[3][0],
		matrix[3][1],
		matrix[3][2],
		matrix[3][3]);
      }
      fclose(E_X); fclose(E_Y);
      fclose(mueller);
      
      /* compress data */
      sprintz(stringbuffer,"gzip -9 %s/EgridX",directory);
      systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/EgridY",directory); 
      systemz(stringbuffer);
      sprintz(stringbuffer,"gzip -9 %s/mueller-AllDir",directory);
      systemz(stringbuffer);
    }
}

void CalculateE(char which,int type) /* x or y polarized incident light */
{ 
  int i, j, jjj, k, l;
  double temp;
  int numiters;		/* number of iterations of the CGNR */
  dcomplex *ebuff;	/* small vector to hold E fields */
  double *robserver;	/* small vector for observer in E calculation */
  dcomplex **fint;	/* small matrix for couple matrix in E calc */
  double theta, dtheta;	/* scattering angle and step in it */
  FILE *filepar, *fileper; /* files for calculated fields */
  FILE *Intfldpntr;	   /* file to store internal fields */
  char *Fname = "CalculateE",stringbuffer[500];
  int base, reminder = nDip%nRing;
  
  unsigned int tstart, tstop, tiostart, tiostop;
  
  extern unsigned long Timing_EField, Timing_FileIO;
  extern double beam_w0,beam_x0,beam_y0,beam_z0;
  extern int beamtype;
  extern doublecomplex cc[10];
  
  if (reminder == 0)
    base = RingId * nlocalDip;
  else if (RingId < reminder)
    base = RingId * nlocalDip;
  else
    base = RingId * nlocalDip + reminder;

  /* initialize robserver and fint */
  if ((ebuff = dCvector (0, 2)) == NULL) {
    LogError (EC_ERROR, Fname,
              "ringID %d, could not malloc ebuff",
	      RingId);
    AbortServer (1);
  }
  
  if ((robserver = dvector (0, 2)) == NULL) {
    LogError (EC_ERROR, Fname,
              "ringID %d, could not malloc robserver",
	      RingId);
    AbortServer (1);
  }
  
  if ((fint = dCmatrix (0, 2, 0, 2)) == NULL) {
    LogError (EC_ERROR, Fname,
              "ringID %d, could not malloc fint",
	      RingId);
    AbortServer (1);
  }
  
  /* calculate the vector b, store it temporarily in x, so that
     x0 is the incident field !
     */
  D("Generating B");
  GenerateB (which, x, DipoleCoord, nlocalDip, nDip, RingId, nRing, WaveNum,
	     beamtype,beam_w0,beam_x0,beam_y0,beam_z0);
  D("done Generating B");
  memcpy(Einc,x,nlocalRows*sizeof(dcomplex));
  
  /* calculate |b|  BLAS1 !!!*/
  
  inprodB = 0.0;
  for (i = 0; i < local_Ndip; ++i) { 
    if (material[i]<Nmat-1) {
      inprodB += x[3*i].r * x[3*i].r + x[3*i].i * x[3*i].i;
      inprodB += x[3*i+1].r * x[3*i+1].r + x[3*i+1].i * x[3*i+1].i;
      inprodB += x[3*i+2].r * x[3*i+2].r + x[3*i+2].i * x[3*i+2].i;
    }
  }
  
  my_inner_product(&inprodB);	/* ACCUMULATE THE INNNERPRODUCT OF B */
  epsB = eps * inprodB;
  
  /* calculate Ax0, r0 (= b - Ax0) and |r| */
  MatVec_nim (x, Avecbuffer,
	      nlocalDip, nDip, RingId, nRing,
	      WaveNum);
  
  inprodR = 0.0;
  
  for (i=0; i < local_Ndip; ++i) {
    if (material[i]<Nmat-1) {
      r[3*i].r = x[3*i].r - Avecbuffer[3*i].r;
      r[3*i].i = x[3*i].i - Avecbuffer[3*i].i;
      inprodR += (r[3*i].r * r[3*i].r + r[3*i].i * r[3*i].i);
      r[3*i+1].r = x[3*i+1].r - Avecbuffer[3*i+1].r;
      r[3*i+1].i = x[3*i+1].i - Avecbuffer[3*i+1].i;
      inprodR += (r[3*i+1].r * r[3*i+1].r + r[3*i+1].i * r[3*i+1].i);
      r[3*i+2].r = x[3*i+2].r - Avecbuffer[3*i+2].r;
      r[3*i+2].i = x[3*i+2].i - Avecbuffer[3*i+2].i;
      inprodR += (r[3*i+2].r * r[3*i+2].r + r[3*i+2].i * r[3*i+2].i);
    }
  }
  
  
  D("inner_product");
  my_inner_product(&inprodR);	/* ACCUMULATE THE INNERPRODUCT OF R */
  D("inner_product done");
  
  /* calculate p0 = Ahr0 */
  D("MatVecHer_nim");
  MatVecHer_nim (r, p,
	         nlocalDip, nDip, RingId, nRing,
	         WaveNum);
  D("MatVecHer_nim done");
  
  if (RingId == 0) {		/* PRINT THE START VALUES */
    fprintz(logfile,"epsilon*b = %1.10e\n",epsB);
    fprintz(logfile,"r0        = %1.10e\n", inprodR);
    printz("r0        = %1.10e\n", inprodR);
    fflush(stdout);
  }
  
  /* calculate solution vector x */
  if ((numiters = CGNR (nDip, 20)) < 0) {
    /* NO CONVERGENCE, DO SOMETHING  */
    LogError (EC_ERROR, Fname,
              "ringID %d, NO CONVERGENCE !\n", RingId);
    AbortServer (1);
  }
  else {
    /* CALCULATE THE ELECTIC FIELD FROM X */
    /* 1: for all angles, summate the condtributions of all dipoles in
     *    this processor;
     * 2: accumulate all local fields and summate
     */
    int nmax;
    int orient,Norient;
    extern int symY,symX;
    char whichfile;
    
    if (yzplane)
      {
	if (type==NORMAL) Norient=1; else Norient=2; /* # orientations */
	tstart = extime ();                  /* Norient is not working !!! */
	
	dtheta = PI / ((double)nTheta);
	
	for(orient=0;orient<Norient;orient++) {
	  /* in case of Rotation symmetry */
	  if (orient==0) whichfile=which;
	  if (orient==1 && which=='X') whichfile='Y'; /* swap X and Y */
	  if (orient==1 && which=='Y') whichfile='X';
	  
	  if ((symY==true && whichfile=='Y') || (symX==true && whichfile=='X'))
	    nmax=nTheta; else nmax=2*nTheta-1;
	  
	  for (i = 0; i <= nmax; ++i) {
	    if ((i%(20*nprocs))==0) test_interrupt();
	    theta = i * dtheta;
	    
	    robserver[0] = 0.0;
	    robserver[1] = sin (theta);
	    robserver[2] = cos (theta);
	    if (orient==0) { /* normal procedure */
	      robserver[0] = 0.0;
	      robserver[1] = sin (theta);
	      robserver[2] = cos (theta);
	    }
	    else { /* Rotation symmetry: calculate per-per from current data */
	    /* CalculateE is called from calculator with 'Y' polarization */
	    /* we now just assume that we have the x-z plane as the scatering plane, */
	    /* rotating in the negative x-direction. This mimics the real case of */
	    /* X polarization with the y-z plane as scattering plane */
	      robserver[0] = -sin (theta);
	      robserver[1] = 0.0;
	      robserver[2] = cos (theta);
	    }
	    
	    ebuff[0].r = ebuff[0].i = 0.0;
	    ebuff[1].r = ebuff[1].i = 0.0;
	    ebuff[2].r = ebuff[2].i = 0.0;
	    calc_field(x,ebuff,DipoleCoord,robserver,WaveNum,fint,nlocalDip);
	    
	    /* convert to (l,r) frame */
	    if (orient==0) {
	      Eper[i].r = ebuff[0].r;
	      Eper[i].i = ebuff[0].i;
	      Epar[i].r = cos (theta) * ebuff[1].r - sin (theta) * ebuff[2].r;
	      Epar[i].i = cos (theta) * ebuff[1].i - sin (theta) * ebuff[2].i;
	    }
	    else {
	      Eper[i].r=ebuff[1].r;
	      Eper[i].i=ebuff[1].i;
	      /* these two lines are not checked for bugs */
	      Epar[i].r = cos (theta) * ebuff[0].r - sin (theta) * ebuff[2].r;
	      Epar[i].i = cos (theta) * ebuff[0].i - sin (theta) * ebuff[2].i;
	    }
	  } /*  end for i */
	  
	  /* ACCUMULATE EPAR AND EPER TO ROOT AND SUMMATE */
	  D("accumulating Epar and Eper");
	  /* accumulate only on processor 0 ! */
	  accumulate((REAL *)Epar,2*(nmax+1));
	  accumulate((REAL *)Eper,2*(nmax+1));
	  D("done accumulating");
	  
	  tstop = extime ();
	  Timing_EField = tstop - tstart;
	  
	  
	  if (ringid == 0) {		 /* WRITE RESULTS TO FILE */
	    tiostart = extime ();
	    switch (whichfile) {
	    case 'X' :
	      strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-par");
	      if(xzplane==true) strcat(stringbuffer,"-yz");
	      filepar = fopen (stringbuffer, "wA");
	      strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-per");
	      if(xzplane==true) strcat(stringbuffer,"-yz");
	      fileper = fopen (stringbuffer, "wA");
	      break;
	    case 'Y' :
	      strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-par");
	      if(xzplane==true) strcat(stringbuffer,"-yz");
	      filepar = fopen (stringbuffer, "wA");
	      strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-per");
	      if(xzplane==true) strcat(stringbuffer,"-yz");
	      fileper = fopen (stringbuffer, "wA");
	      break;
	    }
	    
	    D("writing to files ...\n");
	    for (i = 0; i <= nmax; ++i) {
	      fprintz (filepar,
		       "%.2f % .7E     % .7E\n",
		       180/PI*i*dtheta,Epar[i].r, Epar[i].i);
	      fflush(filepar);
	      fprintz (fileper,
		       "%.2f % .7E     % .7E\n",
		       180/PI*i*dtheta,Eper[i].r, Eper[i].i);
	      fflush(fileper);
	    }
	    fclose (filepar);
	    fclose (fileper);    
	    tiostop = extime ();
	    Timing_FileIO += tiostop - tiostart;
	  }
	} /* end of orient loop */
      }
	
    if(xzplane==true) { /*also calculate fields in xz-plane */
    /* this code was changed and debugged in dec. 2003 */
    /* now it works good */
      if (type==NORMAL) Norient=1; else Norient=2; /* number of orientations */
      tstart = extime ();                  /* Norient is not working !!! */
      
      dtheta = PI / ((double)nTheta);
      
      for(orient=0;orient<Norient;orient++) {
	/* in case of Rotation symmetry */
	if (orient==0) whichfile=which;
	if (orient==1 && which=='X') whichfile='Y'; /* swap X and Y */
	if (orient==1 && which=='Y') whichfile='X';
	
	if ((symY==true && whichfile=='Y') || (symX==true && whichfile=='X'))
	  nmax=nTheta; else nmax=2*nTheta-1;
	
	for (i = 0; i <= nmax; ++i) {
	  if ((i%(20*nprocs))==0) test_interrupt();
	  theta = i * dtheta;
	  
	  robserver[0] = sin (theta);
	  robserver[1] = 0.0;
	  robserver[2] = cos (theta);
	  if (orient==0) { /* normal procedure */
	    robserver[0] = sin (theta);
	    robserver[1] = 0.0;
	    robserver[2] = cos (theta);
	  }
	  else { /* Rotation symmetry: calculate par-par from current data */
	    robserver[0] = 0.0;
	    robserver[1] = sin (theta);
	    robserver[2] = cos (theta);
	  }
	  
	  ebuff[0].r = ebuff[0].i = 0.0;
	  ebuff[1].r = ebuff[1].i = 0.0;
	  ebuff[2].r = ebuff[2].i = 0.0;
	  calc_field(x,ebuff,DipoleCoord, robserver, WaveNum, fint,nlocalDip);
	  
	  /* convert to (l,r) frame */
	  if (orient==0) {
	    if(which=='X'){
	    Eper[i].r = -ebuff[1].r;
	    Eper[i].i = -ebuff[1].i;
	    Epar[i].r = cos (theta) * ebuff[0].r - sin (theta) * ebuff[2].r;
	    Epar[i].i = cos (theta) * ebuff[0].i - sin (theta) * ebuff[2].i;
	    }
	    else {/*which == 'Y', so all fields must be phase shifted 180 degrees */
	    Eper[i].r = ebuff[1].r;
	    Eper[i].i = ebuff[1].i;
	    Epar[i].r = -(cos (theta) * ebuff[0].r - sin (theta) * ebuff[2].r);
	    Epar[i].i = -(cos (theta) * ebuff[0].i - sin (theta) * ebuff[2].i);
	    }
	  }
	  else {
	    Eper[i].r=ebuff[0].r;
	    Eper[i].i=ebuff[0].i;
	    /* these two lines are not checked for bugs */
	    Epar[i].r = cos (theta) * ebuff[1].r - sin (theta) * ebuff[2].r;
	    Epar[i].i = cos (theta) * ebuff[1].i - sin (theta) * ebuff[2].i;
	  }
	} /*  end for i */
	
	/* ACCUMULATE EPAR AND EPER TO ROOT AND SUMMATE */
	D("accumulating Epar and Eper");
	/* accumulate only on processor 0 ! */
	accumulate((REAL *)Epar,2*(nmax+1));
	accumulate((REAL *)Eper,2*(nmax+1));
	D("done accumulating");
	
	tstop = extime ();
	Timing_EField = tstop - tstart;
	
	
	if (ringid == 0) {		/* WRITE RESULTS TO FILE */
	  tiostart = extime ();
	  switch (whichfile) {
	  case 'X' :
	    strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-par-xz");
	    filepar = fopen (stringbuffer, "wA");
	    strcpy(stringbuffer,directory); strcat(stringbuffer,"/par-per-xz");
	    fileper = fopen (stringbuffer, "wA");
	    break;
	  case 'Y' :
	    strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-par-xz");
	    filepar = fopen (stringbuffer, "wA");
	    strcpy(stringbuffer,directory); strcat(stringbuffer,"/per-per-xz");
	    fileper = fopen (stringbuffer, "wA");
	    break;
	  }
	  
	  D("writing to files ...\n");
	  for (i = 0; i <= nmax; ++i) {
	    fprintz (filepar,"%.2f % .7E     % .7E\n",
		     180/PI*i*dtheta,Epar[i].r, Epar[i].i);
	    fflush(filepar);
	    fprintz (fileper,"%.2f % .7E     % .7E\n",
		     180/PI*i*dtheta,Eper[i].r, Eper[i].i);
	    fflush(fileper);
	  }
	  fclose (filepar);
	  fclose (fileper);    
	  tiostop = extime ();
	  Timing_FileIO += tiostop - tiostart;
	}
      } /* end of orient loop */
    }

    if (all_dir)
      /* Calculate the scattered field for the total space-angle
       * in discretized phi and theta */

      calc_alldir(x,DipoleCoord,WaveNum,which);

    if (calc_Cext || calc_Cabs || calc_Csca || calc_asym || calc_mat_force)
      {
	int i,ndips;
	REAL 
	  *Fsca,*Finc,*Frp; /* Scattering force, extinction force and
			       radiation pressure per dipole */
	double 
	  Cext,Cabs,Csca,   /* Cross sections */
	  Csca_diff,        /* formal difference between
			       integrated and direct
			       scattering cross section */
	  g[3],             /* asymmetry paramter */
	  Fsca_tot[3],      /* total scattering force */
	  Finc_tot[3],      /* total extinction force */
	  Frp_tot[3],       /* total radiation pressure */
	  Cnorm,            /* normalizing factor from 
			       Force to cross section */
	  Qnorm,            /* normalizing factor from 
			       Force to efficiency */
	  a_eq,             /* equivalent sphere radius */
	  G;                /* cross section surface */
	FILE 
	  *VisFrp,
	  *CCfile;

	strcpy(stringbuffer,directory);
	if (which == 'X')
	  strcat(stringbuffer,"/CrossSec-X");
	if (which == 'Y')
	  strcat(stringbuffer,"/CrossSec-Y");
	CCfile = fopenz(stringbuffer,"w");

	/* count the number of dipoles i.e. the number of occupied sites */
        ndips = n_non_void_sites(mat_count,Nmat);
	fprintz(logfile,"Total number of dipoles %d\n",ndips);
	a_eq = pow(3*ndips/4/PI,.333333333333)*gridspaceX;
        G = PI*a_eq*a_eq;
        fprintz(CCfile,"x=%lg\n\n",WaveNum*a_eq);

	if (calc_Cext) 
	  {
	    Cext = Ext_cross(x,Einc,WaveNum);
	    printz("Cext\t= %12lg\nQext\t= %12lg\n",Cext,Cext/G);
	    fprintz(CCfile,"Cext\t= %12lg\nQext\t= %12lg\n",Cext,Cext/G);
	  }
	
	if (calc_Cabs) 
	  {
	    Cabs = Abs_cross(x,WaveNum);
	    printz("Cabs\t= %12lg\nQabs\t= %12lg\n",Cabs,Cabs/G);
	    fprintz(CCfile,"Cabs\t= %12lg\nQabs\t= %12lg\n",Cabs,Cabs/G);
	  }

	if (calc_Csca_diff)
	  {
	    MatVec_nim (x, r,
			nlocalDip, nDip, RingId, nRing,
			WaveNum);
  
	    for (i=0; i < local_Ndip; ++i) {
	      if (material[i]<Nmat-1) {
		r[3*i].r -= Einc[3*i].r;
		r[3*i].i -= Einc[3*i].i;
		r[3*i+1].r -= Einc[3*i+1].r;
		r[3*i+1].i -= Einc[3*i+1].i;
		r[3*i+2].r -= Einc[3*i+2].r;
		r[3*i+2].i -= Einc[3*i+2].i;
	      }
	    }

	    Csca_diff = Sca_diff(x,r,WaveNum);
	    printz("Csca_diff\t= %12lg\nQsca_diff\t= %12lg\n",
		   Csca_diff,Csca_diff/G);
	    fprintz(CCfile,"Csca_diff\t= %12lg\nQsca_diff\t= %12lg\n",
		    Csca_diff,Csca_diff/G);
	  }

	if (ringid == 0)
	  {
	    char string[200];
	    
	    if (calc_Csca)
	      {
		fprintf(CCfile,"\nIntegration\n");
                printf("int Csca\n");
		Sca_cross(WaveNum,&Csca);
		printf("Csca\t= %12lg\nQsca\t= %12lg\t  (integration)\n\n", 
		       Csca,Csca/G); 
		fprintf(CCfile,
			"Csca\t= %12lg\nQsca\t= %12lg\n",
			Csca,Csca/G);
	      }

            if (calc_vec)
              {
                double
                  dummy[3];

                printf("\n\nint asym-x\n");
                Asym_parm_x(WaveNum,&dummy[0]);
                printf("g.Csca-x\t= %12lg\n",dummy[0]);
                if (calc_asym)
                  printf("g-x\t\t= %12lg\n",
                         g[0] = dummy[0]/Csca);

                printf("\n\nint asym-y\n");
                Asym_parm_y(WaveNum,&dummy[1]);
                printf("g.Csca-y\t= %12lg\n",dummy[1]);
                if (calc_asym)
                  printf("g-y\t\t= %12lg\n",
                         g[1] = dummy[1]/Csca);

                printf("\n\nint asym-z\n");
                Asym_parm_z(WaveNum,&dummy[2]);
                printf("g.Csca-z\t= %12lg\n",dummy[2]);
                if (calc_asym)
                  printf("g-z\t\t= %12lg\n",
                         g[2] = dummy[2]/Csca);

                fprintf(CCfile,"Csca.g\t=(%12lg,%12lg,%12lg)\n",
                        dummy[0],dummy[1],dummy[2]);
                if (calc_asym)
                  fprintf(CCfile,"g\t=(%12lg,%12lg,%12lg)\n",
                          g[0],g[1],g[2]);
              }
	  }

	if (calc_mat_force)
	  {
	    if ((Fsca = (REAL *) calloc(3*local_Ndip,sizeof(REAL))) == NULL) {
	      LogError (EC_ERROR, Fname,
			"processor %d, ringID %d, could not malloc Fsca",
			MyProcId, RingId);
	      AbortServer (1);
	    }
	    if ((Finc = (REAL *) calloc(3*local_Ndip,sizeof(REAL))) == NULL) {
	      LogError (EC_ERROR, Fname,
			"processor %d, ringID %d, could not malloc Finc",
			MyProcId, RingId);
	      AbortServer (1);
	    }
	    if ((Frp = (REAL *) calloc(3*local_Ndip,sizeof(REAL))) == NULL) {
	      LogError (EC_ERROR, Fname,
			"processor %d, ringID %d, could not malloc Frp",
			MyProcId, RingId);
	      AbortServer (1);
	    }
	    printz("Calculating the force per dipole\n");

	    /* Calculate forces */
	    Frp_mat(Fsca_tot,Fsca,Finc_tot,Finc,Frp_tot,Frp,
		    x,DipoleCoord,WaveNum);

	    /* Write Cross-Sections and Efficiencies to file */
	    Cnorm = 8*PI;
	    Qnorm = 8*PI/G;
	    printz("\nMatrix\n"\
		   "Cext\t= %12lg\nQext\t= %12lg\n"\
		   "Csca.g\t= (%12lg,%12lg,%12lg)\n"\
		   "Cpr\t= (%12lg,%12lg,%12lg)\n"\
		   "Qpr\t= (%12lg,%12lg,%12lg)\n",
		   Cnorm*Finc_tot[2],Qnorm*Finc_tot[2],
		   -Cnorm*Fsca_tot[0],-Cnorm*Fsca_tot[1],-Cnorm*Fsca_tot[2],
		   Cnorm*Frp_tot[0],Cnorm*Frp_tot[1],Cnorm*Frp_tot[2],
		   Qnorm*Frp_tot[0],Qnorm*Frp_tot[1],Qnorm*Frp_tot[2]);	    
	    fprintz(CCfile,"\nMatrix\n"\
		    "Cext\t= %lg\nQext\t= %12lg\n"\
		    "Csca.g\t= (%12lg,%12lg,%12lg)\n"\
		    "Cpr\t= (%12lg,%12lg,%12lg)\n"\
		    "Qpr\t= (%12lg,%12lg,%12lg)\n",
		    Cnorm*Finc_tot[2],Qnorm*Finc_tot[2],
		    -Cnorm*Fsca_tot[0],-Cnorm*Fsca_tot[1],-Cnorm*Fsca_tot[2],
		    Cnorm*Frp_tot[0],Cnorm*Frp_tot[1],Cnorm*Frp_tot[2],
		    Qnorm*Frp_tot[0],Qnorm*Frp_tot[1],Qnorm*Frp_tot[2]);
      
	    if (store_force)
	      {
		/* Write Radiation pressure per dipole to file */
		strcpy(stringbuffer,directory);
		if (which == 'X')
		  strcat(stringbuffer,"/VisFrp-X.dat");
		if (which == 'Y')
		  strcat(stringbuffer,"/VisFrp-Y.dat");
		VisFrp = fopen(stringbuffer,"w");
	  
		fprintf(VisFrp,"#sphere\t\t\tx=%lg\tm=%lg+i%lg\n"\
			"#number of dipoles\t%d\n"\
			"#Forces per dipole\n"\
			"#r.x\t\tr.y\t\tr.z\t\tF.x\t\tF.y\t\tF.z\n",
			WaveNum*a_eq,ref_index[0].r,ref_index[0].i,
			(int) ndips);
		for (j=0;j<local_Ndip;++j)
		  if (material[j]<Nmat-1)
		    {
		      /*if (sizeof(REAL) == sizeof(float))*/
			fprintf(VisFrp,
				"%8g\t%8g\t%8g\t"\
				"%8g\t%8g\t%8g\n",
				DipoleCoord[j][0],DipoleCoord[j][1],
				DipoleCoord[j][2],
				Frp[3*j],Frp[3*j+1],Frp[3*j+2]);
		    }
		
		fclose(VisFrp);
	      }
	    
	    free(Fsca);
	    free(Finc);
	    free(Frp);
	  }

	fclosez(CCfile);
      }
    if (all_dir)
      finish_int();


    if (which=='X') if (store_field==true) pbm_fields(which,0);
    if (which=='Y') if (store_field==true) pbm_fields(which,1);
  }

  /* save effective volume in_mie */
  if (ringid==0) save_in_mie();
  
  /* all processors should write the internal field, stored in */
  /* the vector x to file. Files get the name IntfieldY-procid */
  /* or IntfieldXr-procid */

  if (store_int_field == true) { 
  /* Write internal field on each dipole to file */
     if (which == 'X')
       sprintf(stringbuffer,"%s%s%i",directory,"/IntFieldX_",ringid);
     if (which == 'Y')
       sprintf(stringbuffer,"%s%s%i",directory,"/IntFieldY_",ringid);
     Intfldpntr = fopen(stringbuffer,"w");
     if(ringid==0)
       fprintf(Intfldpntr,
       "x\t\ty\t\tz\t\t|E|^2\t\tEx.r\t\tEx.i\t\tEy.r\t\tEy.i\t\tEz.r\t\tEz.i\t\t\n\n");
     for (i=0;i<local_Ndip;++i)
       if (material[i]<Nmat-1) {
         fprintf(Intfldpntr,
           "%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\t%8g\n",
	   DipoleCoord[i][0],DipoleCoord[i][1], DipoleCoord[i][2],
	   x[3*i].r * x[3*i].r + x[3*i].i * x[3*i].i+
	   x[3*i+1].r * x[3*i+1].r + x[3*i+1].i * x[3*i+1].i+
	   x[3*i+2].r * x[3*i+2].r + x[3*i+2].i * x[3*i+2].i,
	   x[3*i].r,x[3*i].i,x[3*i+1].r,x[3*i+1].i,x[3*i+2].r,x[3*i+2].i);
       }
     }
}
