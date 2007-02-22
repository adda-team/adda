/* FILE: CalculateE.c
 * AUTH: Maxim Yurkin
 * DESCR: The module that will calculate the E field and all scattering quantities.
 *        Routines for most scattering quantities are in crosssec.c
 *        Also saves internal fields to file (optional).
 *        
 *        January 2004 : include module to compute full Mueller Matrix over
 *        full space angle, not very efficient, must be improved (A. Hoekstra)
 *
 *        Previous versions by Alfons Hoekstra
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <time.h>
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "crosssec.h"
#include "Romberg.h"
#include "io.h"
#include "vars.h"

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in calculator.c */
extern double *muel_phi,*muel_phi1;
extern doublecomplex *EplaneX, *EplaneY;
extern double *Eplane_buffer;
extern const double dtheta_deg,dtheta_rad;
extern doublecomplex *ampl_alphaX,*ampl_alphaY;
extern double *muel_alpha;
/* defined and initialized in crosssec.c */
extern const Parms_1D phi_sg;
/* defined and initialized in param.c */
extern const int store_int_field;
extern const int store_scat_grid;
extern const int calc_Cext;
extern const int calc_Cabs;
extern const int calc_Csca;
extern const int calc_vec;
extern const int calc_asym;
extern const int calc_mat_force;
extern const int store_force;
extern const int phi_int_type;
/* defined and initialized in timing.c */
extern clock_t Timing_EFieldPlane,Timing_calc_EField,Timing_comm_EField,
               Timing_IntField,Timing_IntFieldOne,Timing_ScatQuan;
extern unsigned long TotalEFieldPlane;

/* used in iterative.c */
clock_t tstart_CE;

/* EXTERNAL FUNCTIONS */

/* GenerateB.c */
void GenerateB(char which,doublecomplex *x);
/* iterative.c */
int IterativeSolver(int method);

/*============================================================*/

static void ComputeMuellerMatrix(double matrix[4][4], const doublecomplex s1,
   const doublecomplex s2,const doublecomplex s3,const doublecomplex s4)
/* computer mueller matrix from scattering matrix elements s1, s2, s3, s4, accoording
   to formula 3.16 from Bohren and Huffman */
{
	matrix[0][0] = 0.5*(cMultConRe(s1,s1)+cMultConRe(s2,s2)+
	                    cMultConRe(s3,s3)+cMultConRe(s4,s4));
	matrix[0][1] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+
	                    cMultConRe(s4,s4)-cMultConRe(s3,s3));
	matrix[0][2] = cMultConRe(s2,s3)+cMultConRe(s1,s4);
	matrix[0][3] = cMultConIm(s2,s3)-cMultConIm(s1,s4);
	
	matrix[1][0] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+
	                    cMultConRe(s3,s3)-cMultConRe(s4,s4));
	matrix[1][1] = 0.5*(cMultConRe(s2,s2)+cMultConRe(s1,s1)+
	                    -cMultConRe(s3,s3)-cMultConRe(s4,s4));
	matrix[1][2] = cMultConRe(s2,s3)-cMultConRe(s1,s4);
	matrix[1][3] = cMultConIm(s2,s3)+cMultConIm(s1,s4);

	matrix[2][0] = cMultConRe(s2,s4)+cMultConRe(s1,s3);
	matrix[2][1] = cMultConRe(s2,s4)-cMultConRe(s1,s3);
	matrix[2][2] = cMultConRe(s1,s2)+cMultConRe(s3,s4);
	matrix[2][3] = cMultConIm(s2,s1)+cMultConIm(s4,s3);

	matrix[3][0] = cMultConIm(s4,s2)+cMultConIm(s1,s3);
	matrix[3][1] = cMultConIm(s4,s2)-cMultConIm(s1,s3);
	matrix[3][2] = cMultConIm(s1,s2)-cMultConIm(s3,s4);
	matrix[3][3] = cMultConRe(s1,s2)-cMultConRe(s3,s4);
}

/*============================================================*/

static void ComputeMuellerMatrixNorm(double matrix[4][4],const doublecomplex s1,
   const doublecomplex s2,const doublecomplex s3,const doublecomplex s4)
/* computer mueller matrix from scattering matrix elements s1, s2, s3, s4, accoording to 
   formula 3.16 from Bohren and Huffman; normalize all elements to S11 (except itself)*/
{
	matrix[0][0] = 0.5*(cMultConRe(s1,s1)+cMultConRe(s2,s2)+
	                    cMultConRe(s3,s3)+cMultConRe(s4,s4));
	matrix[0][1] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+
	                    cMultConRe(s4,s4)-cMultConRe(s3,s3))/matrix[0][0];
	matrix[0][2] = (cMultConRe(s2,s3)+cMultConRe(s1,s4))/matrix[0][0];
	matrix[0][3] = (cMultConIm(s2,s3)-cMultConIm(s1,s4))/matrix[0][0];

	matrix[1][0] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+
	                    cMultConRe(s3,s3)-cMultConRe(s4,s4))/matrix[0][0];
	matrix[1][1] = 0.5*(cMultConRe(s2,s2)+cMultConRe(s1,s1)+
	                    -cMultConRe(s3,s3)-cMultConRe(s4,s4))/matrix[0][0];
	matrix[1][2] = (cMultConRe(s2,s3)-cMultConRe(s1,s4))/matrix[0][0];
	matrix[1][3] = (cMultConIm(s2,s3)+cMultConIm(s1,s4))/matrix[0][0];

	matrix[2][0] = (cMultConRe(s2,s4)+cMultConRe(s1,s3))/matrix[0][0];
	matrix[2][1] = (cMultConRe(s2,s4)-cMultConRe(s1,s3))/matrix[0][0];
	matrix[2][2] = (cMultConRe(s1,s2)+cMultConRe(s3,s4))/matrix[0][0];
	matrix[2][3] = (cMultConIm(s2,s1)+cMultConIm(s4,s3))/matrix[0][0];

	matrix[3][0] = (cMultConIm(s4,s2)+cMultConIm(s1,s3))/matrix[0][0];
	matrix[3][1] = (cMultConIm(s4,s2)-cMultConIm(s1,s3))/matrix[0][0];
	matrix[3][2] = (cMultConIm(s1,s2)-cMultConIm(s3,s4))/matrix[0][0];
	matrix[3][3] = (cMultConRe(s1,s2)-cMultConRe(s3,s4))/matrix[0][0];
}

/*==============================================================*/

void MuellerMatrix(void)
{
  FILE *mueller,*mueller_int,*mueller_int_c2,*mueller_int_s2,*mueller_int_c4,*mueller_int_s4;
  double *cos2,*sin2,*cos4,*sin4;
  double matrix[4][4];
  double theta,phi,ph,err,
         max_err,max_err_c2,max_err_s2,max_err_c4,max_err_s4;
  doublecomplex s1,s2,s3,s4,s10,s20,s30,s40;
  char fname[MAX_FNAME];
  int index,index1,i,j,k,n,k_or;
  double co,si;
  double alph;
  clock_t tstart;

  if (ringid!=ROOT) return;

  if (orient_avg) {                 /* Amplitude matrix stored in ampl_alplha is */
    index1=index=0;                 /* transformed into Mueller matrix stored in muel_alpha */
    for (k_or=0;k_or<alpha_int.N;k_or++) {
      alph=Deg2Rad(alpha_int.val[k_or]);    /* read current alpha */
      co=cos(alph);
      si=sin(alph);
      for (i=0;i<nTheta;i++) {
        /* read amplitude matrix from memory */
        memcpy(s10,&ampl_alphaX[index],sizeof(doublecomplex));
        memcpy(s30,&ampl_alphaX[index+1],sizeof(doublecomplex));
        memcpy(s40,&ampl_alphaY[index],sizeof(doublecomplex));
        memcpy(s20,&ampl_alphaY[index+1],sizeof(doublecomplex));
        /* transform it, multiplying by rotation matrix (-alpha) */
        cLinComb(s20,s30,co,si,s2);         /* s2 =  co*s20 + si*s30  */
        cLinComb(s20,s30,-si,co,s3);        /* s3 = -si*s20 + co*s30  */
        cLinComb(s40,s10,co,si,s4);         /* s4 =  co*s40 + si*s10  */
        cLinComb(s40,s10,-si,co,s1);        /* s1 = -si*s40 + co*s10  */

        ComputeMuellerMatrix((double (*)[4])(muel_alpha+index1),s1,s2,s3,s4);
        index+=2;
        index1+=16;
      }
    }
  }
  else {                                /* here amplitude matrix is read from file */
    tstart=clock();                     /* and Mueller matrix is saved to file     */
    if (yzplane) {
      strcpy(fname,directory);
      strcat(fname,"/" F_MUEL);
      mueller=FOpenErr(fname,"w",ONE_POS);
      fprintf(mueller,"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");
      for (i=0;i<nTheta;i++) {
        theta=i*dtheta_deg;
	ComputeMuellerMatrix(matrix,EplaneX[2*i],EplaneY[2*i+1],EplaneX[2*i+1],EplaneY[2*i]);
        fprintf(mueller,
		"%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"\
                " %.10E %.10E %.10E %.10E %.10E %.10E %.10E\n",
		theta,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
                matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
		matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
		matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3]);
      }
      FCloseErr(mueller,F_MUEL,ONE_POS);
    }

    if (scat_grid) {
      /* compute Mueller Matrix in full space angle.
       * E-fields are stored in arrays EgridX and EgridY
       * for incoming X and Y polarized light. From this compute the Mueller matrix.
       * It is converted to the scattering matrix elements (see e.g Bohren and Huffman) :
       * s2 = cos(phi)E'X'par + sin(phi)E'Y'par
       * s3 = sin(phi)E'X'par - cos(phi)E'Y'par
       * s4 = cos(phi)E'X'per + sin(phi)E'Y'per
       * s1 = sin(phi)E'X'per - cos(phi)E'Y'per
       * from this the mueller matrix elements are computed */

      /* open files for writing */
      if (store_scat_grid) {
        strcpy(fname,directory);
        strcat(fname,"/" F_MUEL_SG);
        mueller=FOpenErr(fname,"w",ONE_POS);
        fprintf(mueller,"theta phi s11 s12 s13 s14 s21 s22 s23 s24"\
                        " s31 s32 s33 s34 s41 s42 s43 s44\n");
      }
      if (phi_integr) {   /* also initializes arrays of multipliers */
        if (phi_int_type & PHI_UNITY) {
          strcpy(fname,directory);
          strcat(fname,"/" F_MUEL_INT);
          mueller_int=FOpenErr(fname,"w",ONE_POS);
          fprintf(mueller_int,"theta s11 s12 s13 s14 s21 s22 s23 s24"\
                              " s31 s32 s33 s34 s41 s42 s43 s44 RMSE(integr)\n");
        }
        if (phi_int_type & PHI_COS2) {
          strcpy(fname,directory);
          strcat(fname,"/" F_MUEL_C2);
          mueller_int_c2=FOpenErr(fname,"w",ONE_POS);
          fprintf(mueller_int_c2,"theta s11 s12 s13 s14 s21 s22 s23 s24"\
                                 " s31 s32 s33 s34 s41 s42 s43 s44 RMSE(integr)\n");
          if ((cos2 = (double *) malloc(angles.phi.N*sizeof(double))) == NULL)
	    LogError(EC_ERROR,ALL_POS,"Could not malloc cos2");
        }
        if (phi_int_type & PHI_SIN2) {
          strcpy(fname,directory);
          strcat(fname,"/" F_MUEL_S2);
          mueller_int_s2=FOpenErr(fname,"w",ONE_POS);
          fprintf(mueller_int_s2,"theta s11 s12 s13 s14 s21 s22 s23 s24"\
                                 " s31 s32 s33 s34 s41 s42 s43 s44 RMSE(integr)\n");
          if ((sin2 = (double *) malloc(angles.phi.N*sizeof(double))) == NULL)
	    LogError(EC_ERROR,ALL_POS,"Could not malloc sin2");
        }
        if (phi_int_type & PHI_COS4) {
          strcpy(fname,directory);
          strcat(fname,"/" F_MUEL_C4);
          mueller_int_c4=FOpenErr(fname,"w",ONE_POS);
          fprintf(mueller_int_c4,"theta s11 s12 s13 s14 s21 s22 s23 s24"\
                                 " s31 s32 s33 s34 s41 s42 s43 s44 RMSE(integr)\n");
          if ((cos4 = (double *) malloc(angles.phi.N*sizeof(double))) == NULL)
	    LogError(EC_ERROR,ALL_POS,"Could not malloc cos4");
        }
        if (phi_int_type & PHI_SIN4) {
          strcpy(fname,directory);
          strcat(fname,"/" F_MUEL_S4);
          mueller_int_s4=FOpenErr(fname,"w",ONE_POS);
          fprintf(mueller_int_s4,"theta s11 s12 s13 s14 s21 s22 s23 s24"\
                                 " s31 s32 s33 s34 s41 s42 s43 s44 RMSE(integr)\n");
          if ((sin4 = (double *) malloc(angles.phi.N*sizeof(double))) == NULL)
	    LogError(EC_ERROR,ALL_POS,"Could not malloc sin4");
        }
        /* fills arrays with multipliers (optimized) */
        for (j=0;j<angles.phi.N;j++) {
          /* prepare */
          ph=2*Deg2Rad(angles.phi.val[j]);
          if (phi_int_type & (PHI_COS2|PHI_COS4|PHI_SIN4)) co=cos(ph);
          if (phi_int_type & (PHI_SIN2|PHI_SIN4)) si=sin(ph);
          /* fill */
          if (phi_int_type & PHI_COS2) cos2[j]=co;
          if (phi_int_type & PHI_SIN2) sin2[j]=si;
          if (phi_int_type & PHI_COS4) cos4[j]=2*co*co-1;
          if (phi_int_type & PHI_SIN4) sin4[j]=2*si*co;
        }
      }
      /* set type of cycling through angles */
      if (angles.type==SG_GRID) n=angles.phi.N;
      else if (angles.type==SG_PAIRS) n=1;
      /* main cycle */
      index=0;
      max_err=max_err_c2=max_err_s2=max_err_c4=max_err_s4=0;
      for (i=0;i<angles.theta.N;++i) {
        index1=0;
        theta=angles.theta.val[i];
        for (j=0;j<n;++j) {
          if (angles.type==SG_GRID) phi=angles.phi.val[j];
          else if (angles.type==SG_PAIRS) phi=angles.phi.val[i];
          ph=Deg2Rad(phi);
          co=cos(ph);
          si=sin(ph);
          /* read amplitude matrix from memory */
          memcpy(s10,EgridY+index,sizeof(doublecomplex));
          memcpy(s30,EgridY+index+1,sizeof(doublecomplex));
          memcpy(s40,EgridX+index,sizeof(doublecomplex));
          memcpy(s20,EgridX+index+1,sizeof(doublecomplex));
          /* transform it, multiplying by rotation matrix from per-par to X-Y */
          cLinComb(s20,s30,co,si,s2);         /* s2 =  co*s20 + si*s30  */
          cLinComb(s20,s30,si,-co,s3);        /* s3 =  si*s20 - co*s30  */
          cLinComb(s40,s10,co,si,s4);         /* s4 =  co*s40 + si*s10  */
          cLinComb(s40,s10,si,-co,s1);        /* s1 =  si*s40 - co*s10  */

          ComputeMuellerMatrix(matrix,s1,s2,s3,s4);
          index+=2;
          if (phi_integr) {
            memcpy(muel_phi+index1,matrix[0],16*sizeof(double));
            index1+=16;
          }
          if (store_scat_grid)
            fprintf(mueller,
	      "%.2f %.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"\
              " %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E\n",
	      theta,phi,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
              matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
              matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
              matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3]);
        }
        if (phi_integr) {
          if (phi_int_type & PHI_UNITY) {
            err=Romberg1D(phi_sg,16,muel_phi,matrix[0]);
            if (err>max_err) max_err=err;
            fprintf(mueller_int,
	      "%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"\
              " %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.3E\n",
              theta,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
              matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
              matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
              matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3],err);
          }
          if (phi_int_type & PHI_COS2) {
            for (j=0;j<angles.phi.N;j++) for(k=0;k<16;k++)
              muel_phi1[16*j+k]=muel_phi[16*j+k]*cos2[j];
            err=Romberg1D(phi_sg,16,muel_phi1,matrix[0]);
            if (err>max_err_c2) max_err_c2=err;
            fprintf(mueller_int_c2,
	      "%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"\
              " %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.3E\n",
              theta,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
              matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
              matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
              matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3],err);
          }
          if (phi_int_type & PHI_SIN2) {
            for (j=0;j<angles.phi.N;j++) for(k=0;k<16;k++)
              muel_phi1[16*j+k]=muel_phi[16*j+k]*sin2[j];
            err=Romberg1D(phi_sg,16,muel_phi1,matrix[0]);
            if (err>max_err_s2) max_err_s2=err;
            fprintf(mueller_int_s2,
	      "%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"\
              " %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.3E\n",
              theta,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
              matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
              matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
              matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3],err);
          }
          if (phi_int_type & PHI_COS4) {
            for (j=0;j<angles.phi.N;j++) for(k=0;k<16;k++)
              muel_phi1[16*j+k]=muel_phi[16*j+k]*cos4[j];
            err=Romberg1D(phi_sg,16,muel_phi1,matrix[0]);
            if (err>max_err_c4) max_err_c4=err;
            fprintf(mueller_int_c4,
	      "%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"\
              " %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.3E\n",
              theta,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
              matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
              matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
              matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3],err);
          }
          if (phi_int_type & PHI_SIN4) {
            for (j=0;j<angles.phi.N;j++) for(k=0;k<16;k++)
              muel_phi1[16*j+k]=muel_phi[16*j+k]*sin4[j];
            err=Romberg1D(phi_sg,16,muel_phi1,matrix[0]);
            if (err>max_err_s4) max_err_s4=err;
            fprintf(mueller_int_s4,
	      "%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"\
              " %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.3E\n",
              theta,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
              matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
              matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
              matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3],err);
          }
        }
      }
      if (phi_integr) {
        fprintf(logfile,"\nMaximum relative mean-square error of Mueller integration:\n");
        if (phi_int_type & PHI_UNITY) fprintf(logfile,"  1          -> %.3E\n",max_err);
        if (phi_int_type & PHI_COS2) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
        if (phi_int_type & PHI_SIN2) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
        if (phi_int_type & PHI_COS4) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
        if (phi_int_type & PHI_SIN4) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
      }
      /* close files; free arrays */
      if (store_scat_grid) FCloseErr(mueller,F_MUEL_SG,ONE_POS);
      if (phi_integr) {
        if (phi_int_type & PHI_UNITY) FCloseErr(mueller_int,F_MUEL_INT,ONE_POS);
        if (phi_int_type & PHI_COS2) {
          FCloseErr(mueller_int_c2,F_MUEL_C2,ONE_POS);
          free(cos2);
        }
        if (phi_int_type & PHI_SIN2) {
          FCloseErr(mueller_int_s2,F_MUEL_S2,ONE_POS);
          free(sin2);
        }
        if (phi_int_type & PHI_COS4) {
          FCloseErr(mueller_int_c4,F_MUEL_C4,ONE_POS);
          free(cos4);
        }
        if (phi_int_type & PHI_SIN4) {
          FCloseErr(mueller_int_s4,F_MUEL_S4,ONE_POS);
          free(sin4);
        }
      }
    }
    Timing_FileIO += clock() - tstart;
  }
}

/*============================================================*/

static void CalcEplane(const char which,const int type)
  /* calculates scattered electric field in a plane */
{
  double *incPol,*incPolper,*incPolpar;
  doublecomplex *Eplane; /* where to store calculated field for one plane
                            (actually points to different other arrays) */
  int i;
  doublecomplex ebuff[3];	/* small vector to hold E fields */
  double robserver[3];	/* small vector for observer in E calculation */
  double epar[3];       /* unit vector in direction of Epar */
  double theta; 	/* scattering angle */
  double co,si;          /* temporary, cos and sin of some angle */
  double incPol_tmp1[3],incPol_tmp2[3]; /* just allocateed memory for incPolper, incPolpar */
  double alph;
  clock_t tstart, tstart2;
  int k_or;
  int orient,Norient;
  char choice;

  incPolper=incPol_tmp1;  /* initialization of per and par polarizations */
  incPolpar=incPol_tmp2;

  if (type==CE_NORMAL) Norient=1; /* initialize # orientations */
  else if (type==CE_PARPER) Norient=2;

  for (k_or=0;k_or<alpha_int.N;k_or++) {
    /* cycle over alpha - for orientation averaging */
    if (orient_avg) {
      alph=Deg2Rad(alpha_int.val[k_or]);   /* rotate polarization basis vectors by -alpha */
      co=cos(alph);
      si=sin(alph);
      LinComb(incPolX,incPolY,co,-si,incPolper);  /* incPolper = co*incPolX - si*incPolY; */
      LinComb(incPolX,incPolY,si,co,incPolpar);   /* incPolpar = si*incPolX + co*incPolY; */
    }
    else {
      memcpy(incPolper,incPolX,3*sizeof(double));   /* per <=> X */
      memcpy(incPolpar,incPolY,3*sizeof(double));   /* par <=> Y */
    }

    for(orient=0;orient<Norient;orient++) {
      /* in case of Rotation symmetry */
      tstart = clock ();
      if (orient==0) choice=which;
      else if (orient==1) {
        /* Rotation symmetry: calculate per-per from current data. CalculateE is called
         * from calculator with 'Y' polarization - we now just assume that we have
         * the x-z plane as the scatering plane, rotating in the negative x-direction.
         * This mimics the real case of X polarization with the y-z plane as scattering plane
         * Then IncPolY -> -IncPolX; incPolX -> IncPolY */
        if (which=='X') choice='Y';
        else if (which=='Y') choice='X';
        incPol=incPolper;
        incPolper=incPolpar;
        incPolpar=incPol;
        MultScal(-1,incPolpar,incPolpar);
      }
      /* initialize Eplane */
      if (orient_avg) {
        if (choice=='X') Eplane=ampl_alphaX + 2*nTheta*k_or;
        else if (choice=='Y') Eplane=ampl_alphaY + 2*nTheta*k_or;
      }
      else {
        if (choice=='X') Eplane=EplaneX;
        else if (choice=='Y') Eplane=EplaneY;
      }

      for (i=0;i<nTheta;i++) {
        theta = i * dtheta_rad;
        co=cos(theta);
        si=sin(theta);
        LinComb(prop,incPolpar,co,si,robserver);  /* robserver = co*prop + si*incPolpar; */

        CalcField(ebuff,robserver);
         /* convert to (l,r) frame */
        crDotProd(ebuff,incPolper,Eplane[2*i]);     /* Eper[i]=Esca.incPolper */
        LinComb(prop,incPolpar,-si,co,epar);        /* epar=-si*prop+co*incPolpar */
        crDotProd(ebuff,epar,Eplane[2*i+1]); 	      /* Epar[i]=Esca.epar */
      } /*  end for i */

      /* Accumulate Eplane to root and summate */
      D("accumulating Eplane");
      /* accumulate only on processor 0 !, done in one operation */
      tstart2=clock();
      Accumulate((double *)Eplane,4*nTheta,Eplane_buffer);
      Timing_comm_EField = clock() - tstart2;
      D("done accumulating");

      Timing_EFieldPlane = clock() - tstart;
      Timing_calc_EField = Timing_EFieldPlane - Timing_comm_EField;
      Timing_EField += Timing_EFieldPlane;
      TotalEFieldPlane++;
    } /* end of orient loop */
  } /* end of alpha loop */
}

/*============================================================*/

static void CalcIntegralScatQuantities(const char which)
 /* calculates all the scattering crosssections, normalized and unnormalized
    asymmetry parameter, and force on the particle and each dipole.
    Cext and Cabs are overaged over orientation, if needed. */
{
  double *Fsca,*Finc,*Frp; /* Scattering force, extinction force and
  		             radiation pressure per dipole */
  double Cext,Cabs,Csca,   /* Cross sections */
         dummy[3],         /* asymmetry paramter*Csca */
         Fsca_tot[3],      /* total scattering force */
         Finc_tot[3],      /* total extinction force */
         Frp_tot[3],       /* total radiation pressure */
         Cnorm,            /* normalizing factor from force to cross section */
         Qnorm;            /* normalizing factor from force to efficiency */
  FILE *VisFrp,*CCfile;
  clock_t tstart;
  char fname_cs[MAX_FNAME],fname_frp[MAX_FNAME];
  int j;
  double const *incPol;

  D("Calculation of cross sections started");
  tstart = clock();

  strcpy(fname_cs,directory);
  strcat(fname_cs,"/" F_CS);
  if (which == 'X') {
    strcat(fname_cs,F_XSUF);
    incPol=incPolX;
  }
  if (which == 'Y') {
    strcat(fname_cs,F_YSUF);
    incPol=incPolY;
  }
  if (calc_Cext) Cext = ExtCross(incPol);
  if (calc_Cabs) Cabs = AbsCross();
  D("Cext and Cabs calculated");
  if (orient_avg) {
    if (ringid==ROOT) {
      if (which == 'Y') {       /* assumed that first call of CalculateE is with 'Y' flag */
        muel_alpha[-2]=Cext;
        muel_alpha[-1]=Cabs;
      }
      else if (which == 'X') {
        muel_alpha[-2]=(muel_alpha[-2]+Cext)/2;
        muel_alpha[-1]=(muel_alpha[-1]+Cabs)/2;
      }
    }
  }
  else {  /* not orient_avg */
    if (ringid==ROOT) {
      CCfile=FOpenErr(fname_cs,"w",ONE_POS);
      if (calc_Cext) PrintBoth(CCfile,"Cext\t= %.10g\nQext\t= %.10g\n",Cext,Cext*inv_G);
      if (calc_Cabs) PrintBoth(CCfile,"Cabs\t= %.10g\nQabs\t= %.10g\n",Cabs,Cabs*inv_G);
      if (all_dir) fprintf(CCfile,"\nIntegration\n\n");
      if (calc_Csca) {
        Csca=ScaCross();
        PrintBoth(CCfile,"Csca\t= %.10g\nQsca\t= %.10g\n",Csca,Csca*inv_G);
      }
      if (calc_vec) {
        AsymParm_x(dummy);
        AsymParm_y(dummy+1);
        AsymParm_z(dummy+2);
        PrintBoth(CCfile,"Csca.g\t=(%.10g,%.10g,%.10g)\n",dummy[0],dummy[1],dummy[2]);
        if (calc_asym) PrintBoth(CCfile,"g\t=(%.10g,%.10g,%.10g)\n",
                                 dummy[0]/Csca,dummy[1]/Csca,dummy[2]/Csca);
      }
    } /* end of ROOT */
    if (calc_mat_force) {
      if ((Fsca = (double *) calloc(3*local_nvoid_Ndip,sizeof(double))) == NULL)
        LogError(EC_ERROR,ALL_POS,"Could not malloc Fsca");
      if ((Finc = (double *) calloc(3*local_nvoid_Ndip,sizeof(double))) == NULL)
        LogError(EC_ERROR,ALL_POS,"Could not malloc Finc");
      if ((Frp = (double *) calloc(3*local_nvoid_Ndip,sizeof(double))) == NULL)
        LogError(EC_ERROR,ALL_POS,"Could not malloc Frp");

      PRINTZ("Calculating the force per dipole\n");

      /* Calculate forces */
      Frp_mat(Fsca_tot,Fsca,Finc_tot,Finc,Frp_tot,Frp);

      /* Write Cross-Sections and Efficiencies to file */
      if (ringid==ROOT) {
        Cnorm = 8*PI;
        Qnorm = 8*PI*inv_G;
        PrintBoth(CCfile,"\nMatrix\n"\
                  "Cext\t= %.10g\nQext\t= %.10g\n"\
                  "Csca.g\t= (%.10g,%.10g,%.10g)\n"\
                  "Cpr\t= (%.10g,%.10g,%.10g)\n"\
                  "Qpr\t= (%.10g,%.10g,%.10g)\n",
                  Cnorm*Finc_tot[2],Qnorm*Finc_tot[2],
                  -Cnorm*Fsca_tot[0],-Cnorm*Fsca_tot[1],-Cnorm*Fsca_tot[2],
                  Cnorm*Frp_tot[0],Cnorm*Frp_tot[1],Cnorm*Frp_tot[2],
                  Qnorm*Frp_tot[0],Qnorm*Frp_tot[1],Qnorm*Frp_tot[2]);
        if (store_force) {
          /* Write Radiation pressure per dipole to file */
          strcpy(fname_frp,directory);
          strcat(fname_frp,"/" F_FRP);
          if (which == 'X') strcat(fname_frp,F_XSUF);
          else if (which == 'Y') strcat(fname_frp,F_YSUF);
          strcat(fname_frp,".dat"); /* should be removed in the future */
          VisFrp=FOpenErr(fname_frp,"w",ONE_POS);
          fprintf(VisFrp,"#sphere  x=%.10g  m=%.10g%+.10gi\n"\
                  "#number of dipoles  %d\n"\
                  "#Forces per dipole\n"\
                  "#r.x r.y r.z F.x F.y F.z\n",
                  ka_eq,ref_index[0][RE],ref_index[0][IM],
                  nvoid_Ndip);
          for (j=0;j<local_nvoid_Ndip;++j) fprintf(VisFrp,
                "%.10g %.10g %.10g %.10g %.10g %.10g\n",
                DipoleCoord[3*j],DipoleCoord[3*j+1],
                DipoleCoord[3*j+2],
                Frp[3*j],Frp[3*j+1],Frp[3*j+2]);
          FCloseErr(VisFrp,fname_frp,ONE_POS);
        }
      }
      free(Fsca);
      free(Finc);
      free(Frp);
    }
    if (ringid==ROOT) FCloseErr(CCfile,fname_cs,ONE_POS);
  }
  D("Calculation of cross sections complete");
  Timing_ScatQuan += clock() - tstart;
}

/*============================================================*/

static void StoreIntFields(const char which)
    /* Write actual internal fields (not exciting) on each dipole to file
       all processors should write the internal field, stored in the vector x to file.
       Files get the name IntfieldY-procid or IntfieldXr-procid.
       If CE_PARPER is employed then naturally saves only once;
       use '-sym no' if needed */
{
  FILE *intfld;	   /* file to store internal fields */
  double V;
  doublecomplex hi,hi_inv[MAX_NMAT],fld[3];
  char mat;
  int i,j,buf_size;
  clock_t tstart;
  char fname[MAX_FNAME],intfld_fname[MAX_FNAME_SH];

  tstart=clock();
  /* calculate multipliers */
  V=gridspace*gridspace*gridspace;
  for (i=0;i<Nmat;i++) {
    /* hi_inv=1/(V*hi)=4*PI/(V(m^2-1)) */
    cSquare(ref_index[i],hi);
    hi[RE]-=1;
    cMultReal(V,hi,hi);
    cInv(hi,hi_inv[i]);
    cMultReal(4*PI,hi_inv[i],hi_inv[i]);
  }
  /* choose inffld_fname */
  strcpy(intfld_fname,F_INTFLD);
  if (which=='X') strcat(intfld_fname,F_XSUF);
  else if (which=='Y') strcat(intfld_fname,F_YSUF);
  /* choose filename for direct saving */
#ifdef PARALLEL
  sprintf(fname,"%s/" F_INTFLD_TMP,directory,ringid);
#else
  sprintf(fname,"%s/%s",directory,intfld_fname);
#endif
  intfld=FOpenErr(fname,"w",ALL_POS);
  /* print head of file */
#ifdef PARALLEL
  if (ringid==0) {  /* this condition can be different from being ROOT */
#endif
    fprintf(intfld,"x y z |E|^2 Ex.r Ex.i Ey.r Ey.i Ez.r Ez.i\n");
#ifdef PARALLEL
  }     /* end of if */
#endif
  /* saves fields to file */
  for (i=0;i<local_nvoid_Ndip;++i) {
    mat=material[i];
    /* field=P/(V*hi) */
    for (j=0;j<3;j++) cMult(hi_inv[mat],pvec[3*i+j],fld[j]);
    fprintf(intfld,
      "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
      DipoleCoord[3*i],DipoleCoord[3*i+1], DipoleCoord[3*i+2],
      cvNorm2(fld),
      fld[0][RE],fld[0][IM],fld[1][RE],fld[1][IM],fld[2][RE],fld[2][IM]);
  }
  FCloseErr(intfld,fname,ALL_POS);
#ifdef PARALLEL
  /* wait for all processes to save their part of geometry */
  Synchronize();
  if (ringid==ROOT) CatNFiles(directory,F_INTFLD_TMP,intfld_fname);
#endif
  PRINTZ("Internal fields saved to file\n");
  Timing_FileIO += clock() - tstart;
}

/*============================================================*/

int CalculateE(const char which,const int type)
  /* Calculate everything for x or y polarized incident light; or one and use symmetry
     to determine the rest (determined by type) */
{
  int exit_status;

  tstart_CE=clock();
  /* calculate the incident field Einc; vector b=Einc*cc_sqrt */
  D("Generating B");
  GenerateB (which, Einc);
  /* calculate solution vector x */
  D("Iterative solver started");
  exit_status=IterativeSolver(IterMethod);
  Timing_IntFieldOne = clock() - tstart_CE;
  Timing_IntField += Timing_IntFieldOne;
  /* return if checkpoint (normal) occured */
  if (exit_status==CHP_EXIT) return CHP_EXIT;

  if (yzplane) CalcEplane(which,type);  /*generally plane of incPolY and prop*/
  /* Calculate the scattered field for the whole solid-angle */
  if (all_dir) CalcAlldir();
  /* Calculate the scattered field on the given grid of angles */
  if (scat_grid) CalcScatGrid(which);
  /* Calculate integral scattering quantities
     (crosssections, asymmetry parameter, electric forces) */
  if (calc_Cext || calc_Cabs || calc_Csca || calc_asym || calc_mat_force)
    CalcIntegralScatQuantities(which);
  /* saves internal fields to text file */
  if (store_int_field) StoreIntFields(which);
  return 0;
}

