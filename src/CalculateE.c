 /* File: CalculateE.c
 * $Date::                            $
 * Descr: the module to calculate the E field and all scattering quantities
 *
 *        Routines for most scattering quantities are in crosssec.c. Also saves internal fields to
 *        file (optional).
 *
 * Copyright (C) 2006-2013 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#include "const.h" // keep this first
// project headers
#include "cmplx.h"
#include "comm.h"
#include "crosssec.h"
#include "debug.h"
#include "function.h"
#include "io.h"
#include "memory.h"
#include "linalg.h" // for nMult_mat
#include "Romberg.h"
#include "timing.h"
#include "vars.h"
// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>

// SEMI-GLOBAL VARIABLES

// defined and initialized in calculator.c
extern double * restrict muel_phi,* restrict muel_phi_buf;
extern doublecomplex * restrict EplaneX, * restrict EplaneY;
extern double * restrict Eplane_buffer;
extern const double dtheta_deg,dtheta_rad;
extern doublecomplex * restrict ampl_alphaX,* restrict ampl_alphaY;
extern double * restrict muel_alpha;
// defined and initialized in crosssec.c
extern const Parms_1D phi_sg;
// defined and initialized in param.c
extern const bool store_int_field,store_dip_pol,store_beam,store_scat_grid,calc_Cext,calc_Cabs,
calc_Csca,calc_vec,calc_asym,calc_mat_force,store_force,store_mueller,store_ampl;
extern const int phi_int_type;
// defined and initialized in timing.c
extern TIME_TYPE Timing_EPlane,Timing_EPlaneComm,
Timing_IntField,Timing_IntFieldOne,Timing_ScatQuan;
extern size_t TotalEFieldPlane;

// used in iterative.c
TIME_TYPE tstart_CE;

// LOCAL VARIABLES

#define MUEL_HEADER "s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44"
#define AMPL_HEADER "s1.r s1.i s2.r s2.i s3.r s3.i s4.r s4.i"
#define THETA_HEADER "theta"
#define PHI_HEADER "phi"
#define RMSE_HEADER "RMSE(integr)"
#define MUEL_FORMAT EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "\
	EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM
#define AMPL_FORMAT EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM" "EFORM
#define ANGLE_FORMAT "%.2f"
#define RMSE_FORMAT "%.3E"

// EXTERNAL FUNCTIONS

// GenerateB.c
void GenerateB(enum incpol which,doublecomplex *x);
// iterative.c
int IterativeSolver(enum iter method);

//============================================================

static void ComputeMuellerMatrix(double matrix[4][4], const doublecomplex s1,const doublecomplex s2,
	const doublecomplex s3,const doublecomplex s4)
/* computer mueller matrix from scattering matrix elements s1, s2, s3, s4, according to formula
 * 3.16 from Bohren and Huffman
 */
{
	matrix[0][0] = 0.5*(cMultConRe(s1,s1)+cMultConRe(s2,s2)+cMultConRe(s3,s3)+cMultConRe(s4,s4));
	matrix[0][1] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+cMultConRe(s4,s4)-cMultConRe(s3,s3));
	matrix[0][2] = cMultConRe(s2,s3)+cMultConRe(s1,s4);
	matrix[0][3] = cMultConIm(s2,s3)-cMultConIm(s1,s4);

	matrix[1][0] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+cMultConRe(s3,s3)-cMultConRe(s4,s4));
	matrix[1][1] = 0.5*(cMultConRe(s2,s2)+cMultConRe(s1,s1)-cMultConRe(s3,s3)-cMultConRe(s4,s4));
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

//============================================================

INLINE void InitMuellerIntegrFile(const int type,const char * restrict fname,FILE * restrict * file,
	char * restrict buf,const size_t buf_size,double * restrict *mult)
/* If 'phi_int_type' matches 'type', appropriate file (name given by 'fname') is created (with
 * handle '*file'), and heading line is put into it. String buffer 'buf' is used. Vector of
 * multipliers '*mult' is allocated if its pointer is specified.
 */
{
	if (phi_int_type & type) {
		SnprintfErr(ONE_POS,buf,buf_size,"%s/%s",directory,fname);
		(*file)=FOpenErr(buf,"w",ONE_POS);
		fprintf(*file,THETA_HEADER" "MUEL_HEADER" "RMSE_HEADER"\n");
		if (mult!=NULL) MALLOC_VECTOR(*mult,double,angles.phi.N,ALL);
	}
}

//==============================================================

INLINE void PrintToIntegrFile(const int type,FILE * restrict file,double *maxerr,
	const double * restrict muel,double *  restrict muel_buf,const double * restrict mult,
	double matrix[4][4],const double theta)
/* If 'phi_int_type' matches 'type', array 'muel' is integrated over phi (possibly using multiplier
 * 'mult' and buffer 'muel_buf') and saved to 'file' together with 'theta'. Maximum error '*maxerr'
 * is updated, 'matrix' buffer is used.
 */
{
	int k;
	size_t j;
	double err;

	if (phi_int_type & type) {
		if (mult==NULL) err=Romberg1D(phi_sg,16,muel,matrix[0]);
		else {
			for (j=0;j<angles.phi.N;j++) for(k=0;k<16;k++) muel_buf[16*j+k]=muel[16*j+k]*mult[j];
			err=Romberg1D(phi_sg,16,muel_buf,matrix[0]);
		}
		if (err>*maxerr) *maxerr=err;
		fprintf(file,ANGLE_FORMAT" "MUEL_FORMAT" "RMSE_FORMAT"\n",theta,matrix[0][0],
			matrix[0][1],matrix[0][2],matrix[0][3],matrix[1][0],matrix[1][1],matrix[1][2],
			matrix[1][3],matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],matrix[3][0],
			matrix[3][1],matrix[3][2],matrix[3][3],err);
	}
}

//==============================================================

INLINE void CloseIntegrFile(const int type,FILE * restrict file,const char * restrict fname,
	double * restrict mult)
/* If 'phi_int_type' matches 'type', appropriate 'file' (named 'fname') is closed and array 'mult'
 * is freed.
 */
{
	if (phi_int_type & type) {
		FCloseErr(file,fname,ONE_POS);
		Free_general(mult);
	}
}
//==============================================================

void MuellerMatrix(void)
{
	// redundant initializations in the following lines are against warnings
	FILE * restrict mueller,* restrict ampl,* restrict mueller_int=NULL,
		* restrict mueller_int_c2=NULL,* restrict mueller_int_s2=NULL,
		* restrict mueller_int_c4=NULL,* restrict mueller_int_s4=NULL;
	double * restrict cos2=NULL,* restrict sin2=NULL,* restrict cos4=NULL,* restrict sin4=NULL;
	double matrix[4][4];
	double theta,phi,ph,max_err,max_err_c2,max_err_s2,max_err_c4,max_err_s4;
	doublecomplex s1,s2,s3,s4,s10,s20,s30,s40;
	char fname[MAX_FNAME];
	int i;
	size_t index,index1,k_or,j,n,ind;
	double co,si;
	double alph;
	TIME_TYPE tstart;

	// redundant initialization to remove warnings
	mueller=ampl=NULL;
	co=si=0;

	// Everything is done on ROOT only
	if (!IFROOT) return;

	if (orient_avg) { // Amplitude matrix (ampl_alplha) => Mueller matrix (muel_alpha)
		/* amplitude matrix is not integrated (hence not used here). We do check store_mueller
		 * because orient_avg may have sense (though very little) without it (e.g. to compute only
		 * averaged cross sections).
		 */
		if (store_mueller) {
			index1=index=0;
			for (k_or=0;k_or<alpha_int.N;k_or++) {
				alph=Deg2Rad(alpha_int.val[k_or]); // read current alpha
				co=cos(alph);
				si=sin(alph);
				for (i=0;i<nTheta;i++) {
					// read amplitude matrix from memory
					cEqual(ampl_alphaX[index],s10);
					cEqual(ampl_alphaX[index+1],s30);
					cEqual(ampl_alphaY[index],s40);
					cEqual(ampl_alphaY[index+1],s20);
					// transform it, multiplying by rotation matrix (-alpha)
					cLinComb(s20,s30,co,si,s2);  // s2 =  co*s20 + si*s30
					cLinComb(s20,s30,-si,co,s3); // s3 = -si*s20 + co*s30
					cLinComb(s40,s10,co,si,s4);  // s4 =  co*s40 + si*s10
					cLinComb(s40,s10,-si,co,s1); // s1 = -si*s40 + co*s10

					ComputeMuellerMatrix((double (*)[4])(muel_alpha+index1),s1,s2,s3,s4);
					index+=2;
					index1+=16;
				}
			}
		}
	}
	else {
		tstart=GET_TIME(); // here Mueller matrix is saved to file
		if (yzplane) {
			if (store_ampl) {
				SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_AMPL,directory);
				ampl=FOpenErr(fname,"w",ONE_POS);
				fprintf(ampl,THETA_HEADER" "AMPL_HEADER"\n");
				for (i=0;i<nTheta;i++) {
					theta=i*dtheta_deg;
					fprintf(ampl,ANGLE_FORMAT" "AMPL_FORMAT"\n",theta,EplaneX[2*i][RE],
						EplaneX[2*i][IM],EplaneY[2*i+1][RE],EplaneY[2*i+1][IM],EplaneX[2*i+1][RE],
						EplaneX[2*i+1][IM],EplaneY[2*i][RE],EplaneY[2*i][IM]);
				}
				FCloseErr(ampl,F_AMPL,ONE_POS);
			}
			if (store_mueller) {
				SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_MUEL,directory);
				mueller=FOpenErr(fname,"w",ONE_POS);
				fprintf(mueller,THETA_HEADER" "MUEL_HEADER"\n");
				for (i=0;i<nTheta;i++) {
					theta=i*dtheta_deg;
					ComputeMuellerMatrix(matrix,EplaneX[2*i],EplaneY[2*i+1],EplaneX[2*i+1],
						EplaneY[2*i]);
					fprintf(mueller,ANGLE_FORMAT" "MUEL_FORMAT"\n",theta,matrix[0][0],matrix[0][1],
						matrix[0][2],matrix[0][3],matrix[1][0],matrix[1][1],matrix[1][2],
						matrix[1][3],matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
						matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3]);
				}
				FCloseErr(mueller,F_MUEL,ONE_POS);
			}
		}
		if (scat_grid) {
			/* compute Mueller Matrix in full space angle.
			 * E-fields are stored in arrays EgridX and EgridY for incoming X and Y polarized light.
			 * It is converted to the scattering matrix elements (see e.g Bohren and Huffman) :
			 * s2 = cos(phi)E'X'par + sin(phi)E'Y'par
			 * s3 = sin(phi)E'X'par - cos(phi)E'Y'par
			 * s4 = cos(phi)E'X'per + sin(phi)E'Y'per
			 * s1 = sin(phi)E'X'per - cos(phi)E'Y'per
			 * from these the mueller matrix elements are computed
			 */
			// open files for writing
			if (store_scat_grid) {
				if (store_ampl) {
					SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_AMPL_SG,directory);
					ampl=FOpenErr(fname,"w",ONE_POS);
					fprintf(ampl,THETA_HEADER" "PHI_HEADER" "AMPL_HEADER"\n");
				}
				if (store_mueller) {
					SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_MUEL_SG,directory);
					mueller=FOpenErr(fname,"w",ONE_POS);
					fprintf(mueller,THETA_HEADER" "PHI_HEADER" "MUEL_HEADER"\n");
				}
			}
			if (phi_integr) { // also initializes arrays of multipliers
				/* amplitude matrix is not integrated. Moreover, if mueller matrix is not needed,
				 * then phi_integr has no action at all. So we assume that phi_integr implies
				 * store_mueller (to be checked in param.c).
				 */
				InitMuellerIntegrFile(PHI_UNITY,F_MUEL_INT,&mueller_int,fname,MAX_FNAME,NULL);
				InitMuellerIntegrFile(PHI_COS2,F_MUEL_C2,&mueller_int_c2,fname,MAX_FNAME,&cos2);
				InitMuellerIntegrFile(PHI_SIN2,F_MUEL_S2,&mueller_int_s2,fname,MAX_FNAME,&sin2);
				InitMuellerIntegrFile(PHI_COS4,F_MUEL_C4,&mueller_int_c4,fname,MAX_FNAME,&cos4);
				InitMuellerIntegrFile(PHI_SIN4,F_MUEL_S4,&mueller_int_s4,fname,MAX_FNAME,&sin4);
				// fills arrays with multipliers (optimized)
				for (j=0;j<angles.phi.N;j++) {
					// prepare
					ph=2*Deg2Rad(angles.phi.val[j]);
					if (phi_int_type & (PHI_COS2|PHI_COS4|PHI_SIN4)) co=cos(ph);
					if (phi_int_type & (PHI_SIN2|PHI_SIN4)) si=sin(ph);
					// fill
					if (phi_int_type & PHI_COS2) cos2[j]=co;
					if (phi_int_type & PHI_SIN2) sin2[j]=si;
					if (phi_int_type & PHI_COS4) cos4[j]=2*co*co-1;
					if (phi_int_type & PHI_SIN4) sin4[j]=2*si*co;
				}
			}
			// set type of cycling through angles
			if (angles.type==SG_GRID) n=angles.phi.N;
			else n=1; // angles.type==SG_PAIRS
			// main cycle
			index=0;
			max_err=max_err_c2=max_err_s2=max_err_c4=max_err_s4=0;
			for (ind=0;ind<angles.theta.N;++ind) {
				index1=0;
				theta=angles.theta.val[ind];
				for (j=0;j<n;++j) {
					if (angles.type==SG_GRID) phi=angles.phi.val[j];
					else phi=angles.phi.val[ind]; // angles.type==SG_PAIRS
					ph=Deg2Rad(phi);
					co=cos(ph);
					si=sin(ph);
					// read amplitude matrix from memory
					cEqual(EgridY[index],s10);
					cEqual(EgridY[index+1],s30);
					cEqual(EgridX[index],s40);
					cEqual(EgridX[index+1],s20);
					index+=2;
					// transform it, multiplying by rotation matrix from per-par to X-Y
					cLinComb(s20,s30,co,si,s2);  // s2 =  co*s20 + si*s30
					cLinComb(s20,s30,si,-co,s3); // s3 =  si*s20 - co*s30
					cLinComb(s40,s10,co,si,s4);  // s4 =  co*s40 + si*s10
					cLinComb(s40,s10,si,-co,s1); // s1 =  si*s40 - co*s10

					if (store_mueller) ComputeMuellerMatrix(matrix,s1,s2,s3,s4);

					if (phi_integr) {
						memcpy(muel_phi+index1,matrix[0],16*sizeof(double));
						index1+=16;
					}
					if (store_scat_grid) {
						if (store_mueller)
							fprintf(mueller,ANGLE_FORMAT" "ANGLE_FORMAT" "MUEL_FORMAT"\n",
								theta,phi,matrix[0][0],matrix[0][1],matrix[0][2],matrix[0][3],
								matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
								matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],
								matrix[3][0],matrix[3][1],matrix[3][2],matrix[3][3]);
						if (store_ampl)
							fprintf(ampl,ANGLE_FORMAT" "ANGLE_FORMAT" "AMPL_FORMAT"\n",
								theta,phi,s1[RE],s1[IM],s2[RE],s2[IM],s3[RE],s3[IM],s4[RE],s4[IM]);
					}
				}
				if (phi_integr) {
					PrintToIntegrFile(PHI_UNITY,mueller_int,&max_err,muel_phi,NULL,
					NULL,matrix,theta);
				PrintToIntegrFile(PHI_COS2,mueller_int_c2,&max_err_c2,muel_phi,muel_phi_buf,
					cos2,matrix,theta);
				PrintToIntegrFile(PHI_SIN2,mueller_int_s2,&max_err_s2,muel_phi,muel_phi_buf,
					sin2,matrix,theta);
				PrintToIntegrFile(PHI_COS4,mueller_int_c4,&max_err_c4,muel_phi,muel_phi_buf,
					cos4,matrix,theta);
				PrintToIntegrFile(PHI_SIN4,mueller_int_s4,&max_err_s4,muel_phi,muel_phi_buf,
					sin4,matrix,theta);
				}
			}
			if (phi_integr) {
				fprintf(logfile,"\nMaximum relative mean-square error of Mueller integration:\n");
				if (phi_int_type & PHI_UNITY)
					fprintf(logfile,"  1          -> "RMSE_FORMAT"\n",max_err);
				if (phi_int_type & PHI_COS2)
					fprintf(logfile,"  cos(2*phi) -> "RMSE_FORMAT"\n",max_err_c2);
				if (phi_int_type & PHI_SIN2)
					fprintf(logfile,"  sin(2*phi) -> "RMSE_FORMAT"\n",max_err_s2);
				if (phi_int_type & PHI_COS4)
					fprintf(logfile,"  cos(4*phi) -> "RMSE_FORMAT"\n",max_err_c4);
				if (phi_int_type & PHI_SIN4)
					fprintf(logfile,"  sin(4*phi) -> "RMSE_FORMAT"\n",max_err_s4);
			}
			// close files; free arrays
			if (store_scat_grid) {
				if (store_mueller) FCloseErr(mueller,F_MUEL_SG,ONE_POS);
				if (store_ampl) FCloseErr(ampl,F_AMPL_SG,ONE_POS);
			}
			if (phi_integr) {
				CloseIntegrFile(PHI_UNITY,mueller_int,F_MUEL_INT,NULL);
				CloseIntegrFile(PHI_COS2,mueller_int_c2,F_MUEL_C2,cos2);
				CloseIntegrFile(PHI_SIN2,mueller_int_s2,F_MUEL_S2,sin2);
				CloseIntegrFile(PHI_COS4,mueller_int_c4,F_MUEL_C4,cos4);
				CloseIntegrFile(PHI_SIN4,mueller_int_s4,F_MUEL_S4,sin4);
			}
		}
		Timing_FileIO += GET_TIME() - tstart;
	}
}

//============================================================

static void CalcEplane(const enum incpol which,const enum Eftype type)
// calculates scattered electric field in a plane
{
	double *incPol,*incPolper,*incPolpar;
	// where to store calculated field for one plane (actually points to different other arrays)
	doublecomplex *Eplane;
	int i;
	doublecomplex ebuff[3]; // small vector to hold E fields
	double robserver[3];    // small vector for observer in E calculation
	double epar[3];         // unit vector in direction of Epar
	double theta;           // scattering angle
	double co,si;           // temporary, cos and sin of some angle
	double incPol_tmp1[3],incPol_tmp2[3]; // just allocated memory for incPolper, incPolpar
	double alph;
	TIME_TYPE tstart;
	size_t k_or;
	int orient,Norient;
	enum incpol choice;

	incPolper=incPol_tmp1; // initialization of per and par polarizations
	incPolpar=incPol_tmp2;

	if (type==CE_NORMAL) Norient=1; // initialize # orientations
	else Norient=2;                 // type==CE_PARPER

	for (k_or=0;k_or<alpha_int.N;k_or++) {
		// cycle over alpha - for orientation averaging
		if (orient_avg) {
			alph=Deg2Rad(alpha_int.val[k_or]); // rotate polarization basis vectors by -alpha
			co=cos(alph);
			si=sin(alph);
			LinComb(incPolX,incPolY,co,-si,incPolper); // incPolper = co*incPolX - si*incPolY;
			LinComb(incPolX,incPolY,si,co,incPolpar);  // incPolpar = si*incPolX + co*incPolY;
		}
		else {
			memcpy(incPolper,incPolX,3*sizeof(double)); // per <=> X
			memcpy(incPolpar,incPolY,3*sizeof(double)); // par <=> Y
		}

		for(orient=0;orient<Norient;orient++) {
			// in case of Rotation symmetry
			tstart = GET_TIME ();
			if (orient==0) choice=which;
			else { // orient==1
				/* Rotation symmetry: calculate per-per from current data. CalculateE is called
				 * from calculator with Y polarization - we now just assume that we have
				 * the x-z plane as the scattering plane, rotating in the negative x-direction.
				 * This mimics the real case of X polarization with the y-z plane as scattering plane
				 * Then IncPolY -> -IncPolX; incPolX -> IncPolY
				 */
				if (which==INCPOL_Y) choice=INCPOL_X;
				else choice=INCPOL_Y; // which==INCPOL_X
				incPol=incPolper;
				incPolper=incPolpar;
				incPolpar=incPol;
				vInvSign(incPolpar);
			}
			// initialize Eplane
			if (orient_avg) {
				if (choice==INCPOL_Y) Eplane=ampl_alphaY + 2*nTheta*k_or;
				else Eplane=ampl_alphaX + 2*nTheta*k_or; // choice==INCPOL_X
			}
			else {
				if (choice==INCPOL_Y) Eplane=EplaneY;
				else Eplane=EplaneX; // choice==INCPOL_X
			}

			for (i=0;i<nTheta;i++) {
				theta = i * dtheta_rad;
				co=cos(theta);
				si=sin(theta);
				LinComb(prop,incPolpar,co,si,robserver); // robserver = co*prop + si*incPolpar;

				CalcField(ebuff,robserver);
				// convert to (l,r) frame
				crDotProd(ebuff,incPolper,Eplane[2*i]); // Eper[i]=Esca.incPolper
				LinComb(prop,incPolpar,-si,co,epar);    // epar=-si*prop+co*incPolpar
				crDotProd(ebuff,epar,Eplane[2*i+1]);    // Epar[i]=Esca.epar
			} //  end for i

			// Accumulate Eplane to root and sum
			D("Accumulating Eplane started");
			// accumulate only on processor 0 !, done in one operation
			Accumulate((double *)Eplane,4*nTheta,Eplane_buffer,&Timing_EPlaneComm);
			D("Accumulating Eplane finished");

			Timing_EPlane = GET_TIME() - tstart;
			Timing_EField += Timing_EPlane;
			TotalEFieldPlane++;
		} // end of orient loop
	} // end of alpha loop
}

//============================================================

static void StoreFields(const enum incpol which,doublecomplex * restrict cmplxF,
	double * restrict realF,const char * restrict fname_preffix,const char * restrict tmpl UOIP,
	const char * restrict field_name,const char * restrict fullname)
/* Write any fields on each dipole to file (internal fields, incident beam, polarization, etc.).
 * Accepts both complex and real fields. Processes the one, which is not null.
 *
 * All processors should write the 'field' to temporary file. These files are named by template
 * 'tmpl' and afterwards are concatenated into the file, which name is build by adding a small
 * suffix to 'fname_preffix'. If CE_PARPER is employed then naturally saves only once; use '-sym no'
 * if needed. 'field_name' is used to build column labels (i.e. there is difference in the first row
 * between different fields). 'fullname' is for standard output.
 */
{
	FILE * restrict file; // file to store the fields
	size_t j;
	TIME_TYPE tstart;
	char fname[MAX_FNAME],fname_sh[MAX_FNAME_SH];
	bool cmplx_mode; // whether complex (true) or real (false) field is processed

	tstart=GET_TIME();
	// choose operational mode
	if ((cmplxF==NULL) ^ (realF==NULL)) cmplx_mode=(realF==NULL);
	else LogError(ONE_POS,"One field (either real or complex) must be given to StoreFields()");
	// build file name (without directory)
	strcpy(fname_sh,fname_preffix);
	if (which==INCPOL_Y) strcat(fname_sh,F_YSUF);
	else strcat(fname_sh,F_XSUF); // which==INCPOL_X
	// choose filename for direct saving
#ifdef PARALLEL
	size_t shift=SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/",directory);
	/* the following will cause warning by GCC if -Wformat-nonliteral (or -Wformat=2) is used,
	 * but we do not know any other convenient way to make this function work for different file
	 * name templates.
	 */
	SnprintfShiftErr(ALL_POS,shift,fname,MAX_FNAME,tmpl,ringid);
#else
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/%s",directory,fname_sh);
#endif
	file=FOpenErr(fname,"w",ALL_POS);
	// print head of file
#ifdef PARALLEL
	if (ringid==0) { // this condition can be different from being root
#endif
		if (cmplx_mode) fprintf(file,"x y z |%s|^2 %sx.r %sx.i %sy.r %sy.i %sz.r %sz.i\n",
			field_name,field_name,field_name,field_name,field_name,field_name,field_name);
		else fprintf(file,"x y z |%s|^2 %sx %sy %sz\n",field_name,field_name,field_name,field_name);
#ifdef PARALLEL
	} // end of if
#endif
	// saves fields to file
	if (cmplx_mode) for (j=0;j<local_nRows;j+=3) fprintf(file,GFORM10L"\n",
		DipoleCoord[j],DipoleCoord[j+1],DipoleCoord[j+2],cvNorm2(cmplxF+j),cmplxF[j][RE],
		cmplxF[j][IM],cmplxF[j+1][RE],cmplxF[j+1][IM],cmplxF[j+2][RE],cmplxF[j+2][IM]);
	else for (j=0;j<local_nRows;j+=3) fprintf(file,GFORM7L"\n",DipoleCoord[j],DipoleCoord[j+1],
		DipoleCoord[j+2],DotProd(realF+j,realF+j),realF[j],realF[j+1],realF[j+2]);
	FCloseErr(file,fname,ALL_POS);
#ifdef PARALLEL
	// wait for all processes to save their part of geometry
	Synchronize();
	if (IFROOT) CatNFiles(directory,tmpl,fname_sh);
#endif
	if (IFROOT) printf("%s saved to file\n",fullname);
	Timing_FileIO += GET_TIME() - tstart;
}

//============================================================

INLINE void ParticleToBeamRF(const double vp[static restrict 3],double vb[static restrict 3])
// transform real vector vp from particle to beam reference frame, vb and vp must not alias!!!
{
	vb[0]=DotProd(vp,incPolX);
	vb[1]=DotProd(vp,incPolY);
	vb[2]=DotProd(vp,prop);
}

//============================================================

static void CalcIntegralScatQuantities(const enum incpol which)
/* calculates all the scattering cross sections, normalized and unnormalized asymmetry parameter,
 * and force on the particle and each dipole. Cext and Cabs are averaged over orientation,
 * if needed.
 */
{
	// Scattering force, extinction force and radiation pressure per dipole
	double * restrict Frp;
	double Cext,Cabs,Csca,   // Cross sections
	dummy[3],                // asymmetry parameter*Csca
	Finc_tot[3],Fsca_tot[3], // total extinction and scattering force in particle reference frame
	Finc_brf[3],Fsca_brf[3],Frp_brf[3], // extinction, scattering, and sum in beam reference frame
	Cnorm,            // normalizing factor from force to cross section
	Qnorm;            // normalizing factor from force to efficiency
	FILE * restrict CCfile;
	TIME_TYPE tstart;
	char fname_cs[MAX_FNAME];
	const double *incPol;
	const char *f_suf;

	// redundant initialization to remove warnings
	Cext=Cabs=Csca=0;
	CCfile=NULL;

	D("Calculation of cross sections started");
	tstart = GET_TIME();

	if (which == INCPOL_Y) {
		f_suf=F_YSUF;
		incPol=incPolY;
	}
	else { // which == INCPOL_X
		f_suf=F_XSUF;
		incPol=incPolX;
	}
	/* order of calculations is important, when using SQ_FINDIP. Then first dCabs should be
	 * calculated, which is further used to correct Cext
	 */
	if (calc_Cabs) Cabs = AbsCross();
	if (calc_Cext) Cext = ExtCross(incPol);
	D("Cext and Cabs calculated");
	if (orient_avg) {
		if (IFROOT) {
			if (which == INCPOL_Y) { // assumed that first call of CalculateE is with INCPOL_Y flag
				muel_alpha[-2]=Cext;
				muel_alpha[-1]=Cabs;
			}
			else { // which == INCPOL_X
				muel_alpha[-2]=(muel_alpha[-2]+Cext)/2;
				muel_alpha[-1]=(muel_alpha[-1]+Cabs)/2;
			}
		}
	}
	else { // not orient_avg
		if (IFROOT) {
			SnprintfErr(ONE_POS,fname_cs,MAX_FNAME,"%s/"F_CS"%s",directory,f_suf);
			CCfile=FOpenErr(fname_cs,"w",ONE_POS);
			if (calc_Cext) PrintBoth(CCfile,"Cext\t= "GFORM"\nQext\t= "GFORM"\n",
				Cext,Cext*inv_G);
			if (calc_Cabs) PrintBoth(CCfile,"Cabs\t= "GFORM"\nQabs\t= "GFORM"\n",
				Cabs,Cabs*inv_G);
			if (all_dir) fprintf(CCfile,"\nIntegration\n\n");
			if (calc_Csca) {
				Csca=ScaCross(f_suf);
				PrintBoth(CCfile,"Csca\t= "GFORM"\nQsca\t= "GFORM"\n",Csca,Csca*inv_G);
			}
			if (calc_vec) {
				AsymParm_x(dummy,f_suf);
				AsymParm_y(dummy+1,f_suf);
				AsymParm_z(dummy+2,f_suf);
				PrintBoth(CCfile,"Csca.g\t= "GFORM3V"\n",dummy[0],dummy[1],dummy[2]);
				if (calc_asym) PrintBoth(CCfile,"g\t= "GFORM3V"\n",
					dummy[0]/Csca,dummy[1]/Csca,dummy[2]/Csca);
			}
		} // end of root
		if (calc_mat_force) {
			if (IFROOT) printf("Calculating the force per dipole\n");
			if (store_force) MALLOC_VECTOR(Frp,double,local_nRows,ALL);
			else Frp=NULL;
			Frp_mat(Finc_tot,Fsca_tot,Frp);
			// Write Cross-Sections and Efficiencies to file
			/* This output contains a number of redundant quantities, like Cext and Csca.g. The main
			 * purpose of this is to use it as independent test (it can be quickly compared against
			 * the same quantities calculated by different means).
			 * It also partly addresses issue 133.
			 */
			if (IFROOT) {
				Cnorm = EIGHT_PI;
				Qnorm = EIGHT_PI*inv_G;
				ParticleToBeamRF(Finc_tot,Finc_brf);
				ParticleToBeamRF(Fsca_tot,Fsca_brf);
				vAdd(Finc_brf,Fsca_brf,Frp_brf);
				PrintBoth(CCfile,"\nBased on radiation forces:\n"
				                 "Cext.a\t= "GFORM3V"\n"
				                 "Csca.g\t= "GFORM3V"\n"
				                 "Cpr\t= "GFORM3V"\n"
				                 "Qpr\t= "GFORM3V"\n",
				                 Cnorm*Finc_brf[0],Cnorm*Finc_brf[1],Cnorm*Finc_brf[2],
				                 -Cnorm*Fsca_brf[0],-Cnorm*Fsca_brf[1],-Cnorm*Fsca_brf[2],
				                 Cnorm*Frp_brf[0],Cnorm*Frp_brf[1],Cnorm*Frp_brf[2],
				                 Qnorm*Frp_brf[0],Qnorm*Frp_brf[1],Qnorm*Frp_brf[2]);
			}
			if (store_force) {
				StoreFields(which,NULL,Frp,F_FRP,F_FRP_TMP,"F","Radiation forces");
				Free_general(Frp);
			}
		}
		if (IFROOT) FCloseErr(CCfile,fname_cs,ONE_POS);
	}
	D("Calculation of cross sections finished");
	Timing_ScatQuan += GET_TIME() - tstart;
}

//============================================================

static void StoreIntFields(const enum incpol which)
// Write actual internal fields (not exciting) on each dipole to file
{
	// calculate fields; e_field=P/(V*chi)=chi_inv*P; for anisotropic - by components
	nMult_mat(xvec,pvec,chi_inv);
	// save fields to file
	StoreFields(which,xvec,NULL,F_INTFLD,F_INTFLD_TMP,"E","Internal fields");
}

//============================================================

int CalculateE(const enum incpol which,const enum Eftype type)
/* Calculate everything for x or y polarized incident light; or one and use symmetry to determine
 * the rest (determined by type)
 */
{
	int exit_status;

	tstart_CE=GET_TIME();
	// calculate the incident field Einc; vector b=Einc*cc_sqrt
	D("Generating B");
	GenerateB (which, Einc);
	if (store_beam) StoreFields(which,Einc,NULL,F_BEAM,F_BEAM_TMP,"Einc","Incident beam");
	// calculate solution vector x
	D("Iterative solver started");
	exit_status=IterativeSolver(IterMethod);
	D("Iterative solver finished");
	Timing_IntFieldOne = GET_TIME() - tstart_CE;
	Timing_IntField += Timing_IntFieldOne;
	// return if checkpoint (normal) occurred
	if (exit_status==CHP_EXIT) return CHP_EXIT;

	if (yzplane) CalcEplane(which,type); //generally plane of incPolY and prop
	// Calculate the scattered field for the whole solid-angle
	if (all_dir) CalcAlldir();
	// Calculate the scattered field on the given grid of angles
	if (scat_grid) CalcScatGrid(which);
	/* Calculate integral scattering quantities (cross sections, asymmetry parameter,
	 * electric forces)
	 */
	if (calc_Cext || calc_Cabs || calc_Csca || calc_asym || calc_mat_force)
		CalcIntegralScatQuantities(which);
	// saves internal fields and/or dipole polarizations to text file
	if (store_int_field) StoreIntFields(which);
	if (store_dip_pol)
		StoreFields(which,pvec,NULL,F_DIPPOL,F_DIPPOL_TMP,"P","Dipole polarizations");
	return 0;
}

//============================================================

void SaveMuellerAndCS(double * restrict in)
/* saves Mueller matrix and cross sections (averaged) to files; vector in contains values of cross
 * sections and then array of Mueller matrix elements; designed to be called from ROOT only
 */
{
	FILE * restrict mueller,* restrict CCfile;
	char fname[MAX_FNAME];
	int i;
	double Cext,Cabs,*muel;
	TIME_TYPE tstart;

	tstart=GET_TIME();
	// distribute input values
	Cext=in[0];
	Cabs=in[1];
	muel=in+2;

	if (store_mueller) { // save Mueller matrix
		SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_MUEL,directory);
		mueller=FOpenErr(fname,"w",ONE_POS);
		fprintf(mueller,THETA_HEADER" "MUEL_HEADER"\n");
		for (i=0;i<nTheta;i++) {
			fprintf(mueller,ANGLE_FORMAT" "MUEL_FORMAT"\n",
				i*dtheta_deg,muel[0],muel[1],muel[2],muel[3],muel[4],muel[5],muel[6],muel[7],
				muel[8],muel[9],muel[10],muel[11],muel[12],muel[13],muel[14],muel[15]);
			muel+=16;
		}
		FCloseErr(mueller,F_MUEL,ONE_POS);
	}
	// save cross sections
	SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_CS,directory);
	CCfile=FOpenErr(fname,"w",ONE_POS);
	PrintBoth(CCfile,"Cext\t= "GFORM"\nQext\t= "GFORM"\n",Cext,Cext*inv_G);
	PrintBoth(CCfile,"Cabs\t= "GFORM"\nQabs\t= "GFORM"\n",Cabs,Cabs*inv_G);
	FCloseErr(CCfile,F_CS,ONE_POS);

	Timing_FileIO += GET_TIME() - tstart;
}
