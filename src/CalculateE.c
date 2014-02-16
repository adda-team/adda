 /* File: CalculateE.c
 * $Date::                            $
 * Descr: the module to calculate the E field and all scattering quantities
 *
 *        Routines for most scattering quantities are in crosssec.c. Also saves internal fields to
 *        file (optional).
 *
 * Copyright (C) 2006-2014 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
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
extern doublecomplex * restrict EplaneX, * restrict EplaneY, * restrict EyzplX, * restrict EyzplY;
extern const double dtheta_deg,dtheta_rad;
extern doublecomplex * restrict ampl_alphaX,* restrict ampl_alphaY;
extern double * restrict muel_alpha;
// defined and initialized in crosssec.c
extern const Parms_1D phi_sg;
extern const double ezLab[3],exSP[3];
// defined and initialized in GenerateB.c
extern const double C0dipole,C0dipole_refl;
// defined and initialized in param.c
extern const bool store_int_field,store_dip_pol,store_beam,store_scat_grid,calc_Cext,calc_Cabs,
	calc_Csca,calc_vec,calc_asym,calc_mat_force,store_force,store_ampl;
extern const int phi_int_type;
// defined and initialized in timing.c
extern TIME_TYPE Timing_EPlane,Timing_EPlaneComm,Timing_IntField,Timing_IntFieldOne,Timing_ScatQuan,Timing_IncBeam;
extern size_t TotalEFieldPlane;

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
#define COMP44M(a) (a)[0][0],(a)[0][1],(a)[0][2],(a)[0][3],(a)[1][0],(a)[1][1],(a)[1][2],(a)[1][3],(a)[2][0],\
	(a)[2][1],(a)[2][2],(a)[2][3],(a)[3][0],(a)[3][1],(a)[3][2],(a)[3][3]

// EXTERNAL FUNCTIONS

// GenerateB.c
void GenerateB(enum incpol which,doublecomplex *x);
// iterative.c
int IterativeSolver(enum iter method,enum incpol which);

//======================================================================================================================

static void ComputeMuellerMatrix(double matrix[4][4], const doublecomplex s1,const doublecomplex s2,
	const doublecomplex s3,const doublecomplex s4,const double theta)
/* compute Mueller matrix from scattering matrix elements s1, s2, s3, s4, according to formula 3.16 from Bohren and
 * Huffman
 * In surface mode the result is additionally multiplied by factor Re(m_sca)/Re(m_inc) to account for corresponding
 * factor relating squared amplitude of the field and Stokes vector (irradiance). For that theta (in deg) is used.
 */
{
	matrix[0][0] = 0.5*(cAbs2(s1) + cAbs2(s2) + cAbs2(s3) + cAbs2(s4));
	matrix[0][1] = 0.5*(cAbs2(s2) - cAbs2(s1) + cAbs2(s4) - cAbs2(s3));
	matrix[0][2] = creal(s2*conj(s3) + s1*conj(s4));
	matrix[0][3] = cimag(s2*conj(s3) - s1*conj(s4));

	matrix[1][0] = 0.5*(cAbs2(s2) - cAbs2(s1) + cAbs2(s3) - cAbs2(s4));
	matrix[1][1] = 0.5*(cAbs2(s2) + cAbs2(s1) - cAbs2(s3) - cAbs2(s4));
	matrix[1][2] = creal(s2*conj(s3) - s1*conj(s4));
	matrix[1][3] = cimag(s2*conj(s3) + s1*conj(s4));

	matrix[2][0] = creal(s2*conj(s4) + s1*conj(s3));
	matrix[2][1] = creal(s2*conj(s4) - s1*conj(s3));
	matrix[2][2] = creal(s1*conj(s2) + s3*conj(s4));
	matrix[2][3] = cimag(s2*conj(s1) + s4*conj(s3));

	matrix[3][0] = cimag(s4*conj(s2) + s1*conj(s3));
	matrix[3][1] = cimag(s4*conj(s2) - s1*conj(s3));
	matrix[3][2] = cimag(s1*conj(s2) - s3*conj(s4));
	matrix[3][3] = creal(s1*conj(s2) - s3*conj(s4));

	if (surface) {
		double scale=inc_scale;
		if (TestBelowDeg(theta)) scale*=creal(msub);
		if (fabs(scale-1)>ROUND_ERR) for (int i=0;i<4;i++) for (int j=0;j<4;j++) matrix[i][j]*=scale;
	}
}

//======================================================================================================================

static void InitMuellerIntegrFile(const int type,const char * restrict fname,FILE * restrict * file,
	char * restrict buf,const size_t buf_size,double * restrict *mult)
/* If 'phi_int_type' matches 'type', appropriate file (name given by 'fname') is created (with handle '*file'), and
 * heading line is put into it. String buffer 'buf' is used. Vector of multipliers '*mult' is allocated if its pointer
 * is specified.
 */
{
	if (phi_int_type & type) {
		SnprintfErr(ONE_POS,buf,buf_size,"%s/%s",directory,fname);
		(*file)=FOpenErr(buf,"w",ONE_POS);
		fprintf(*file,THETA_HEADER" "MUEL_HEADER" "RMSE_HEADER"\n");
		if (mult!=NULL) MALLOC_VECTOR(*mult,double,angles.phi.N,ALL);
	}
}

//======================================================================================================================

static inline void PrintToIntegrFile(const int type,FILE * restrict file,double *maxerr,const double * restrict muel,
	double *  restrict muel_buf,const double * restrict mult,double matrix[4][4],const double theta)
/* If 'phi_int_type' matches 'type', array 'muel' is integrated over phi (possibly using multiplier 'mult' and buffer
 * 'muel_buf') and saved to 'file' together with 'theta'. Maximum error '*maxerr' is updated, 'matrix' buffer is used.
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
		fprintf(file,ANGLE_FORMAT" "MUEL_FORMAT" "RMSE_FORMAT"\n",theta,COMP44M(matrix),err);
	}
}

//======================================================================================================================

static void CloseIntegrFile(const int type,FILE * restrict file,const char * restrict fname,double * restrict mult)
// If 'phi_int_type' matches 'type', appropriate 'file' (named 'fname') is closed and array 'mult' is freed.
{
	if (phi_int_type & type) {
		FCloseErr(file,fname,ONE_POS);
		Free_general(mult);
	}
}
//======================================================================================================================

void MuellerMatrix(void)
{
	// redundant initializations in the following lines are against warnings
	FILE * restrict mueller,* restrict ampl,* restrict mueller_int=NULL,* restrict mueller_int_c2=NULL,
		* restrict mueller_int_s2=NULL,* restrict mueller_int_c4=NULL,* restrict mueller_int_s4=NULL;
	double * restrict cos2=NULL,* restrict sin2=NULL,* restrict cos4=NULL,* restrict sin4=NULL;
	double matrix[4][4];
	double theta,phi,ph,max_err,max_err_c2,max_err_s2,max_err_c4,max_err_s4;
	doublecomplex s1,s2,s3,s4;
	char fname[MAX_FNAME];
	int i;
	size_t index,index1,k_or,j,n,ind;
	double co,si,alph;
	TIME_TYPE tstart;

	// redundant initialization to remove warnings
	mueller=ampl=NULL;
	co=si=0;

	// Everything is done on ROOT only
	if (!IFROOT) return;

	if (orient_avg) { // Amplitude matrix (ampl_alplha) => Mueller matrix (muel_alpha)
		/* amplitude matrix is not integrated (hence not used here). We do check store_mueller because orient_avg may
		 * have sense (though very little) without it (e.g. to compute only averaged cross sections).
		 */
		if (store_mueller) {
			index1=index=0;
			for (k_or=0;k_or<alpha_int.N;k_or++) {
				alph=Deg2Rad(alpha_int.val[k_or]); // read current alpha
				co=cos(alph);
				si=sin(alph);
				for (i=0;i<nTheta;i++) {
					// transform amplitude matrix, multiplying by rotation matrix (-alpha)
					if (yzplane) { // here the default (alpha=0) is yz-plane, so par=Y, per=X
						s2 =  co*ampl_alphaY[index+1] + si*ampl_alphaX[index+1]; // s2 =  co*s20 + si*s30
						s3 = -si*ampl_alphaY[index+1] + co*ampl_alphaX[index+1]; // s3 = -si*s20 + co*s30
						s4 =  co*ampl_alphaY[index] + si*ampl_alphaX[index];     // s4 =  co*s40 + si*s10
						s1 = -si*ampl_alphaY[index] + co*ampl_alphaX[index];     // s1 = -si*s40 + co*s10
					}
					else { // scat_plane; here the default (alpha=0) is xz-plane, so par=X, per=-Y
						s2 =  co*ampl_alphaX[index+1] - si*ampl_alphaY[index+1]; // s2 =  co*s20 + si*s30
						s3 = -si*ampl_alphaX[index+1] - co*ampl_alphaY[index+1]; // s3 = -si*s20 + co*s30
						s4 =  co*ampl_alphaX[index] - si*ampl_alphaY[index];     // s4 =  co*s40 + si*s10
						s1 = -si*ampl_alphaX[index] - co*ampl_alphaY[index];     // s1 = -si*s40 + co*s10
					}
					theta=i*dtheta_deg;
					ComputeMuellerMatrix((double (*)[4])(muel_alpha+index1),s1,s2,s3,s4,theta);
					index+=2;
					index1+=16;
				}
			}
		}
	}
	else {
		tstart=GET_TIME(); // here Mueller matrix is saved to file
		if (yzplane) { // par=Y, per=X
			if (store_ampl) {
				SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_AMPL,directory);
				ampl=FOpenErr(fname,"w",ONE_POS);
				fprintf(ampl,THETA_HEADER" "AMPL_HEADER"\n");
				for (i=0;i<nTheta;i++) {
					theta=i*dtheta_deg;
					fprintf(ampl,ANGLE_FORMAT" "AMPL_FORMAT"\n",theta,REIM(EyzplX[2*i]),REIM(EyzplY[2*i+1]),
						REIM(EyzplX[2*i+1]),REIM(EyzplY[2*i]));
				}
				FCloseErr(ampl,F_AMPL,ONE_POS);
			}
			if (store_mueller) {
				SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_MUEL,directory);
				mueller=FOpenErr(fname,"w",ONE_POS);
				fprintf(mueller,THETA_HEADER" "MUEL_HEADER"\n");
				for (i=0;i<nTheta;i++) {
					theta=i*dtheta_deg;
					ComputeMuellerMatrix(matrix,EyzplX[2*i],EyzplY[2*i+1],EyzplX[2*i+1],EyzplY[2*i],theta);
					fprintf(mueller,ANGLE_FORMAT" "MUEL_FORMAT"\n",theta,COMP44M(matrix));
				}
				FCloseErr(mueller,F_MUEL,ONE_POS);
			}
		}
		if (scat_plane) { // par=X, per=-Y
			if (store_ampl) {
				SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_AMPL,directory);
				ampl=FOpenErr(fname,"w",ONE_POS);
				fprintf(ampl,THETA_HEADER" "AMPL_HEADER"\n");
				for (i=0;i<nTheta;i++) {
					theta=i*dtheta_deg;
					fprintf(ampl,ANGLE_FORMAT" "AMPL_FORMAT"\n",theta,REIM(-EplaneY[2*i]),REIM(EplaneX[2*i+1]),
						REIM(-EplaneY[2*i+1]),REIM(EplaneX[2*i]));
				}
				FCloseErr(ampl,F_AMPL,ONE_POS);
			}
			if (store_mueller) {
				SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_MUEL,directory);
				mueller=FOpenErr(fname,"w",ONE_POS);
				fprintf(mueller,THETA_HEADER" "MUEL_HEADER"\n");
				for (i=0;i<nTheta;i++) {
					theta=i*dtheta_deg;
					ComputeMuellerMatrix(matrix,-EplaneY[2*i],EplaneX[2*i+1],-EplaneY[2*i+1],EplaneX[2*i],theta);
					fprintf(mueller,ANGLE_FORMAT" "MUEL_FORMAT"\n",theta,COMP44M(matrix));
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
				/* amplitude matrix is not integrated. Moreover, if mueller matrix is not needed, then phi_integr has no
				 * action at all. So we assume that phi_integr implies store_mueller (to be checked in param.c).
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
			double tmp3[3],incPolper[3],th,cthet,sthet;
			for (ind=0;ind<angles.theta.N;++ind) {
				index1=0;
				theta=angles.theta.val[ind];
				th=Deg2Rad(theta);
				cthet=cos(th);
				sthet=sin(th);
				for (j=0;j<n;++j) {
					if (angles.type==SG_GRID) phi=angles.phi.val[j];
					else phi=angles.phi.val[ind]; // angles.type==SG_PAIRS
					ph=Deg2Rad(phi);
					// rather complicated (but general) approach to determine rotation angle from XY to scattering plane
					SetScatPlane(cthet,sthet,ph,tmp3,incPolper);
					co=-DotProd(incPolper,incPolY);
					si=DotProd(incPolper,incPolX);
					// transform the amplitude matrix, multiplying by rotation matrix from per-par to X-Y
					s2 = co*EgridX[index+1] + si*EgridY[index+1]; // s2 =  co*s20 + si*s30
					s3 = si*EgridX[index+1] - co*EgridY[index+1]; // s3 =  si*s20 - co*s30
					s4 = co*EgridX[index] + si*EgridY[index]; // s4 =  co*s40 + si*s10
					s1 = si*EgridX[index] - co*EgridY[index]; // s1 =  si*s40 - co*s10
					index+=2;

					if (store_mueller) ComputeMuellerMatrix(matrix,s1,s2,s3,s4,theta);

					if (phi_integr) {
						memcpy(muel_phi+index1,matrix[0],16*sizeof(double));
						index1+=16;
					}
					if (store_scat_grid) {
						if (store_mueller) fprintf(mueller,ANGLE_FORMAT" "ANGLE_FORMAT" "MUEL_FORMAT"\n",theta,phi,
							COMP44M(matrix));
						if (store_ampl) fprintf(ampl,ANGLE_FORMAT" "ANGLE_FORMAT" "AMPL_FORMAT"\n",theta,phi,
							REIM(s1),REIM(s2),REIM(s3),REIM(s4));
					}
				}
				if (phi_integr) {
					PrintToIntegrFile(PHI_UNITY,mueller_int,&max_err,muel_phi,NULL,NULL,matrix,theta);
					PrintToIntegrFile(PHI_COS2,mueller_int_c2,&max_err_c2,muel_phi,muel_phi_buf,cos2,matrix,theta);
					PrintToIntegrFile(PHI_SIN2,mueller_int_s2,&max_err_s2,muel_phi,muel_phi_buf,sin2,matrix,theta);
					PrintToIntegrFile(PHI_COS4,mueller_int_c4,&max_err_c4,muel_phi,muel_phi_buf,cos4,matrix,theta);
					PrintToIntegrFile(PHI_SIN4,mueller_int_s4,&max_err_s4,muel_phi,muel_phi_buf,sin4,matrix,theta);
				}
			}
			if (phi_integr) {
				fprintf(logfile,"\nMaximum relative mean-square error of Mueller integration:\n");
				if (phi_int_type & PHI_UNITY) fprintf(logfile,"  1          -> "RMSE_FORMAT"\n",max_err);
				if (phi_int_type & PHI_COS2) fprintf(logfile,"  cos(2*phi) -> "RMSE_FORMAT"\n",max_err_c2);
				if (phi_int_type & PHI_SIN2) fprintf(logfile,"  sin(2*phi) -> "RMSE_FORMAT"\n",max_err_s2);
				if (phi_int_type & PHI_COS4) fprintf(logfile,"  cos(4*phi) -> "RMSE_FORMAT"\n",max_err_c4);
				if (phi_int_type & PHI_SIN4) fprintf(logfile,"  sin(4*phi) -> "RMSE_FORMAT"\n",max_err_s4);
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

//======================================================================================================================

static bool TestSymVec(const double a[static 3])
/* tests whether a and -a are equivalent under existing symmetries of the scattering problem, i.e. if there exist a
 * combination of symmetries that transforms a into -a. In particular, symR is sufficient for any vector in xy-plane,
 * since double such rotation is equivalent to the in-plane inversion.
 */
{
	return ( ( ((symX||fabs(a[0])<ROUND_ERR) && (symY||fabs(a[1])<ROUND_ERR)) || symR )
		     && (symZ||fabs(a[2])<ROUND_ERR) );
}

//======================================================================================================================

bool TestExtendThetaRange(void)
/* Decides whether the range of [0,180] deg should be extended to 360 deg (for a scattering plane). The test is based on
 * orient_avg and the symmetry of the used scattering plane (either yz or scat_plane). The latter corresponds to logic
 * of calling CalcEplaneYZ() and CalcScatPlane() from CalculateE() below
 */
{
	/* we test the scattering planes by incPolY and exSP assuming one-run-one-SP regime. If two SPs are used for one run
	 * then symR is already valid, which implies that symmetry condition is the same for both scattering planes
	 */
	return !( orient_avg || (yzplane&&TestSymVec(incPolY)) || (scat_plane&&TestSymVec(exSP)) );

}

//======================================================================================================================

static void CalcEplaneYZ(const enum incpol which,const enum Eftype type)
// calculates scattered electric field in a yz-plane (through prop and incPolY)
{
	double incPol[3],incPolper[3],incPolpar[3];
	// where to store calculated field for one plane (actually points to different other arrays)
	doublecomplex *Eplane;
	int i;
	doublecomplex ebuff[3]; // small vector to hold E fields
	double robserver[3];    // small vector for observer in E calculation
	double epar[3];         // unit vector in direction of Epar
	double theta;           // scattering angle
	double co,si;           // temporary, cos and sin of some angle
	double alph;
	TIME_TYPE tstart;
	size_t k_or;
	int orient,Norient;
	enum incpol choice;

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
		else { // special case of the above for alpha=0
			vCopy(incPolX,incPolper); // per <=> X
			vCopy(incPolY,incPolpar); // par <=> Y
		}

		for(orient=0;orient<Norient;orient++) {
			// in case of Rotation symmetry
			tstart = GET_TIME ();
			if (orient==0) choice=which;
			else { // orient==1
				/* Rotation symmetry: calculate per-per from current data. CalculateE is called from calculator with Y
				 * polarization - we now just assume that we have the x-z plane as the scattering plane, rotating in the
				 * negative x-direction. This mimics the real case of X polarization with the y-z plane as scattering
				 * plane. Then IncPolY (par) -> -IncPolX; incPolX (per) -> IncPolY
				 */
				if (which==INCPOL_Y) choice=INCPOL_X;
				else choice=INCPOL_Y; // which==INCPOL_X
				vCopy(incPolper,incPol);
				vCopy(incPolpar,incPolper);
				vMultScal(-1,incPol,incPolpar);
			}
			// initialize Eplane
			if (orient_avg) {
				if (choice==INCPOL_Y) Eplane=ampl_alphaY + 2*nTheta*k_or;
				else Eplane=ampl_alphaX + 2*nTheta*k_or; // choice==INCPOL_X
			}
			else {
				if (choice==INCPOL_Y) Eplane=EyzplY;
				else Eplane=EyzplX; // choice==INCPOL_X
			}

			for (i=0;i<nTheta;i++) {
				theta = i * dtheta_rad;
				co=cos(theta);
				si=sin(theta);
				LinComb(prop,incPolpar,co,si,robserver); // robserver = co*prop + si*incPolpar;
				CalcField(ebuff,robserver);
				// convert to (l,r) frame
				Eplane[2*i]=crDotProd(ebuff,incPolper); // Eper[i]=Esca.incPolper
				LinComb(prop,incPolpar,-si,co,epar);    // epar=-si*prop+co*incPolpar
				Eplane[2*i+1]=crDotProd(ebuff,epar);    // Epar[i]=Esca.epar
			} //  end for i

			// Accumulate Eplane to root and sum
			D("Accumulating Eplane started");
			// accumulate only on processor 0 !, done in one operation
			Accumulate(Eplane,cmplx_type,2*nTheta,&Timing_EPlaneComm);
			D("Accumulating Eplane finished");

			Timing_EPlane = GET_TIME() - tstart;
			Timing_EField += Timing_EPlane;
			TotalEFieldPlane++;
		} // end of orient loop
	} // end of alpha loop
}

//======================================================================================================================

static void CalcScatPlane(const enum incpol which,const enum Eftype type)
// calculates scattered electric field in the plane through ez, prop, incPolX - xz-plane by default
{
	double tmp3[3],incPolper[3],incPolpar[3];
	// where to store calculated field for one plane (actually points to different other arrays)
	doublecomplex *Eplane;
	int i;
	doublecomplex ebuff[3]; // small vector to hold E fields
	double robserver[3];    // small vector for observer in E calculation
	double epar[3];         // unit vector in direction of Epar (for scattered field)
	double theta;           // scattering angle
	double co,si;           // temporary, cos and sin of some angle
	double unitSP[3];       // unit vector (perpendicular to ezLab), which determines the scattering plane
	double alph;
	TIME_TYPE tstart;
	size_t k_or;
	int orient,Norient;
	enum incpol choice;

	if (type==CE_NORMAL) Norient=1; // initialize # orientations
	else Norient=2;                 // type==CE_PARPER

	for (k_or=0;k_or<alpha_int.N;k_or++) {
		// cycle over alpha - for orientation averaging
		if (orient_avg) {
			// rotate polarization basis vectors by -alpha; so per and par are relative to rotated (by -alpha) xz-plane
			alph=Deg2Rad(alpha_int.val[k_or]);
			co=cos(alph);
			si=sin(alph);
			LinComb(incPolX,incPolY,-si,-co,incPolper);  // incPolper = - si*incPolX - co*incPolY;
			LinComb(incPolX,incPolY,co,-si,incPolpar);   // incPolpar = co*incPolX - si*incPolY;
			vCopy(incPolpar,unitSP); // for orientation averaging prop=ez, and epar determines the scattering plane
		}
		else {
			/* the scattering plane is (ez,prop), which contains X-pol;
			 * for the default prop this is xz-plane, a special case of the above for alpha=0
			 */
			vMultScal(-1,incPolY,incPolper); // per <=> -Y
			vCopy(incPolX,incPolpar);        // par <=> X
			vCopy(exSP,unitSP);              // use exSP defined before
		}

		for(orient=0;orient<Norient;orient++) {
			// in case of Rotation symmetry
			tstart = GET_TIME ();
			if (orient==0) choice=which;
			else { // orient==1
				 /* Rotation symmetry: calculate per-per from current data. CalculateE is called from calculator with Y
				 * polarization. To mimic the result for X-polarization, we rotate everything (scattering plane and per,
				 * par directions) by 90 degrees over prop (since that transforms eX into eY).
				 * Since, per,par,prop is RHS orthonormal basis, the rotation is equivalent to
				 * IncPolper=incPolpar(old), incPolpar=-IncPolper(old)
				 */
				if (which==INCPOL_Y) choice=INCPOL_X;
				else choice=INCPOL_Y; // which==INCPOL_X
				// redefine the scattering plane (relative to the code above)
				vCopy(incPolper,tmp3);
				vCopy(incPolpar,incPolper);
				vMultScal(-1,tmp3,incPolpar);
				/* This part can be executed only if prop is along eZ (+ or -).
				 * TODO: remove this limitation, but that requires considering a completely different scattering plane
				 * (rotated with respect to prop), which no more contains ez. That is similar to calculating scat_grid
				 * from a single incident polarization.
				 */
				if (!propAlongZ) LogError(ONE_POS,"Incompatibility error in CalcScatPlane");
				// a general formula, implementing rotation of scattering plane for both prop=+-ez: v->prop x v
				CrossProd(prop,unitSP,tmp3);
				vCopy(tmp3,unitSP);
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
				LinComb(ezLab,unitSP,co,si,robserver); // robserver = co*ezLab + si*unitSP;
				CalcField(ebuff,robserver);
				// convert to (l,r) frame
				Eplane[2*i]=crDotProd(ebuff,incPolper); // Eper[i]=Esca.incPolper
				LinComb(ezLab,unitSP,-si,co,epar);        // epar=-si*ezLab+co*unitSP
				Eplane[2*i+1]=crDotProd(ebuff,epar);    // Epar[i]=Esca.epar
			} //  end for i

			// Accumulate Eplane to root and sum
			D("Accumulating Eplane started");
			// accumulate only on processor 0 !, done in one operation
			Accumulate(Eplane,cmplx_type,2*nTheta,&Timing_EPlaneComm);
			D("Accumulating Eplane finished");

			Timing_EPlane = GET_TIME() - tstart;
			Timing_EField += Timing_EPlane;
			TotalEFieldPlane++;
		} // end of orient loop
	} // end of alpha loop
}

//======================================================================================================================

static void StoreFields(const enum incpol which,doublecomplex * restrict cmplxF,
	double * restrict realF,const char * restrict fname_preffix,const char * restrict tmpl UOIP,
	const char * restrict field_name,const char * restrict fullname)
/* Write any fields on each dipole to file (internal fields, incident beam, polarization, etc.). Accepts both complex
 * and real fields. Processes the one, which is not null.
 *
 * All processors should write the 'field' to temporary file. These files are named by template 'tmpl' and afterwards
 * are concatenated into the file, which name is build by adding a small suffix to 'fname_preffix'. If CE_PARPER is
 * employed then naturally saves only once; use '-sym no' if needed. 'field_name' is used to build column labels (i.e.
 * there is difference in the first row between different fields). 'fullname' is for standard output.
 *
 * This (parallel) algorithm is far from being optimal due to the (redundant) concatenation step. However, this is
 * mainly the limitation of the text file. The only feasible way to improve it is to use binary format like NetCDF
 * (which may work on top of MPI_IO).
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
	/* the following will cause warning by GCC if -Wformat-nonliteral (or -Wformat=2) is used, but we do not know any
	 * other convenient way to make this function work for different file name templates.
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
	if (cmplx_mode) for (j=0;j<local_nRows;j+=3) fprintf(file,GFORM10L"\n",COMP3V(DipoleCoord+j),cvNorm2(cmplxF+j),
		REIM3V(cmplxF+j));
	else for (j=0;j<local_nRows;j+=3) fprintf(file,GFORM7L"\n",COMP3V(DipoleCoord+j),DotProd(realF+j,realF+j),
		COMP3V(realF+j));
	FCloseErr(file,fname,ALL_POS);
#ifdef PARALLEL
	// wait for all processes to save their part of geometry
	Synchronize();
	if (IFROOT) CatNFiles(directory,tmpl,fname_sh);
#endif
	if (IFROOT) printf("%s saved to file\n",fullname);
	Timing_FileIO += GET_TIME() - tstart;
}

//======================================================================================================================

static void ParticleToBeamRF(double vec[static restrict 3])
// transform real vector from particle to beam reference frame
{
	double tmp[3];

	tmp[0]=DotProd(vec,incPolX);
	tmp[1]=DotProd(vec,incPolY);
	tmp[2]=DotProd(vec,prop);
	vCopy(tmp,vec);
}

//======================================================================================================================

static void CalcIntegralScatQuantities(const enum incpol which)
/* calculates all the scattering cross sections, normalized and unnormalized asymmetry parameter, and force on the'
 * particle and each dipole. Cext and Cabs are averaged over orientation, if needed.
 */
{
	// Scattering force, extinction force and radiation pressure per dipole
	double * restrict Frp;
	double Cext,Cabs,Csca,Cdec, // Cross sections
	dummy[3],                // asymmetry parameter*Csca
	Finc_tot[3],Fsca_tot[3],Frp_tot[3], // total extinction and scattering forces, and their sum (radiation pressure)
	Cnorm,            // normalizing factor from force to cross section
	Qnorm;            // normalizing factor from force to efficiency
	FILE * restrict CCfile;
	TIME_TYPE tstart;
	char fname_cs[MAX_FNAME];
	const double *incPol;
	const char *f_suf;

	// redundant initialization to remove warnings
	Cext=Cabs=Csca=Cdec=0;
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
	/* order of calculations is important, when using SQ_FINDIP. Then first dCabs should be calculated, which is further
	 * used to correct Cext
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
				/* These formulae assume that cross sections need to be calculated for unpolarized incident light. It
				 * doesn't matter (equal to the result for any linear polarization), when full averaging over alpha is
				 * performed. However, the difference occur for non-full range of alpha (which is hardly ever used).
				 */
				muel_alpha[-2]=(muel_alpha[-2]+Cext)/2;
				muel_alpha[-1]=(muel_alpha[-1]+Cabs)/2;
			}
		}
	}
	else { // not orient_avg
		if (beamtype==B_DIPOLE) Cdec=DecayCross(); // this is here to be run by all processors
		if (IFROOT) {
			SnprintfErr(ONE_POS,fname_cs,MAX_FNAME,"%s/"F_CS"%s",directory,f_suf);
			CCfile=FOpenErr(fname_cs,"w",ONE_POS);
			if (calc_Cext) PrintBoth(CCfile,"Cext\t= "GFORM"\nQext\t= "GFORM"\n",Cext,Cext*inv_G);
			if (calc_Cabs) PrintBoth(CCfile,"Cabs\t= "GFORM"\nQabs\t= "GFORM"\n",Cabs,Cabs*inv_G);
			if (beamtype==B_DIPOLE) {
				double self=1;
				if (surface) self+=C0dipole_refl/C0dipole;
				double tot=self+Cdec/C0dipole;
				fprintf(CCfile,"\nDecay-rate enhancement\n\n");
				printf("\nDecay-rate enhancement:\n");
				PrintBoth(CCfile,"Total\t= "GFORM"\n",tot);
				if (calc_Cabs) { // for simplicity we keep a single condition here
					double nonrad=Cabs/C0dipole;
					double rad=tot-nonrad;
					// TODO: !!! how nonrad should account for absorption inside the surface???
					PrintBoth(CCfile,"Radiat\t= "GFORM"\n",rad);
					PrintBoth(CCfile,"Nonrad\t= "GFORM"\n",nonrad);
				}
				if (surface) PrintBoth(CCfile,"Surface\t= "GFORM"\n",self);
			}
			if (all_dir) fprintf(CCfile,"\nIntegration\n\n");
			if (calc_Csca) {
				Csca=ScaCross(f_suf);
				PrintBoth(CCfile,"Csca\t= "GFORM"\nQsca\t= "GFORM"\n",Csca,Csca*inv_G);
			}
			if (calc_vec) {
				AsymParm_x(dummy,f_suf);
				AsymParm_y(dummy+1,f_suf);
				AsymParm_z(dummy+2,f_suf);
				PrintBoth(CCfile,"Csca.g\t= "GFORM3V"\n",COMP3V(dummy));
				if (calc_asym) PrintBoth(CCfile,"g\t= "GFORM3V"\n",dummy[0]/Csca,dummy[1]/Csca,dummy[2]/Csca);
			}
		} // end of root
		if (calc_mat_force) {
			if (IFROOT) printf("Calculating the force per dipole\n");
			if (store_force) MALLOC_VECTOR(Frp,double,local_nRows,ALL);
			else Frp=NULL;
			Frp_mat(Finc_tot,Fsca_tot,Frp);
			// Write Cross-Sections and Efficiencies to file
			/* This output contains a number of redundant quantities, like Cext and Csca.g. The main purpose of this is
			 * to use it as independent test (it can be quickly compared against the same quantities calculated by
			 * different means). It also partly addresses issue 133.
			 */
			if (IFROOT) {
				Cnorm = EIGHT_PI;
				Qnorm = EIGHT_PI*inv_G;
				/* in surface mode we keep the result in particle RF that is equal to laboratory one
				 * in free-space mode, the force is transformed into the beam reference frame
				 */
				if (!surface) {
					ParticleToBeamRF(Finc_tot);
					ParticleToBeamRF(Fsca_tot);
				}
				vAdd(Finc_tot,Fsca_tot,Frp_tot);
				PrintBoth(CCfile,"\nBased on radiation forces:\n"
				                 "Cext.a\t= "GFORM3V"\n"
				                 "Csca.g\t= "GFORM3V"\n"
				                 "Cpr\t= "GFORM3V"\n"
				                 "Qpr\t= "GFORM3V"\n",
				                 Cnorm*Finc_tot[0],Cnorm*Finc_tot[1],Cnorm*Finc_tot[2],
				                 -Cnorm*Fsca_tot[0],-Cnorm*Fsca_tot[1],-Cnorm*Fsca_tot[2],
				                 Cnorm*Frp_tot[0],Cnorm*Frp_tot[1],Cnorm*Frp_tot[2],
				                 Qnorm*Frp_tot[0],Qnorm*Frp_tot[1],Qnorm*Frp_tot[2]);
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

//======================================================================================================================

static void StoreIntFields(const enum incpol which)
// Write actual internal fields (not exciting) on each dipole to file
{
	// calculate fields; e_field=P/(V*chi)=chi_inv*P; for anisotropic - by components
	nMult_mat(xvec,pvec,chi_inv);
	// save fields to file
	StoreFields(which,xvec,NULL,F_INTFLD,F_INTFLD_TMP,"E","Internal fields");
}

//======================================================================================================================

int CalculateE(const enum incpol which,const enum Eftype type)
/* Calculate everything for x or y polarized incident light; or one and use symmetry to determine the rest (determined
 * by type)
 */
{
	int exit_status;
	TIME_TYPE tstart;

	tstart=GET_TIME();
	// calculate the incident field Einc; vector b=Einc*cc_sqrt
	D("Generating B");
	GenerateB (which,Einc);
	if (store_beam) StoreFields(which,Einc,NULL,F_BEAM,F_BEAM_TMP,"Einc","Incident beam");
	Timing_IncBeam = GET_TIME() - tstart;
	// calculate solution vector x
	D("Iterative solver started");
	exit_status=IterativeSolver(IterMethod,which);
	D("Iterative solver finished");
	Timing_IntFieldOne = GET_TIME() - tstart;
	Timing_IntField += Timing_IntFieldOne;
	// return if checkpoint (normal) occurred
	if (exit_status==CHP_EXIT) return CHP_EXIT;

	if (yzplane) CalcEplaneYZ(which,type);     // generally plane of incPolY and prop
	if (scat_plane) CalcScatPlane(which,type); // the scattering plane through ez,prop,incPolX - xz by default
	// Calculate the scattered field for the whole solid-angle
	if (all_dir) CalcAlldir();
	// Calculate the scattered field on the given grid of angles
	if (scat_grid) CalcScatGrid(which);
	// Calculate integral scattering quantities (cross sections, asymmetry parameter, electric forces)
	if (calc_Cext || calc_Cabs || calc_Csca || calc_asym || calc_mat_force) CalcIntegralScatQuantities(which);
	// saves internal fields and/or dipole polarizations to text file
	if (store_int_field) StoreIntFields(which);
	if (store_dip_pol) StoreFields(which,pvec,NULL,F_DIPPOL,F_DIPPOL_TMP,"P","Dipole polarizations");
	return 0;
}

//======================================================================================================================

void SaveMuellerAndCS(double * restrict in)
/* saves Mueller matrix and cross sections (averaged) to files; vector in contains values of cross sections and then
 * array of Mueller matrix elements; designed to be called from ROOT only
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
			fprintf(mueller,ANGLE_FORMAT" "MUEL_FORMAT"\n",i*dtheta_deg,COMP16V(muel));
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
