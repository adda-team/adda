/* File: CalculateE.c
 * $Author$
 * $Date::                            $
 * Descr: the module to calculate the E field and all scattering quantities
 *
 *        Routines for most scattering quantities are in crosssec.c. Also saves internal fields to
 *        file (optional).
 *
 *        January 2004 : include module to compute full Mueller Matrix over full space angle, not
 *        very efficient, must be improved (A. Hoekstra)
 *
 *        Previous versions by Alfons Hoekstra
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "crosssec.h"
#include "Romberg.h"
#include "io.h"
#include "vars.h"
#include "memory.h"
#include "timing.h"
#include "function.h"

// SEMI-GLOBAL VARIABLES

// defined and initialized in calculator.c
extern double *muel_phi,*muel_phi_buf;
extern doublecomplex *EplaneX, *EplaneY;
extern double *Eplane_buffer;
extern const double dtheta_deg,dtheta_rad;
extern doublecomplex *ampl_alphaX,*ampl_alphaY;
extern double *muel_alpha;
// defined and initialized in crosssec.c
extern const Parms_1D phi_sg;
// defined and initialized in param.c
extern const int store_int_field,store_dip_pol,store_beam,store_scat_grid,calc_Cext,calc_Cabs,
calc_Csca,calc_vec,calc_asym,calc_mat_force,store_force,phi_int_type;
// defined and initialized in timing.c
extern TIME_TYPE Timing_EFieldPlane,Timing_comm_EField,
Timing_IntField,Timing_IntFieldOne,Timing_ScatQuan;
extern unsigned long TotalEFieldPlane;

// used in iterative.c
TIME_TYPE tstart_CE;

// EXTERNAL FUNCTIONS

// GenerateB.c
void GenerateB(char which,doublecomplex *x);
// iterative.c
int IterativeSolver(int method);

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
// this function is currently not used
static void ComputeMuellerMatrixNorm(double [4][4],const doublecomplex,const doublecomplex,
	const doublecomplex,const doublecomplex) ATT_UNUSED;

static void ComputeMuellerMatrixNorm(double matrix[4][4],const doublecomplex s1,
	const doublecomplex s2,const doublecomplex s3,const doublecomplex s4)
/* computer mueller matrix from scattering matrix elements s1, s2, s3, s4, according to formula
 * 3.16 from Bohren and Huffman; normalize all elements to S11 (except itself)
 */
{
	matrix[0][0] = 0.5*(cMultConRe(s1,s1)+cMultConRe(s2,s2)+cMultConRe(s3,s3)+cMultConRe(s4,s4));
	matrix[0][1] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+cMultConRe(s4,s4)-cMultConRe(s3,s3))
	             / matrix[0][0];
	matrix[0][2] = (cMultConRe(s2,s3)+cMultConRe(s1,s4))/matrix[0][0];
	matrix[0][3] = (cMultConIm(s2,s3)-cMultConIm(s1,s4))/matrix[0][0];

	matrix[1][0] = 0.5*(cMultConRe(s2,s2)-cMultConRe(s1,s1)+cMultConRe(s3,s3)-cMultConRe(s4,s4))
	             / matrix[0][0];
	matrix[1][1] = 0.5*(cMultConRe(s2,s2)+cMultConRe(s1,s1)-cMultConRe(s3,s3)-cMultConRe(s4,s4))
	             / matrix[0][0];
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

//==============================================================
INLINE void InitMuellerIntegrFile(const int type,const char *fname,FILE **file,char *buf,
	double **mult)
/* If 'phi_int_type' matches 'type', appropriate file (name given by 'fname') is created (with
 * handle '*file'), and heading line is put into it. String buffer 'buf' is used. Vector of
 * multipliers '*mult' is allocated if its pointer is specified.
 */
{
	if (phi_int_type & type) {
		sprintf(buf,"%s/%s",directory,fname);
		(*file)=FOpenErr(buf,"w",ONE_POS);
		fprintf(*file,"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44 "
			"RMSE(integr)\n");
		if (mult!=NULL) MALLOC_VECTOR(*mult,double,angles.phi.N,ALL);
	}
}

//==============================================================

INLINE void PrintToIntegrFile(const int type,FILE *file,double *maxerr,const double *muel,
	double *muel_buf,const double *mult,double matrix[4][4],const double theta)
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
		fprintf(file,"%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E"
			" %.10E %.10E %.10E %.10E %.3E\n",theta,matrix[0][0],matrix[0][1],matrix[0][2],
			matrix[0][3],matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],matrix[2][0],
			matrix[2][1],matrix[2][2],matrix[2][3],matrix[3][0],matrix[3][1],matrix[3][2],
			matrix[3][3],err);
	}
}

//==============================================================

INLINE void CloseIntegrFile(const int type,FILE *file,const char *fname,double *mult)
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
	FILE *mueller,*mueller_int,*mueller_int_c2,*mueller_int_s2,*mueller_int_c4,*mueller_int_s4;
	double *cos2,*sin2,*cos4,*sin4;
	double matrix[4][4];
	double theta,phi,ph,
	max_err,max_err_c2,max_err_s2,max_err_c4,max_err_s4;
	doublecomplex s1,s2,s3,s4,s10,s20,s30,s40;
	char fname[MAX_FNAME];
	int i;
	size_t index,index1,k_or,j,n,ind;
	double co,si;
	double alph;
	TIME_TYPE tstart;

	if (ringid!=ROOT) return;

	if (orient_avg) { // Amplitude matrix stored in ampl_alplha is
		index1=index=0; // transformed into Mueller matrix stored in muel_alpha
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
	else {
		tstart=GET_TIME(); // here Mueller matrix is saved to file
		if (yzplane) {
			strcpy(fname,directory);
			strcat(fname,"/" F_MUEL);
			mueller=FOpenErr(fname,"w",ONE_POS);
			fprintf(mueller,
				"theta s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");
			for (i=0;i<nTheta;i++) {
				theta=i*dtheta_deg;
				ComputeMuellerMatrix(matrix,EplaneX[2*i],EplaneY[2*i+1],EplaneX[2*i+1],
					EplaneY[2*i]);
				fprintf(mueller,"%.2f %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E %.10E "
					"%.10E %.10E %.10E %.10E %.10E %.10E\n",theta,matrix[0][0],matrix[0][1],
					matrix[0][2],matrix[0][3],matrix[1][0],matrix[1][1],matrix[1][2],matrix[1][3],
					matrix[2][0],matrix[2][1],matrix[2][2],matrix[2][3],matrix[3][0],matrix[3][1],
					matrix[3][2],matrix[3][3]);
			}
			FCloseErr(mueller,F_MUEL,ONE_POS);
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
				strcpy(fname,directory);
				strcat(fname,"/" F_MUEL_SG);
				mueller=FOpenErr(fname,"w",ONE_POS);
				fprintf(mueller,
					"theta phi s11 s12 s13 s14 s21 s22 s23 s24 s31 s32 s33 s34 s41 s42 s43 s44\n");
			}
			if (phi_integr) { // also initializes arrays of multipliers
				InitMuellerIntegrFile(PHI_UNITY,F_MUEL_INT,&mueller_int,fname,NULL);
				InitMuellerIntegrFile(PHI_COS2,F_MUEL_C2,&mueller_int_c2,fname,&cos2);
				InitMuellerIntegrFile(PHI_SIN2,F_MUEL_S2,&mueller_int_s2,fname,&sin2);
				InitMuellerIntegrFile(PHI_COS4,F_MUEL_C4,&mueller_int_c4,fname,&cos4);
				InitMuellerIntegrFile(PHI_SIN4,F_MUEL_S4,&mueller_int_s4,fname,&sin4);
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
			else if (angles.type==SG_PAIRS) n=1;
			// main cycle
			index=0;
			max_err=max_err_c2=max_err_s2=max_err_c4=max_err_s4=0;
			for (ind=0;ind<angles.theta.N;++ind) {
				index1=0;
				theta=angles.theta.val[ind];
				for (j=0;j<n;++j) {
					if (angles.type==SG_GRID) phi=angles.phi.val[j];
					else if (angles.type==SG_PAIRS) phi=angles.phi.val[ind];
					ph=Deg2Rad(phi);
					co=cos(ph);
					si=sin(ph);
					// read amplitude matrix from memory
					cEqual(EgridY[index],s10);
					cEqual(EgridY[index+1],s30);
					cEqual(EgridX[index],s40);
					cEqual(EgridX[index+1],s20);
					// transform it, multiplying by rotation matrix from per-par to X-Y
					cLinComb(s20,s30,co,si,s2);  // s2 =  co*s20 + si*s30
					cLinComb(s20,s30,si,-co,s3); // s3 =  si*s20 - co*s30
					cLinComb(s40,s10,co,si,s4);  // s4 =  co*s40 + si*s10
					cLinComb(s40,s10,si,-co,s1); // s1 =  si*s40 - co*s10

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
				if (phi_int_type & PHI_UNITY) fprintf(logfile,"  1          -> %.3E\n",max_err);
				if (phi_int_type & PHI_COS2) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
				if (phi_int_type & PHI_SIN2) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
				if (phi_int_type & PHI_COS4) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
				if (phi_int_type & PHI_SIN4) fprintf(logfile,"  cos(2*phi) -> %.3E\n",max_err_c2);
			}
			// close files; free arrays
			if (store_scat_grid) FCloseErr(mueller,F_MUEL_SG,ONE_POS);
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

static void CalcEplane(const char which,const int type)
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
	char choice;

	incPolper=incPol_tmp1; // initialization of per and par polarizations
	incPolpar=incPol_tmp2;

	if (type==CE_NORMAL) Norient=1; // initialize # orientations
	else if (type==CE_PARPER) Norient=2;

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
			else if (orient==1) {
				/* Rotation symmetry: calculate per-per from current data. CalculateE is called
				 * from calculator with 'Y' polarization - we now just assume that we have
				 * the x-z plane as the scattering plane, rotating in the negative x-direction.
				 * This mimics the real case of X polarization with the y-z plane as scattering plane
				 * Then IncPolY -> -IncPolX; incPolX -> IncPolY
				 * */
				if (which=='X') choice='Y';
				else if (which=='Y') choice='X';
				incPol=incPolper;
				incPolper=incPolpar;
				incPolpar=incPol;
				MultScal(-1,incPolpar,incPolpar);
			}
			// initialize Eplane
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
			Accumulate((double *)Eplane,4*nTheta,Eplane_buffer,&Timing_comm_EField);
			D("Accumulating Eplane finished");

			Timing_EFieldPlane = GET_TIME() - tstart;
			Timing_EField += Timing_EFieldPlane;
			TotalEFieldPlane++;
		} // end of orient loop
	} // end of alpha loop
}

//============================================================

static void CalcIntegralScatQuantities(const char which)
/* calculates all the scattering cross sections, normalized and unnormalized asymmetry parameter,
 * and force on the particle and each dipole. Cext and Cabs are averaged over orientation,
 * if needed.
 */
{
	double *Fsca,*Finc,*Frp; // Scattering force, extinction force and radiation pressure per dipole
	double Cext,Cabs,Csca,   // Cross sections
	dummy[3],         // asymmetry parameter*Csca
	Fsca_tot[3],      // total scattering force
	Finc_tot[3],      // total extinction force
	Frp_tot[3],       // total radiation pressure
	Cnorm,            // normalizing factor from force to cross section
	Qnorm;            // normalizing factor from force to efficiency
	FILE *VisFrp,*CCfile;
	TIME_TYPE tstart;
	char fname_cs[MAX_FNAME],fname_frp[MAX_FNAME];
	size_t j;
	double const *incPol;
	char f_suf[MAX_WORD];

	D("Calculation of cross sections started");
	tstart = GET_TIME();

	strcpy(fname_cs,directory);
	strcat(fname_cs,"/" F_CS);
	if (which == 'X') {
		strcpy(f_suf,F_XSUF);
		incPol=incPolX;
	}
	if (which == 'Y') {
		strcpy(f_suf,F_YSUF);
		incPol=incPolY;
	}
	strcat(fname_cs,f_suf);
	if (calc_Cext) Cext = ExtCross(incPol);
	if (calc_Cabs) Cabs = AbsCross();
	D("Cext and Cabs calculated");
	if (orient_avg) {
		if (ringid==ROOT) {
			if (which == 'Y') { // assumed that first call of CalculateE is with 'Y' flag
				muel_alpha[-2]=Cext;
				muel_alpha[-1]=Cabs;
			}
			else if (which == 'X') {
				muel_alpha[-2]=(muel_alpha[-2]+Cext)/2;
				muel_alpha[-1]=(muel_alpha[-1]+Cabs)/2;
			}
		}
	}
	else { // not orient_avg
		if (ringid==ROOT) {
			CCfile=FOpenErr(fname_cs,"w",ONE_POS);
			if (calc_Cext) PrintBoth(CCfile,"Cext\t= %.10g\nQext\t= %.10g\n",Cext,Cext*inv_G);
			if (calc_Cabs) PrintBoth(CCfile,"Cabs\t= %.10g\nQabs\t= %.10g\n",Cabs,Cabs*inv_G);
			if (all_dir) fprintf(CCfile,"\nIntegration\n\n");
			if (calc_Csca) {
				Csca=ScaCross(f_suf);
				PrintBoth(CCfile,"Csca\t= %.10g\nQsca\t= %.10g\n",Csca,Csca*inv_G);
			}
			if (calc_vec) {
				AsymParm_x(dummy,f_suf);
				AsymParm_y(dummy+1,f_suf);
				AsymParm_z(dummy+2,f_suf);
				PrintBoth(CCfile,"Csca.g\t= (%.10g,%.10g,%.10g)\n",dummy[0],dummy[1],dummy[2]);
				if (calc_asym) PrintBoth(CCfile,"g\t= (%.10g,%.10g,%.10g)\n",
					dummy[0]/Csca,dummy[1]/Csca,dummy[2]/Csca);
			}
		} // end of ROOT
		if (calc_mat_force) {
			MALLOC_VECTOR(Fsca,double,3*local_nvoid_Ndip,ALL);
			MALLOC_VECTOR(Finc,double,3*local_nvoid_Ndip,ALL);
			MALLOC_VECTOR(Frp,double,3*local_nvoid_Ndip,ALL);
			for (j=0;j<3*local_nvoid_Ndip;j++) Fsca[j]=Finc[j]=Frp[j]=0;
			PRINTZ("Calculating the force per dipole\n");
			// Calculate forces
			Frp_mat(Fsca_tot,Fsca,Finc_tot,Finc,Frp_tot,Frp);
			// Write Cross-Sections and Efficiencies to file
			if (ringid==ROOT) {
				Cnorm = EIGHT_PI;
				Qnorm = EIGHT_PI*inv_G;
				PrintBoth(CCfile,"\nMatrix\n"\
				                 "Cext\t= %.10g\nQext\t= %.10g\n"\
				                 "Csca.g\t= (%.10g,%.10g,%.10g)\n"\
				                 "Cpr\t= (%.10g,%.10g,%.10g)\n"\
				                 "Qpr\t= (%.10g,%.10g,%.10g)\n",Cnorm*Finc_tot[2],Qnorm*Finc_tot[2],
				                 -Cnorm*Fsca_tot[0],-Cnorm*Fsca_tot[1],-Cnorm*Fsca_tot[2],
				                 Cnorm*Frp_tot[0],Cnorm*Frp_tot[1],Cnorm*Frp_tot[2],
				                 Qnorm*Frp_tot[0],Qnorm*Frp_tot[1],Qnorm*Frp_tot[2]);
				if (store_force) {
					// Write Radiation pressure per dipole to file
					strcpy(fname_frp,directory);
					strcat(fname_frp,"/" F_FRP);
					strcat(fname_frp,f_suf);
					strcat(fname_frp,".dat"); // TODO: should be removed in the future
					VisFrp=FOpenErr(fname_frp,"w",ONE_POS);
					fprintf(VisFrp,"#sphere  x=%.10g  m=%.10g%+.10gi\n"\
						"#number of real dipoles  %.0f\n"\
						"#Forces per dipole\n"\
						"#r.x r.y r.z F.x F.y F.z\n",
						ka_eq,ref_index[0][RE],ref_index[0][IM],nvoid_Ndip);
					for (j=0;j<local_nvoid_Ndip;++j) fprintf(VisFrp,
						"%.10g %.10g %.10g %.10g %.10g %.10g\n",
						DipoleCoord[3*j],DipoleCoord[3*j+1],
						DipoleCoord[3*j+2],
						Frp[3*j],Frp[3*j+1],Frp[3*j+2]);
					FCloseErr(VisFrp,fname_frp,ONE_POS);
				}
			}
			Free_general(Fsca);
			Free_general(Finc);
			Free_general(Frp);
		}
		if (ringid==ROOT) FCloseErr(CCfile,fname_cs,ONE_POS);
	}
	D("Calculation of cross sections finished");
	Timing_ScatQuan += GET_TIME() - tstart;
}

//============================================================

static void StoreFields(const char which,doublecomplex *field,const char *fname_preffix,
	const char *tmpl,const char *field_name,const char *fullname)
/* Write any fields on each dipole to file (internal fields, incident beam, polarization, etc.).
 * All processors should write the 'field' to temporary file. These files are named by template
 * 'tmpl' and afterwards are concatenated into the file, which name is build by adding a small
 * suffix to 'fname_preffix'. If CE_PARPER is employed then naturally saves only once; use '-sym no'
 * if needed. 'field_name' is used to build column labels (i.e. there is difference in the first row
 * between different fields). 'fullname' is for standard output.
 */
{
	FILE *file; // file to store the fields
	size_t i,j;
	TIME_TYPE tstart;
	char fname[MAX_FNAME],fname_sh[MAX_FNAME_SH];

	tstart=GET_TIME();
	// build file name (without directory)
	strcpy(fname_sh,fname_preffix);
	if (which=='X') strcat(fname_sh,F_XSUF);
	else if (which=='Y') strcat(fname_sh,F_YSUF);
	// choose filename for direct saving
#ifdef PARALLEL
	sprintf(fname,"%s/",directory);
	sprintf(fname+strlen(fname),tmpl,ringid);
#else
	sprintf(fname,"%s/%s",directory,fname_sh);
#endif
	file=FOpenErr(fname,"w",ALL_POS);
	// print head of file
#ifdef PARALLEL
	if (ringid==0) { // this condition can be different from being ROOT
#endif
		fprintf(file,"x y z |%s|^2 %sx.r %sx.i %sy.r %sy.i %sz.r %sz.i\n",field_name,field_name,
			field_name,field_name,field_name,field_name,field_name);
#ifdef PARALLEL
	} // end of if
#endif
	// saves fields to file
	for (i=0;i<local_nvoid_Ndip;++i) {
		j=3*i;
		fprintf(file,
			"%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g\n",
			DipoleCoord[j],DipoleCoord[j+1],DipoleCoord[j+2],cvNorm2(field+j),
			field[j][RE],field[j][IM],field[j+1][RE],field[j+1][IM],field[j+2][RE],field[j+2][IM]);
	}
	FCloseErr(file,fname,ALL_POS);
#ifdef PARALLEL
	// wait for all processes to save their part of geometry
	Synchronize();
	if (ringid==ROOT) CatNFiles(directory,tmpl,fname_sh);
#endif
	PRINTZ("%s saved to file\n",fullname);
	Timing_FileIO += GET_TIME() - tstart;
}

//============================================================

static void StoreIntFields(const char which)
// Write actual internal fields (not exciting) on each dipole to file
{
	double V;
	doublecomplex hi,hi_inv[MAX_NMAT];
	unsigned char mat;
	size_t i;
	int j;

	// calculate multipliers
	V=gridspace*gridspace*gridspace;
	for (j=0;j<Ncomp*Nmat;j++) {
		// hi_inv=1/(V*hi)=4*PI/(V(m^2-1)); for anisotropic - by components
		cSquare(ref_index[j],hi);
		hi[RE]-=1;
		cMultReal(V,hi,hi);
		cInv(hi,hi_inv[j]);
		cMultReal(FOUR_PI,hi_inv[j],hi_inv[j]);
	}
	// calculate fields
	for (i=0;i<local_nvoid_Ndip;++i) {
		mat=(unsigned char)(material[i]*Ncomp);
		// e_field=P/(V*hi); for anisotropic - by components
		for (j=0;j<3;j++) {
			cMult(hi_inv[mat],pvec[3*i+j],xvec[3*i+j]);
			if (anisotropy) mat++;
		}
	}
	// save fields to file
	StoreFields(which,xvec,F_INTFLD,F_INTFLD_TMP,"E","Internal fields");
}

//============================================================

int CalculateE(const char which,const int type)
/* Calculate everything for x or y polarized incident light; or one and use symmetry to determine
 * the rest (determined by type)
 */
{
	int exit_status;

	tstart_CE=GET_TIME();
	// calculate the incident field Einc; vector b=Einc*cc_sqrt
	D("Generating B");
	GenerateB (which, Einc);
	if (store_beam) StoreFields(which,Einc,F_BEAM,F_BEAM_TMP,"Einc","Incident beam");
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
	if (store_dip_pol) StoreFields(which,pvec,F_DIPPOL,F_DIPPOL_TMP,"P","Dipole polarizations");
	return 0;
}
