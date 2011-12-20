/* FILE : crosssec.c
 * $Date::                            $
 * Descr: all the functions to calculate scattering quantities (except Mueller matrix); to read
 *        different parameters from files; and initialize orientation of the particle
 *
 * Copyright (C) 2006-2011 ADDA contributors
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <limits.h>
#include "vars.h"
#include "cmplx.h"
#include "const.h"
#include "Romberg.h"
#include "crosssec.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"
#include "io.h"
#include "timing.h"
#include "function.h"

// SEMI-GLOBAL VARIABLES

// defined and initialized in calculator.c
extern double * restrict E2_alldir,* restrict E2_alldir_buffer;
extern const doublecomplex cc[][3];
extern doublecomplex * restrict expsX,* restrict expsY,* restrict expsZ;
// defined and initialized in GenerateB.c
extern const double beam_center_0[3];
// defined and initialized in param.c
extern const double prop_0[3],incPolX_0[3],incPolY_0[3];
extern const enum scat ScatRelation;
// defined and initialized in timing.c
extern TIME_TYPE Timing_EFieldAD,Timing_EFieldADComm,Timing_EFieldSG,Timing_EFieldSGComm,
Timing_ScatQuanComm;

// used in CalculateE.c
Parms_1D phi_sg;
// used in calculator.c
Parms_1D parms_alpha; // parameters of integration over alpha
Parms_1D parms[2];    // parameters for integration over theta,phi or beta,gamma
angle_set beta_int,gamma_int,theta_int,phi_int; // sets of angles
// used in param.c
char avg_string[MAX_PARAGRAPH]; // string for output of function that reads averaging parameters
// used in Romberg.c
bool full_al_range; // whether full range of alpha angle is used

// LOCAL VARIABLES

double dCabs;           // difference between Cabs calculated by 'dr' and 'fin' formulations
bool dCabs_ready=false; // whether dCabs is already calculated

//=====================================================================

INLINE int AlldirIndex(const int theta,const int phi)
// Convert the (theta,phi) couple into a linear array index
{
	return (theta*phi_int.N + phi);
}

//=====================================================================

void InitRotation (void)
/* initialize matrices used for reference frame transformation; based on Mishchenko M.I.
 * "Calculation of the amplitude matrix for a nonspherical particle in a fixed orientation",
 * Applied Optics 39(6):1026-1031. This is so-called zyz-notation or y-convention.
 */
{
	double ca,sa,cb,sb,cg,sg;
	double beta_matr[3][3];
	double alph,bet,gam; // in radians

	// initialization of angle values in radians
	alph=Deg2Rad(alph_deg);
	bet=Deg2Rad(bet_deg);
	gam=Deg2Rad(gam_deg);
	// calculation of rotation matrix
	ca=cos(alph);
	sa=sin(alph);
	cb=cos(bet);
	sb=sin(bet);
	cg=cos(gam);
	sg=sin(gam);

	beta_matr[0][0]=ca*cb*cg-sa*sg;
	beta_matr[0][1]=sa*cb*cg+ca*sg;
	beta_matr[0][2]=-sb*cg;
	beta_matr[1][0]=-ca*cb*sg-sa*cg;
	beta_matr[1][1]=-sa*cb*sg+ca*cg;
	beta_matr[1][2]=sb*sg;
	beta_matr[2][0]=ca*sb;
	beta_matr[2][1]=sa*sb;
	beta_matr[2][2]=cb;
	// rotation of incident field
	MatrVec(beta_matr,prop_0,prop);
	MatrVec(beta_matr,incPolY_0,incPolY);
	MatrVec(beta_matr,incPolX_0,incPolX);
	if (beam_asym) MatrVec(beta_matr,beam_center_0,beam_center);
}

//=====================================================================

static int ATT_UNUSED ReadLine(FILE * restrict file, // opened file
	const char * restrict fname,                     // ... its filename
	char * restrict buf,const int buf_size)          // buffer for line and its size
// reads the first uncommented line; returns 1 if EOF reached
{
	while (!feof(file)) {
		fgets(buf,buf_size,file);
		if (*buf!='#') { // if uncommented
			if (strstr(buf,"\n")==NULL && !feof(file)) LogError(ONE_POS,
				"Buffer overflow while reading '%s' (size of uncommented line > %d)",
				fname,buf_size-1);
			else return 0; // complete line is read
		} // finish reading the commented line
		else while (strstr(buf,"\n")==NULL  && !feof(file)) fgets(buf,buf_size,file);
	}
	return 1;
}

//=====================================================================

static void ReadLineStart(FILE *  restrict file,                  // opened file
	                      const char * restrict fname,            // ... its filename
                          char * restrict buf,const int buf_size, // buffer for line and its size
                          const char * restrict start)            // beginning of the line to search
// reads the first line that starts with 'start'
{
	while (!feof(file)) {
		fgets(buf,buf_size,file);
		if (strstr(buf,start)==buf) { // if correct beginning
			if (strstr(buf,"\n")==NULL && !feof(file)) LogError(ONE_POS,
				"Buffer overflow while reading '%s' (size of essential line > %d)",
				fname,buf_size-1);
			else return; // line found and fits into buffer
		} // finish reading unmatched line
		else while (strstr(buf,"\n")==NULL && !feof(file)) fgets(buf,buf_size,file);
	}
	LogError(ONE_POS,"String '%s' is not found (in correct place) in file '%s'",start,fname);
}

//=====================================================================

INLINE void ScanDouble(FILE * restrict file,const char * restrict fname,char * restrict buf,
	const int buf_size,const char * restrict start,double *res)
/* scans double value from a line starting with exactly 'start'; contains the same arguments as
 * ReadLineStart function, plus pointer to where the result should be placed
 */
{
	ReadLineStart(file,fname,buf,buf_size,start);
	if (sscanf(buf+strlen(start),"%lf",res)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
}

//=====================================================================

INLINE void ScanInt(FILE * restrict file,const char * restrict fname,char * restrict buf,
	const int buf_size,const char * restrict start,int *res)
/* scans integer value from a line starting with exactly 'start'; contains the same arguments as
 * ReadLineStart function, plus pointer to where the result should be placed
 */
{
	double tmp;

	ReadLineStart(file,fname,buf,buf_size,start);
	if (sscanf(buf+strlen(start),"%lf",&tmp)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
	if (tmp<INT_MIN || tmp>INT_MAX)
		LogError(ONE_POS,"Value after '%s' in file '%s' is out of integer bounds",start,fname);
	if (sscanf(buf+strlen(start),"%d",res)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
}

//=====================================================================

INLINE void ScanSizet(FILE * restrict file,const char * restrict fname,char * restrict buf,
	const int buf_size,const char * restrict start,size_t *res)
/* scans large integer value from a line starting with exactly 'start'; contains the same arguments
 * as ReadLineStart function, plus pointer to where the result should be placed.
 * MinGW already provides C99-compliant printf-style functions, but not yet scanf-style. So we have
 * to use workarounds instead of straightforward "%zu" format specifier.
 * TODO: change to %zu, when libmingwex will include scanf.
 */
{
	double tmp;
	unsigned long res_tmp;

	ReadLineStart(file,fname,buf,buf_size,start);
	if (sscanf(buf+strlen(start),"%lf",&tmp)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
	if (tmp<0 || tmp>SIZE_MAX)
		LogError(ONE_POS,"Value after '%s' in file '%s' is out of size_t bounds",start,fname);
	if (sscanf(buf+strlen(start),"%lu",&res_tmp)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
	*res=(size_t)res_tmp;
}

//=====================================================================

INLINE void ScanString(FILE * restrict file,const char * restrict fname,char * restrict buf,
	const int buf_size,const char * restrict start,char * restrict res)
/* scans string value from a line starting with exactly 'start'; contains the same arguments as
 * ReadLineStart function, plus pointer to where the result should be placed; the memory allocated
 * to 'res' should be at least buf_size and independent of buf.
 */
{
	ReadLineStart(file,fname,buf,buf_size,start);
	if (sscanf(buf+strlen(start),"%s",res)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
	/* More secure would be to put field width in format string above (like "%Ns"), however
	 * this field width is defined by the variable buf_size. The latter can only be implemented by
	 * a preliminary printf to get a format string.
	 */
}

//=====================================================================

static void ScanIntegrParms(
	FILE * restrict file,const char * restrict fname, // opened file and filename
	angle_set *a,                                     // pointer to angle set
	Parms_1D *b,                                      // pointer to parameters of integration
	const bool ifcos,                                 // if space angles equally in cos
	char * restrict buf,char*  restrict temp,         // 2 independent buffers
	const int buf_size)                               // and their size
// scan integration parameters for angles from file
{
	size_t i;
	double unit;

	// scan file
	ScanDouble(file,fname,buf,buf_size,"min=",&(a->min));
	ScanDouble(file,fname,buf,buf_size,"max=",&(a->max));
	ScanInt(file,fname,buf,buf_size,"Jmin=",&(b->Jmin));
	ScanInt(file,fname,buf,buf_size,"Jmax=",&(b->Jmax));
	ScanDouble(file,fname,buf,buf_size,"eps=",&(b->eps));
	ScanString(file,fname,buf,buf_size,"equiv=",temp);
	if (strcmp(temp,"true")==0) b->equival=true;
	else if (strcmp(temp,"false")==0) b->equival=false;
	else LogError(ONE_POS,"Wrong argument of 'equiv' option in file %s",fname);
	ScanString(file,fname,buf,buf_size,"periodic=",temp);
	if (strcmp(temp,"true")==0) b->periodic=true;
	else if (strcmp(temp,"false")==0) b->periodic=false;
	else LogError(ONE_POS,"Wrong argument of 'periodic' option in file %s",fname);

	// fill all parameters
	if (a->min==a->max) {
		a->N=b->Grid_size=1;
		b->Jmax=1;
	}
	else {
		// consistency check
		if (a->min>a->max) LogError(ONE_POS,"Wrong range (min="GFORMDEF", max="GFORMDEF") in file "
			"%s (max must be >= min)",a->min,a->max,fname);
		if (b->Jmax<b->Jmin) LogError(ONE_POS,
			"Wrong Jmax (%d) in file %s; it must be >= Jmin (%d)",b->Jmax,fname,b->Jmin);
		if (b->Jmin<1)
			LogError(ONE_POS,"Wrong Jmin (%d) in file %s (must be >=1)",b->Jmin,fname);
		if (b->eps<0)
			LogError(ONE_POS,"Wrong eps ("GFORMDEF") in file %s (must be >=0)",b->eps,fname);
		if (b->Jmax >= (int)(8*sizeof(int))) LogError(ONE_POS,
			"Too large Jmax(%d) in file %s, it will cause integer overflow",b->Jmax,fname);

		a->N=b->Grid_size=(1 << b->Jmax) + 1;
		if (b->equival && a->N>1) (a->N)--;
	}
	// initialize points of integration
	MALLOC_VECTOR(a->val,double,a->N,ALL);
	memory += a->N*sizeof(double);

	if (ifcos) { // make equal intervals in cos(angle)
		// consistency check
		if (a->min<0) LogError(ONE_POS,
			"Wrong min ("GFORMDEF") in file %s (must be >=0 for this angle)",a->min,fname);
		if (a->max>180) LogError(ONE_POS,
			"Wrong max ("GFORMDEF") in file %s (must be <=180 for this angle)",a->max,fname);
		b->min=cos(Deg2Rad(a->max));
		b->max=cos(Deg2Rad(a->min));
		if (fabs(b->min)<ROUND_ERR) b->min=0; // just for convenience of display in log file
		if (fabs(b->max)<ROUND_ERR) b->max=0;
		if (b->Grid_size==1) a->val[0]=a->min;
		else {
			unit = (b->max - b->min)/(b->Grid_size-1);
			for (i=0;i<a->N;i++) a->val[i] = Rad2Deg(acos(b->min+unit*i));
		}
	}
	else { // make equal intervals in angle
		b->min=Deg2Rad(a->min);
		b->max=Deg2Rad(a->max);
		if (b->Grid_size==1) a->val[0]=a->min;
		else {
			unit = (a->max - a->min)/(b->Grid_size-1);
			for (i=0;i<a->N;i++) a->val[i] = a->min + unit*i;
		}
	}
}

//=====================================================================

static enum angleset ScanAngleSet(
	FILE * restrict file,const char * restrict fname, // opened file and filename
	angle_set *a,                                     // pointers to angle set
	char * restrict buf,char * restrict temp,         // 2 independent buffers
	const int buf_size)                               // and their size
// scan range or set of angles (theta or phi) from file (used for scat_grid)
{
	size_t i;
	double unit;
	enum angleset out;

	ScanString(file,fname,buf,buf_size,"type=",temp);
	ScanSizet(file,fname,buf,buf_size,"N=",&(a->N));
	// initialize angle array
	MALLOC_VECTOR(a->val,double,a->N,ALL);
	memory += a->N*sizeof(double);

	if (strcmp(temp,"range")==0) {
		ScanDouble(file,fname,buf,buf_size,"min=",&(a->min));
		ScanDouble(file,fname,buf,buf_size,"max=",&(a->max));
		if (a->min>a->max) LogError(ONE_POS,"Wrong range (min="GFORMDEF", max="GFORMDEF") in file "
			"%s (max must be >= min)",a->min,a->max,fname);
		if (a->N==1) a->val[0]=(a->max + a->min)/2;
		else {
			unit = (a->max - a->min)/(a->N - 1);
			for (i=0;i<a->N;i++) a->val[i] = a->min + unit*i;
		}
		out=AS_RANGE;
	}
	else if (strcmp(temp,"values")==0) {
		ReadLineStart(file,fname,buf,buf_size,"values=");
		for (i=0;i<a->N;i++) {
			fgets(buf,buf_size,file);
			if (strstr(buf,"\n")==NULL  && !feof(file)) LogError(ONE_POS,
				"Buffer overflow while scanning lines in file '%s' (line size > %d)",
				fname,buf_size-1);
			if (sscanf(buf,"%lf\n",a->val+i)!=1)
				LogError(ONE_POS,"Failed scanning values from line '%s' in file '%s'",buf,fname);
		}
		out=AS_VALUES;
	}
	else LogError(ONE_POS,"Unknown type '%s' in file '%s'",temp,fname);
	return out;
}

//=====================================================================

void ReadAvgParms(const char * restrict fname)
// read parameters of orientation averaging from a file
{
	FILE * restrict input;
	char buf[BUF_LINE],temp[BUF_LINE];

	// open file
	input=FOpenErr(fname,"r",ALL_POS);
	//scan file
	ReadLineStart(input,fname,buf,BUF_LINE,"alpha:");
	ScanIntegrParms(input,fname,&alpha_int,&parms_alpha,false,buf,temp,BUF_LINE);
	full_al_range=fabs(alpha_int.max-alpha_int.min-FULL_ANGLE)<FULL_ANGLE*ROUND_ERR;
	ReadLineStart(input,fname,buf,BUF_LINE,"beta:");
	ScanIntegrParms(input,fname,&beta_int,&parms[THETA],true,buf,temp,BUF_LINE);
	ReadLineStart(input,fname,buf,BUF_LINE,"gamma:");
	ScanIntegrParms(input,fname,&gamma_int,&parms[PHI],false,buf,temp,BUF_LINE);
	// close file
	FCloseErr(input,fname,ALL_POS);
	// print info to string
	if (IFROOT) SnprintfErr(ONE_POS,avg_string,MAX_PARAGRAPH,
		"alpha: from "GFORMDEF" to "GFORMDEF" in %zu steps\n"
		"beta: from "GFORMDEF" to "GFORMDEF" in (up to) %zu steps (equally spaced in cosine "
			"values)\n"
		"gamma: from "GFORMDEF" to "GFORMDEF" in (up to) %zu steps\n"
		"see file 'log_orient_avg' for details\n",
		alpha_int.min,alpha_int.max,alpha_int.N,beta_int.min,beta_int.max,beta_int.N,gamma_int.min,
		gamma_int.max,gamma_int.N);
	D("ReadAvgParms finished");
}

//=====================================================================

void ReadAlldirParms(const char * restrict fname)
/* read integration parameters for asymmetry-parameter & C_sca; should not be used together with
 * orientation averaging because they use the same storage space - parms
 */
{
	FILE * restrict input;
	char buf[BUF_LINE],temp[BUF_LINE];

	// open file
	input=FOpenErr(fname,"r",ALL_POS);
	//scan file
	ReadLineStart(input,fname,buf,BUF_LINE,"theta:");
	ScanIntegrParms(input,fname,&theta_int,&parms[THETA],true,buf,temp,BUF_LINE);
	ReadLineStart(input,fname,buf,BUF_LINE,"phi:");
	ScanIntegrParms(input,fname,&phi_int,&parms[PHI],false,buf,temp,BUF_LINE);
	// close file
	FCloseErr(input,fname,ALL_POS);
	// print info
	if (IFROOT) fprintf(logfile,"\n"
		"Scattered field is calculated for all directions (for integrated scattering quantities)\n"
		"theta: from "GFORMDEF" to "GFORMDEF" in (up to) %zu steps (equally spaced in cosine "
			"values)\n"
		"phi: from "GFORMDEF" to "GFORMDEF" in (up to) %zu steps\n"
		"see files 'log_int_***' for details\n\n",
		theta_int.min,theta_int.max,theta_int.N,phi_int.min,phi_int.max,phi_int.N);
	D("ReadAlldirParms finished");
}

//=====================================================================

void ReadScatGridParms(const char * restrict fname)
// read parameters of the grid on which to calculate scattered field
{
	FILE * restrict input;
	char buf[BUF_LINE],temp[BUF_LINE];
	enum angleset theta_type,phi_type;
	size_t i;

	// redundant initialization to remove warnings
	theta_type=phi_type=AS_RANGE;

	// open file
	input=FOpenErr(fname,"r",ALL_POS);
	// scan file
	ScanString(input,fname,buf,BUF_LINE,"global_type=",temp);
	if (strcmp(temp,"grid")==0) {
		angles.type = SG_GRID;
		ReadLineStart(input,fname,buf,BUF_LINE,"theta:");
		theta_type=ScanAngleSet(input,fname,&(angles.theta),buf,temp,BUF_LINE);
		if (phi_integr) {
			ReadLineStart(input,fname,buf,BUF_LINE,"phi_integr:");
			ScanIntegrParms(input,fname,&(angles.phi),&phi_sg,false,buf,temp,BUF_LINE);
			phi_type = AS_RANGE;
		}
		else {
			ReadLineStart(input,fname,buf,BUF_LINE,"phi:");
			phi_type=ScanAngleSet(input,fname,&(angles.phi),buf,temp,BUF_LINE);
		}
		angles.N=MultOverflow(angles.theta.N,angles.phi.N,ONE_POS,"angles.N");;
	}
	else if (strcmp(temp,"pairs")==0) {
		if (phi_integr)
			LogError(ONE_POS,"Integration over phi can't be done with 'global_type=pairs'");
		angles.type = SG_PAIRS;
		ScanSizet(input,fname,buf,BUF_LINE,"N=",&(angles.N));
		angles.theta.N=angles.phi.N=angles.N;
		// malloc angle arrays
		MALLOC_VECTOR(angles.theta.val,double,angles.N,ALL);
		MALLOC_VECTOR(angles.phi.val,double,angles.N,ALL);
		memory += 2*angles.N*sizeof(double);

		ReadLineStart(input,fname,buf,BUF_LINE,"pairs=");
		for (i=0;i<angles.N;i++) {
			fgets(buf,BUF_LINE,input);
			if (strstr(buf,"\n")==NULL && !feof(input)) LogError(ONE_POS,
				"Buffer overflow while scanning lines in file '%s' (line size > %d)",
				fname,BUF_LINE-1);
			if (sscanf(buf,"%lf %lf\n",angles.theta.val+i,angles.phi.val+i)!=2) LogError(
				ONE_POS,"Failed scanning values from line '%s' in file '%s'",buf,fname);
		}
	}
	else LogError(ONE_POS,"Unknown global_type '%s' in file '%s'",temp,fname);
	// close file
	FCloseErr(input,fname,ALL_POS);
	// print info
	if (IFROOT) {
		fprintf(logfile,"\nScattered field is calculated for multiple directions\n");
		if (angles.type==SG_GRID) {
			if (theta_type==AS_RANGE)
				fprintf(logfile,"theta: from "GFORMDEF" to "GFORMDEF" in %zu steps\n",
					angles.theta.min,angles.theta.max,angles.theta.N);
			else if (theta_type==AS_VALUES)
				fprintf(logfile,"theta: %zu given values\n",angles.theta.N);
			if (phi_type==AS_RANGE) {
				fprintf(logfile,"phi: from "GFORMDEF" to "GFORMDEF" in %zu steps\n",
					angles.phi.min,angles.phi.max,angles.phi.N);
				if (phi_integr) fprintf(logfile,"(Mueller matrix is integrated over phi)\n");
			}
			else if (phi_type==AS_VALUES)
				fprintf(logfile,"phi: %zu given values\n",angles.phi.N);
		}
		else if (angles.type==SG_PAIRS)
			fprintf(logfile,"Total %zu given (theta,phi) pairs\n",angles.N);
		fprintf(logfile,"\n");
	}
	D("ReadScatGridParms finished");
}

//=====================================================================*/

void CalcField (doublecomplex * restrict ebuff, // where to write calculated scattering amplitude
                const double * restrict n)      // scattering direction
/* Near-optimal routine to compute the scattered fields at one specific angle (more exactly -
 * scattering amplitude); Specific optimization are possible when e.g. n[0]=0 for scattering in
 * yz-plane, however in this case it is very improbable that the routine will become a bottleneck.
 * The latter happens mostly for cases, when  grid of scattering angles is used with only small
 * fraction of n, allowing simplifications.
 */
{
	double kkk;
	doublecomplex a,m2,dpr;
	doublecomplex sum[3],tbuff[3],tmp={0,0}; // redundant initialization to remove warnings
	int i;
	unsigned short ix,iy1,iy2,iz1,iz2;
	size_t j,jjj;
	double temp, na;
	doublecomplex mult_mat[MAX_NMAT];
	const bool scat_avg=true; // temporary fixed option for SO formulation

	if (ScatRelation==SQ_SO) {
		// !!! this should never happen
		if (anisotropy) LogError(ONE_POS,"Incompatibility error in CalcField");
		// calculate correction coefficient
		if (scat_avg) na=0;
		else na=DotProd(n,prop);
		temp=kd*kd/24;
		for(i=0;i<Nmat;i++) {
			cSquare(ref_index[i],m2);
			// mult_mat=1-(kd^2/24)(m^2-2(n.a)m+1)
			mult_mat[i][RE]=1-temp*(m2[RE]-2*na*ref_index[i][RE]+1);
			mult_mat[i][IM]=temp*(2*na*ref_index[i][IM]-m2[IM]);
		}
	}
	for(i=0;i<3;i++) sum[i][RE]=sum[i][IM]=0.0;
	// prepare values of exponents, along each of the coordinates
	imExp_arr(-kd*n[0],boxX,expsX);
	imExp_arr(-kd*n[1],boxY,expsY);
	imExp_arr(-kd*n[2],local_Nz_unif,expsZ);
	/* not to double the code in the source we use two temporary defines,since the following 'if'
	 * cases differ only by one line of code; (taking 'if' inside the cycle will affect performance)
	 */
	/* this piece of code tries to use that usually only x position changes from dipole to dipole,
	 * saving a complex multiplication seems to be beneficial, even considering bookkeeping
	 * overhead; it may not be as good for very porous particles though, but for them this part of
	 * code is anyway fast relative to the FFT on a large grid; Further optimization is possible
	 * using some kind of plans, i.e. by preliminary analyzing the position of the real dipoles on
	 * the grid.
	 */
#define PART1\
	iy1=iz1=UNDEF;\
	for (j=0;j<local_nvoid_Ndip;++j) {\
		jjj=3*j;\
		/* a=exp(-ikr.n), but r is taken relative to the first dipole of the local box */\
		ix=position[jjj];\
		iy2=position[jjj+1];\
		iz2=position[jjj+2];\
		/* the second part is very improbable, but needed for robustness */\
		if (iy2!=iy1 || iz2!=iz1) {\
			iy1=iy2;\
			iz1=iz2;\
			cMult(expsY[iy2],expsZ[iz2],tmp);\
		}\
		cMult(tmp,expsX[ix],a);
#define PART2\
	/* sum(P*exp(-ik*r.n)) */\
		for(i=0;i<3;i++) {\
			sum[i][RE]+=pvec[jjj+i][RE]*a[RE]-pvec[jjj+i][IM]*a[IM];\
			sum[i][IM]+=pvec[jjj+i][RE]*a[IM]+pvec[jjj+i][IM]*a[RE];\
		}\
	} /* end for j */
	if (ScatRelation==SQ_SO) {
		PART1
		cMultSelf(a,mult_mat[material[j]]);
		PART2
	}
	else if (ScatRelation==SQ_DRAINE || ScatRelation==SQ_FINDIP || ScatRelation==SQ_IGT_SO) {
		PART1
		PART2
	}
#undef PART1
#undef PART2
	// tbuff=(I-nxn).sum=sum-n*(n.sum)
	crDotProd(sum,n,dpr);
	cScalMultRVec(n,dpr,tbuff);
	cvSubtrSelf(sum,tbuff);
	// ebuff=(-i*k^3)*exp(-ikr0.n)*tbuff, where r0=box_origin_unif
	imExp(-WaveNum*DotProd(box_origin_unif,n),a); // a=exp(-ikr0.n)
	kkk=WaveNum*WaveNum*WaveNum;
	// the following additional multiplier implements IGT_SO
	if (ScatRelation==SQ_IGT_SO) kkk*=(1-kd*kd/24);
	tmp[RE]=a[IM]*kkk; // tmp=(-i*k^3)*exp(-ikr0.n)
	tmp[IM]=-a[RE]*kkk;
	cvMultScal_cmplx(tmp,tbuff,ebuff);
}

//=====================================================================

double ExtCross(const double * restrict incPol)
// Calculate the Extinction cross-section
{
	doublecomplex ebuff[3];
	double sum;
	size_t i;

	if (beamtype==B_PLANE) {
		CalcField (ebuff,prop);
		sum=crDotProd_Re(ebuff,incPol); // incPol is real, so no conjugate is needed
		MyInnerProduct(&sum,double_type,1,&Timing_ScatQuanComm);
		sum*=FOUR_PI/(WaveNum*WaveNum);
	}
	else { /* more general formula; normalization is done assuming the unity amplitude of the
	        * electric field in the focal point of the beam; It does not comply with
	        * ScatRelation: SO and IGT_SO (since the corrections used in these formulations are
	        * derived assuming plane incident field. So, effectively, SO and IGT_SO are replaced by
	        * DRAINE when calculating Cext for non-plane beams
	        */
		sum=0;
		for (i=0;i<local_nvoid_Ndip;++i) sum+=cDotProd_Im(pvec+3*i,Einc+3*i); // sum{Im(P.E_inc*)}
		MyInnerProduct(&sum,double_type,1,&Timing_ScatQuanComm);
		sum*=FOUR_PI*WaveNum;
	}
	/* TO ADD NEW BEAM
	 * The formulae above works only if the amplitude of the beam is unity at the focal point.
	 * Either make sure that new beam satisfies this condition or add another case here with
	 * different formulae.
	*/
	if (ScatRelation==SQ_FINDIP) {
		if (dCabs_ready) sum+=dCabs;
		else LogError(ONE_POS,"When using 'fin' scattering quantities formulation, Cabs should be "
			"calculated before Cext");
	}
	return sum;
}

//=====================================================================

double AbsCross(void)
// Calculate the Absorption cross-section for process 0
{
	size_t dip,index;
	int i,j;
	unsigned char mat;
	double sum,temp1,temp2;
	doublecomplex m2;
	double *m; // not doublecomplex=double[2] to allow assignment to it
	double multdr[MAX_NMAT][3];  // multiplier for draine formulation
	double multfin[MAX_NMAT][3]; // multiplier for finite dipoles
	double mult1[MAX_NMAT];    // multiplier, which is always isotropic

	// Cabs = 4*pi*sum
	// In this function IGT_SO is equivalent to DRAINE
	if (ScatRelation==SQ_DRAINE || ScatRelation==SQ_FINDIP || ScatRelation==SQ_IGT_SO) {
		/* code below is applicable only for diagonal (possibly anisotropic) polarizability and
		 * should be rewritten otherwise
		 */

		/* based on Eq.(35) from Yurkin and Hoekstra, "The discrete dipole approximation: an
		 * overview and recent developments," JQSRT 106:558-589 (2007).
		 * summand: Im(P.Eexc(*))-(2/3)k^3*|P|^2=|P|^2*(-Im(1/cc)-(2/3)k^3)
		 */
		temp1 = 2*WaveNum*WaveNum*WaveNum/3;
		for (i=0;i<Nmat;i++) for (j=0;j<3;j++) multdr[i][j]=-cInvIm(cc[i][j])-temp1;
		if (ScatRelation==SQ_FINDIP) {
			/* based on Eq.(31) or equivalently Eq.(58) from the same paper (ref. above)
			 * summand: Im(P.E(*))=-|P|^2*Im(chi_inv), chi_inv=1/(V*chi)
			 * Difference between this formulation and the classical one is also calculated, which
			 * is further used to correct Cext.
			 */
			for (i=0;i<Nmat;i++) for (j=0;j<3;j++) multfin[i][j]=-chi_inv[i][j][IM];
		}
		// main cycle
		if (ScatRelation==SQ_DRAINE || ScatRelation==SQ_IGT_SO) {
			for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip) {
				mat=material[dip];
				index=3*dip;
				for(i=0;i<3;i++) sum+=multdr[mat][i]*cAbs2(pvec[index+i]);
			}
		}
		else if (ScatRelation==SQ_FINDIP) {
			for (dip=0,sum=0,dCabs=0;dip<local_nvoid_Ndip;++dip) {
				mat=material[dip];
				index=3*dip;
				for(i=0;i<3;i++) {
					temp1=cAbs2(pvec[index+i]);
					sum+=multfin[mat][i]*temp1;
					dCabs+=(multfin[mat][i]-multdr[mat][i])*temp1;
				}
			}
		}
	}
	else if (ScatRelation==SQ_SO) {
		// !!! this should never happen
		if (anisotropy) LogError(ONE_POS,"Incompatibility error in AbsCross");
		// calculate mult1
		temp1=kd*kd/6;
		temp2=FOUR_PI/dipvol;
		for (i=0;i<Nmat;i++) {
			m=ref_index[i];
			cSquare(m,m2);
			m2[RE]-=1;
			// mult1=-Im(1/chi)*(1+(kd*Im(m))^2)/d^3;  chi=(m^2-1)/(4*PI)
			mult1[i]=temp2*m2[IM]*(1+temp1*m[IM]*m[IM])/cAbs2(m2);
		}
		// main cycle
		for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip)
			sum+=mult1[material[dip]]*cvNorm2(pvec+3*dip);
	}
	if (ScatRelation==SQ_FINDIP) {
		MyInnerProduct(&dCabs,double_type,1,&Timing_ScatQuanComm);
		dCabs*=(FOUR_PI*WaveNum);
		dCabs_ready=true;
	}
	MyInnerProduct(&sum,double_type,1,&Timing_ScatQuanComm);
	return FOUR_PI*WaveNum*sum;
}

//=====================================================================

void CalcAlldir(void)
// calculate scattered field in many directions
{
	int index,npoints,point;
	size_t i,j;
	TIME_TYPE tstart;
	double robserver[3],incPolpar[3],incPolper[3],cthet,sthet,cphi,sphi,th,ph;
	doublecomplex ebuff[3];

	// Calculate field
	tstart = GET_TIME();
	npoints = theta_int.N*phi_int.N;
	if (IFROOT) printf("Calculating scattered field for the whole solid angle:\n");
	for (i=0,point=0;i<theta_int.N;++i) {
		th=Deg2Rad(theta_int.val[i]);
		cthet=cos(th);
		sthet=sin(th);
		for (j=0;j<phi_int.N;++j) {
			ph=Deg2Rad(phi_int.val[j]);
			cphi=cos(ph);
			sphi=sin(ph);
			// robserver = cos(theta)*prop + sin(theta)*[cos(phi)*incPolX + sin(phi)*incPolY];
			LinComb(incPolX,incPolY,cphi,sphi,robserver);
			LinComb(prop,robserver,cthet,sthet,robserver);
			// calculate scattered field - main bottleneck
			CalcField(ebuff,robserver);
			/* set Epar and Eper - use E2_alldir array to store them this is done to decrease
			 * communications in 1.5 times
			 */
			// incPolper = sin(phi)*incPolX - cos(phi)*incPolY;
			LinComb(incPolX,incPolY,sphi,-cphi,incPolper);
			// incPolpar = -sin(theta)*prop + cos(theta)*[cos(phi)*incPolX + sin(phi)*incPolY];
			LinComb(incPolX,incPolY,cphi,sphi,incPolpar);
			LinComb(prop,incPolpar,-sthet,cthet,incPolpar);
			index=2*point;
			crDotProd(ebuff,incPolper,((doublecomplex*)E2_alldir)[index]);
			crDotProd(ebuff,incPolpar,((doublecomplex*)E2_alldir)[index+1]);
			point++;
			// show progress
			if (((10*point)%npoints)<10 && IFROOT) printf(" %d%%",100*point/npoints);
		}
	}
	// accumulate fields
	Accumulate(E2_alldir,4*npoints,E2_alldir_buffer,&Timing_EFieldADComm);
	// calculate square of the field
	for (point=0;point<npoints;point++)
		E2_alldir[point] = cAbs2(((doublecomplex*)E2_alldir)[2*point]) +
		cAbs2(((doublecomplex*)E2_alldir)[2*point+1]);
	if (IFROOT) printf("  done\n");
	// timing
	Timing_EFieldAD = GET_TIME() - tstart;
	Timing_EField += Timing_EFieldAD;
}

//=====================================================================

void CalcScatGrid(const enum incpol which)
// calculate scattered field in many directions
{
	size_t i,j,n,point,index;
	TIME_TYPE tstart;
	double robserver[3],incPolpar[3],incPolper[3],cthet,sthet,cphi,sphi,th,ph;
	doublecomplex ebuff[3];
	doublecomplex *Egrid; // either EgridX or EgridY

	// Calculate field
	tstart = GET_TIME();
	// choose which array to fill
	if (which==INCPOL_Y) Egrid=EgridY;
	else Egrid=EgridX; // which==INCPOL_X
	// set type of cycling through angles
	if (angles.type==SG_GRID) n=angles.phi.N;
	else n=1; // angles.type==SG_PAIRS
	if (IFROOT) printf("Calculating grid of scattered field:\n");
	// main cycle
	for (i=0,point=0;i<angles.theta.N;++i) {
		th=Deg2Rad(angles.theta.val[i]);
		cthet=cos(th);
		sthet=sin(th);
		for (j=0;j<n;++j) {
			if (angles.type==SG_GRID) ph=Deg2Rad(angles.phi.val[j]);
			else ph=Deg2Rad(angles.phi.val[i]); // angles.type==SG_PAIRS
			cphi=cos(ph);
			sphi=sin(ph);
			// robserver = cos(theta)*prop + sin(theta)*[cos(phi)*incPolX + sin(phi)*incPolY];
			LinComb(incPolX,incPolY,cphi,sphi,robserver);
			LinComb(prop,robserver,cthet,sthet,robserver);
			// calculate scattered field - main bottleneck
			CalcField(ebuff,robserver);
			/* set Epar and Eper - use Egrid array to store them this is done to decrease
			 * communications in 1.5 times
			 */
			// incPolper = sin(phi)*incPolX - cos(phi)*incPolY;
			LinComb(incPolX,incPolY,sphi,-cphi,incPolper);
			// incPolpar = -sin(theta)*prop + cos(theta)*[cos(phi)*incPolX + sin(phi)*incPolY];
			LinComb(incPolX,incPolY,cphi,sphi,incPolpar);
			LinComb(prop,incPolpar,-sthet,cthet,incPolpar);
			index=2*point;
			crDotProd(ebuff,incPolper,Egrid[index]);
			crDotProd(ebuff,incPolpar,Egrid[index+1]);
			point++;
			// show progress; the value is always from 0 to 100, so conversion to int is safe
			if (((10*point)%angles.N)<10 && IFROOT) printf(" %d%%",(int)(100*point/angles.N));
		}
	}
	// accumulate fields; timing
	Accumulate((double *)Egrid,4*angles.N,Egrid_buffer,&Timing_EFieldSGComm);
	if (IFROOT) printf("  done\n");
	Timing_EFieldSG = GET_TIME() - tstart;
	Timing_EField += Timing_EFieldSG;
}

//=====================================================================

static double CscaIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating Csca
{
	res[0]=E2_alldir[AlldirIndex(theta,phi)];
	return 0;
}

//=====================================================================

double ScaCross(char * restrict f_suf)
// Calculate the scattering cross section from the integral
{
	TIME_TYPE tstart;
	char fname[MAX_FNAME];
	double res;

	SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_LOG_INT_CSCA "%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,CscaIntegrand,1,&res,fname);
	res*=FOUR_PI/(WaveNum*WaveNum);
	Timing_Integration += GET_TIME() - tstart;
	return res;
}

//=====================================================================

static double gIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g
{
	double E_square,th,ph;
	th=Deg2Rad(theta_int.val[theta]);
	ph=Deg2Rad(phi_int.val[phi]);

	E_square=E2_alldir[AlldirIndex(theta,phi)];
	res[0] = E_square*sin(th)*cos(ph);
	res[1] = E_square*sin(th)*sin(ph);
	res[2] = E_square*cos(th);
	return 0;
}

//=====================================================================

void AsymParm(double *vec,char * restrict f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca
{
	int comp;
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM "%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gIntegrand,3,vec,log_int);
	for (comp=0;comp<3;++comp) vec[comp]*=FOUR_PI/(WaveNum*WaveNum);
	Timing_Integration += GET_TIME() - tstart;
}

//=====================================================================

static double gxIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g_x
{
	double th=theta_int.val[theta];

	/* a separate case is used to avoid negligibly small non-zero results, which further cause
	 * unusually large relative errors in integration log
	 */
	if (th==180) res[0]=0;
	else res[0]=E2_alldir[AlldirIndex(theta,phi)]*sin(Deg2Rad(th))*cos(Deg2Rad(phi_int.val[phi]));
	return 0;
}

//=====================================================================

void AsymParm_x(double *vec,char * restrict f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca
{
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM F_LOG_X"%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gxIntegrand,1,vec,log_int);
	vec[0] *= FOUR_PI/(WaveNum*WaveNum);
	Timing_Integration += GET_TIME() - tstart;
}

//=====================================================================

static double gyIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g_y
{
	double th=theta_int.val[theta];

	/* a separate case is used to avoid negligibly small non-zero results, which further cause
	 * unusually large relative errors in integration log
	 */
	if (th==180) res[0]=0;
	else res[0]=E2_alldir[AlldirIndex(theta,phi)]*sin(Deg2Rad(th))*sin(Deg2Rad(phi_int.val[phi]));
	return 0;
}

//=====================================================================

void AsymParm_y(double *vec,char * restrict f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca
{
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM F_LOG_Y"%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gyIntegrand,1,vec,log_int);
	vec[0] *= FOUR_PI/(WaveNum*WaveNum);
	Timing_Integration += GET_TIME() - tstart;
}

//=====================================================================

static double gzIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g_z
{
	res[0]=E2_alldir[AlldirIndex(theta,phi)]*cos(Deg2Rad(theta_int.val[theta]));
	return 0;
}

//=====================================================================

void AsymParm_z(double *vec,char * restrict f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca

{
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM F_LOG_Z"%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gzIntegrand,1,vec,log_int);
	vec[0] *= FOUR_PI/(WaveNum*WaveNum);
	Timing_Integration += GET_TIME() - tstart;
}

//=====================================================================

void Frp_mat(double Fsca_tot[3],double * restrict Fsca,double Finc_tot[3],double * restrict Finc,
	double Frp_tot[3],double * restrict Frp)
/* Calculate the Radiation Pressure by direct calculation of the scattering force. Per dipole the
 * force of the incoming photons, the scattering force and the radiation pressure are calculated
 * as intermediate results
 *
 * This should be completely rewritten to work through FFT. Moreover, it should comply with
 * '-scat ...' command line option.
 */
{
	size_t j,l,lll,jjj,jjj_loc,index,comp;
	double * restrict rdipT;
	doublecomplex * restrict pT;
	doublecomplex temp;
	doublecomplex dummy,_E_inc;
	double r,r2; // (squared) absolute distance
	doublecomplex
	n[3],                  // unit vector in the direction of r_{jl}; complex part is always zero
	a,ab1,ab2,c1[3],c2[3], // see chapter ...
	x_cg[3], // complex conjugate P*_j
	Pn_j,    // n_jl.P_l
	Pn_l,    // P*_j.n_jl
	inp;     // P*_j.P_l

	// initialize
	for (comp=0;comp<3;++comp) Fsca_tot[comp]=Finc_tot[comp]=Frp_tot[comp]=0.0;
	// Calculate incoming force per dipole
	for (j=0;j<local_nvoid_Ndip;++j) {
		dummy[RE]=dummy[IM]=0.0;
		for (comp=0;comp<3;++comp) {
			index = 3*j+comp;
			// Im(P.E*inc)
			_E_inc[RE] = Einc[index][RE];
			_E_inc[IM] = -Einc[index][IM];
			cMult(pvec[index],_E_inc,temp);
			cAdd(dummy,temp,dummy);
		}
		Finc[3*j+2] = WaveNum*dummy[IM]/2;
		Finc_tot[2] += Finc[3*j+2];
	}
#ifdef PARALLEL
	/* Because of the parallelization by row-block decomposition the distributed arrays involved
	 * need to be gathered on each node a) DipoleCoord -> rdipT; b) pvec -> pT.
	 * Actually this routine is usually called for two polarizations and rdipT does not change
	 * between the calls. So one AllGather of rdipT can be removed. Number of memory allocations can
	 * also be reduced. But this should be replaced by Fourier anyway.
	 */
	// check if it can work at all
	size_t nRows=MultOverflow(3,nvoid_Ndip,ONE_POS_FUNC);
	// allocates a lot of additional memory
	MALLOC_VECTOR(rdipT,double,nRows,ALL);
	MALLOC_VECTOR(pT,complex,nRows,ALL);
	// gathers everything
	AllGather(DipoleCoord,rdipT,double3_type,&Timing_ScatQuanComm);
	AllGather(pvec,pT,cmplx3_type,&Timing_ScatQuanComm);
#else
	pT=pvec;
	rdipT=DipoleCoord;
#endif
	// Calculate scattering force per dipole
	/* Currently, testing the correctness of the following is very hard because the original code
	 * lacks comments. So the best we can do before rewriting it completely is to test that it
	 * produces reasonable results for a number of test cases.
	 */
	for (j=local_nvoid_d0;j<local_nvoid_d1;++j) {
		jjj = 3*j;
		jjj_loc=3*(j-local_nvoid_d0);

		for (l=0;l<nvoid_Ndip;++l) if (j!=l) {
			lll = 3*l;
			r2 = 0;
			Pn_j[RE]=Pn_j[IM]=Pn_l[RE]=Pn_l[IM]=inp[RE]=inp[IM]=0.0;
			// Set distance related variables
			for (comp=0;comp<3;++comp) {
				n[comp][IM] = 0;
				n[comp][RE] = rdipT[jjj+comp] - rdipT[lll+comp];
				r2 += n[comp][RE]*n[comp][RE];
			}
			r = sqrt(r2);
			n[0][RE]/=r; n[1][RE]/=r; n[2][RE]/=r;
			// Set the scalar products a.b1 and a.b2
			a[RE] = cos(WaveNum*r);
			a[IM] = sin(WaveNum*r);
			ab1[RE] = 3/(r2*r2) - WaveNum*WaveNum/r2;
			ab2[RE] = -WaveNum*WaveNum/r2;
			ab1[IM] = -3*WaveNum/(r*r2);
			ab2[IM] = WaveNum*WaveNum*WaveNum/r;
			cMultSelf(ab1,a);
			cMultSelf(ab2,a);
			// Prepare c1 and c2
			for (comp=0;comp<3;++comp) {
				x_cg[comp][RE] = pT[jjj+comp][RE];
				x_cg[comp][IM] = -pT[jjj+comp][IM];
				cMult(x_cg[comp],n[comp],temp);
				cAdd(Pn_j,temp,Pn_j);
				cMult(n[comp],pT[lll+comp],temp);
				cAdd(Pn_l,temp,Pn_l);
				cMult(x_cg[comp],pT[lll+comp],temp);
				cAdd(inp,temp,inp);
			}
			for (comp=0;comp<3;++comp) {
				// Set c1
				cMult(Pn_j,Pn_l,temp);
				cMult(n[comp],temp,c1[comp]);
				c1[comp][RE] *= -5;
				c1[comp][IM] *= -5;
				cMult(inp,n[comp],temp);
				cAdd(c1[comp],temp,c1[comp]);
				cMult(Pn_j,pT[lll+comp],temp);
				cAdd(c1[comp],temp,c1[comp]);
				cMult(x_cg[comp],Pn_l,temp);
				cAdd(c1[comp],temp,c1[comp]);
				// Set c2
				cMult(Pn_j,Pn_l,temp);
				cMult(n[comp],temp,c2[comp]);
				c2[comp][RE] *= -1;
				c2[comp][IM] *= -1;
				cMult(inp,n[comp],temp);
				cAdd(c2[comp],temp,c2[comp]);
				// Fsca_{jl} = ...
				cMultSelf(c1[comp],ab1);
				cMultSelf(c2[comp],ab2);
				Fsca[jjj_loc+comp] += (c1[comp][RE] + c2[comp][RE])/2;
			}
		} // end l-loop
		// Concluding
		for (comp=0;comp<3;++comp) {
			Fsca_tot[comp] += Fsca[jjj_loc+comp];
			Frp[jjj_loc+comp] = Finc[jjj_loc+comp] + Fsca[jjj_loc+comp];
			Frp_tot[comp] += Frp[jjj_loc+comp];
		}
	} // end j-loop

	// Accumulate the total forces on all nodes
	MyInnerProduct(Finc_tot,double_type,3,&Timing_ScatQuanComm);
	MyInnerProduct(Fsca_tot,double_type,3,&Timing_ScatQuanComm);
	MyInnerProduct(Frp_tot,double_type,3,&Timing_ScatQuanComm);
#ifdef PARALLEL
	Free_general(rdipT);
	Free_cVector(pT);
#endif
}
