/* FILE : crosssec.c
 * $Date::                            $
 * Descr: all the functions to calculate scattering quantities (except Mueller matrix); to read different parameters
 *        from files; and initialize orientation of the particle
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
#include "crosssec.h" // corresponding header
// project headers
#include "cmplx.h"
#include "comm.h"
#include "debug.h"
#include "io.h"
#include "memory.h"
#include "Romberg.h"
#include "timing.h"
#include "vars.h"
// system headers
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

// SEMI-GLOBAL VARIABLES

// defined and initialized in calculator.c
extern doublecomplex * restrict E_ad;
extern double * restrict E2_alldir;
extern const doublecomplex cc[][3];
#ifndef SPARSE
extern doublecomplex * restrict expsX,* restrict expsY,* restrict expsZ;
#endif
// defined and initialized in GenerateB.c
extern const double beam_center_0[3];
//extern doublecomplex eIncRefl[3],eIncTran[3];
// defined and initialized in param.c
extern const double incPolX_0[3],incPolY_0[3];
extern const enum scat ScatRelation;
// defined and initialized in timing.c
extern TIME_TYPE Timing_EFieldAD,Timing_EFieldADComm,Timing_EFieldSG,Timing_EFieldSGComm,
Timing_ScatQuanComm;

// used in CalculateE.c
Parms_1D phi_sg;
double ezLab[3]; // basis vector ez of laboratory RF transformed into the RF of particle
double exSP[3]; // second vector, which determines the scattering plane passing through prop and ez (in particle RF)
// used in calculator.c
Parms_1D parms_alpha; // parameters of integration over alpha
Parms_1D parms[2];    // parameters for integration over theta,phi or beta,gamma
angle_set beta_int,gamma_int,theta_int,phi_int; // sets of angles
// used in param.c
char avg_string[MAX_PARAGRAPH]; // string for output of function that reads averaging parameters
// used in Romberg.c
bool full_al_range; // whether full range of alpha angle is used

// LOCAL VARIABLES

static double exLab[3],eyLab[3]; // basis vectors of laboratory RF transformed into the RF of particle

//======================================================================================================================

static inline int AlldirIndex(const int theta,const int phi)
// Convert the (theta,phi) couple into a linear array index
{
	return (theta*phi_int.N + phi);
}

//======================================================================================================================

void InitRotation (void)
/* initialize matrices used for reference frame transformation; based on Mishchenko M.I. "Calculation of the amplitude
 * matrix for a nonspherical particle in a fixed orientation", Applied Optics 39(6):1026-1031. This is so-called
 * zyz-notation or y-convention.
 */
{
	double ca,sa,cb,sb,cg,sg;
	double beta_matr[3][3];
	double alph,bet,gam; // in radians
	double ex[3]; // second vector of scattering plane in the laboratory RF

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
	// rotation of incident field and basis vectors of laboratory RF
	MatrVec(beta_matr,prop_0,prop);
	MatrVec(beta_matr,incPolY_0,incPolY);
	MatrVec(beta_matr,incPolX_0,incPolX);
	MatrColumn(beta_matr,0,exLab);
	MatrColumn(beta_matr,1,eyLab);
	MatrColumn(beta_matr,2,ezLab);
	/* define the second vector (x) of the scattering plane (through prop and ez)
	 * if prop is along ez, we choose the original ex (that is - we generally assume that beam propagates from the left
	 * side and positive x-direction is in the right side (left and right is with respect to the possible rotation of
	 * ex,ey over the z-axis
	 */
	if (prop_0[2]==1) vCopy(incPolX_0,ex);
	else if (prop_0[2]==-1) vMultScal(-1,incPolX_0,ex);
	else {
		ex[0]=prop_0[0];
		ex[1]=prop_0[1];
		ex[2]=0;
		vNormalize(ex);
	}
	MatrVec(beta_matr,ex,exSP);
	// if needed rotate beam center
	if (beam_asym) MatrVec(beta_matr,beam_center_0,beam_center);
}

//======================================================================================================================

static void ReadLineStart(FILE *  restrict file,                  // opened file
	                      const char * restrict fname,            // ... its filename
                          char * restrict buf,const int buf_size, // buffer for line and its size
                          const char * restrict start)            // beginning of the line to search
// reads the first line that starts with 'start'
{
	while (!feof(file)) {
		fgets(buf,buf_size,file);
		if (strstr(buf,start)==buf) { // if correct beginning
			if (strstr(buf,"\n")==NULL && !feof(file))
				LogError(ONE_POS,"Buffer overflow while reading '%s' (size of essential line > %d)",fname,buf_size-1);
			else return; // line found and fits into buffer
		} // finish reading unmatched line
		else while (strstr(buf,"\n")==NULL && !feof(file)) fgets(buf,buf_size,file);
	}
	LogError(ONE_POS,"String '%s' is not found (in correct place) in file '%s'",start,fname);
}

//======================================================================================================================

static inline void ScanDouble(FILE * restrict file,const char * restrict fname,char * restrict buf,const int buf_size,
	const char * restrict start,double *res)
/* scans double value from a line starting with exactly 'start'; contains the same arguments as ReadLineStart function,
 * plus pointer to where the result should be placed
 */
{
	ReadLineStart(file,fname,buf,buf_size,start);
	if (sscanf(buf+strlen(start),"%lf",res)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
}

//======================================================================================================================

static inline void ScanInt(FILE * restrict file,const char * restrict fname,char * restrict buf,const int buf_size,
	const char * restrict start,int *res)
/* scans integer value from a line starting with exactly 'start'; contains the same arguments as ReadLineStart function,
 * plus pointer to where the result should be placed
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

//======================================================================================================================

static inline void ScanSizet(FILE * restrict file,const char * restrict fname,char * restrict buf,const int buf_size,
	const char * restrict start,size_t *res)
/* scans large integer value from a line starting with exactly 'start'; contains the same arguments as ReadLineStart
 * function, plus pointer to where the result should be placed. MinGW already provides C99-compliant printf-style
 * functions, but not yet scanf-style. So we have to use workarounds instead of straightforward "%zu" format specifier.
 * TODO: change to %zu, when libmingwex will include scanf.
 */
{
	double tmp;
	unsigned long res_tmp;

	ReadLineStart(file,fname,buf,buf_size,start);
	if (sscanf(buf+strlen(start),"%lf",&tmp)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
	if (tmp<0 || tmp>SIZE_MAX) LogError(ONE_POS,"Value after '%s' in file '%s' is out of size_t bounds",start,fname);
	if (sscanf(buf+strlen(start),"%lu",&res_tmp)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
	*res=(size_t)res_tmp;
}

//======================================================================================================================

static inline void ScanString(FILE * restrict file,const char * restrict fname,char * restrict buf,const int buf_size,
	const char * restrict start,char * restrict res)
/* scans string value from a line starting with exactly 'start'; contains the same arguments as ReadLineStart function,
 * plus pointer to where the result should be placed; the memory allocated to 'res' should be at least buf_size and
 * independent of buf.
 */
{
	ReadLineStart(file,fname,buf,buf_size,start);
	if (sscanf(buf+strlen(start),"%s",res)!=1)
		LogError(ONE_POS,"Error reading value after '%s' in file '%s'",start,fname);
	/* More secure would be to put field width in format string above (like "%.Ns"), however this field width is
	 * defined by the variable buf_size. The latter can only be implemented by a preliminary printf to get a format
	 * string.
	 */
}

//======================================================================================================================

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
		if (a->min>a->max) LogError(ONE_POS,
			"Wrong range (min="GFORMDEF", max="GFORMDEF") in file %s (max must be >= min)",a->min,a->max,fname);
		if (b->Jmax<b->Jmin)
			LogError(ONE_POS,"Wrong Jmax (%d) in file %s; it must be >= Jmin (%d)",b->Jmax,fname,b->Jmin);
		if (b->Jmin<1) LogError(ONE_POS,"Wrong Jmin (%d) in file %s (must be >=1)",b->Jmin,fname);
		if (b->eps<0) LogError(ONE_POS,"Wrong eps ("GFORMDEF") in file %s (must be >=0)",b->eps,fname);
		if (b->Jmax >= (int)(8*sizeof(int)))
			LogError(ONE_POS,"Too large Jmax(%d) in file %s, it will cause integer overflow",b->Jmax,fname);

		a->N=b->Grid_size=(1 << b->Jmax) + 1;
		if (b->equival && a->N>1) (a->N)--;
	}
	// initialize points of integration
	MALLOC_VECTOR(a->val,double,a->N,ALL);
	memory += a->N*sizeof(double);

	if (ifcos) { // make equal intervals in cos(angle)
		// consistency check
		if (a->min<0) LogError(ONE_POS,"Wrong min ("GFORMDEF") in file %s (must be >=0 for this angle)",a->min,fname);
		if (a->max>180)
			LogError(ONE_POS,"Wrong max ("GFORMDEF") in file %s (must be <=180 for this angle)",a->max,fname);
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

//======================================================================================================================

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
		if (a->min>a->max) LogError(ONE_POS,
			"Wrong range (min="GFORMDEF", max="GFORMDEF") in file %s (max must be >= min)",a->min,a->max,fname);
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
			if (strstr(buf,"\n")==NULL  && !feof(file))
				LogError(ONE_POS,"Buffer overflow while scanning lines in file '%s' (line size > %d)",fname,buf_size-1);
			if (sscanf(buf,"%lf\n",a->val+i)!=1)
				LogError(ONE_POS,"Failed scanning values from line '%s' in file '%s'",buf,fname);
		}
		out=AS_VALUES;
	}
	else LogError(ONE_POS,"Unknown type '%s' in file '%s'",temp,fname);
	return out;
}

//======================================================================================================================

void ReadAvgParms(const char * restrict fname)
// read parameters of orientation averaging from a file
{
	FILE * restrict input;
	char buf[BUF_LINE],temp[BUF_LINE];

	TIME_TYPE tstart=GET_TIME();
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
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================

void ReadAlldirParms(const char * restrict fname)
/* read integration parameters for asymmetry-parameter & C_sca; should not be used together with orientation averaging
 * because they use the same storage space - parms
 */
{
	FILE * restrict input;
	char buf[BUF_LINE],temp[BUF_LINE];

	TIME_TYPE tstart=GET_TIME();
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
		"theta: from "GFORMDEF" to "GFORMDEF" in (up to) %zu steps (equally spaced in cosine values)\n"
		"phi: from "GFORMDEF" to "GFORMDEF" in (up to) %zu steps\n"
		"see files 'log_int_***' for details\n\n",
		theta_int.min,theta_int.max,theta_int.N,phi_int.min,phi_int.max,phi_int.N);
	D("ReadAlldirParms finished");
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================

void ReadScatGridParms(const char * restrict fname)
// read parameters of the grid on which to calculate scattered field
{
	FILE * restrict input;
	char buf[BUF_LINE],temp[BUF_LINE];
	enum angleset theta_type,phi_type;
	size_t i;

	TIME_TYPE tstart=GET_TIME();
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
		if (phi_integr) LogError(ONE_POS,"Integration over phi can't be done with 'global_type=pairs'");
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
			if (strstr(buf,"\n")==NULL && !feof(input))
				LogError(ONE_POS,"Buffer overflow while scanning lines in file '%s' (line size > %d)",fname,BUF_LINE-1);
			if (sscanf(buf,"%lf %lf\n",angles.theta.val+i,angles.phi.val+i)!=2)
				LogError(ONE_POS,"Failed scanning values from line '%s' in file '%s'",buf,fname);
		}
	}
	else LogError(ONE_POS,"Unknown global_type '%s' in file '%s'",temp,fname);
	// close file
	FCloseErr(input,fname,ALL_POS);
	// print info
	if (IFROOT) {
		fprintf(logfile,"\nScattered field is calculated for multiple directions\n");
		switch (angles.type) {
			case SG_GRID:
				switch (theta_type) {
					case AS_RANGE:
						fprintf(logfile,"theta: from "GFORMDEF" to "GFORMDEF" in %zu steps\n",angles.theta.min,
							angles.theta.max,angles.theta.N);
						break;
					case AS_VALUES: fprintf(logfile,"theta: %zu given values\n",angles.theta.N); break;
				}
				switch (phi_type) {
					case AS_RANGE:
						fprintf(logfile,"phi: from "GFORMDEF" to "GFORMDEF" in %zu steps\n",angles.phi.min,
							angles.phi.max,angles.phi.N);
						if (phi_integr) fprintf(logfile,"(Mueller matrix is integrated over phi)\n");
						break;
					case AS_VALUES: fprintf(logfile,"phi: %zu given values\n",angles.phi.N); break;
				}
				break;
			case SG_PAIRS: fprintf(logfile,"Total %zu given (theta,phi) pairs\n",angles.N); break;
		}
		fprintf(logfile,"\n");
	}
	D("ReadScatGridParms finished");
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================

static void CalcFieldFree(doublecomplex ebuff[static restrict 3], // where to write calculated scattering amplitude
                          const double n[static restrict 3])      // scattering direction
/* Near-optimal routine to compute the scattered fields at one specific angle (more exactly - scattering amplitude);
 * Specific optimization are possible when e.g. n[0]=0 for scattering in yz-plane, however in this case it is very
 * improbable that the routine will become a bottleneck. The latter happens mostly for cases, when grid of scattering
 * angles is used with only small fraction of n, allowing simplifications.
 */
{
	double kkk;
	doublecomplex a,dpr;
	doublecomplex sum[3],tbuff[3],tmp=0; // redundant initialization to remove warnings
	int i;
	unsigned short ix,iy1,iy2,iz1,iz2;
	size_t j,jjj;
	double temp, na;
	doublecomplex mult_mat[MAX_NMAT];
	const bool scat_avg=true; // temporary fixed option for SO formulation
#ifdef SPARSE
	doublecomplex expX, expY, expZ;
#endif

	if (ScatRelation==SQ_SO) {
		// !!! this should never happen
		if (anisotropy) LogError(ONE_POS,"Incompatibility error in CalcField");
		// calculate correction coefficient
		if (scat_avg) na=0;
		else na=DotProd(n,prop);
		temp=kd*kd/24;
		// mult_mat=1-(kd^2/24)(m^2-2(n.a)m+1)
		for(i=0;i<Nmat;i++) mult_mat[i]=1-temp*(ref_index[i]*ref_index[i]-2*na*ref_index[i]+1);
	}
	cvInit(sum);
#ifndef SPARSE
	// prepare values of exponents, along each of the coordinates
	imExp_arr(-kd*n[0],boxX,expsX);
	imExp_arr(-kd*n[1],boxY,expsY);
	imExp_arr(-kd*n[2],local_Nz_unif,expsZ);
#endif // !SPARSE
	/* this piece of code tries to use that usually only x position changes from dipole to dipole, saving a complex
	 * multiplication seems to be beneficial, even considering bookkeeping overhead; it may not be as good for very
	 * porous particles though, but for them this part of code is anyway fast relative to the FFT on a large grid;
	 * Further optimization is possible using some kind of plans, i.e. by preliminary analyzing the position of the
	 * real dipoles on the grid.
	 */
	iy1=iz1=UNDEF;
	for (j=0;j<local_nvoid_Ndip;++j) {
		jjj=3*j;
		// a=exp(-ikr.n), but r is taken relative to the first dipole of the local box
		ix=position[jjj];
		iy2=position[jjj+1];
		iz2=position[jjj+2];
		// the second part is very improbable, but needed for robustness
		if (iy2!=iy1 || iz2!=iz1) {
			iy1=iy2;
			iz1=iz2;
#ifndef SPARSE // FFT mode
			tmp=expsY[iy2]*expsZ[iz2];
		}
		a=tmp*expsX[ix];
#else // sparse mode - the difference is that exponents are not precomputed
			expY=imExp(-kd*n[1]*iy2);
			expZ=imExp(-kd*n[2]*iz2);
			tmp=expY*expZ;
		}
		expX=imExp(-kd*n[0]*ix);
		a=tmp*expX;
#endif // SPARSE
		/* the following line may incur certain overhead (from 0% to 5% depending on tests).
		 * It is possible to remove this overhead by separating the complete loop for SQ_SO in a separate case (and it
		 * was like that at r1209). However, the code was much harder to read and maintain. Since there are several
		 * ideas that may speed up this calculation by a factor of a few times, we should not worry about 5%.
		 */
		if (ScatRelation==SQ_SO) a*=mult_mat[material[j]];
		// sum(P*exp(-ik*r.n))
		for(i=0;i<3;i++) sum[i]+=pvec[jjj+i]*a;
	} /* end for j */
	

	// tbuff=(I-nxn).sum=sum-n*(n.sum)
	dpr=crDotProd(sum,n);
	cvMultScal_RVec(dpr,n,tbuff);
	cvSubtr(sum,tbuff,tbuff);
	// ebuff=(-i*k^3)*exp(-ikr0.n)*tbuff, where r0=box_origin_unif
	a=imExp(-WaveNum*DotProd(box_origin_unif,n)); // a=exp(-ikr0.n)
	kkk=WaveNum*WaveNum*WaveNum;
	// the following additional multiplier implements IGT_SO
	if (ScatRelation==SQ_IGT_SO) kkk*=(1-kd*kd/24);
	tmp=-I*a*kkk; // tmp=(-i*k^3)*exp(-ikr0.n)
	cvMultScal_cmplx(tmp,tbuff,ebuff);
}

//======================================================================================================================

static void CalcFieldSurf(doublecomplex ebuff[static restrict 3], // where to write calculated scattering amplitude
                          const double nF[static restrict 3])     // scattering direction (at infinity)
/* Same as CalcFieldFree but for particle near surface.
 * For scattering into the substrate we employ the reciprocity principle. The scattered field is obtained from field of
 * the plane wave incoming from the scattered direction at the dipole position. In particular,
 * E_sca(s,p) = eF(s,p)*(k_0^2/r)*exp(ikr)*t'(s,p) * Sum[P_j.eN_(s,P)*exp(-i*k_0*nN.r_j)],
 * where eF,eN are unit [e.e=1] vectors at far and near-field, nN is the normalized transmitted k-vector (also nN.nN=1).
 * t' is transmittance coefficient from substrate into the vacuum, k is wavevector in the substrate.
 * Total scattered field is obtained by summing s and p components. The actual computed quantity is scattering
 * amplitude F, defined as E_sca = F*exp(ikr)/(-ikr).
 * Reciprocity should be valid even for absorbing substrate (with any symmetric tensor), not depending on the symmetry
 * of the refractive index of the particle itself. In principle, the same formula can be obtained by transmitting
 * cylindrical waves emitted by dipole, but that needs additional coefficient due to stretching of wavefront during
 * transmission (similar to the difference between amplitude and intensity transmission coefficients).
 * Moreover, the consideration of specific scattering angle (and wavenumber) implies that we completely ignore the
 * surface plasmon polaritons (SPP), which can, in principle, be considered as scattering at 90 degrees. These SPPs may
 * be important for energy balance for metallic substrates. However, for such substrates energy balance is not perfect
 * anyway due to absorption.
 */
{
	doublecomplex aF,aN,phSh;
	doublecomplex cs,cp; // coefficients (reflectance or transmittance) for s- and p-polarizations
	doublecomplex ki,kt; // normal component of wavevector above and below the surface
	doublecomplex sumF[3],sumN[3],t3[3],tmpF=0,tmpN=0; // redundant initialization to remove warnings
	doublecomplex nN[3]; // scattering direction (n.n=1) at near field (corresponds to ktVec in GenerateB.c)
	double epF[3],es[3]; // unit vectors of s- and p-polarization, ep differs for near- and far-field
	doublecomplex epN[3]; // ep at near-field can be complex
	int i;
	unsigned short ix,iy1,iy2,iz1,iz2;
	size_t j,jjj;
#ifdef SPARSE
	doublecomplex expX, expY, expZ;
#endif

	const bool above=(nF[2]>-ROUND_ERR); // we assume above-the-surface scattering for all boundary cases (like 90 deg)
	// Using SQ_SO for particles near surface seems even beyond "under development"
	if (ScatRelation==SQ_SO) LogError(ONE_POS,"Incompatibility error in CalcFieldSurf");
	cvInit(sumN);
	if (above) cvInit(sumF); //additional storage for directly propagated scattering

	/* There is an inherent discontinuity for msub approaching 1 and scattering angle 90 degrees (nF[2]=0). The problem
	 * is that for m=1+-0, but |m-1|>>(nF[2])^2, ki<<kt<<1 => rs=rp=-1
	 * while for m=1 (exactly) the limit of nF[2]->0 results in kt=ki => rs=rp=0
	 * Therefore, below is a certain logic, which behaves in an intuitively expected way, for common special cases.
	 * However, it is still not expected to be continuous for fine-changing parameters (like msub approaching 1).
	 * In particular, the jump occurs when msub crosses 1+-ROUND_ERR boundary.
	 * Still, the discontinuity should apply only to scattering at exactly 90 degrees, but not to, e.g., integral
	 * quantities, like Csca (if sufficient large number of integration points is chosen).
	 */
	// calculate nN, ki, kt, cs, cp, and phSh
	if (above) { // simple reflection
		/* No scattering at exactly 90 degrees for non-trivial surface (to avoid randomness for this case).
		 * See A. Small, J. Fung, and V.N. Manoharan, “Generalization of the optical theorem for light scattering from
		 * a particle at a planar interface,” J. Opt. Soc. Am. A 30, 2519–2525 (2013) for theoretical discussion of
		 * this fact.
		 */
		if (fabs(nF[2])<ROUND_ERR && cabs(msub-1)>ROUND_ERR) {
			cvInit(ebuff);
			return;
		}
		cvBuildRe(nF,nN);
		nN[2]*=-1;
		ki=nF[2];
		if (msubInf) {
			cs=-1;
			cp=1;
		}
		  // since kt is not further needed, we directly calculate cs and cp (equivalent to kt=ki)
		else if (cabs(msub-1)<ROUND_ERR && fabs(ki)<SQRT_RND_ERR) cs=cp=0;
		else { // no special treatment here, since other cases, including 90deg-scattering, are taken care above.
			kt=cSqrtCut(msub*msub - (nN[0]*nN[0]+nN[1]*nN[1]));
			cs=FresnelRS(ki,kt);
			cp=FresnelRP(ki,kt,msub);
		}
		phSh=imExp(2*WaveNum*hsub*ki);
	}
	else { // transmission; here nF[2] is negative
		// formulae correspond to plane wave incoming from below, but with change ki<->kt
		if (msubInf) { // no transmission for perfectly reflecting substrate => zero result
			cvInit(ebuff);
			return;
		}
		kt=-msub*nF[2];
		if (cabs(msub-1)<ROUND_ERR && fabs(kt)<SQRT_RND_ERR) ki=kt;
		else ki=cSqrtCut(1 - msub*msub*(nF[0]*nF[0]+nF[1]*nF[1]));
		// here nN may be complex, but normalized to n.n=1
		nN[0]=msub*nF[0];
		nN[1]=msub*nF[1];
		nN[2]=-ki;
		// these formulae works fine for ki=kt (even very small), and ki=kt=0 is impossible here
		cs=FresnelTS(kt,ki);
		cp=FresnelTP(kt,ki,1/msub);
		// coefficient comes from  k0->k in definition of F(n) (in denominator)
		phSh=msub*cexp(I*WaveNum*hsub*(ki-kt));
	}
#ifndef SPARSE
	// prepare values of exponents, along each of the coordinates
	imExp_arr(-kd*nN[0],boxX,expsX);
	imExp_arr(-kd*nN[1],boxY,expsY);
	imExp_arr(-kd*nN[2],local_Nz_unif,expsZ);
#endif // !SPARSE
	/* this piece of code tries to use that usually only x position changes from dipole to dipole, saving a complex
	 * multiplication seems to be beneficial, even considering bookkeeping overhead; it may not be as good for very
	 * porous particles though, but for them this part of code is anyway fast relative to the FFT on a large grid;
	 * Further optimization is possible using some kind of plans, i.e. by preliminary analyzing the position of the
	 * real dipoles on the grid.
	 */
	iy1=iz1=UNDEF;
	if (above) for (j=0;j<local_nvoid_Ndip;++j) { // two sums need to be calculated
		jjj=3*j;
		// a=exp(-ikr.n), but r is taken relative to the first dipole of the local box
		ix=position[jjj];
		iy2=position[jjj+1];
		iz2=position[jjj+2];
		// the second part is very improbable, but needed for robustness
		if (iy2!=iy1 || iz2!=iz1) {
			iy1=iy2;
			iz1=iz2;
#ifndef SPARSE // FFT mode
			tmpN=expsY[iy2]*expsZ[iz2];
			tmpF=expsY[iy2]*conj(expsZ[iz2]);
		}
		aN=tmpN*expsX[ix];
		aF=tmpF*expsX[ix];
#else // sparse mode - the difference is that exponents are not precomputed; cexp is used since argument can be complex
			expY=cexp(-I*kd*nN[1]*iy2);
			expZ=cexp(-I*kd*nN[2]*iz2);
			tmpN=expY*expZ;
			tmpF=expY*conj(expZ);
		}
		expX=cexp(-I*kd*nN[0]*ix);
		aN=tmpN*expX;
		aF=tmpF*expX;
#endif // SPARSE
		// sum(P*exp(-ik*r.nN,F))
		for(i=0;i<3;i++) {
			sumN[i]+=pvec[jjj+i]*aN;
			sumF[i]+=pvec[jjj+i]*aF;
		}
	} /* end for j above surface */
	else for (j=0;j<local_nvoid_Ndip;++j) { // below surface, single sum - similar to free-space scattering
		jjj=3*j;
		// a=exp(-ikr.n), but r is taken relative to the first dipole of the local box
		ix=position[jjj];
		iy2=position[jjj+1];
		iz2=position[jjj+2];
		// the second part is very improbable, but needed for robustness
		if (iy2!=iy1 || iz2!=iz1) {
			iy1=iy2;
			iz1=iz2;
#ifndef SPARSE // FFT mode
			tmpN=expsY[iy2]*expsZ[iz2];
		}
		aN=tmpN*expsX[ix];
#else // sparse mode - the difference is that exponents are not precomputed; cexp is used since argument can be complex
			expY=cexp(-I*kd*nN[1]*iy2);
			expZ=cexp(-I*kd*nN[2]*iz2);
			tmpN=expY*expZ;
		}
		expX=cexp(-I*kd*nN[0]*ix);
		aN=tmpN*expX;
#endif // SPARSE
		// sum(P*exp(-ik*r.nN))
		for(i=0;i<3;i++) sumN[i]+=pvec[jjj+i]*aN;
	} /* end for j below surface */
	// Reflected or transmitted light phSh*(Rs*es(es.sumN) + Rp*epF(epN.sumN)), [dot product w/o conjugation]
	/* If reciprocal configuration is rigorously considered signs of vectors nF and nN should be changed, along with
	 * either es or ep. However, such sign change would not change the final result.
	 */
		// set unit vectors for s- and p-polarizations; es is the same for all cases
	if (vAlongZ(nF)) { // special case - es=ey
		es[0]=0;
		es[1]=1;
		es[2]=0;
	}
	else { // general case: es = ez x nF /||...||; here we implicitly use that lab RF coincides with particle RF
		CrossProd(ezLab,nF,es);
		vNormalize(es);
	}
	CrossProd(es,nF,epF); // epF = es x nF
	crCrossProd(es,nN,epN); // epN = es x nN (complex)
	  // finalize fields
	cvMultScal_RVec(phSh*cs*crDotProd(sumN,es),es,ebuff);
	cvMultScal_RVec(phSh*cp*cDotProd_conj(sumN,epN),epF,t3);
	cvAdd(t3,ebuff,ebuff);
	// add directly scattered light, when above the surface (phase shift due to main direction being nN)
	if (above) { // ebuff+= [(I-nxn).sum=sum-nF*(nF.sum)] * exp(-2ik*r0*nz), where r0=box_origin_unif
		cvMultScal_RVec(crDotProd(sumF,nF),nF,t3);
		cvSubtr(sumF,t3,t3);
		cvMultScal_cmplx(imExp(-2*WaveNum*ki*box_origin_unif[2]),t3,t3);
		cvAdd(t3,ebuff,ebuff);
	}
	// ebuff=(-i*k^3)*exp(-ikr0.n)*tbuff, where r0=box_origin_unif
	// All m-scaling for substrate has been accounted in phSh above
	doublecomplex sc=-I*WaveNum*WaveNum*WaveNum*cexp(-I*WaveNum*crDotProd(nN,box_origin_unif));
	// the following additional multiplier implements IGT_SO
	if (ScatRelation==SQ_IGT_SO) sc*=(1-kd*kd/24);
	cvMultScal_cmplx(sc,ebuff,ebuff);
}

//======================================================================================================================*/

void CalcField(doublecomplex ebuff[static restrict 3], // where to write calculated scattering amplitude
               const double n[static restrict 3])      // scattering direction
// wrapper, which redirects the calculation of the field to one of two functions
{
	if (surface) CalcFieldSurf(ebuff,n);
	else CalcFieldFree(ebuff,n);
}

//======================================================================================================================

double ExtCross(const double * restrict incPol)
// Calculate the Extinction cross-section
{
	doublecomplex ebuff[3];
	double sum;
	size_t i;

	if (beamtype==B_PLANE && !surface) {
		CalcField (ebuff,prop);
		sum=crDotProd_Re(ebuff,incPol); // incPol is real, so no conjugate is needed
		MyInnerProduct(&sum,double_type,1,&Timing_ScatQuanComm);
		sum*=FOUR_PI/(WaveNum*WaveNum);
	}
	/* more general formula; normalization is done assuming the unity amplitude of the electric field in the focal point
	 * of the beam; It does not comply with ScatRelation SO. So SO is, effectively, replaced by DRAINE when calculating
	 * Cext for non-plane beams.
	 */
	else {
		sum=0;
		for (i=0;i<local_nvoid_Ndip;++i) sum+=cDotProd_Im(pvec+3*i,Einc+3*i); // sum{Im(P.E_inc*)}
		MyInnerProduct(&sum,double_type,1,&Timing_ScatQuanComm);
		sum*=FOUR_PI*WaveNum;
		/* Surprisingly, this little trick is enough to satisfy IGT_SO, because this factor is applied in CalcField()
		 * and is independent of propagation or scattering direction. Thus it can be applied to any linear combination
		 * of plane waves, i.e. any field.
		 *
		 * Unfortunately, the same reasoning fails for SO of full IGT, because there the correction factor does
		 * (slightly) depend on the propagation direction.
		 */
		if (ScatRelation==SQ_IGT_SO) sum*=(1-kd*kd/24);
	}
	/* TO ADD NEW BEAM
	 * The formulae above works only if the amplitude of the beam is unity at the focal point. Either make sure that new
	 * beam satisfies this condition or add another case here with different formulae.
	*/
	if (surface) sum*=inc_scale;
	return sum;
}

//======================================================================================================================

double AbsCross(void)
// Calculate the Absorption cross-section for process 0
{
	size_t dip,index;
	int i,j;
	unsigned char mat;
	double sum,temp1,temp2;
	doublecomplex m,m2m1;
	double mult[MAX_NMAT][3]; // multiplier (possibly anisotropic)
	double mult1[MAX_NMAT];   // multiplier, which is always isotropic

	// Cabs = 4*pi*sum
	/* In this function IGT_SO is equivalent to DRAINE. It may seem more logical to make IGT_SO same as FINDIP. However,
	 * the result is different only for LDR (and similar), for which using IGT does not make a lot of sense anyway.
	 * Overall, peculiar details related to optical theorem warrant a further study.
	 */
	switch (ScatRelation) {
		/* code below is applicable only for diagonal (for some cases - possibly anisotropic) polarizability and should
		 * be rewritten otherwise
		 */
		case SQ_IGT_SO:
		case SQ_DRAINE:
			/* based on Eq.(35) from Yurkin and Hoekstra, "The discrete dipole approximation: an overview and recent
			 * developments," JQSRT 106:558-589 (2007).
			 * summand: Im(P.Eexc(*))-(2/3)k^3*|P|^2=|P|^2*(-Im(1/cc)-(2/3)k^3)
			 */
			temp1 = 2*WaveNum*WaveNum*WaveNum/3;
			for (i=0;i<Nmat;i++) for (j=0;j<3;j++) mult[i][j]=-cimag(1/cc[i][j])-temp1;
			for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip) {
				mat=material[dip];
				index=3*dip;
				for(i=0;i<3;i++) sum+=mult[mat][i]*cAbs2(pvec[index+i]);
			}
			break;
		case SQ_FINDIP:
			/* based on Eq.(31) or equivalently Eq.(58) from the same paper (ref. above)
			 * summand: Im(P.E(*))=-|P|^2*Im(chi_inv), chi_inv=1/(V*chi)
			 */
			temp1 = 2*WaveNum*WaveNum*WaveNum/3;
			for (i=0;i<Nmat;i++) for (j=0;j<3;j++) mult[i][j]=-cimag(chi_inv[i][j]);
			for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip) {
				mat=material[dip];
				index=3*dip;
				for(i=0;i<3;i++) sum+=mult[mat][i]*cAbs2(pvec[index+i]);
			}
			break;
		case SQ_SO:
			// !!! the following should never happen
			if (anisotropy) LogError(ONE_POS,"Incompatibility error in AbsCross");
			// calculate mult1
			temp1=kd*kd/6;
			temp2=FOUR_PI/dipvol;
			for (i=0;i<Nmat;i++) {
				m=ref_index[i];
				m2m1=m*m-1;
				// mult1=-Im(1/chi)*(1+(kd*Im(m))^2)/d^3;  chi=(m^2-1)/(4*PI)
				mult1[i]=temp2*cimag(m2m1)*(1+temp1*cimag(m)*cimag(m))/cAbs2(m2m1);
			}
			// main cycle
			for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip) sum+=mult1[material[dip]]*cvNorm2(pvec+3*dip);
			break;
	}
	MyInnerProduct(&sum,double_type,1,&Timing_ScatQuanComm);
	if (surface) sum*=inc_scale;
	return FOUR_PI*WaveNum*sum;
}

//======================================================================================================================

double DecayCross(void)
// computes total cross section for the dipole incident field; similar to Cext
// 4pi*k*Im[p0(*).Escat(r0)]
{
	double sum;
	size_t i;

	/* This is a correct expression only _if_ exciting p0 is real, then
	 * (using G(r1,r2) = G(r2,r1)^T, valid also for surface)
	 * p0(*).Escat_i(r0) = p0(*).G_0i.p_i = p_i.G_i0.p0(*) = p_i.G_i0.p0 = p_i.Einc_i
	 * => Im(p0(*).Escat(r0)) = sum{Im(P.E_inc)}
	 *
	 * For complex p0 an efficient calculation strategy (not to waste evaluations of interaction) is to compute an array
	 * of G_i0.p0(*) together with Einc and use it here afterwards.
	 */
	sum=0;
	for (i=0;i<local_nvoid_Ndip;++i) sum+=cimag(cDotProd_conj(pvec+3*i,Einc+3*i)); // sum{Im(P.E_inc)}
	MyInnerProduct(&sum,double_type,1,&Timing_ScatQuanComm);
	return FOUR_PI*WaveNum*sum;
}

//======================================================================================================================

void SetScatPlane(const double ct,const double st,const double phi,double robs[static restrict 3],
	double polPer[static restrict 3])
/* Given theta (cos,sin) and phi, calculates scattering direction and unit vector perpendicular to the scattering plane.
 * Currently, two alternative definitions are used: theta and phi are either angles in the laboratory reference frame,
 * or in the beam reference frame. Generally, polPer = robs x prop (normalized to unit vector), but cases when the two
 * latter vectors are aligned requires special treatment.
 *
 * For special cases with prop||ez we assume that the result is continuous for theta changing from 0 to pi. In other
 * words th=0 and pi is effective replaced by th=+0 and pi-0 respectively. Then phi determines the scattering plane.
 *
 * For other special cases (prop not aligned with ez) the scattering plane is taken to include ez (and hence IncPolX).
 *
 * Overall, we spend a lot of effort to make consistent definition of polPer. However, its sign is irrelevant because
 * change of sign (simultaneously for polPer and both incident and scattering polPar) doesn't change the amplitude
 * (and Mueller) matrix.
 */
{
	double cphi,sphi;

	cphi=cos(phi);
	sphi=sin(phi);
	if (surface) { // definitions are in the laboratory reference frame
		// robs = cos(theta)*ezLab + sin(theta)*[cos(phi)*exLab + sin(phi)*eyLab];
		LinComb(exLab,eyLab,cphi,sphi,robs);
		LinComb(ezLab,robs,ct,st,robs);
	}
	else { // definitions are with respect to the incident beam
		// robs = cos(theta)*prop + sin(theta)*[cos(phi)*incPolX + sin(phi)*incPolY];
		LinComb(incPolX,incPolY,cphi,sphi,robs);
		LinComb(prop,robs,ct,st,robs);
	}
	// set unit vectors which defines Eper
	// two special cases; testing ct is more accurate than st==0
	// 1) robs has +0 component along exLab (for phi=0): incPolper = (+-)(sin(phi)*exLab - cos(phi)*eyLab)
	if (surface && propAlongZ && fabs(st)<ROUND_ERR) LinComb(exLab,eyLab,sphi*propAlongZ,-cphi*propAlongZ,polPer);
	// 2) robs has +0 component along incPolX (for phi=0): incPolper = sin(phi)*incPolX - cos(phi)*incPolY
	else if (!surface && fabs(st)<ROUND_ERR) LinComb(incPolX,incPolY,sphi,-cphi,polPer);
	else { // general case: incPolPer = robserver x prop /||...||
		CrossProd(robs,prop,polPer);
		// when robs || prop scattering plane contains ez, then polPer = prop x ez/||...|| = -incPolY
		// currently this may only occur for surface, since otherwise it is handled by a special case above
		if (DotProd(polPer,polPer)==0) vMultScal(-1,incPolY,polPer);
		else vNormalize(polPer);
	}
}
//======================================================================================================================

void CalcAlldir(void)
// calculate scattered field in many directions
{
	int index,npoints,point;
	size_t i,j;
	TIME_TYPE tstart;
	double robserver[3],incPolpar[3],incPolper[3],cthet,sthet,th,ph;
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
			/* set robserver and unit vector for Eper (determines scattering plane); actually Eper and Epar are
			 * irrelevant, since we are interested only in |E|^2. But projecting the vector on two axes helps to
			 * somewhat decrease communication time. We also may need these components in the future.
			 */
			SetScatPlane(cthet,sthet,ph,robserver,incPolper);
			// set unit vector for Epar
			CrossProd(robserver,incPolper,incPolpar);
			// calculate scattered field - main bottleneck
			CalcField(ebuff,robserver);
			/* Set Epar and Eper - use separate E_ad array to store them (to decrease communications in 1.5 times).
			 * Writing a special case for sequential mode can eliminate the need of E_ad altogether. Moreover,
			 * E2_alldir can be stored in 1/4 of memory allocated for E_ad. However, we do not do it, because it doesn't
			 * seem so significant. And, more importantly, complex fields may also be useful in the future, e.g.
			 * for radiation force calculation through integration of the far-field
			 */
			index=2*point;
			E_ad[index]=crDotProd(ebuff,incPolper);
			E_ad[index+1]=crDotProd(ebuff,incPolpar);
			point++;
			// show progress
			if (((10*point)%npoints)<10 && IFROOT) printf(" %d%%",100*point/npoints);
		}
	}
	// accumulate fields
	Accumulate(E_ad,cmplx_type,2*npoints,&Timing_EFieldADComm);
	// calculate square of the field
	for (point=0;point<npoints;point++) E2_alldir[point] = cAbs2(E_ad[2*point]) + cAbs2(E_ad[2*point+1]);
	/* when below surface we scale E2 by Re(1/msub) in accordance with formula for the Poynting vector (and factor of
	 * k_sca^2). After that Csca (and g) computed using the standard formula should correctly describe the energy and
	 * momentum balance for any (even complex) msub (since the scattered wave is homogeneous at far-field), but doesn't
	 * include energy or momentum absorbed (obtained) by the medium.
	 */
	if (surface && !msubInf) {
		double scale=creal(1/msub); // == Re(msub)/|msub|^2
		for (i=0;i<theta_int.N;++i) if (TestBelowDeg(theta_int.val[i]))
			for (j=0,point=phi_int.N*i;j<phi_int.N;j++,point++) E2_alldir[point]*=scale;
	}
	if (IFROOT) printf("  done\n");
	// timing
	Timing_EFieldAD = GET_TIME() - tstart;
	Timing_EField += Timing_EFieldAD;
}

//======================================================================================================================

void CalcScatGrid(const enum incpol which)
// calculate scattered field in many directions
{
	size_t i,j,n,point,index;
	TIME_TYPE tstart;
	double robserver[3],incPolpar[3],incPolper[3],cthet,sthet,th,ph;
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
			// set robserver and unit vector for Eper (determines scattering plane)
			SetScatPlane(cthet,sthet,ph,robserver,incPolper);
			// set unit vector for Epar
			CrossProd(robserver,incPolper,incPolpar);
			// calculate scattered field - main bottleneck
			CalcField(ebuff,robserver);
			// set Epar and Eper - use Egrid array to store them (to decrease communications in 1.5 times)
			index=2*point;
			Egrid[index]=crDotProd(ebuff,incPolper);
			Egrid[index+1]=crDotProd(ebuff,incPolpar);
			point++;
			// show progress; the value is always from 0 to 100, so conversion to int is safe
			if (((10*point)%angles.N)<10 && IFROOT) printf(" %d%%",(int)(100*point/angles.N));
		}
	}
	// accumulate fields; timing
	Accumulate(Egrid,cmplx_type,2*angles.N,&Timing_EFieldSGComm);
	if (IFROOT) printf("  done\n");
	Timing_EFieldSG = GET_TIME() - tstart;
	Timing_EField += Timing_EFieldSG;
}

//======================================================================================================================

static double CscaIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating Csca
{
	res[0]=E2_alldir[AlldirIndex(theta,phi)];

	return 0;
}

//======================================================================================================================

double ScaCross(const char *f_suf)
// Calculate the scattering cross section from the integral
{
	TIME_TYPE tstart;
	char fname[MAX_FNAME];
	double res;

	SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_LOG_INT_CSCA "%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,CscaIntegrand,1,&res,fname);
	res*=FOUR_PI/(WaveNum*WaveNum);
	if (surface) res*=inc_scale;
	Timing_Integration += GET_TIME() - tstart;
	return res;
}

//======================================================================================================================

static double gIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g
/* We use definition that g is n, averaged with weight dC/dOmega. For free-space it is proportional to momentum of
 * scattered field, but not in surface mode. For the latter to get momentum, it should be additionally weighted by
 * msca(n) (and probably normalized differently). But it doesn't make a lot of sense, because part of the momentum is
 * transferred to the surface (and hard to calculate)
 */
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

//======================================================================================================================

static void AsymParm(double *vec,const char *f_suf) ATT_UNUSED;
static void AsymParm(double *vec,const char *f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca
{
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM "%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gIntegrand,3,vec,log_int);
	vMultScal(FOUR_PI/(WaveNum*WaveNum),vec,vec);
	if (surface) vMultScal(inc_scale,vec,vec);
	Timing_Integration += GET_TIME() - tstart;
}

//======================================================================================================================

static double gxIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g_x, see also gIntegrand()
{
	double th=theta_int.val[theta];

	/* a separate case is used to avoid negligibly small non-zero results, which further cause unusually large relative
	 * errors in integration log
	 */
	if (th==180) res[0]=0;
	else res[0]=E2_alldir[AlldirIndex(theta,phi)]*sin(Deg2Rad(th))*cos(Deg2Rad(phi_int.val[phi]));
	return 0;
}

//======================================================================================================================

void AsymParm_x(double *vec,const char *f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca
{
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM F_LOG_X"%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gxIntegrand,1,vec,log_int);
	vec[0] *= FOUR_PI/(WaveNum*WaveNum);
	if (surface) vec[0]*=inc_scale;
	Timing_Integration += GET_TIME() - tstart;
}

//======================================================================================================================

static double gyIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g_y, see also gIntegrand()
{
	double th=theta_int.val[theta];

	/* a separate case is used to avoid negligibly small non-zero results, which further cause unusually large relative
	 * errors in integration log
	 */
	if (th==180) res[0]=0;
	else res[0]=E2_alldir[AlldirIndex(theta,phi)]*sin(Deg2Rad(th))*sin(Deg2Rad(phi_int.val[phi]));
	return 0;
}

//======================================================================================================================

void AsymParm_y(double *vec,const char *f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca
{
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM F_LOG_Y"%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gyIntegrand,1,vec,log_int);
	vec[0] *= FOUR_PI/(WaveNum*WaveNum);
	if (surface) vec[0]*=inc_scale;
	Timing_Integration += GET_TIME() - tstart;
}

//======================================================================================================================

static double gzIntegrand(const int theta,const int phi,double * restrict res)
// function that is transferred to integration module when calculating g_z, see also gIntegrand()
{
	double th=Deg2Rad(theta_int.val[theta]);
	res[0]=E2_alldir[AlldirIndex(theta,phi)]*cos(th);
	return 0;
}

//======================================================================================================================

void AsymParm_z(double *vec,const char *f_suf)
// Calculate the unnormalized asymmetry parameter, i.e. not yet normalized by Csca

{
	TIME_TYPE tstart;
	char log_int[MAX_FNAME];

	SnprintfErr(ONE_POS,log_int,MAX_FNAME,"%s/"F_LOG_INT_ASYM F_LOG_Z"%s",directory,f_suf);

	tstart = GET_TIME();
	Romberg2D(parms,gzIntegrand,1,vec,log_int);
	vec[0] *= FOUR_PI/(WaveNum*WaveNum);
	if (surface) vec[0]*=inc_scale;
	Timing_Integration += GET_TIME() - tstart;
}

//======================================================================================================================

void Frp_mat(double Finc_tot[static restrict 3],double Fsca_tot[static restrict 3],
	double * restrict Frp)
/* Calculate the Radiation Pressure (separately incident and scattering part by direct calculation of the scattering
 * force. The total force per dipole is calculated as intermediate results. It is saved to Frp, if the latter is not
 * NULL. mem denotes the specific memory allocated before function call
 *
 * This should be completely rewritten to work through FFT. Moreover, it should comply with '-scat ...' command line
 * option.
 */
{
	size_t j,k,jg,comp;
	double * restrict rdipT;
	doublecomplex * restrict pT;
	doublecomplex temp;
	double Fsca[3],Finc[3];
	double *vec;
	double r,r2; // (squared) absolute distance
	double n[3]; // unit vector in the direction of r_{jl}
	doublecomplex
	a,ab1,ab2,c1[3],c2[3], // see chapter ...
	x_cg[3], // complex conjugate P*_j
	Pn_j,    // n_jk.P_k
	Pn_k,    // P*_j.n_jk
	inp;     // P*_j.P_k
	size_t mem=0; // memory count

	// initialize
	vInit(Fsca_tot);
	vInit(Finc_tot);
	// Calculate incoming force per dipole
	if (Frp==NULL) vec=Finc;
	else mem+=sizeof(double)*local_nRows; // memory allocated before for Frp
	/* The following expression F_inc=k(v)*0.5*Sum(P.Einc(*)) is valid only for the plane wave
	 * TODO: Implement formulae for arbitrary Gaussian beams
	 */
	for (j=0;j<local_nRows;j+=3) {
		if (Frp!=NULL) vec=Frp+j;
		vMultScal(WaveNum*cDotProd_Im(pvec+j,Einc+j)/2,prop,vec);
		vAdd(vec,Finc_tot,Finc_tot);
	}
	// check if it can work at all; check is redundant for sequential mode
	size_t nRows=MultOverflow(3,nvoid_Ndip,ONE_POS_FUNC);
#ifdef PARALLEL
	/* Because of the parallelization by row-block decomposition the distributed arrays involved need to be gathered on
	 * each node a) DipoleCoord -> rdipT; b) pvec -> pT. Actually this routine is usually called for two polarizations
	 * and rdipT does not change between the calls. So one AllGather of rdipT can be removed. Number of memory
	 * allocations can also be reduced. But this should be replaced by Fourier anyway.
	 */
	/* The following is somewhat redundant in sparse mode, since "full" (containing information about all dipoles)
	 * vectors are already present in that mode. However, we do not optimize it now, since in standard mode radiation
	 * forces should be computed by FFT anyway. Moreover, there are certain ideas to optimize sparse mode, so it will
	 * not use full vectors - if done, this improvement can be also adjusted to the code below.
	 */
	// allocates a lot of additional memory
	MALLOC_VECTOR(rdipT,double,nRows,ALL);
	MALLOC_VECTOR(pT,complex,nRows,ALL);
	mem+=nRows*(sizeof(double)+sizeof(doublecomplex));
	// this is approximate value, but not far
	if (IFROOT)
		PrintBoth(logfile,"Additional memory usage for radiation forces (per processor): "FFORMM" MB\n",mem/MBYTE);
	// gathers everything
	AllGather(DipoleCoord,rdipT,double3_type,&Timing_ScatQuanComm);
	AllGather(pvec,pT,cmplx3_type,&Timing_ScatQuanComm);
#else
	pT=pvec;
	rdipT=DipoleCoord;
	if (mem!=0) PrintBoth(logfile,"Additional memory usage for radiation forces: "FFORMM" MB\n",mem/MBYTE);
#endif
	// Calculate scattering force per dipole
	/* Currently, testing the correctness of the following is very hard because the original code lacks comments. So the
	 * best we can do before rewriting it completely is to test that it produces reasonable results for a number of test
	 * cases.
	 */
	for (j=0,jg=3*local_nvoid_d0;j<local_nRows;j+=3,jg+=3) {
		vInit(Fsca);
		for (k=0;k<nRows;k+=3) if (jg!=k) {
			// Set distance related variables
			vSubtr(rdipT+jg,rdipT+k,n);
			r2=DotProd(n,n);
			r = sqrt(r2);
			vMultScal(1/r,n,n);
			// Set the scalar products a.b1 and a.b2
			a = imExp(WaveNum*r);
			ab1 = (3/(r2*r2) - I*3*WaveNum/(r*r2) - WaveNum*WaveNum/r2)*a;
			ab2 = (-WaveNum*WaveNum/r2 + I*WaveNum*WaveNum*WaveNum/r)*a;
			// Prepare c1 and c2
			vConj(pT+jg,x_cg);
			Pn_j=crDotProd(x_cg,n);
			Pn_k=crDotProd(pT+k,n);
			inp=cDotProd(pT+k,pT+jg);
			temp=-Pn_j*Pn_k;
			// Set c1
			cvMultScal_RVec(5*temp+inp,n,c1);
			cvLinComb1_cmplx(pT+k,c1,Pn_j,c1);
			cvLinComb1_cmplx(x_cg,c1,Pn_k,c1);
			// Set c2
			cvMultScal_RVec(temp+inp,n,c2);
			// Fsca_{jk} = ...
			for (comp=0;comp<3;++comp) Fsca[comp] += creal(c1[comp]*ab1 + c2[comp]*ab2)/2;
		} // end k-loop
		// Concluding
		vAdd(Fsca,Fsca_tot,Fsca_tot);
		if (Frp!=NULL) vAdd(Fsca,Frp+j,Frp+j);
	} // end j-loop

	// Accumulate the total forces on all nodes
	MyInnerProduct(Finc_tot,double_type,3,&Timing_ScatQuanComm);
	MyInnerProduct(Fsca_tot,double_type,3,&Timing_ScatQuanComm);
#ifdef PARALLEL
	Free_general(rdipT);
	Free_cVector(pT);
#endif
}
