/* File: iterative.c
 * $Date::                            $
 * Descr: a few iterative techniques to solve DDA equations
 *
 *        The linear system is composed so that diagonal terms are equal to 1, therefore use of Jacobi preconditioners
 *        does not have any effect.
 *
 *        CS methods still converge to the right result even when matrix is slightly non-symmetric (e.g. -int so),
 *        however they do it much slowly than usually. It is recommended then to use BiCGStab or BCGS2.
 *
 * Copyright (C) 2006-2015 ADDA contributors
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
#include "debug.h"
#include "io.h"
#include "linalg.h"
#include "memory.h"
#include "timing.h"
#include "vars.h"
// system headers
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // for time_t & time

#ifdef OCL_BLAS
#	include "oclcore.h"
#	include <clBLAS.h> //external library
#endif

// SEMI-GLOBAL VARIABLES

// defined and initialized in calculator.c
extern doublecomplex *rvec; // can't be declared restrict due to SwapPointers
extern doublecomplex * restrict vec1,* restrict vec2,* restrict vec3,* restrict vec4,* restrict Avecbuffer;
// defined and initialized in fft.c
#if !defined(OPENCL) && !defined(SPARSE)
extern doublecomplex * restrict Xmatrix; // used as storage for arrays in WKB init field
#endif
// defined and initialized in param.c
extern const double iter_eps;
extern const enum init_field InitField;
extern const char *infi_fnameY,*infi_fnameX;
extern const bool recalc_resid;
extern const enum chpoint chp_type;
extern const time_t chp_time;
extern const char *chp_dir;
// defined and initialized in timing.c
extern time_t last_chp_wt;
extern TIME_TYPE Timing_OneIter,Timing_OneIterComm,Timing_InitIter,Timing_InitIterComm,Timing_IntFieldOneComm,
	Timing_MVP,Timing_MVPComm,Timing_OneIterMVP,Timing_OneIterMVPComm;
extern size_t TotalIter;

// LOCAL VARIABLES

#define RESID_STRING "RE_%03d = "EFORM // string containing residual value
#define FFORM_PROG "% .6f"  // format for progress value

static double inprodR;     // used as |r_0|^2 and best squared norm of residual up to some iteration
static double inprodRp1;   // used as |r_k+1|^2 and squared norm of current residual
static double epsB;        // stopping criterion
static double resid_scale; // scale to get square of relative error
static double prev_err;    // previous relative error; used in ProgressReport, initialized in IterativeSolver
static int ind_m;          // index of iterative method
static int niter;          // iteration count
static int counter;        // number of successive iterations without residual decrease
static bool chp_exit;      // checkpoint occurred - exit
static bool complete;      // complete iteration was performed (not stopped in the middle)
	// whether matrix-vector product computed during initialization can be reused at first iteration
static bool matvec_ready;
typedef struct // data for checkpoints
{
	void *ptr; // pointer to the data
	int size;  // size of one element
} chp_data;
chp_data * restrict scalars,* restrict vectors;
enum phase {
	PHASE_VARS, // Initialization of variables, and linking them to scalars and vectors
	PHASE_INIT, // Initialization of iterations
	PHASE_ITER  // Each iteration
};
struct iter_params_struct {
	enum iter meth;   // identifier
	int mc;           // maximum allowed number of iterations without residual decrease
	int sc_N;         // number of additional scalars to describe the state
	int vec_N;        // number of additional vectors to describe the state
	void (*func)(const enum phase); // pointer to implementation of the iterative solver
};
static doublecomplex dumb ATT_UNUSED; // dumb variable, used in workaround for issue 146

#define ITER_FUNC(name) static void name(const enum phase ph)

ITER_FUNC(BCGS2);
ITER_FUNC(BiCG_CS);
ITER_FUNC(BiCGStab);
ITER_FUNC(CGNR);
ITER_FUNC(CSYM);
ITER_FUNC(QMR_CS);
ITER_FUNC(QMR_CS_2);
/* TO ADD NEW ITERATIVE SOLVER
 * Add the line to this list in the alphabetical order, analogous to the ones already present. The variable part is the
 * name of the function, implementing the method. The macros expands to a function prototype.
 */

static const struct iter_params_struct params[]={
	{IT_BCGS2,15000,2,1,BCGS2},
	{IT_BICG_CS,50000,1,0,BiCG_CS},
	{IT_BICGSTAB,30000,3,3,BiCGStab},
	{IT_CGNR,10,1,0,CGNR},
	{IT_CSYM,10,6,2,CSYM},
	{IT_QMR_CS,50000,8,3,QMR_CS},
	{IT_QMR_CS_2,50000,5,2,QMR_CS_2}
	/* TO ADD NEW ITERATIVE SOLVER
	 * Add its parameters to this list in the alphabetical order. The parameters, in order of appearance, are identifier
	 * (specified in const.h), maximum allowed number of iterations without the residual decrease, numbers of additional
	 * scalars and vectors to describe the state of the iterative solver (see comment before function SaveCheckpoint),
	 * and name of a function, implementing the method.
	 */
};

// EXTERNAL FUNCTIONS

// matvec.c
void MatVec(doublecomplex * restrict in,doublecomplex * restrict out,double * inprod,bool her,TIME_TYPE *timing,
	TIME_TYPE *comm_timing);

//======================================================================================================================

static void MatVec_wrapper(doublecomplex * restrict in,doublecomplex * restrict out,double * inprod,bool her,
	TIME_TYPE *timing,TIME_TYPE *comm_timing)
/* fuction wrapper for MatVec to be called within the iterative solver if the solver is able to use clBLAS, i.e.
 * the host and GPU memory does not have to be synchronized. Currently it is only used in the BiCG solver.
 */
{
#ifdef OCL_BLAS
	bufupload=false;
#endif
	MatVec(in,out,inprod,her,timing,comm_timing);
#ifdef OCL_BLAS
	bufupload=true;
#endif
}

//======================================================================================================================

static inline void SwapPointers(doublecomplex **a,doublecomplex **b)
/* swap two pointers of (doublecomplex *) type; should work for others but will give "Suspicious pointer conversion"
 * warning. While this is a convenient function that can save some copying between memory blocks, it doesn't allow
 * consistent usage of 'restrict' keyword for affected pointers. This may hamper some optimizations. Hopefully, the most
 * important optimizations are those in the linalg.c, which can be improved by using 'restrict' keyword in the functions
 * themselves.
 */
{
	doublecomplex *tmp;

	tmp=*a;
	*a=*b;
	*b=tmp;
}

//======================================================================================================================

/* Checkpoint systems saves the current state of the iterative solver to the file. By default (for every iterative
 * solver) a number of scalars and vectors are saved. The scalars include, among others, inprodR. There are 3 default
 * vectors: xvec, rvec, pvec (Avecbuffer is _not_ saved). If the iterative solver requires any other scalars or vectors
 * to describe its state, this information should be specified in structure arrays 'scalars' and 'vectors'.
 */

static void SaveIterChpoint(void)
/* save a binary checkpoint; only limitedly foolproof - user should take care to load checkpoints on the same machine
 * (number of processors) and with the same command line.
 */
{
	int i;
	char fname[MAX_FNAME];
	FILE * restrict chp_file;
	TIME_TYPE tstart;

	tstart=GET_TIME();
	if (IFROOT) {
		// create directory "chp_dir" if needed and open info file
		SnprintfErr(ONE_POS,fname,MAX_FNAME,"%s/"F_CHP_LOG,chp_dir);
		if ((chp_file=fopen(fname,"w"))==NULL) {
			MkDirErr(chp_dir,ONE_POS);
			chp_file=FOpenErr(fname,"w",ONE_POS);
		}
		// write info and close file
		fprintf(chp_file,"Info about the run, which produced the checkpoint, can be found in ../%s",directory);
		FCloseErr(chp_file,fname,ONE_POS);
	}
	// wait to ensure that directory exists
	Synchronize();
	// open output file; writing errors are checked only for vectors
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/"F_CHP,chp_dir,ringid);
	chp_file=FOpenErr(fname,"wb",ALL_POS);
	// write common scalars
	fwrite(&ind_m,sizeof(int),1,chp_file);
	fwrite(&local_nRows,sizeof(size_t),1,chp_file);
	fwrite(&niter,sizeof(int),1,chp_file);
	fwrite(&counter,sizeof(int),1,chp_file);
	fwrite(&inprodR,sizeof(double),1,chp_file);
	fwrite(&prev_err,sizeof(double),1,chp_file); // written on ALL processors but used only on root
	fwrite(&resid_scale,sizeof(double),1,chp_file);
	// write specific scalars
	for (i=0;i<params[ind_m].sc_N;i++) fwrite(scalars[i].ptr,scalars[i].size,1,chp_file);
	// write common vectors
	if (fwrite(xvec,sizeof(doublecomplex),local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed writing to file '%s'",fname);
	if (fwrite(rvec,sizeof(doublecomplex),local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed writing to file '%s'",fname);
	if (fwrite(pvec,sizeof(doublecomplex),local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed writing to file '%s'",fname);
	// write specific vectors
	for (i=0;i<params[ind_m].vec_N;i++) if (fwrite(vectors[i].ptr,vectors[i].size,local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed writing to file '%s'",fname);
	// close file
	FCloseErr(chp_file,fname,ALL_POS);
	// write info to logfile after everyone is finished
	Synchronize();
	if (IFROOT) PrintBoth(logfile,"Checkpoint (iteration) saved\n");
	Timing_FileIO+=GET_TIME()-tstart;
	Synchronize(); // this is to ensure that message above appears if and only if OK
}

//======================================================================================================================

static void LoadIterChpoint(void)
/* load a binary checkpoint; only limitedly foolproof - user should take care to load checkpoints on the same machine
 * (number of processors) and with the same command line.
 */
{
	int i;
	int ind_m_new;
	size_t local_nRows_new;
	char fname[MAX_FNAME],ch;
	FILE * restrict chp_file;
	TIME_TYPE tstart;

	tstart=GET_TIME();
	// open input file; reading errors are checked only for vectors
	SnprintfErr(ALL_POS,fname,MAX_FNAME,"%s/"F_CHP,chp_dir,ringid);
	chp_file=FOpenErr(fname,"rb",ALL_POS);
	/* check for consistency. This implies that the same index corresponds to the same iterative solver in list params.
	 * So if the ADDA executable was changed, e.g. by adding a new iterative solver, between writing and reading
	 * checkpoint, this test may fail.
	 */
	fread(&ind_m_new,sizeof(int),1,chp_file);
	if (ind_m_new!=ind_m) LogError(ALL_POS,"File '%s' is for different iterative method",fname);
	fread(&local_nRows_new,sizeof(size_t),1,chp_file);
	if (local_nRows_new!=local_nRows) LogError(ALL_POS,"File '%s' is for different vector size",fname);
	// read common scalars
	fread(&niter,sizeof(int),1,chp_file);
	fread(&counter,sizeof(int),1,chp_file);
	fread(&inprodR,sizeof(double),1,chp_file);
	fread(&prev_err,sizeof(double),1,chp_file); // read on ALL processors but used only on root
	fread(&resid_scale,sizeof(double),1,chp_file);
	// read specific scalars
	for (i=0;i<params[ind_m].sc_N;i++) fread(scalars[i].ptr,scalars[i].size,1,chp_file);
	// read common vectors
	if (fread(xvec,sizeof(doublecomplex),local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed reading from file '%s'",fname);
	if (fread(rvec,sizeof(doublecomplex),local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed reading from file '%s'",fname);
	if (fread(pvec,sizeof(doublecomplex),local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed reading from file '%s'",fname);
	// read specific vectors
	for (i=0;i<params[ind_m].vec_N;i++) if (fread(vectors[i].ptr,vectors[i].size,local_nRows,chp_file)!=local_nRows)
		LogError(ALL_POS,"Failed reading from file '%s'",fname);
	// check if EOF reached and close file
	if(fread(&ch,1,1,chp_file)!=0) LogError(ALL_POS,"File '%s' is too long",fname);
	FCloseErr(chp_file,fname,ALL_POS);
	// initialize auxiliary variables
	epsB=iter_eps*iter_eps/resid_scale;
	// print info
	if (IFROOT) {
		PrintBoth(logfile,"Checkpoint (iteration) loaded\n");
		// if residual is stagnating print info about last minimum
		if (counter!=0) fprintf(logfile,"Residual has been stagnating already for %d iterations since:\n"
			RESID_STRING"\n...\n",counter,niter-counter-1,sqrt(resid_scale*inprodR));
	}
	Timing_FileIO+=GET_TIME()-tstart;
}

//======================================================================================================================

static void ProgressReport(void)
// Do common procedures; show progress in logfile and stdout; also check for checkpoint condition
{
	double err,progr,elapsed;
	char progr_string[MAX_LINE];
	const char *temp;
	time_t wt;

	if (inprodRp1<=inprodR) {
		inprodR=inprodRp1;
		counter=0;
	}
	else counter++;
	if (IFROOT) {
		err=sqrt(resid_scale*inprodRp1);
		progr=1-err/prev_err;
		if (counter==0) temp="+ ";
		else if (progr>0) temp="-+";
		else temp="- ";
		SnprintfErr(ONE_POS,progr_string,MAX_LINE,RESID_STRING"  %s",niter,err,temp);
		if (!orient_avg) fprintf(logfile,"%s  progress ="FFORM_PROG"\n",progr_string,progr);
		printf("%s\n",progr_string);
		prev_err=err;
	}
	niter++;
	TotalIter++;
	// check condition for checkpoint; checkpoint is saved at first time
	if (chp_type!=CHP_NONE && chp_time!=UNDEF && complete) {
		time(&wt);
		elapsed=difftime(wt,last_chp_wt);
		if (chp_time<elapsed) {
			SaveIterChpoint();
			time(&last_chp_wt);
			if (chp_type!=CHP_REGULAR) chp_exit=true;
		}
	}
}

//======================================================================================================================

static double ResidualNorm2(doublecomplex * restrict x,doublecomplex * restrict r,doublecomplex * restrict buffer,
	TIME_TYPE *mvp_timing,TIME_TYPE *mvp_comm_timing,TIME_TYPE *comm_timing)
/* Computes ||Ax-b||^2, where b=sqrt(C).Einc; buffer is used for Ax, r contains Ax-b at the end; comm_timing is
 * incremented with communication time. If only the norm is required, the calculation can be done without using vector
 * r, but this does not make a lot of sense, since memory is allocated anyway.
 */
{
	double res;

	TIME_TYPE mc_time=0;
	MatVec(x,buffer,NULL,false,mvp_timing,&mc_time);
	(*mvp_comm_timing) += mc_time;
	(*comm_timing) += mc_time;
	nMult_mat(r,Einc,cc_sqrt);
	nDecrem(r,buffer,&res,comm_timing);
	return res;
}

//======================================================================================================================

ITER_FUNC(BCGS2)
/* Enhanced Bi-CGStab(2) method.
 * Based on the code by M.A. Botchev and D.R. Fokkema - http://www.math.uu.nl/people/vorst/zbcg2.f90 and
 * D. R. Fokkema, "Enhanced implementation of BiCGstab(l) for solving linear systems of equations," Preprint 976,
 * Department of Mathematics, Utrecht University (1996).
 *
 * "Reliable update part" was removed, since tests using '-recalc_resid' show that it is (almost) never needed.
 *
 * For l=1, the method is equivalent to BiCGStab, rewritten through 2-term recurrences (as QMR2 is equivalent to QMR),
 * so we use l=2 here. In many cases one iteration of this method is similar to two iterations of BiCGStab, but overall
 * convergence is slightly better.
 * Breakdown tests were made to coincide with that for BiCGStab for l=1.
 *
 * !!! This iterative solver produces segmentation fault when compiled with icc 11.1. Probably that is related to
 * issue 146. But we leave it be (assume that this is a compiler bug). Even if someone uses this compiler, he can
 * live fine without this iterative solver.
 */
{
#define LL 2 // potentially the method will also work for l=1 (but memory allocation and freeing need to be adjusted)
#define EPS1 1E-10 // for 1/|beta|
#define EPS2 1E-10 // for |u_j+1.r~|/|r_j.r~|
	static doublecomplex * restrict r[LL+1],* restrict u[LL+1];
	static doublecomplex matrix_z[LL+1][LL+1],y0[LL+1],yl[LL+1],zy0[LL+1],zyl[LL+1];
	static int i,j;
	static doublecomplex alpha,beta,omega,rho0,rho1,sigma,varrho,hatgamma,temp1;
	static double kappa0,kappal,dtmp;

	switch (ph) {
		case PHASE_VARS:
			/* rename some vectors; this doesn't contradict with 'restrict' keyword, since new names are not used
			 * together with old names
			 */
			r[0]=rvec;
			r[1]=vec1;
			u[0]=vec2;
			u[1]=Avecbuffer;
			if (LL==2) {
				r[2]=vec3;
				u[2]=vec4;
			}
			// initialize data structure for checkpoints
			scalars[0].ptr=&rho0;
			scalars[1].ptr=&alpha;
			scalars[0].size=scalars[1].size=sizeof(doublecomplex);
			vectors[0].ptr=vec2; // u[0]
			vectors[0].size=sizeof(doublecomplex);
			return;
		case PHASE_INIT:
			if (!load_chpoint) {
				nCopy(pvec,rvec); // (pvec = r~0) = r0
				rho0=-1;
			}
			return;
		case PHASE_ITER:
			// --- The BiCG part ---
			for (j=0;j<LL;j++) {
				rho1=nDotProd(r[j],pvec,&Timing_OneIterComm); // rho1 = r_j.r~0
				// u_i = r_i - beta*u_i
				if (niter==1 && j==0) nCopy(u[0],r[0]);
				else {
					// test for zero rho0 (1/beta)
					dtmp=cabs(rho0)/(cabs(rho1)*cabs(alpha)); // assume that rho1 is not exactly zero
					Dz("1/|beta|="GFORM_DEBUG,dtmp);
					if (dtmp<EPS1) LogError(ONE_POS,"BCGS2 fails: 1/|beta| is too small ("GFORM_DEBUG").",dtmp);
					beta=alpha*rho1/rho0;
					// u_i = r_i - beta*u_i
					temp1=-beta;
					for (i=0;i<=j;i++) nIncrem10_cmplx(u[i],r[i],temp1,NULL,NULL);
				}
				rho0=rho1;
				// u_j+1 = A.u_j
				if (niter==1 && j==0 && matvec_ready) {} // do nothing; u[1]<=>Avecbuffer already contains matvec result
				else MatVec(u[j],u[j+1],NULL,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
				sigma=nDotProd(u[j+1],pvec,&Timing_OneIterComm); // sigma = u_j+1.r~0
				// test for zero sigma (1/alpha)
				dtmp=cabs(sigma)/cabs(rho1); // assume that rho1 is not exactly zero
				Dz("|u_%d.r~|/|r_%d.r~|="GFORM_DEBUG,j+1,j,dtmp);
				if (dtmp<EPS1)
					LogError(ONE_POS,"BCGS2 fails: |u_%d.r~|/|r_%d.r~| is too small ("GFORM_DEBUG").",j+1,j,dtmp);
				alpha = rho1/sigma;
				nIncrem01_cmplx(xvec,u[0],alpha,NULL,NULL); // x = x + alpha*u_0
				// r_i = r_i - alpha*u_i+1
				temp1=-alpha;
				for (i=0;i<=j;i++) nIncrem01_cmplx(r[i],u[i+1],temp1,NULL,NULL);
				MatVec(r[j],r[j+1],NULL,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			}
			// --- The convex polynomial part ---
			// Z = R'R
			for(i=0;i<=LL;i++) for (j=0;j<=i;j++) {
				matrix_z[i][j]=nDotProd(r[j],r[i],&Timing_OneIterComm);
				if (i!=j) matrix_z[j][i]=conj(matrix_z[i][j]);
			}
			// small vectors y0 and yl
			y0[0]=-1;
			if (LL==2) y0[1]=matrix_z[1][0]/matrix_z[1][1]; // works only for l<=2
			y0[LL]=0;
			yl[0]=0;
			if (LL==2) yl[1]=matrix_z[1][2]/matrix_z[1][1]; // works only for l<=2
			yl[LL]=-1;
			//  Convex combination
			// compute Zy0 and Zyl
			for (i=0;i<=LL;i++) {
				zy0[i]=zyl[i]=0;
				for (j=0;j<=LL;j++) {
					zy0[i]+=matrix_z[i][j]*y0[j];
					zyl[i]+=matrix_z[i][j]*yl[j];
				}
			}
			// kappa0 = sqrt(y0.Zy0); kappal = sqrt(yl.Zyl); employs that dot products are always real
			dtmp=0;
			for (i=0;i<=LL;i++) dtmp+=creal(zy0[i]*conj(y0[i]));
			kappa0=sqrt(dtmp);
			dtmp=0;
			for (i=0;i<=LL;i++) dtmp+=creal(zyl[i]*conj(yl[i]));
			kappal=sqrt(dtmp);
			// varrho = Zy0.yl/(kappa0*kappal)
			varrho=0;
			for (i=0;i<=LL;i++) varrho+=zy0[i]*conj(yl[i]);
			varrho/=kappa0*kappal;
			// hatgamma = varrho/abs(varrho) * max( abs(varrho),0.7)
			dtmp=cabs(varrho);
			hatgamma=varrho*MAX(dtmp,0.7)/dtmp;
			// y0 = y0 - (hatgamma*kappa0/kappal)*yl
			temp1=-hatgamma*kappa0/kappal;
			for (i=0;i<=LL;i++) y0[i]+=temp1*yl[i];
			// Update
			omega = y0[LL];
			for (i=1;i<=LL;i++) {
				temp1=-y0[i];
				nIncrem01_cmplx(u[0],u[i],temp1,NULL,NULL);   // u_0 = u_0 - y0[i]*u_i
				nIncrem01_cmplx(xvec,r[i-1],y0[i],NULL,NULL); // x = x + y0[i]*r_i-1
				nIncrem01_cmplx(r[0],r[i],temp1,NULL,NULL);   // r_0 = r_0 - y0[i]*r_i
			}
			// y0 has changed; compute Zy0 once more
			for (i=0;i<=LL;i++) {
				zy0[i]=0;
				for (j=0;j<=LL;j++) zy0[i]+=matrix_z[i][j]*y0[j];
			}
			// |r|^2 = y0.Zy0
			inprodRp1=0;
			for (i=0;i<=LL;i++) inprodRp1+=creal(zy0[i]*conj(y0[i]));
			// rho0 = -omega*rho0; moved from the beginning of the iteration
			rho0*=-omega;
			return; // end of PHASE_ITER
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
}
#undef LL
#undef EPS1
#undef EPS2

//======================================================================================================================

ITER_FUNC(BiCG_CS)
/* Bi-Conjugate Gradient for Complex Symmetric systems, based on:
 * Freund R.W. "Conjugate gradient-type methods for linear systems with complex symmetric coefficient matrices",
 * SIAM Journal of Scientific Statistics and Computation, 13(1):425-448,1992.
 *
 * it is also identical to COCG, described in:
 * van der Vorst H.A., Melissen J.B.M. "A Petrov-Galerkin type method for solving Ax=b, where A is symmetric complex",
 * IEEE Transactions on Magnetics, 26(2):706-708, 1990.
 */
{
#define EPS1 1E-10 // for (rT.r)/(r.r)
#define EPS2 1E-10 // for (pT.A.p)/(rT.r)
	static doublecomplex alpha, mu;
	static doublecomplex beta,ro_new,ro_old,temp;
	static double dtmp,abs_ro_new;
#ifdef OCL_BLAS
	cl_mem bufro_new;
	cl_mem bufmu;
	cl_mem bufinprodRp1;
	cl_mem bufpvec=bufargvec;
	cl_mem bufAvecbuffer=bufresultvec;
#endif

	switch (ph) {
		case PHASE_VARS:
			scalars[0].ptr=&ro_old;
			scalars[0].size=sizeof(doublecomplex);
			return;
		case PHASE_INIT: 
#ifdef OCL_BLAS
			D("clblasSetup started");
			CL_CH_ERR(clblasSetup());
#	ifdef DEBUGFULL
			cl_uint major,minor,patch;
			CL_CH_ERR(clblasGetVersion(&major,&minor,&patch));
			D("clBLAS library version - %u.%u.%u",major,minor,patch);
#	endif
			CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufpvec,CL_FALSE,0,sizeof(doublecomplex)*local_nRows,pvec,0,
				NULL,NULL));
			CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufrvec,CL_FALSE,0,sizeof(doublecomplex)*local_nRows,rvec,0,
				NULL,NULL));
			CL_CH_ERR(clEnqueueWriteBuffer(command_queue,bufxvec,CL_FALSE,0,sizeof(doublecomplex)*local_nRows,xvec,0,
				NULL,NULL));
#endif
			return; // no specific initialization required
		case PHASE_ITER:
#ifdef OCL_BLAS
			/* TODO: Initialization of this two and one other scalar buffers (and then their release) happens at each
			 * iteration. This is not a bottleneck, but still seems redundant. But to solve this problem we probably
			 * need to add additional phase of the iterative solver, like PHASE_RELEASE, where all such buffers can be
			 * released.
			 */
			CREATE_CL_BUFFER(bufro_new,CL_MEM_READ_WRITE,sizeof(doublecomplex),NULL);
			CREATE_CL_BUFFER(bufmu,CL_MEM_READ_WRITE,sizeof(doublecomplex),NULL);
			CL_CH_ERR(clblasZdotu(local_nRows,bufro_new,0,bufrvec,0,1,bufrvec,0,1,buftmp,1,&command_queue,0,NULL,
				NULL));
			CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufro_new,CL_TRUE,0,sizeof(doublecomplex),&ro_new,0,NULL,NULL));
#else
			ro_new=nDotProdSelf_conj(rvec,&Timing_OneIterComm);
#endif
			// ro_k-1=r_k-1(*).r_k-1; check for ro_k-1!=0
			abs_ro_new=cabs(ro_new);
			dtmp=abs_ro_new/inprodR;
			Dz("|rT.r|/(r.r)="GFORM_DEBUG,dtmp);
			if (dtmp<EPS1) LogError(ONE_POS,"BiCG_CS fails: |rT.r|/(r.r) is too small ("GFORM_DEBUG").",dtmp);
			if (niter==1) {
#ifdef OCL_BLAS
				clEnqueueCopyBuffer(command_queue,bufrvec,bufpvec,0,0,sizeof(doublecomplex)*local_nRows,0,NULL,NULL);
#else
				nCopy(pvec,rvec); // p_1=r_0
#endif 
			}
			else {
				// beta_k-1=ro_k-1/ro_k-2
				beta=ro_new/ro_old;
				// p_k=beta_k-1*p_k-1+r_k-1
#ifdef OCL_BLAS
				cl_double2 clbeta = {.s={creal(beta),cimag(beta)}};
				CL_CH_ERR(clblasZscal(local_nRows,clbeta,bufpvec,0,1,1,&command_queue,0,NULL,NULL));
				cl_double2 clunit = {.s={1,0}};
				CL_CH_ERR(clblasZaxpy(local_nRows,clunit,bufrvec,0,1,bufpvec,0,1,1,&command_queue,0,NULL,NULL));
#else
				nIncrem10_cmplx(pvec,rvec,beta,NULL,NULL);
#endif
			}
			// q_k=Avecbuffer=A.p_k
			if (niter==1 && matvec_ready) {} // do nothing, Avecbuffer is ready to use
			else MatVec_wrapper(pvec,Avecbuffer,NULL,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			// mu_k=p_k.q_k; check for mu_k!=0
#ifdef OCL_BLAS
			CL_CH_ERR(clblasZdotu(local_nRows,bufmu,0,bufpvec,0,1,bufAvecbuffer,0,1,buftmp,1,&command_queue,0,NULL,
				NULL));
			CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufmu,CL_TRUE,0,sizeof(doublecomplex),&mu,0,NULL,NULL));
#else
			mu=nDotProd_conj(pvec,Avecbuffer,&Timing_OneIterComm);
#endif
			dtmp=cabs(mu)/abs_ro_new;
			Dz("|pT.A.p|/(rT.r)="GFORM_DEBUG,dtmp);
			if (dtmp<EPS2) LogError(ONE_POS,"BiCG_CS fails: |pT.A.p|/(rT.r) is too small ("GFORM_DEBUG").",dtmp);
			// alpha_k=ro_k/mu_k
			alpha=ro_new/mu;
			// x_k=x_k-1+alpha_k*p_k
#ifdef OCL_BLAS
			cl_double2 clalpha = {.s={creal(alpha),cimag(alpha)}};
			CL_CH_ERR(clblasZaxpy(local_nRows,clalpha,bufpvec,0,1,bufxvec,0,1,1,&command_queue,0,NULL,NULL));
#else
			nIncrem01_cmplx(xvec,pvec,alpha,NULL,NULL);
#endif
			// r_k=r_k-1-alpha_k*A.p_k and |r_k|^2
			temp=-alpha;
#ifdef OCL_BLAS
			cl_double2 cltemp = {.s={creal(temp),cimag(temp)}};
			CREATE_CL_BUFFER(bufinprodRp1,CL_MEM_READ_WRITE,sizeof(double),NULL);
			CL_CH_ERR(clblasZaxpy(local_nRows,cltemp,bufAvecbuffer,0,1,bufrvec,0,1,1,&command_queue,0,NULL,NULL));
			CL_CH_ERR(clblasDznrm2(local_nRows,bufinprodRp1,0,bufrvec,0,1,buftmp,1,&command_queue,0,NULL,NULL));
			CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufinprodRp1,CL_TRUE,0,sizeof(double),&inprodRp1,0,NULL,NULL));
			inprodRp1=inprodRp1*inprodRp1;
#else
			nIncrem01_cmplx(rvec,Avecbuffer,temp,&inprodRp1,&Timing_OneIterComm);
#endif
#ifdef OCL_BLAS
			my_clReleaseBuffer(bufinprodRp1);
			my_clReleaseBuffer(bufro_new);
			my_clReleaseBuffer(bufmu);
			if (inprodRp1<epsB){
				CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufpvec,CL_FALSE,0,sizeof(doublecomplex)*local_nRows,pvec,0,
					NULL,NULL));
				CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufrvec,CL_FALSE,0,sizeof(doublecomplex)*local_nRows,rvec,0,
					NULL,NULL));
				CL_CH_ERR(clEnqueueReadBuffer(command_queue,bufxvec,CL_TRUE,0,sizeof(doublecomplex)*local_nRows,xvec,0,
					NULL,NULL));
			}
#endif
			// initialize ro_old -> ro_k-2 for next iteration
			ro_old=ro_new;
			return; // end of PHASE_ITER
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
}
#undef EPS1
#undef EPS2

//======================================================================================================================

ITER_FUNC(BiCGStab)
/* Bi-Conjugate Gradient Stabilized, based on
 * "Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods",
 * http://www.netlib.org/templates/Templates.html .
 */
{
#define EPS1 1E-10 // for 1/|beta|
#define EPS2 1E-10 // for |v.r~|/|r.r~|
	static double denumOmega,dtmp;
	static doublecomplex beta,ro_new,ro_old,omega,alpha,temp1,temp2;
	static doublecomplex * restrict v,* restrict s,* restrict rtilda;

	switch (ph) {
		case PHASE_VARS:
			/* rename some vectors; this doesn't contradict with 'restrict' keyword, since new names are not used
			 * together with old names
			 */
			v=vec1;
			s=vec2;
			rtilda=vec3;
			// initialize data structure for checkpoints
			scalars[0].ptr=&ro_old;
			scalars[1].ptr=&omega;
			scalars[2].ptr=&alpha;
			scalars[0].size=scalars[1].size=scalars[2].size=sizeof(doublecomplex);
			vectors[0].ptr=vec1; // v
			vectors[1].ptr=vec2; // s
			vectors[2].ptr=vec3; // rtilda
			vectors[0].size=vectors[1].size=vectors[2].size=sizeof(doublecomplex);
			return;
		case PHASE_INIT:
			if (!load_chpoint) nCopy(rtilda,rvec); // r~=r_0
			return;
		case PHASE_ITER:
			// ro_k-1=r_k-1.r~ ; check for ro_k-1!=0
			ro_new=nDotProd(rvec,rtilda,&Timing_OneIterComm);
			if (niter==1) nCopy(pvec,rvec); // p_1=r_0
			else {
				// beta_k-1=(ro_k-1/ro_k-2)*(alpha_k-1/omega_k-1)
				temp1=ro_new*alpha;
				temp2=ro_old*omega;
				// check that omega_k-1!=0; assume that ro_new is not exactly zero
				dtmp=cabs(temp2)/cabs(temp1);
				Dz("1/|beta|="GFORM_DEBUG,dtmp);
				if (dtmp<EPS1) LogError(ONE_POS,"BiCGStab fails: 1/|beta| is too small ("GFORM_DEBUG").",dtmp);
				beta=temp1/temp2;
				// p_k=beta_k-1*(p_k-1-omega_k-1*v_k-1)+r_k-1
				temp1=-beta*omega;
				nIncrem110_cmplx(pvec,v,rvec,beta,temp1);
			}
			// calculate v_k=A.p_k
			if (niter==1 && matvec_ready) nCopy(v,Avecbuffer);
			else MatVec(pvec,v,NULL,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			// alpha_k=ro_new/(v_k.r~)
			temp1=nDotProd(v,rtilda,&Timing_OneIterComm);
			dtmp=cabs(temp1)/cabs(ro_new); // assume that ro_new is not exactly zero
			Dz("|v.r~|/|r.r~|="GFORM_DEBUG,dtmp);
			if (dtmp<EPS2) LogError(ONE_POS,"BiCGStab fails: |v.r~|/|r.r~| is too small ("GFORM_DEBUG").",dtmp);
			alpha=ro_new/temp1;
			// s=r_k-1-alpha*v_k-1
			temp1=-alpha;
			nLinComb1_cmplx(s,v,rvec,temp1,&inprodRp1,&Timing_OneIterComm);
			// check convergence at this step; if yes, checkpoint should not be saved afterwards
			if (inprodRp1<epsB && chp_type!=CHP_ALWAYS) {
				// x_k=x_k-1+alpha_k*p_k
				nIncrem01_cmplx(xvec,pvec,alpha,NULL,NULL);
				complete=false;
			}
			else {
				// t=Avecbuffer=A.s
				MatVec(s,Avecbuffer,&denumOmega,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
				// omega_k=s.t/|t|^2
				omega=nDotProd(s,Avecbuffer,&Timing_OneIterComm)/denumOmega;
				// x_k=x_k-1+alpha_k*p_k+omega_k*s
				nIncrem011_cmplx(xvec,pvec,s,alpha,omega);
				// r_k=s-omega_k*t and |r_k|^2
				temp1=-omega;
				nLinComb1_cmplx(rvec,Avecbuffer,s,temp1,&inprodRp1,&Timing_OneIterComm);
				// initialize ro_old -> ro_k-2 for next iteration
				ro_old=ro_new;
			}
			return; // end of PHASE_ITER
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
}
#undef EPS1
#undef EPS2

//======================================================================================================================

ITER_FUNC(CGNR)
/* Conjugate Gradient applied to Normalized Equations with minimization of Residual Norm, based on
 * "Templates for the Solution of Linear Systems: Building Blocks for Iterative Methods",
 * http://www.netlib.org/templates/Templates.html .
 */
{
	static double alpha, denumeratorAlpha;
	static double beta,ro_new,ro_old;

	switch (ph) {
		case PHASE_VARS:
			scalars[0].ptr=&ro_old;
			scalars[0].size=sizeof(double);
			return;
		case PHASE_INIT: return; // no specific initialization required
		case PHASE_ITER:
			// p_1=Ah.r_0 and ro_new=ro_0=|Ah.r_0|^2
			// since first product is with Ah , matvec_ready can't be employed
			if (niter==1) MatVec(rvec,pvec,&ro_new,true,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			else {
				// Avecbuffer=AH.r_k-1, ro_new=ro_k-1=|AH.r_k-1|^2
				MatVec(rvec,Avecbuffer,&ro_new,true,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
				// beta_k-1=ro_k-1/ro_k-2
				beta=ro_new/ro_old;
				// p_k=beta_k-1*p_k-1+AH.r_k-1
				nIncrem10(pvec,Avecbuffer,beta,NULL,NULL);
			}
			// alpha_k=ro_k-1/|A.p_k|^2
			// Avecbuffer=A.p_k
			MatVec(pvec,Avecbuffer,&denumeratorAlpha,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			alpha=ro_new/denumeratorAlpha;
			// x_k=x_k-1+alpha_k*p_k
			nIncrem01(xvec,pvec,alpha,NULL,NULL);
			// r_k=r_k-1-alpha_k*A.p_k and |r_k|^2
			nIncrem01(rvec,Avecbuffer,-alpha,&inprodRp1,&Timing_OneIterComm);
			// initialize ro_old -> ro_k-2 for next iteration
			ro_old=ro_new;
			return; // end of PHASE_ITER
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
}

//======================================================================================================================

ITER_FUNC(CSYM)
/* Bi-Conjugate Gradient for Complex Symmetric systems, based on
 * A. Bunse-Gerstner and R. Stover, "On a conjugate gradient-type method for solving complex symmetric linear systems,"
 * Lin. Alg. Appl. 287, 105-123 (1999). with rearrangement of operations (MatVec is now calculated in the beginning of
 * the iteration)
 *
 * Consumes one less vector than QMR-CS, because rvec does not need to be explicitly computed. The residual should
 * always decrease and always be smaller than that of CGNR (for the same number of matrix-vector products).
 */
{
	static doublecomplex alpha,gamma,invksi,theta,eta,tau,temp1,temp2,s_old,s_new;
	static double dtmp,beta,c_old,c_new;
	static doublecomplex *q_new,*q_old,*p_new,*p_old; // can't be declared restrict due to SwapPointers

	switch (ph) {
		case PHASE_VARS:
			// rename some vectors
			q_new=rvec;  // q_k
			q_old=vec1;  // q_k-1
			p_new=pvec;  // p_k-1
			p_old=vec2;  // p_k-2
			// initialize data structure for checkpoints
			scalars[0].ptr=&beta;
			scalars[1].ptr=&c_old;
			scalars[2].ptr=&c_new;
			scalars[3].ptr=&tau;
			scalars[4].ptr=&s_old;
			scalars[5].ptr=&s_new;
			scalars[0].size=scalars[1].size=scalars[2].size=sizeof(double);
			scalars[3].size=scalars[4].size=scalars[5].size=sizeof(doublecomplex);
			vectors[0].ptr=vec1; // now it is q_old, but can be changed further by swapping
			vectors[1].ptr=vec2; // now it is p_old, but can be changed further by swapping
			vectors[0].size=vectors[1].size=sizeof(doublecomplex);
			return;
		case PHASE_INIT:
			if (load_chpoint) { // change pointers names according to count parity
				if (IS_EVEN(niter)) SwapPointers(&q_old,&q_new);
				else SwapPointers(&p_old,&p_new);
			}
			else {
				// tau_1 = ||r_0||; q_1 = r_0(*)/||r_0||; here r_0 is already stored in q_old
				tau=sqrt(inprodR);
				nMultSelf_conj(q_new,1/creal(tau));
				// c_0=1; c_-1=0; s_0=s_-1=0
				c_new=1;
				c_old=0;
				s_new=s_old=0;
			}
			return;
		case PHASE_ITER:
			/* Avecbuffer = A.q_k. Since q_1 is r_0(*), mat-vec product for niter==1 is equivalent to Ah.r_0 (as in
			 * CGNR). Thus, matvec_ready can't be employed.
			 */
			MatVec(q_new,Avecbuffer,NULL,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			// alpha_k = q_k(T).A.q_k
			alpha=nDotProd_conj(q_new,Avecbuffer,&Timing_OneIterComm);
			// eta_k = c_k-2*c_k-1*beta_k + s_k-1(*)*alpha_k
			eta = c_old*c_new*beta + alpha*conj(s_new);
			// gamma_k = c_k-1*alpha_k - c_k-2*s_k-1*beta_k
			gamma = c_new*alpha - c_old*s_new*beta;
			// theta_k = s_k-2(*)*beta_k
			theta=beta*conj(s_old);
			// w = Aq_k - alpha_k*q_k(*) - beta_k*q_k-1(*); w is stored in q_old
			temp1=-alpha; // temp1 = -alpha_k
			// use explicitly that q_0=0
			if (niter==1) nLinComb1_cmplx_conj(q_old,q_new,Avecbuffer,temp1,&dtmp,&Timing_OneIterComm);
			else nIncrem110_d_c_conj(q_old,q_new,Avecbuffer,-beta,temp1,&dtmp,&Timing_OneIterComm);
			// beta_k+1 = ||w|| (after that beta is beta_k+1)
			// if beta=0 this is the last iteration, following formulae work fine in this case
			beta=sqrt(dtmp);
			// (k-1)-th values are moved to "old", while new (k-th) values are calculated next
			c_old=c_new;
			s_old=s_new;
			dtmp=cabs(gamma);
			if (dtmp==0) {
				// the following condition should never occur
				if (beta==0) LogError(ONE_POS,"Fatal error in CSYM iterative solver. Interaction matrix is singular");
				c_new=0;
				s_new=1;
				invksi=1/beta;
			}
			else {
				// c_k = |gamma_k| / sqrt(|gamma_k|^2 + beta_k+1^2), computed to avoid overflows
				if (dtmp<beta) {
					dtmp=dtmp/beta;
					c_new=dtmp/sqrt(1+dtmp*dtmp);
				}
				else {
					dtmp=beta/dtmp;
					c_new=1/sqrt(1+dtmp*dtmp);
				}
				// 1/ksi_k = c_k/gamma_k
				invksi=c_new/gamma;
				// s_k = beta_k+1/ksi_k = beta_k+1*c_k/gamma_k
				s_new=beta*invksi;
			}
			// p_k=(-theta_k*p_k-2-eta_k*p_k-1+q_k)/ksi_k
			if (niter==1) nMult_cmplx(p_new,q_new,invksi); // use implicitly that p_0=p_-1=0
			else {
				temp1=-eta*invksi;
				if (niter==2) nLinComb_cmplx(p_old,p_new,q_new,temp1,invksi,NULL,NULL); // use explicitly that p_0=0
				else {
					temp2=-theta*invksi;
					nIncrem111_cmplx(p_old,p_new,q_new,temp2,temp1,invksi);
				}
				SwapPointers(&p_old,&p_new);
			}
			// x_k=x_k-1+tau_k*c_k*p_k
			temp1=c_new*tau;
			nIncrem01_cmplx(xvec,p_new,temp1,NULL,NULL);
			// q_k+1 = w(*)/beta_k+1; it is first stored into q_old and then swapped
			nMultSelf_conj(q_old,1/beta);
			SwapPointers(&q_old,&q_new);
			// tau_k+1 = -s_k*tau_k; ||r_k|| = |tau_k+1|
			tau*=-s_new;
			inprodRp1=cAbs2(tau);
			dumb=tau; // dumb statement to workaround issue 146
			return; // end of PHASE_ITER
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
}

//======================================================================================================================

ITER_FUNC(QMR_CS)
/* Quasi Minimum Residual for Complex Symmetric systems, based on:
 * Freund R.W. "Conjugate gradient-type methods for linear systems with complex symmetric coefficient matrices",
 * SIAM Journal of Scientific Statistics and Computation, 13(1):425-448,1992.
 */
{
#define EPS1 1E-10 // for (vT.v)/(v.v)
#define EPS2 1E-40 // for overflow of exponent number
	static double c_old,c_new,omega_old,omega_new,zetaabs,dtmp1,dtmp2;
	static doublecomplex alpha,beta,theta,eta,zeta,zetatilda,tau,tautilda;
	static doublecomplex s_new,s_old,temp1,temp2,temp4;
	static doublecomplex *v,*vtilda,*p_new,*p_old; // can't be declared restrict due to SwapPointers

	switch (ph) {
		case PHASE_VARS:
			// rename some vectors
			v=vec1;      // v_k
			vtilda=vec2; // also v_k-1
			p_new=pvec;  // p_k
			p_old=vec3;  // p_k-1
			// initialize data structure for checkpoints
			scalars[0].ptr=&omega_old;
			scalars[1].ptr=&omega_new;
			scalars[2].ptr=&c_old;
			scalars[3].ptr=&c_new;
			scalars[4].ptr=&beta;
			scalars[5].ptr=&tautilda;
			scalars[6].ptr=&s_old;
			scalars[7].ptr=&s_new;
			scalars[0].size=scalars[1].size=scalars[2].size=scalars[3].size=sizeof(double);
			scalars[4].size=scalars[5].size=scalars[6].size=scalars[7].size=sizeof(doublecomplex);
			vectors[0].ptr=vec1; // now it is v, but can be changed further by swapping
			vectors[1].ptr=vec2; // now it is vtilda, but can be changed further by swapping
			vectors[2].ptr=vec3; // now it is p_old, but can be changed further by swapping
			vectors[0].size=vectors[1].size=vectors[2].size=sizeof(doublecomplex);
			return;
		case PHASE_INIT:
			if (load_chpoint) { // change pointers names according to count parity
				if (IS_EVEN(niter)) SwapPointers(&v,&vtilda);
				else SwapPointers(&p_old,&p_new);
			}
			else {
				// omega_0=||v_0||=0
				omega_old=0.0;
				// beta_1=sqrt(v~_1(*).v~_1); omega_1=||v~_1||/|beta_1|; (v~_1=r_0)
				beta=csqrt(nDotProdSelf_conj(rvec,&Timing_InitIterComm));
				omega_new=sqrt(inprodR)/cabs(beta); // inprodR=nNorm2(r_0)
				// v_1=v~_1/beta_1
				temp1=1/beta;
				nMult_cmplx(v,rvec,temp1);
				// tau~_1=omega_1*beta_1
				tautilda=omega_new*beta;
				// c_0=c_-1=1; s_0=s_-1=0
				c_new=c_old=1.0;
				s_new=s_old=0.0;
				dumb=beta; // dumb statement to workaround issue 146
			}
			return;
		case PHASE_ITER:
			// check for very high omega (very small beta/||v||)
			dtmp1=1/(omega_new*omega_new);
			Dz("|vT.v|/(v.v)="GFORM_DEBUG,dtmp1);
			if (dtmp1<EPS1) LogError(ONE_POS,"QMR_CS fails: |vT.v|/(v.v) is too small ("GFORM_DEBUG").",dtmp1);
			// A.v_k; alpha_k=v_k(*).(A.v_k)
			if (niter==1 && matvec_ready) { // uses that v_1=r_0/beta
				temp1=1/beta;
				nMultSelf_cmplx(Avecbuffer,temp1);
			}
			else MatVec(v,Avecbuffer,NULL,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			alpha=nDotProd_conj(v,Avecbuffer,&Timing_OneIterComm);
			// v~_k+1=-beta_k*v_k-1-alpha_k*v_k+A.v_k
			temp2=-alpha;
			if (niter==1) nLinComb1_cmplx(vtilda,v,Avecbuffer,temp2,NULL,NULL); // use explicitly that v_0=0
			else {
				temp1=-beta;
				nIncrem110_cmplx(vtilda,v,Avecbuffer,temp1,temp2);
			}
			// theta_k=s_k-2(*)*omega_k-1*beta_k
			theta=conj(s_old)*omega_old*beta;
			// eta_k=c_k-1*c_k-2*omega_k-1*beta_k+s_k-1(*)*omega_k*alpha_k
			eta = c_old*c_new*omega_old*beta;
			eta += alpha*conj(s_new)*omega_new;
			// zeta~_k=c_k-1*omega_k*alpha_k-s_k-1*c_k-2*omega_k-1*beta_k
			zetatilda = c_new*omega_new*alpha - s_new*c_old*omega_old*beta;
			// beta_k+1=sqrt(v~_k+1(*).v~_k+1); omega_k+1=||v~_k+1||/|beta_k+1|
			omega_old=omega_new;
			temp1=nDotProdSelf_conj_Norm2(vtilda,&dtmp1,&Timing_OneIterComm); // dtmp1=||v~||^2
			beta=csqrt(temp1);
			/* Here we do not check for zero beta, since exact zero is very improbable and the following code (until the
			 * end of iteration) employs only the product omega_k+1*beta_k+1. So the (almost) breakdown is instead
			 * checked at the beginning of the iteration.
			 */
			omega_new=sqrt(dtmp1)/cabs(beta);
			// |zeta_k|=sqrt(|zeta~_k|^2+omega_k+1^2*|beta_k+1|^2)
			dtmp2=cAbs2(zetatilda); // dtmp2=|zeta~_k|^2
			zetaabs=sqrt(dtmp2+dtmp1);
			dtmp1=sqrt(dtmp2); // dtmp1=|zeta~_k|
			// if (|zeta~_k|==0) zeta_k=|zeta_k|; else zeta=|zeta_k|*zeta~_k/|zeta~_k|
			if (dtmp1<EPS2) zeta=zetaabs;
			else zeta=(zetaabs/dtmp1)*zetatilda;
			// c_k=zeta~_k/zeta_k = |zeta~_k|/|zeta_k|
			c_old=c_new;
			c_new=dtmp1/zetaabs;
			// s_k+1=omega_k+1*beta_k+1/zeta_k
			s_old=s_new;
			s_new=omega_new*beta/zeta;
			// p_k=(-theta_k*p_k-2-eta_k*p_k-1+v_k)/zeta_k
			temp4=1/zeta; // temp4=1/zeta_k;
			if (niter==1) nMult_cmplx(p_new,v,temp4); // use implicitly that p_0=p_-1=0
			else {
				temp2=-eta*temp4; // temp2=-eta_k/zeta_k
				if (niter==2) nLinComb_cmplx(p_old,p_new,v,temp2,temp4,NULL,NULL);
				else {
					temp1=-theta*temp4; // temp1=-theta_k/zeta_k
					nIncrem111_cmplx(p_old,p_new,v,temp1,temp2,temp4);
				}
				SwapPointers(&p_old,&p_new);
			}
			// tau_k=c_k*tau~_k
			tau=c_new*tautilda;
			// tau~_k+1=-s_k*tau~_k
			tautilda=-s_new*tautilda;
			// x_k=x_k-1+tau_k*p_k
			nIncrem01_cmplx(xvec,p_new,tau,NULL,NULL);
			// v_k+1=v~_k+1/beta_k+1
			temp1=1/beta;
			nMultSelf_cmplx(vtilda,temp1);
			SwapPointers(&v,&vtilda); // v~ is as v_k-1 at next iteration
			// r_k = |s_k|^2*r_k-1 + (c_k*tau~_k+1/omega_k+1)*v_k+1
			temp1=(c_new/omega_new)*tautilda;
			nIncrem11_d_c(rvec,v,cAbs2(s_new),temp1,&inprodRp1,&Timing_OneIterComm);
			return; // end of PHASE_ITER
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
}
#undef EPS1L
#undef EPS1H
#undef EPS2

//======================================================================================================================

ITER_FUNC(QMR_CS_2)
/* Quasi Minimum Residual for Complex Symmetric systems based on:
 * R.W. Freund and N.M. Nachtigal, "An implementation of the qmr method based on coupled 2-term recurrences,"
 * SIAM J. Sci. Comp. 15, 313-337 (1994).
 *
 * We use recommended values of omega_k=1, which correspond to omega_k=||v_k|| used in QMR_CS above
 */
{
// breakdown tests are similar to that of BiCG_CS
#define EPS1  1E-10 // for vT.v
#define EPS2  1E-10 // for pT.A.p
	static double c_old,c_new,theta_old,theta_new,ro_old,ro_new,sabs2,dtmp1;
	static doublecomplex eps,beta,delta,eta,temp1;
	static doublecomplex * restrict v,* restrict d;

	switch (ph) {
		case PHASE_VARS:
			// rename some vectors
			v=vec1;      // v_k and v~_k+1
			d=vec2;      // d_k
			// initialize data structure for checkpoints
			scalars[0].ptr=&c_old;
			scalars[1].ptr=&theta_old;
			scalars[2].ptr=&ro_old;
			scalars[3].ptr=&eps;
			scalars[4].ptr=&eta;
			scalars[0].size=scalars[1].size=scalars[2].size=sizeof(double);
			scalars[3].size=scalars[4].size=sizeof(doublecomplex);
			vectors[0].ptr=vec1; // v
			vectors[1].ptr=vec2; // d
			vectors[0].size=vectors[1].size=sizeof(doublecomplex);
			return;
		case PHASE_INIT:
			if (!load_chpoint) {
				// ro_1=||r_0||; v~_1=r_0
				ro_old=sqrt(inprodR);
				nCopy(v,rvec);
				// c_0=eps_0=1; theta_0=0; eta_0=-1
				c_old=1;
				eps=1;
				theta_old=0;
				eta=-1;
				dumb=eps; // dumb statement to workaround issue 146
			}
			return;
		case PHASE_ITER:
			// v_k = v~_k/ro_k; this is rearranged as compared to the original algorithm
			// if ro_k=0 then c_k=1,v~_k=0, hence r_k-1=0 and iteration should have stopped by now
			nMultSelf(v,1/ro_old);
			// delta_k = v(*).v ; test it to be non zero
			delta=nDotProdSelf_conj(v,&Timing_OneIterComm);
			dtmp1=cabs(delta);
			Dz("|vT.v|="GFORM_DEBUG,dtmp1);
			if (dtmp1<EPS1) LogError(ONE_POS,"QMR_CS_2 fails: |vT.v| is too small ("GFORM_DEBUG").",dtmp1);
			// p_k = v_k - p_k-1*ro_k*delta_k/eps_k-1
			if (niter==1) nCopy(pvec,v); // use explicitly that p_0=0
			else {
				temp1=-ro_old*delta/eps;
				nIncrem10_cmplx(pvec,v,temp1,NULL,NULL);
			}
			// A.p_k
			if (niter==1 && matvec_ready) { // uses that p_1=v_1=r_0/ro_1
				nMultSelf(Avecbuffer,1/ro_old);
			}
			else MatVec(pvec,Avecbuffer,NULL,false,&Timing_OneIterMVP,&Timing_OneIterMVPComm);
			// eps_k = p_k(*).(A.p_k); beta_k = eps_k/delta_k
			eps=nDotProd_conj(pvec,Avecbuffer,&Timing_OneIterComm);
			beta=eps/delta;
			dtmp1=cabs(eps);
			Dz("|pT.A.p|="GFORM_DEBUG,dtmp1);
			if (dtmp1<EPS1) LogError(ONE_POS,"QMR_CS_2 fails: |pT.A.p| is too small ("GFORM_DEBUG").",dtmp1);
			// v~_k+1 = A.p_k - beta_k*v_k; stored in the same vector v
			temp1=-beta;
			nIncrem10_cmplx(v,Avecbuffer,temp1,&dtmp1,&Timing_OneIterComm);
			ro_new=sqrt(dtmp1); // ro_k+1 = ||v~_k+1||
			// theta_k = ro_k+1/(c_k-1*|beta_k|);
			theta_new=ro_new/(c_old*cabs(beta));
			// c_k = 1/sqrt(1+theta_k^2), |s_k|^2 = 1-c_k^2
			dtmp1=theta_new*theta_new;
			c_new=1/sqrt(1+dtmp1);
			sabs2=dtmp1/(1+dtmp1);
			// eta_k = -eta_k-1*ro_k*c_k^2/(beta_k*c_k-1^2)
			dtmp1=c_new/c_old;
			eta=-ro_old*dtmp1*dtmp1*eta/beta;
			// d_k = p_k*eta_k + d_k-1*(theta_k-1*c_k)^2
			if (niter==1) nMult_cmplx(d,pvec,eta); // use explicitly that d_0=0
			else {
				dtmp1=theta_old*c_new;
				nIncrem11_d_c(d,pvec,dtmp1*dtmp1,eta,NULL,NULL);
			}
			// x_k = x_k-1 + d_k
			nIncrem(xvec,d,NULL,NULL);
			/* The following formula to update residual was not given in the original publication, we derived it
			 * ourselves; r_k = (1-c_k^2)*r_k-1 - eta_k*v~_k+1
			 */
			temp1=-eta;
			nIncrem11_d_c(rvec,v,sabs2,temp1,&inprodRp1,&Timing_OneIterComm);
			// update variables for next iteration
			ro_old=ro_new;
			theta_old=theta_new;
			c_old=c_new;
			return; // end of PHASE_ITER
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
}
#undef EPS1
#undef EPS2

//======================================================================================================================

/* TO ADD NEW ITERATIVE SOLVER
 * Add the function implementing the iterative method to the list above in the alphabetical order. The template for the
 * function is provided below together with additional comments. Please also look at the iterative solvers, already
 * present, for examples. For operations on complex numbers you are advised to use functions from cmplx.h, for switching
 * vectors - SwapPointers (above), for linear algebra - functions from linalg.c, for multiplication of vector with
 * matrix of the linear system - MatVec function from matvec.c. Some of these functions take account of the time spent
 * on communication between different processors (in parallel mode), and increment their last argument by the
 * corresponding amount. You may also use values of variables, defined in the beginning of this source file, especially
 * niter, resid_scale, and epsB.
 */
#if 0
ITER_FUNC(_name_) // only '_name_' should be changed, the macro expansion will do the rest
// Short comment, providing full name of the iterative solver
{
/* It is recommended to define all nontrivial constants here. Do not forget to undef them at the end of this function to
 * avoid conflicts with other iterative solvers.
 */
#define EPS1 1E-30
	/* all internal variables should be defined here as static, since the function will be called many times (once per
	 * iteration).
	 */
	static double xxx;

	/* The function accepts a single argument 'ph' describing a phase, which it should perform at a particular run. This
	 * is done to move all common parts to the function IterativeSolver. Possible phases are defined and briefly
	 * explained in the definition of 'enum phase' in the beginning of this source file.
	 */
	switch (ph) {
		case PHASE_VARS:
			/* Here variables are linked to structure arrays 'scalars' and 'vectors' to initialize checkpoint system
			 * (see comment before function SaveCheckpoint). For example:
			 */
			scalars[0].ptr=&xxx;
			scalars[0].size=sizeof(double);
			/* Also, if auxiliary vectors vec1,... are used, their names may be changed to a more meaningful ones (using
			 * pointer assignments)
			 */
			return;
		case PHASE_INIT:
			/* Initialization of the iterative solver. You may use 'load_chpoint' to distinguish between the plain run
			 * and the one restarted from a checkpoint. Actual loading of checkpoint happens just before this phase. For
			 * gathering communication time use variable Timing_InitIterComm.
			 */
			return;
		case PHASE_ITER:
			/* Performs a general iteration. As a result, inprodRp1 (current residual) should be calculated. For
			 * gathering communication time use variable Timing_OneIterComm.
			 */

			// an example for checking of convergence failure (optional)
			if (xxx<EPS1) LogError(ONE_POS,"_name_ fails: xxx is too small ("GFORM_DEBUG").",xxx);

			/* _Some_ iterative solvers contain extra checks for convergence in the _middle_ of an iteration, designed
			 * to save time of, e.g., one matrix-vector product in some cases. They should be performed as follows. In
			 * particular, the intermediate test will be skipped if final checkpoint is required (checkpoint of type
			 * 'always').
			 */
			if (inprodRp1<epsB && chp_type!=CHP_ALWAYS) {
				// Additional code, e.g. to set xvec to the final value
				complete=false; // this is required to skip saving checkpoint and some timing
			}
			else {
				// Continue the iteration normally
			}
			/* Common check for convergence at the end of an iteration should not be done here, because it is performed in
			 * the function IterativeSolver.
			 */
			return;
	}
	LogError(ONE_POS,"Unknown phase (%d) of the iterative solver",(int)ph);
#undef EPS1
}
#endif

#undef ITER_FUNC

//======================================================================================================================

static void CalcFieldWKB(doublecomplex * restrict Efield)
// Calculate internal electric field in the WKB approximation; uses Einc and stores the result in Efield
{
#ifndef SPARSE	//currently no support for WKB in sparse mode
	if (prop[2]!=1) LogError(ONE_POS,"WKB initial field currently works only with default incident direction "
		"of the incoming wave (along z-axis)");
	doublecomplex vals[Nmat+1],tmpc;
	int i,k; // for traversing single-axis dimensions
	size_t dip,ind,dip_sl; // for traversing slices or up to local_nRows
	size_t boxX_l=(size_t)boxX; // to remove type conversion in indexing
#define INDEX_GRID(i) (position[(i)+2]*boxXY+position[(i)+1]*boxX_l+position[i])
	/* can be optimized by reusing material_tmp from make_particle.c or keeping the values between the calls. But
	 * this will require usage of extra memory. So the current option can be considered as corresponding to
	 * '-opt mem'
	 */
	unsigned char *mat; // same as material, but defined on whole grid (local_Ndip)
	doublecomplex *arg; // argument of exponent for corrections of incident field
#ifdef PARALLEL
	doublecomplex *bottom; // value of arg at bottom of current processor
#endif
	doublecomplex *top; // propagating value of arg at planes between the dipoles

#ifdef OPENCL // Xmatrix is not used in OpenCL, hence a complicated logic to save memory if possible
	bool a_arg=false;
	bool a_mat=false;
	bool a_top=false;
	MAXIMIZE(memPeak,memory);
	if (local_Ndip<=local_nRows) arg=Avecbuffer;
	else {
		MALLOC_VECTOR(arg,complex,local_Ndip,ALL);
		memPeak+=local_Ndip*sizeof(doublecomplex);
		a_arg=true;
	}
	if (local_Ndip*sizeof(char)<=sizeof(doublecomplex)*local_nRows) mat=(unsigned char *)Efield;
	else {
		MALLOC_VECTOR(mat,uchar,local_Ndip,ALL);
		memPeak+=local_Ndip*sizeof(char);
		a_mat=true;
	}
	if (boxXY<local_nRows) top=rvec;
	else {
		MALLOC_VECTOR(top,complex,boxXY,ALL);
		memPeak+=boxXY*sizeof(doublecomplex);
		a_top=true;
	}
#else // define all vectors using memory assigned to Xmatrix; kind of weird but should be OK
	arg=Xmatrix;
#	ifdef PARALLEL
	bottom=Xmatrix+local_Ndip;
	top=bottom+boxXY;
#	else
	top=Xmatrix+local_Ndip;
#	endif
	mat=(unsigned char *)(top + boxXY);
#endif
	// calculate function of refractive index
	for (i=0;i<Nmat;i++) vals[i]=I*(ref_index[i]-1)*kd/2;
	vals[Nmat]=0;
	// calculate values of mat (the same algorithm as in matvec), for void dipoles mat=Nmat
	for (dip=0;dip<local_Ndip;dip++) mat[dip]=(unsigned char)Nmat;
	for (dip=0,ind=0;dip<local_nvoid_Ndip;dip++,ind+=3) mat[INDEX_GRID(ind)]=material[dip];
	/* main part responsible for calculation of arg; arg[i,j,k+1]=arg[i,j,k]+vals[i,j,k]+vals[i,j,k+1]
	 * but that is done with temporary variables (not to index both k and k+1 simultaneously
	 * 'ind' traverses one slice, and 'dip' - all dipoles
	 */
	// First, calculate shifts relative to the bottom of current processor
	for(ind=0;ind<boxXY;ind++) top[ind]=0;
	for(k=local_z0,dip_sl=0;k<local_z1_coer;k++,dip_sl+=boxXY) for(ind=0,dip=dip_sl;ind<boxXY;ind++,dip++) {
		arg[dip]=top[ind]+vals[mat[dip]];
		top[ind]=arg[dip]+vals[mat[dip]];
	}
#ifdef PARALLEL
	// Second, fulfill boundary by exchanging shift values at top and bottom
	if (ExchangePhaseShifts(bottom,top,&Timing_InitIterComm))
		// Third (if required) update shift from the obtained values on the bottom
		for(k=local_z0,dip_sl=0;k<local_z1_coer;k++,dip_sl+=boxXY) for(ind=0,dip=dip_sl;ind<boxXY;ind++,dip++)
			arg[dip]+=bottom[ind];
#endif
	// E=Einc*Exp(arg), but arg is defined on a set of all (including void) dipoles
	for (ind=0;ind<local_nRows;ind+=3) {
		tmpc=cexp(arg[INDEX_GRID(ind)]);
		cvMultScal_cmplx(tmpc,Einc+ind,Efield+ind);
	}
#ifdef OPENCL // free those buffers that were allocated
	if (a_arg) Free_cVector(arg);
	if (a_mat) Free_general(mat);
	if (a_top) Free_cVector(top);
#endif
// end of !SPARSE
#else
	// should never come to this
	LogError(ONE_POS,"WKB initial field is not supported in sparse mode");
	Efield[0]=0; // redundant, to eliminate unused parameter warning
#endif
}

//======================================================================================================================

static void InitFieldfromE(void)
/* sets starting vector for linear system x_0, as well as A.x_0, r_0=b-A.x_0, and |r_0|^2 from given electric field;
 * assumes that xvec contains initial electric field, it is then replaced by x_0
 */
{
	// calculate x = (1/cc_sqrt)*V*chi*E (both x and E are stored in xvec)
	doublecomplex mult[MAX_NMAT][3];
	int i,j;
	for (i=0;i<Nmat;i++) for (j=0;j<3;j++) mult[i][j]=1/(cc_sqrt[i][j]*chi_inv[i][j]);
	nMultSelf_mat(xvec,mult);
	// calculate A.x_0, r_0=b-A.x_0, and |r_0|^2
	MatVec(xvec,Avecbuffer,NULL,false,&Timing_MVP,&Timing_MVPComm);
	nSubtr(rvec,pvec,Avecbuffer,&inprodR,&Timing_InitIterComm);
}

//======================================================================================================================

static const char *CalcInitField(double zero_resid,const enum incpol which)
/* Initializes the field as the starting point of the iterative solver. Assumes that pvec contains the right-hand side
 * of equations (b). At the end of this function xvec should contain initial vector for the iterative solver (x_0), rvec
 * - corresponding residual r_0, and inprodR - the norm of the latter residual. Returns string containing description of
 * the initial field used.
 */
{
	switch (InitField) {
		case IF_AUTO:
			/* This code is somewhat inelegant, but there seem to be no easy way to completely reuse code for other
			 * cases. Moreover, this option will probably be changed afterwards.
			 */
			// calculate A.(x_0=b), r_0=b-A.(x_0=b) and |r_0|^2
			MatVec(pvec,Avecbuffer,NULL,false,&Timing_MVP,&Timing_MVPComm);
			nSubtr(rvec,pvec,Avecbuffer,&inprodR,&Timing_InitIterComm);
			// check which x_0 is better
			if (zero_resid<inprodR) { // use x_0=0
				nInit(xvec);
				nCopy(rvec,pvec);
				inprodR=zero_resid;
				matvec_ready=true; // here Avecbuffer = A.r_0
				return "x_0 = 0\n";
			}
			else { // use x_0=Einc
				nCopy(xvec,pvec);
				return "x_0 = E_inc\n";
			}
		case IF_ZERO:
			nInit(xvec); // x_0=0
			nCopy(rvec,pvec); // r_0=b
			inprodR=zero_resid;
			return "x_0 = 0\n";
		case IF_INC:
			nCopy(xvec,pvec); // x_0=b, i.e. E_exc=E_inc
			// calculate A.(x_0=b), r_0=b-A.(x_0=b) and |r_0|^2
			MatVec(xvec,Avecbuffer,NULL,false,&Timing_MVP,&Timing_MVPComm);
			nSubtr(rvec,pvec,Avecbuffer,&inprodR,&Timing_InitIterComm);
			return "x_0 = E_inc\n";
		case IF_WKB:
			CalcFieldWKB(xvec); // calculate WKB electric field
			InitFieldfromE(); // transform it into starting vector
			return "x_0 = result of WKB\n";
		case IF_READ: {
			const char *fname;
			if (which==INCPOL_Y) fname=infi_fnameY;
			else fname=infi_fnameX; // which==INCPOL_X
			ReadField(fname,xvec); // read electric field
			InitFieldfromE(); // transform it into starting vector
			return dyn_sprintf("x_0 = from file %s\n",fname);
		}
	}
	LogError(ONE_POS,"Unknown method to calculate initial field (%d)",(int)InitField);
}

//======================================================================================================================

int IterativeSolver(const enum iter method_in,const enum incpol which)
/* choose required iterative method; do common initialization part;
 * 'which' is used only if the initial field is read from file
 */
{
	double temp;
	char tmp_str[MAX_LINE];
	TIME_TYPE tstart,time_tmp,time_tmp2,time_tmp3;

	// redundant initialization to remove warnings
	time_tmp=time_tmp2=time_tmp3=0;

	/* Instead of solving system (I+D.C).x=b , C - diagonal matrix with couple constants
	 *                                         D - symmetric interaction matrix of Green's tensor
	 * we solve system (I+S.D.S).(S.x)=(S.b), S=sqrt(C), then total interaction matrix is symmetric and
	 * Jacobi-preconditioned for any distribution of refractive index.
	 */
	/* p=b=(S.Einc) is right part of the linear system; used only here. In iteration methods themselves p is completely
	 * different vector. To avoid confusion this is done before any other initializations, specific to iterative solvers
	 */
	Timing_InitIterComm=Timing_MVP=Timing_MVPComm=0;
	tstart=GET_TIME();
	matvec_ready=false; // can be set to true only in CalcInitField (if !load_chpoint)
	if (!load_chpoint) {
		nMult_mat(pvec,Einc,cc_sqrt);
		temp=nNorm2(pvec,&Timing_InitIterComm); // |r_0|^2 when x_0=0
		resid_scale=1/temp;
		epsB=iter_eps*iter_eps*temp;
		// Calculate initial field
		const char *descr=CalcInitField(temp,which);
		// print start values
		if (IFROOT) {
			prev_err=sqrt(resid_scale*inprodR);
			sprintf(tmp_str,"%s"RESID_STRING"\n",descr,0,prev_err);
			if (!orient_avg) fprintf(logfile,"%s",tmp_str);
			printf("%s",tmp_str);
		}
		// initialize counters
		niter=1;
		counter=0;
	}
	/* determine index of the iterative solver, which is further used to get its parameters from list 'params'. This way
	 * it should be resistant to inconsistencies in orders of iterative solvers inside the list of identifiers in
	 * const.h and in the list 'params' above.
	 */
	ind_m=0;
	while (params[ind_m].meth!=method_in) {
		ind_m++;
		if (ind_m>=LENGTH(params))
			LogError(ONE_POS,"Parameters for the given iterative solver are not found in list 'params'");
	}
	// initialize data required for checkpoints and specific variables
	chp_exit=false;
	complete=true;
	if (params[ind_m].sc_N>0)
		scalars=(chp_data *)voidVector(params[ind_m].sc_N*sizeof(chp_data),ALL_POS,"list of scalars");
	else scalars=NULL;
	if (params[ind_m].vec_N>0)
		vectors=(chp_data *)voidVector(params[ind_m].vec_N*sizeof(chp_data),ALL_POS,"list of scalars");
	else vectors=NULL;
	(*params[ind_m].func)(PHASE_VARS);
	// load checkpoint, if needed, and finish initialization of the iterative solver
	if (load_chpoint) LoadIterChpoint();
	(*params[ind_m].func)(PHASE_INIT);
	// Initialization time includes generating the incident beam
	Timing_InitIter = GET_TIME() - tstart;
	Timing_InitIterComm += Timing_MVPComm; // Timing_MVPComm should (by here) include only iteration initialization
	Timing_IntFieldOneComm=Timing_InitIterComm;
	// main iteration cycle
	while (inprodR>epsB && niter<=maxiter && counter<=params[ind_m].mc && !chp_exit) {
		// initialize time
		Timing_OneIterComm=Timing_OneIterMVP=Timing_OneIterMVPComm=0;
		tstart=GET_TIME();
		// main execution
		(*params[ind_m].func)(PHASE_ITER);
		// finalize time; time for incomplete iteration may be inadequate
		Timing_OneIterComm+=Timing_OneIterMVPComm;
		Timing_IntFieldOneComm+=Timing_OneIterComm;
		Timing_MVP+=Timing_OneIterMVP;
		Timing_MVPComm+=Timing_OneIterMVPComm;
		if (complete) {
			Timing_OneIter=GET_TIME()-tstart;
			time_tmp=Timing_OneIterComm;
			time_tmp2=Timing_OneIterMVP;
			time_tmp3=Timing_OneIterMVPComm;
		}
		// use result from the previous iteration (assumed to be available by this time)
		else {
			Timing_OneIterComm=time_tmp;
			Timing_OneIterMVP=time_tmp2;
			Timing_OneIterMVPComm=time_tmp3;
		}
		/* check progress; it takes negligible time by itself (O(1) operations), but may lead to saving checkpoint.
		 * Since the latter is not relevant to the iteration itself, the ProgressReport is called after finalizing the
		 * time of a single iteration.
		 */
		ProgressReport();
	}
	// Save checkpoint of type always
	if (chp_type==CHP_ALWAYS && !chp_exit) SaveIterChpoint();
	/* process incomplete convergence
	 * Since maxiter can be used in several reasonable ways, e.g. to control execution time, we allow calculation of
	 * (potentially inaccurate) scattering quantities, when it is reached. We leave the warning although it may be
	 * redundant in e.g. scattering-order-formulation controlled by maxiter.
	 * By contrast, stagnation (especially, e.g. for CGNR) is rarely something that should be tolerated. In principle,
	 * this may happen in the end of a long run at already low residual. However, to account for such cases one should
	 * better use maxiter.
	 */
	if (inprodR>epsB) {
		if (niter>maxiter) LogWarning(EC_WARN,ONE_POS,"Iterations haven't converged in %d iterations. Further "
			"calculated scattering quantities may be less accurate.",maxiter);
		else if (counter>params[ind_m].mc) LogError(ONE_POS,"Residual norm haven't decreased for maximum allowed "
			"number of iterations (%d)",params[ind_m].mc);
	}
	if (recalc_resid) { // compute and print final residual norm
		inprodR=ResidualNorm2(xvec,rvec,Avecbuffer,&Timing_MVP,&Timing_MVPComm,&Timing_IntFieldOneComm);
		if (IFROOT) {
			temp=sqrt(resid_scale*inprodR);
			SnprintfErr(ONE_POS,tmp_str,MAX_LINE,"Final (recalculated) residual norm: "EFORM"\n",temp);
			if (!orient_avg) fprintf(logfile,"%s",tmp_str);
			printf("%s",tmp_str);
		}
	}
	// post-processing
	if (params[ind_m].sc_N>0) Free_general(scalars);
	if (params[ind_m].vec_N>0) Free_general(vectors);
	/* x is a solution of a modified system, not exactly internal field; should not be used further except for adaptive
	 * technique (as starting vector for next system)
	 */
	nMult_mat(pvec,xvec,cc_sqrt); // p now contains polarizations. Can be used to calculate e.g. scattered field faster.
	if (chp_exit) return CHP_EXIT; // check if exiting after checkpoint
	return (niter-1); // the number of iterations elapsed
}
