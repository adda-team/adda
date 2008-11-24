/* FILE: iterative.c
 * AUTH: Maxim Yurkin
 * DESCR: Few iterative techniques to solve DDA equations
 *        Currently CGNR,BiCGStab,BiCG-CS,QMR-CS implemented
 *
 *        CGNR and BiCGStab are based on "Templates for the Solution of Linear Systems:
 *          Building Blocks for Iterative Methods"
 *          http://www.netlib.org/templates/Templates.html
 *
 *        BiCG-CS and QMR-CS are based on: Freund,R.W. "Conjugate gradient-type methods for linear
 *          systems with complex symmetric coefficient matrices",
 *          SIAM Journal of Scientific Statistics and Computation, 13(1):425-448, 1992.
 *
 *        BiCG-CS is identical to COCG, described in:
 *          van der Vorst H.A., Melissen J.B.M. "A Petrov-Galerkin type method for solving Ax=b,
 *          where A is symmetric complex", IEEE Transactions on Magnetics, 26(2):706-708, 1990.
 *
 *        CGNR was first implemented by Alfons Hoekstra
 *
 *        The linear system is composed so that diagonal terms are equal to 1, therefore
 *          use of Jacobi preconditioners does not have any effect
 *
 *        CS methods still converge to the right result even when matrix is slightly
 *           non-symmetric (e.g. -int so), however they do it much slowly than usually.
 *           It is recommended then to use BiCGStab
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <time.h>
#include "vars.h"
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "linalg.h"
#include "io.h"

/* maximum allowed iterations without residual decrease */
#define MAXCOUNT_CGNR     10
#define MAXCOUNT_BICGSTAB 10000
#define MAXCOUNT_BICG_CS  10000
#define MAXCOUNT_QMR_CS   10000
/* zero value for checks */
#define EPS_BICGSTAB      1E-30
#define EPS_BICG_CS       1E-30
#define EPS_QMR_CS        1E-30
#define EPS_QMR_CS_1      1E-40   /* problem can only occur if overflow of exponent number */

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in CalculateE.c */
extern const clock_t tstart_CE;
/* defined and initialized in calculator.c */
extern doublecomplex *xvec,*rvec,*vec1,*vec2,*vec3,*Avecbuffer;
/* defined and initialized in param.c */
extern const double eps;
/* defined and initialized in timing.c */
extern clock_t Timing_OneIter,Timing_OneIterCalc,Timing_InitIter;
extern unsigned long TotalIter;

/* LOCAL VARIABLES */

static double inprodR;     /* uses as r_0 (and main residual) in each method */
static double inprodR_0;   /* r_0; saved and load on checkpoint */
static double epsB;        /* stopping criterion */
static double resid_scale; /* scale to get square of relative error */
static double prev_err;    /* previous Rel.Error; used in ProgressReport,
                              initilized in IterativeSolver */
static int method;         /* iteration method */
static int count;          /* iteration count */
static int counter;        /* number of successive iterations without residual decrease */
static int max_count;      /* maximum allowed value of counter */
static int chp_exit;       /* checkpoint occured - exit */
static int chp_skip;       /* skip checkpoint, even if it is time to do */
typedef struct      /* data for checkpoints */
{
  void *ptr;        /* pointer to the data */
  int size;         /* size of one element */
} chp_data;
typedef struct      /* structure to hold information about different scalars and vectors */
{
  chp_data *sc;     /* array of scalar data */
  int sc_N;         /* number of scalars */
  chp_data *vec;    /* array of vector data */
  int vec_N;        /* number of vectors */
} iter_data_type;
static iter_data_type iter_data;  /* actually the structure */

/* EXTERNAL FUNCTIONS */

/* matvec.c */
void MatVec(doublecomplex *in,doublecomplex *out,double *inprod,int her);

/*============================================================*/

INLINE void SwapPointers(doublecomplex **a,doublecomplex **b)
  /* swap two pointers of (doublecomplex *) type;
     should work for others but will give "Suspisious pointer conversion" warning */
{
  doublecomplex *tmp;

  tmp=*a;
  *a=*b;
  *b=tmp;
}

/*============================================================*/

static void SaveIterChpoint(void)
  /* save a binary checkpoint;
     only limitedly foolproof - user should take care to load checkpoints
     on the same machine (number of processors) and with the same command line   */
{
  int i;
  char fname[MAX_FNAME];
  FILE *chp_file;
  clock_t tstart;

  tstart=clock();
  if (ringid==ROOT) {
    /* create directory "chp_dir" if needed and open info file */
    sprintf(fname,"%s/" F_CHP_LOG,chp_dir);
    if ((chp_file=fopen(fname,"w"))==NULL) {
      MkDirErr(chp_dir,ONE_POS);
      chp_file=FOpenErr(fname,"w",ONE_POS);
    }
    /* write info and close file */
    fprintf(chp_file,
      "Info about the run, which produced the checkpoint, can be found in ../%s",directory);
    FCloseErr(chp_file,fname,ONE_POS);
  }
  /* wait to ensure that directory exists */
  Synchronize();
  /* open output file; writing errors are checked only for vectors */
  sprintf(fname,"%s/" F_CHP,chp_dir,ringid);
  chp_file=FOpenErr(fname,"wb",ALL_POS);
  /* write commmon scalars */
  fwrite(&method,sizeof(int),1,chp_file);
  fwrite(&nlocalRows,sizeof(int),1,chp_file);
  fwrite(&count,sizeof(int),1,chp_file);
  fwrite(&counter,sizeof(int),1,chp_file);
  fwrite(&inprodR,sizeof(double),1,chp_file);
  fwrite(&inprodR_0,sizeof(double),1,chp_file);
  fwrite(&resid_scale,sizeof(double),1,chp_file);
  /* write specific scalars*/
  for (i=0;i<iter_data.sc_N;i++)
    fwrite(iter_data.sc[i].ptr,iter_data.sc[i].size,1,chp_file);
  /* write commmon vectors */
  if (fwrite(xvec,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed writing to file '%s'",fname);
  if (fwrite(rvec,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed writing to file '%s'",fname);
  if (fwrite(pvec,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed writing to file '%s'",fname);
  if (fwrite(Avecbuffer,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed writing to file '%s'",fname);
  /* write specific vectors*/
  for (i=0;i<iter_data.vec_N;i++)
    if (fwrite(iter_data.vec[i].ptr,iter_data.vec[i].size,nlocalRows,chp_file)!=nlocalRows)
      LogError(EC_ERROR,ALL_POS,"Failed writing to file '%s'",fname);
  /* close file */
  FCloseErr(chp_file,fname,ALL_POS);
  /* write info to logfile after everyone is finished */
  Synchronize();
  PRINTBOTHZ(logfile,"Checkpoint (iteration) saved\n");

  Timing_FileIO+=clock()-tstart;
}

/*============================================================*/

static void LoadIterChpoint(void)
  /* load a binary checkpoint;
     only limitedly foolproof - user should take care to load checkpoints
     on the same machine (number of processors) and with the same command line  */
{
  int i, method_new;
  size_t nlocalRows_new;
  char fname[MAX_FNAME],ch;
  FILE *chp_file;
  clock_t tstart;

  tstart=clock();
  /* open input file; reading errors are checked only for vectors */
  sprintf(fname,"%s/" F_CHP,chp_dir,ringid);
  chp_file=FOpenErr(fname,"rb",ALL_POS);
  /* check for consistency */
  fread(&method_new,sizeof(int),1,chp_file);
  if (method_new!=method)
    LogError(EC_ERROR,ALL_POS,"File '%s' is for different iterative method",fname);
  fread(&nlocalRows_new,sizeof(int),1,chp_file);
  if (nlocalRows_new!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"File '%s' is for different vector size",fname);
  /* read commmon scalars */
  fread(&count,sizeof(int),1,chp_file);
  fread(&counter,sizeof(int),1,chp_file);
  fread(&inprodR,sizeof(double),1,chp_file);
  fread(&inprodR_0,sizeof(double),1,chp_file);
  fread(&resid_scale,sizeof(double),1,chp_file);
  /* read specific scalars*/
  for (i=0;i<iter_data.sc_N;i++)
    fread(iter_data.sc[i].ptr,iter_data.sc[i].size,1,chp_file);
  /* read commmon vectors */
  if (fread(xvec,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed reading from file '%s'",fname);
  if (fread(rvec,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed reading from file '%s'",fname);
  if (fread(pvec,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed reading from file '%s'",fname);
  if (fread(Avecbuffer,sizeof(doublecomplex),nlocalRows,chp_file)!=nlocalRows)
    LogError(EC_ERROR,ALL_POS,"Failed reading from file '%s'",fname);
  /* read specific vectors*/
  for (i=0;i<iter_data.vec_N;i++)
    if (fread(iter_data.vec[i].ptr,iter_data.vec[i].size,nlocalRows,chp_file)!=nlocalRows)
      LogError(EC_ERROR,ALL_POS,"Failed reading from file '%s'",fname);
  /* check if EOF reached and close file */
  if(fread(&ch,1,1,chp_file)!=0) LogError(EC_ERROR,ALL_POS,"File '%s' is too long",fname);
  FCloseErr(chp_file,fname,ALL_POS);
  /* initialize auxiliary variables */
  epsB=eps*eps*inprodR_0;
  if (ringid==ROOT) prev_err=sqrt(resid_scale*inprodR);
  /* print info */
  PRINTBOTHZ(logfile,"Checkpoint (iteration) loaded\n");

  Timing_FileIO+=clock()-tstart;
}

/*============================================================*/

static void ProgressReport(const double inprod)
  /* Do common procedures; show progress in logfile and stdout
     also check for checkpoint condition */
{
  double err,progr,elapsed;
  char progr_string[MAX_LINE];
  char temp[5];
  time_t wt;

  if (inprod<=inprodR) {
    inprodR=inprod;
    counter=0;
  }
  else counter++;

  if (ringid==ROOT) {
    err=sqrt(resid_scale*inprod);
    progr=1-err/prev_err;
    if (counter==0) strcpy(temp,"+ ");
    else if (progr>0) strcpy(temp,"-+");
    else strcpy(temp,"- ");
    sprintf(progr_string,"RE_%03d = %.10E  %s",count,err,temp);
    if (!orient_avg) {
      fprintf(logfile,"%s  progress = %.6f\n",progr_string,progr);
      fflush(logfile);
    }
    printf("%s\n",progr_string);
    fflush(stdout);

    prev_err=err;
  }
  count++;
  TotalIter++;

  /* check condition for checkpoint; checkpoint is saved at first time  */
  if (chp_type!=CHP_NONE && chp_time!=UNDEF && !chp_skip) {
    time(&wt);
    elapsed=difftime(wt,last_chp_wt);
    if (chp_time<elapsed) {
      SaveIterChpoint();
      time(&last_chp_wt);
      if (chp_type!=CHP_REGULAR) chp_exit=TRUE;
    }
  }
}

/*============================================================*/

static void AfterIterFinished(void)
  /* Do common procedures after the iterations has finished */
{
  Timing_OneIterCalc = Timing_OneIter - Timing_OneIterComm;
  if (chp_type==CHP_ALWAYS && !chp_exit) SaveIterChpoint();
}

/*============================================================*/

static void CGNR(const int mc)
   /* Conjugate Gradient applied to Normalized Equations with minimization of Residual Norm */
{
  double inprodRplus1;		 /* inner product of rk+1 */
  double alpha, denumeratorAlpha;
  double beta,ro_new,ro_old=0;   /* initialization to remove compiler warning */
  clock_t tstart;
  chp_data scalars[1];

  max_count=mc;
  /* initialize data structure for checkpoints */
  scalars[0].ptr=&ro_old;
  scalars[0].size=sizeof(double);
  iter_data.sc=scalars;
  iter_data.sc_N=1;
  iter_data.vec=NULL;
  iter_data.vec_N=0;
  /* initialization of constants and vectors */
  if (load_chpoint) LoadIterChpoint();
  Timing_InitIter = clock() - tstart_CE;  /* initialization complete */
  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count && !chp_exit) {
    Timing_OneIterComm=0;  /* initialize time */
    tstart=clock();
            /* p_1=Ah.r_0 and ro_new=ro_0=|Ah.r_0|^2 */
    if (count==1) MatVec(rvec,pvec,&ro_new,TRUE);
    else {
      /* Avecbuffer=AH.r_k-1, ro_new=ro_k-1=|AH.r_k-1|^2 */
      MatVec(rvec,Avecbuffer,&ro_new,TRUE);
      /* beta_k-1=ro_k-1/ro_k-2 */
      beta=ro_new/ro_old;
      /* p_k=beta_k-1*p_k-1+AH.r_k-1 */
      nIncrem10(pvec,Avecbuffer,beta,NULL);
    }
    /* alpha_k=ro_k-1/|A.p_k|^2 */
    /* Avecbuffer=A.p_k */
    MatVec(pvec,Avecbuffer,&denumeratorAlpha,FALSE);
    alpha=ro_new/denumeratorAlpha;
    /* x_k=x_k-1+alpha_k*p_k */
    nIncrem01(xvec,pvec,alpha,NULL);
    /* r_k=r_k-1-alpha_k*A.p_k and |r_k|^2 */
    nIncrem01(rvec,Avecbuffer,-alpha,&inprodRplus1);
    /* initialize ro_old -> ro_k-2 for next iteration */
    ro_old=ro_new;

    Timing_OneIter=clock()-tstart;
    /* check progress */
    ProgressReport(inprodRplus1);
  } /* end of the big while loop */
  AfterIterFinished();
}

/*============================================================*/

static void BiCGStab(const int mc)
   /* Bi-Conjugate Gradient Stabilized */
{
  double inprodRplus1;		/* inner product of rk+1 */
  double denumOmega;
  doublecomplex beta,ro_new,ro_old,omega,alpha,temp1,temp2;
  doublecomplex *v,*s,*rtilda;
  clock_t tstart;
  chp_data scalars[3],vectors[3];

  max_count=mc;
  /* rename some vectors */
  v=vec1;
  s=vec2;
  rtilda=vec3;
  /* initialize data structure for checkpoints */
  scalars[0].ptr=&ro_old;
  scalars[1].ptr=&omega;
  scalars[2].ptr=&alpha;
  scalars[0].size=scalars[1].size=scalars[2].size=sizeof(doublecomplex);
  vectors[0].ptr=v;
  vectors[1].ptr=s;
  vectors[2].ptr=rtilda;
  vectors[0].size=vectors[1].size=vectors[2].size=sizeof(doublecomplex);
  iter_data.sc=scalars;
  iter_data.sc_N=3;
  iter_data.vec=vectors;
  iter_data.vec_N=3;
  /* initialization of constants and vectors */
  if (load_chpoint) LoadIterChpoint();
  else nCopy(rtilda,rvec); /* r~=r_0 */
  Timing_InitIter=clock()-tstart_CE;  /* initialization complete */
  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count && !chp_exit) {
    Timing_OneIterComm = 0;  /* initialize time */
    tstart = clock();
    /* ro_k-1=r_k-1.r~ ; check for ro_k-1!=0 */
    nDotProd(rvec,rtilda,ro_new);
    if (sqrt(cAbs2(ro_new))<EPS_BICGSTAB)
      LogError(EC_ERROR,ONE_POS,"Bi-CGStab fails: zero ro_new (%.2g).",ro_new);

    if (count==1) nCopy(pvec,rvec); /* p_1=r_0 */
    else {
      /* beta_k-1=(ro_k-1/ro_k-2)*(alpha_k-1/omega_k-1) */
      cMult(ro_new,alpha,temp1);
      cMult(ro_old,omega,temp2);
      cDiv(temp1,temp2,beta);
      /* p_k=beta_k-1*(p_k-1-omega_k-1*v_k-1)+r_k-1 */
      cMult(beta,omega,temp1);
      cInvSign(temp1);
      nIncrem110_cmplx(pvec,v,rvec,beta,temp1);
    }
    /* calculate v_k=A.p_k */
    MatVec(pvec,v,NULL,FALSE);
    /* alpha_k=ro_new/(v_k.r~) */
    nDotProd(v,rtilda,temp1);
    cDiv(ro_new,temp1,alpha);
    /* s=r_k-1-alpha*v_k-1 */
    cInvSign2(alpha,temp1);
    nLinComb1_cmplx(s,v,rvec,temp1,&inprodRplus1);
    /* check convergence at this step; if yes, checkpoint should not be saved afterwards */
    if (inprodRplus1<epsB && chp_type!=CHP_ALWAYS) {
      inprodR=inprodRplus1;
      /* x_k=x_k-1+alpha_k*p_k */
      nIncrem01_cmplx(xvec,pvec,alpha,NULL);
      chp_skip=TRUE;
    }
    else {
      /* t=Avecbuffer=A.s */
      MatVec(s,Avecbuffer,&denumOmega,FALSE);
      /* omega_k=s.t/|t|^2 ; check that omega_k!=0 */
      nDotProd(s,Avecbuffer,temp1);
      cMultReal(1/denumOmega,temp1,omega);
      if (sqrt(cAbs2(omega))<EPS_BICGSTAB)
        LogError(EC_ERROR,ONE_POS,"Bi-CGStab fails: zero omega (%.2g).",omega);
      /* x_k=x_k-1+alpha_k*p_k+omega_k*s */
      nIncrem011_cmplx(xvec,pvec,s,alpha,omega);
      /* r_k=s-omega_k*t and |r_k|^2 */
      cInvSign2(omega,temp1);
      nLinComb1_cmplx(rvec,Avecbuffer,s,temp1,&inprodRplus1);
      /* initialize ro_old -> ro_k-2 for next iteration */
      memcpy(ro_old,ro_new,sizeof(doublecomplex));
      /* take time stamp here, not to measure time of incomplete iteration
         (interrupted at the check above */
      Timing_OneIter=clock()-tstart;
    }
    /* check progress */
    ProgressReport(inprodRplus1);
  } /* end of the big while loop */
  AfterIterFinished();
}
/*============================================================*/

static void BiCG_CS(const int mc)
    /* Bi-Conjugate Gradient for Complex Symmetric systems */
{
  double inprodRplus1;		/* inner product of rk+1 */
  doublecomplex alpha, mu;
  doublecomplex beta,ro_new,ro_old,temp;
  clock_t tstart;
  chp_data scalars[1];

  max_count=mc;
  /* initialize data structure for checkpoints */
  scalars[0].ptr=&ro_old;
  scalars[0].size=sizeof(doublecomplex);
  iter_data.sc=scalars;
  iter_data.sc_N=1;
  iter_data.vec=NULL;
  iter_data.vec_N=0;
  /* initialization of constants and vectors */
  if (load_chpoint) LoadIterChpoint();
  Timing_InitIter = clock() - tstart_CE;  /* initialization complete */
  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count && !chp_exit) {
    Timing_OneIterComm=0;  /* initialize time */
    tstart=clock();
    /* ro_k-1=r_k-1(*).r_k-1; check for ro_k-1!=0 */
    nDotProdSelf_conj(rvec,ro_new);
    if (sqrt(cAbs2(ro_new))<EPS_BICG_CS) LogError(EC_ERROR,ONE_POS,
        "BiCG_CS fails: zero ro_new (%.2g%+.2gi).",ro_new[RE],ro_new[IM]);

    if (count==1) nCopy(pvec,rvec); /* p_1=r_0 */
    else {
      /* beta_k-1=ro_k-1/ro_k-2 */
      cDiv(ro_new,ro_old,beta);
      /* p_k=beta_k-1*p_k-1+r_k-1 */
      nIncrem10_cmplx(pvec,rvec,beta,NULL);
    }
    /* q_k=Avecbuffer=A.p_k */
    MatVec(pvec,Avecbuffer,NULL,FALSE);
    /* mu_k=p_k.q_k; check for mu_k!=0 */
    nDotProd_conj(pvec,Avecbuffer,mu);
    if (sqrt(cAbs2(mu))<EPS_BICG_CS) LogError(EC_ERROR,ONE_POS,
        "BiCG_CS fails: zero ro_new (%.2g%+.2gi).",ro_new[RE],ro_new[IM]);
    /* alpha_k=ro_k/mu_k */
    cDiv(ro_new,mu,alpha);
    /* x_k=x_k-1+alpha_k*p_k */
    nIncrem01_cmplx(xvec,pvec,alpha,NULL);
    /* r_k=r_k-1-alpha_k*A.p_k and |r_k|^2 */
    cInvSign2(alpha,temp);
    nIncrem01_cmplx(rvec,Avecbuffer,temp,&inprodRplus1);
    /* initialize ro_old -> ro_k-2 for next iteration */
    memcpy(ro_old,ro_new,sizeof(doublecomplex));

    Timing_OneIter=clock()-tstart;
    /* check progress */
    ProgressReport(inprodRplus1);
  } /* end of the big while loop */
  AfterIterFinished();
}

/*============================================================*/

static void QMR_CS(const int mc)
  /* Quasi Minimum Residual for Complex Symmetric systems */
{
  double inprodRplus1;		/* inner product of rk+1 */
  double c_old,c_new,omega_old,omega_new,zetaabs,dtmp1,dtmp2;
  doublecomplex alpha,beta,theta,eta,zeta,zetatilda,tau,tautilda;
  doublecomplex s_new,s_old,temp1,temp2,temp3,temp4;
  doublecomplex *v,*vtilda,*p_new,*p_old;
  clock_t tstart;
  chp_data scalars[8],vectors[3];

  max_count=mc;
  /* rename some vectors */
  v=vec1;       /* v_k */
  vtilda=vec2;  /* also v_k-1 */
  p_new=pvec;     /* p_k */
  p_old=vec3;  /* p_k-1 */
  /* initialize data structure for checkpoints */
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
  vectors[0].ptr=v;
  vectors[1].ptr=vtilda;
  vectors[2].ptr=p_old;
  vectors[0].size=vectors[1].size=vectors[2].size=sizeof(doublecomplex);
  iter_data.sc=scalars;
  iter_data.sc_N=8;
  iter_data.vec=vectors;
  iter_data.vec_N=3;
  /* initialization of constants and vectors */
  if (load_chpoint) {
    LoadIterChpoint();
    /* change pointers names according to count parity */
    if ((count%2)==0) SwapPointers(&v,&vtilda);
    else SwapPointers(&p_old,&p_new);
  }
  else {
      /* omega_0=||v_0||=0 */
    omega_old=0.0;
      /* beta_1=sqrt(v~_1(*).v~_1); omega_1=||v~_1||/|beta_1|; (v~_1=r_0) */
    nDotProdSelf_conj(rvec,temp1);
    cSqrt(temp1,beta);
    omega_new=sqrt(inprodR/cAbs2(beta));    /* inprodR=nNorm2(r_0) */
      /* v_1=v~_1/beta_1 */
    cInv(beta,temp1);
    nMult_cmplx(v,rvec,temp1);
      /* tau~_1=omega_1*beta_1 */
    cMultReal(omega_new,beta,tautilda);
      /* c_0=c_-1=1; s_0=s_-1=0 */
    c_new=c_old=1.0;
    s_new[RE]=s_new[IM]=s_old[RE]=s_old[IM]=0.0;
  }
  Timing_InitIter = clock() - tstart_CE;  /* initialization complete */
  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count && !chp_exit) {
    Timing_OneIterComm=0;  /* initialize time */
    tstart=clock();
    /* check for zero beta */
    if (sqrt(cAbs2(beta))<EPS_QMR_CS)
      LogError(EC_ERROR,ONE_POS,"QMR_CS fails: zero beta (%.2g%+.2g).",beta[RE],beta[IM]);
    /* A.v_k; alpha_k=v_k(*).(A.v_k) */
    MatVec(v,Avecbuffer,NULL,FALSE);
    nDotProd_conj(v,Avecbuffer,alpha);
    /* v~_k+1=-beta_k*v_k-1-alpha_k*v_k+A.v_k */
    cInvSign2(alpha,temp2);
       /* use explicitly that v_0=0 */
    if (count==1) nLinComb1_cmplx(vtilda,v,Avecbuffer,temp2,NULL);
    else {
      cInvSign2(beta,temp1);
      nIncrem110_cmplx(vtilda,v,Avecbuffer,temp1,temp2);
    }
    /* theta_k=s_k-2(*)*omega_k-1*beta_k */
    cMultReal(omega_old,beta,temp3);  /* temp3=omega_k-1*beta_k */
    s_old[IM]=-s_old[IM];             /* s_old is only used here, hence can be changed */
    cMult(s_old,temp3,theta);
    /* eta_k=c_k-1*c_k-2*omega_k-1*beta_k+s_k-1(*)*omega_k*alpha_k */
    cMultReal(omega_new,alpha,temp4);  /* temp4=omega_k*alpha_k */
    cMultReal(c_old*c_new,temp3,eta);
    cConj(s_new,temp1);
    cMult(temp1,temp4,temp2);
    cAdd(eta,temp2,eta);
    /* zeta~_k=c_k-1*omega_k*alpha_k-s_k-1*c_k-2*omega_k-1*beta_k */
    cMult(s_new,temp3,temp1);
    cLinComb(temp4,temp1,c_new,-c_old,zetatilda);
    /* beta_k+1=sqrt(v~_k+1(*).v~_k+1); omega_k+1=||v~_k+1||/|beta_k+1| */
    omega_old=omega_new;
    nDotProdSelf_conj_Norm2(vtilda,temp1,&dtmp1);  /* dtmp1=||v~||^2 */
    cSqrt(temp1,beta);
    omega_new=sqrt(dtmp1/cAbs2(beta));
    /* |zeta_k|=sqrt(|zeta~_k|^2+omega_k+1^2*|beta_k+1|^2) */
    dtmp2=cAbs2(zetatilda);     /* dtmp2=|zeta~_k|^2 */
    zetaabs=sqrt(dtmp2+dtmp1);
    dtmp1=sqrt(dtmp2);          /* dtmp1=|zeta~_k| */
    /* if (|zeta~_k|==0) zeta_k=|zeta_k|; else zeta=|zeta_k|*zeta~_k/|zeta~_k| */
    if (dtmp1<EPS_QMR_CS_1) {
      zeta[RE]=zetaabs;
      zeta[IM]=0.0;
    }
    else cMultReal(zetaabs/dtmp1,zetatilda,zeta);
    /* c_k=zeta~_k/zeta_k = |zeta~_k|/|zeta_k| */
    c_old=c_new;
    c_new=dtmp1/zetaabs;
    /* s_k+1=omega_k+1*beta_k+1/zeta_k */
    memcpy(s_old,s_new,sizeof(doublecomplex));
    cInv(zeta,temp4);   /* temp4=1/zeta */
    cMult(beta,temp4,temp1);
    cMultReal(omega_new,temp1,s_new);
    /* p_k=(-theta_k*p_k-2-eta_k*p_k-1+v_k)/zeta_k */
       /* use explicitly that p_0=p_-1=0 */
    if (count==1) nMult_cmplx(p_new,v,temp4);
    else {
      cMult(eta,temp4,temp2);
      cInvSign(temp2);       /* temp2=-eta_k/zeta_k */
      if (count==2) nLinComb_cmplx(p_old,p_new,v,temp2,temp4,NULL);
      else {
        cMult(theta,temp4,temp1);
        cInvSign(temp1);     /* temp1=-theta_k/zeta_k */
        nIncrem111_cmplx(p_old,p_new,v,temp1,temp2,temp4);
      }
      SwapPointers(&p_old,&p_new);
    }
    /* tau_k=c_k*tau~_k */
    cMultReal(c_new,tautilda,tau);
    /* tau~_k+1=-s_k*tau~_k */
    cMult(s_new,tautilda,temp1);
    cInvSign2(temp1,tautilda);
    /* x_k=x_k-1+tau_k*p_k */
    nIncrem01_cmplx(xvec,p_new,tau,NULL);
    /* v_k+1=v~_k+1/beta_k+1 */
    cInv(beta,temp1);
    nMultSelf_cmplx(vtilda,temp1);
    SwapPointers(&v,&vtilda);     /* vtilda is as v_k-1 at next iteration */
    /* r_k=r_k-1+(c_k*tau~_k+1/omega_k+1)*v_k+1 */
    cMultReal(c_new/omega_new,tautilda,temp1);
    nIncrem11_d_c(rvec,v,cAbs2(s_new),temp1,&inprodRplus1);

    Timing_OneIter=clock()-tstart;
    /* check progress */
    ProgressReport(inprodRplus1);
  } /* end of the big while loop */
  AfterIterFinished();
}

/*============================================================*/

int IterativeSolver(const int method_in)
    /* choose required iterative method
       do common initialization part */
{
  double temp;
  char tmp_str[MAX_LINE];

  method=method_in;
  chp_exit=FALSE;
  chp_skip=FALSE;
  /* instead of solving system (I+D.C).x=b , C - diagonal matrix with couple constants
   *                                         D - symmetric interaction matrix of Green's tensor
   * we solve system (I+S.D.S).(S.x)=(S.b), S=sqrt(C), them
   * total interaction matrix is symmetric and Jacobi-preconditioned for any discribution of m */

  /* p=b=(S.Einc) is right part of the linear system; used only here,
     in iteration methods themselves p is completely different vector */
  if (!load_chpoint) {
    nMult_mat(pvec,Einc,cc_sqrt);

    temp=nNorm2(pvec); /* |r_0|^2 when x_0=0 */
    resid_scale=1/temp;
    /* calculate A.(x_0=b), r_0=b-A.(x_0=b) and |r_0|^2 */
    MatVec(pvec,Avecbuffer,NULL,FALSE);
    nSubtr(rvec,pvec,Avecbuffer,&inprodR);
    /* check which x_0 is better */
    if (temp<inprodR) {  /* use x_0=0 */
      nInit(xvec);
      /* r=p, but faster than copy, p is not used afterwards */
      SwapPointers(&rvec,&pvec);
      inprodR=temp;
      strcpy(tmp_str,"x_0 = 0\n");
    }
    else {              /* use x_0=Einc */
      /* x=p, but faster than copy, p is not used afterwards */
      SwapPointers(&xvec,&pvec);
      strcpy(tmp_str,"x_0 = E_inc\n");
    }
    inprodR_0=inprodR;
    epsB=eps*eps*inprodR_0;
    /* print start values */
    if (ringid==ROOT) {
      prev_err=sqrt(resid_scale*inprodR);
      sprintf(tmp_str+strlen(tmp_str),"RE_000 = %.10E\n",prev_err);
      if (!orient_avg) {
        fprintf(logfile,"%s",tmp_str);
        fflush(logfile);
      }
      printf("%s",tmp_str);
      fflush(stdout);
    }
    /* initialize counters */
    count=1;
    counter=0;
  }
  /* call appropriate iterative method */
  if (method==IT_CGNR) CGNR(MAXCOUNT_CGNR);
  else if (method==IT_BICGSTAB) BiCGStab(MAXCOUNT_BICGSTAB);
  else if (method==IT_BICG_CS) BiCG_CS(MAXCOUNT_BICG_CS);
  else if (method==IT_QMR_CS) QMR_CS(MAXCOUNT_QMR_CS);
  /* error output */
  if (count>maxiter) LogError(EC_ERROR,ONE_POS,
    "Iterations haven't converged in maximum allowed number of iterations (%d)",maxiter);
  else if (counter>max_count) LogError(EC_ERROR,ONE_POS,
    "Residual norm haven't decreased for maximum allowed number of iterations (%d)",max_count);
  /* postprocessing */
  /* x is a solution of a modified system, not exactly internal field
     should not be used further except fot adaptive technique
        (as starting vector for next system) */
  nMult_mat(pvec,xvec,cc_sqrt);   /* p is now vector of polarizations -  */
                                  /* faster to calculate ,e.g. scattered field */
  /* check if exiting after checkpoint */
  if (chp_exit) return CHP_EXIT;
  return count;
}
