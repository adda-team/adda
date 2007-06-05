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
 */
#include <stdio.h>
#include <time.h>
#include "cmplx.h"
#include "const.h"
#include "comm.h"
#include "debug.h"
#include "linalg.h"

/* macro to swap two doublecomplex pointers */
#define SWAP_PCMPLX(a,b) {pctmp=a;a=b;b=pctmp;}
doublecomplex *pctmp;

/* maximum allowed iterations without residual decrease */
#define MAXCOUNT_CGNR     10
#define MAXCOUNT_BICGSTAB 5000
#define MAXCOUNT_BICG_CS  5000
#define MAXCOUNT_QMR_CS   5000
/* zero value for checks */
#define EPS_BICGSTAB      1E-20
#define EPS_BICG_CS       1E-20
#define EPS_QMR_CS        1E-20
#define EPS_QMR_CS_1      1E-30   /* problem can only occur if overflow of exponent number */

extern doublecomplex *x,*r,*p,*Avecbuffer,*vec1,*vec2,*vec3;
extern doublecomplex *Einc;

extern FILE *logfile;           /* defined in main */
extern double eps;

extern int nDip;		/* defined in calculator */
extern double WaveNum;		/* defined in calcularor */

extern int maxiter;

/* function that performs matvec products */
extern void MatVec(doublecomplex *in,doublecomplex *out,double *inprod,int her);

extern int orient_avg;

extern clock_t tstart_CE;
extern clock_t Timing_OneIter,Timing_OneIterCalc,Timing_OneIterComm,
      Timing_InitIter;
extern unsigned long TotalIter;

double inprodR;     /* uses as r_0 (and main residual) in each method */
double epsB;        /* stopping criterion */
double resid_scale; /* scale to get square of relative error */
double prev_err;    /* previous Rel.Error; used in ProgressReport, initilized in iterative_solver */

/*============================================================*/

void ProgressReport(double inprod, int count, int *Pcounter)
  /* Do common procedures; show progress in logfile and stdout */
{
  double err,progr;
  static char progr_string[30];
  char temp[5];

  if (inprod<=inprodR) {
    inprodR=inprod;
    (*Pcounter)=0;
  }
  else (*Pcounter)++;

  if (ringid==ROOT) {
    err=sqrt(resid_scale*inprod);
    progr=1-err/prev_err;
    if ((*Pcounter)==0) strcpy(temp,"+ ");
    else if (progr>0) strcpy(temp,"-+");
    else strcpy(temp,"- ");
    sprintf(progr_string,"RE_%03d = %1.10e  %s",count,err,temp);
    if (!orient_avg) {
      fprintf(logfile,"%s  progress = %.6f\n",progr_string,progr);
      fflush(logfile);
    }
    printf("%s\n",progr_string);
    fflush(stdout);

    prev_err=err;
  }
  TotalIter++;
}

/*============================================================*/

int CGNR(int max_count)
   /* Conjugate Gradient applied to Normalized Equations with minimization of Residual Norm */
{
  int count,counter;	        /* two counters to control the convergence */
  double inprodRplus1;		/* inner product of rk+1 */
  double alpha, denumeratorAlpha;
  double beta,ro_new,ro_old=0;   /* initialization to remove compiler warning */
  clock_t tstart;

  count=1;          /* initialize counters */
  counter=0;
  Timing_InitIter = clock() - tstart_CE;  /* initialization complete */
  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count) {
    Timing_OneIterComm=0;  /* initialize time */
    tstart=clock();
            /* p_1=Ah.r_0 and ro_new=ro_0=|Ah.r_0|^2 */
    if (count==1) MatVec(r,p,&ro_new,true);
    else {
      /* Avecbuffer=AH.r_k-1, ro_new=ro_k-1=|AH.r_k-1|^2 */
      MatVec(r,Avecbuffer,&ro_new,true);
      /* beta_k-1=ro_k-1/ro_k-2 */      
      beta=ro_new/ro_old;
      /* p_k=beta_k-1*p_k-1+AH.r_k-1 */
      nIncrem10(p,Avecbuffer,beta,NULL);
    }
    /* alpha_k=ro_k-1/|A.p_k|^2 */
    /* Avecbuffer=A.p_k */
    MatVec(p,Avecbuffer,&denumeratorAlpha,false);
    alpha=ro_new/denumeratorAlpha;
    /* x_k=x_k-1+alpha_k*p_k */
    nIncrem01(x,p,alpha,NULL);
    /* r_k=r_k-1-alpha_k*A.p_k and |r_k|^2 */
    nIncrem01(r,Avecbuffer,-alpha,&inprodRplus1);
    /* initialize ro_old -> ro_k-2 for next iteration */
    ro_old=ro_new;
    /* check progress */
    ProgressReport(inprodRplus1,count,&counter);
    count++;
    Timing_OneIter=clock()-tstart;
  } /* end of the big while loop */
  
  Timing_OneIterCalc = Timing_OneIter - Timing_OneIterComm;

  /* so, let us check for convergence */
  if (count>maxiter) return -1;  /* no convergence, too many iterations */
  else if (counter>max_count) return -2; /* no convergence, residuals increase too much */
  else return count; /* convergence */
}

/*============================================================*/

int BiCGStab(int max_count)
   /* Bi-Conjugate Gradient Stabilized */
{
  int count,counter;	        /* two counters to control the convergence */
  double inprodRplus1;		/* inner product of rk+1 */
  double denumOmega;
  doublecomplex beta,ro_new,ro_old,omega,alpha,temp1,temp2;
  doublecomplex *v,*s,*rtilda;
  clock_t tstart;

  /* rename some vectors */
  v=vec1;
  s=vec2;
  rtilda=vec3;

  count=1;          /* initialize counters */
  counter=0;
  nCopy(rtilda,r); /* r~=r_0 */

  Timing_InitIter=clock()-tstart_CE;  /* initialization complete */

  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count) {
    Timing_OneIterComm = 0;  /* initialize time */
    tstart = clock();
    /* ro_k-1=r_k-1.r~ ; check for ro_k-1!=0 */
    nDotProd(r,rtilda,ro_new);
    if (sqrt(cAbs2(ro_new))<EPS_BICGSTAB)
      LogError(EC_ERROR,ONE,POSIT,"Bi-CGStab fails: zero ro_new (%.2g).",ro_new);

    if (count==1) nCopy(p,r); /* p_1=r_0 */
    else {
      /* beta_k-1=(ro_k-1/ro_k-2)*(alpha_k-1/omega_k-1) */
      cMult(ro_new,alpha,temp1);
      cMult(ro_old,omega,temp2);
      cDiv(temp1,temp2,beta);
      /* p_k=beta_k-1*(p_k-1-omega_k-1*v_k-1)+r_k-1 */
      cMult(beta,omega,temp1);
      cInvSign(temp1);
      nIncrem110_cmplx(p,v,r,beta,temp1);
    }
    /* calculate v_k=A.p_k */
    MatVec(p,v,NULL,false);
    /* alpha_k=ro_new/(v_k.r~) */
    nDotProd(v,rtilda,temp1);
    cDiv(ro_new,temp1,alpha);
    /* s=r_k-1-alpha*v_k-1 */
    cInvSign2(alpha,temp1);
    nLinComb1_cmplx(s,v,r,temp1,&inprodRplus1);

    if (inprodRplus1<epsB) {   /* check convergence at this step */
      inprodR=inprodRplus1;
      /* x_k=x_k-1+alpha_k*p_k */
      nIncrem01_cmplx(x,p,alpha,NULL);
    }
    else {
      /* t=Avecbuffer=A.s */
      MatVec(s,Avecbuffer,&denumOmega,false);
      /* omega_k=s.t/|t|^2 ; check that omega_k!=0 */
      nDotProd(s,Avecbuffer,temp1);
      cMultReal(1/denumOmega,temp1,omega);
      if (sqrt(cAbs2(omega))<EPS_BICGSTAB)
        LogError(EC_ERROR,ONE,POSIT,"Bi-CGStab fails: zero omega (%.2g).",omega);
      /* x_k=x_k-1+alpha_k*p_k+omega_k*s */
      nIncrem011_cmplx(x,p,s,alpha,omega);
      /* r_k=s-omega_k*t and |r_k|^2 */
      cInvSign2(omega,temp1);
      nLinComb1_cmplx(r,Avecbuffer,s,temp1,&inprodRplus1);
      /* initialize ro_old -> ro_k-2 for next iteration */
      memcpy(ro_old,ro_new,sizeof(doublecomplex));
      /* take time stamp here, not to measure time of incomplete iteration
         (interrupted at the check above */
      Timing_OneIter=clock()-tstart;
    }
    /* check progress */
    ProgressReport(inprodRplus1,count,&counter);
    count++;
  } /* end of the big while loop */
  Timing_OneIterCalc = Timing_OneIter - Timing_OneIterComm;
  /* so, let us check for convergence */
  if (count>maxiter) return -1;  /* no convergence, too many iterations */
  else if (counter>max_count) return -2; /* no convergence, residuals increase too much */
  else return count; /* convergence */
}

/*============================================================*/

int BiCG_CS(int max_count)
    /* Bi-Conjugate Gradient for Complex Symmetric systems */
{
  int count,counter;	        /* two counters to control the convergence */
  double inprodRplus1;		/* inner product of rk+1 */
  doublecomplex alpha, mu;
  doublecomplex beta,ro_new,ro_old,temp;
  clock_t tstart;

  count=1;          /* initialize counters */
  counter=0;
  Timing_InitIter = clock() - tstart_CE;  /* initialization complete */
  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count) {
    Timing_OneIterComm=0;  /* initialize time */
    tstart=clock();
    /* ro_k-1=r_k-1(*).r_k-1; check for ro_k-1!=0 */
    nDotProdSelf_conj(r,ro_new);
    if (sqrt(cAbs2(ro_new))<EPS_BICG_CS)
      LogError(EC_ERROR,ONE,POSIT,"COCG fails: zero ro_new (%.2g%+.2g).",ro_new[re],ro_new[im]);

    if (count==1) nCopy(p,r); /* p_1=r_0 */
    else {
      /* beta_k-1=ro_k-1/ro_k-2 */
      cDiv(ro_new,ro_old,beta);
      /* p_k=beta_k-1*p_k-1+r_k-1 */
      nIncrem10_cmplx(p,r,beta,NULL);
    }
    /* q_k=Avecbuffer=A.p_k */
    MatVec(p,Avecbuffer,NULL,false);
    /* mu_k=p_k.q_k; check for mu_k!=0 */
    nDotProd_conj(p,Avecbuffer,mu);
    if (sqrt(cAbs2(mu))<EPS_BICG_CS)
      LogError(EC_ERROR,ONE,POSIT,"COCG fails: zero ro_new (%.2g%+.2g).",ro_new[re],ro_new[im]);
    /* alpha_k=ro_k/mu_k */
    cDiv(ro_new,mu,alpha);
    /* x_k=x_k-1+alpha_k*p_k */
    nIncrem01_cmplx(x,p,alpha,NULL);
    /* r_k=r_k-1-alpha_k*A.p_k and |r_k|^2 */
    cInvSign2(alpha,temp);
    nIncrem01_cmplx(r,Avecbuffer,temp,&inprodRplus1);
    /* initialize ro_old -> ro_k-2 for next iteration */
    memcpy(ro_old,ro_new,sizeof(doublecomplex));
    /* check progress */
    ProgressReport(inprodRplus1,count,&counter);
    count++;
    Timing_OneIter=clock()-tstart;
  } /* end of the big while loop */

  Timing_OneIterCalc = Timing_OneIter - Timing_OneIterComm;

  /* so, let us check for convergence */
  if (count>maxiter) return -1;  /* no convergence, too many iterations */
  else if (counter>max_count) return -2; /* no convergence, residuals increase too much */
  else return count; /* convergence */
}

/*============================================================*/

int QMR_CS(int max_count)
  /* Quasi Minimum Residual for Complex Symmetric systems */
{
  int count,counter;	        /* two counters to control the convergence */
  double inprodRplus1;		/* inner product of rk+1 */
  double c_old,c_new,omega_old,omega_new,zetaabs,dtmp1,dtmp2;
  doublecomplex alpha,beta,theta,eta,zeta,zetatilda,tau,tautilda;
  doublecomplex s_new,s_old,temp1,temp2,temp3,temp4;
  doublecomplex *v,*vtilda,*p_new,*p_old;
  clock_t tstart;

  count=1;          /* initialize counters */
  counter=0;
  /* rename some vectors */
  v=vec1;       /* v_k */
  vtilda=vec2;  /* also v_k-1 */
  p_new=p;     /* p_k */
  p_old=vec3;  /* p_k-1 */
  /* initialization of constants and vectors */
    /* omega_0=||v_0||=0 */
  omega_old=0.0;
    /* beta_1=sqrt(v~_1(*).v~_1); omega_1=||v~_1||/|beta_1|; (v~_1=r_0) */
  nDotProdSelf_conj(r,temp1);
  cSqrt(temp1,beta);
  omega_new=sqrt(inprodR/cAbs2(beta));    /* inprodR=nNorm2(r_0) */
    /* v_1=v~_1/beta_1 */
  cInv(beta,temp1);
  nMult_cmplx(v,r,temp1);
    /* tau~_1=omega_1*beta_1 */
  cMultReal(omega_new,beta,tautilda);
    /* c_0=c_-1=1; s_0=s_-1=0 */
  c_new=c_old=1.0;
  s_new[re]=s_new[im]=s_old[re]=s_old[im]=0.0;

  Timing_InitIter = clock() - tstart_CE;  /* initialization complete */
  /* main iteration cycle */
  while (inprodR>=epsB && count<=maxiter && counter<=max_count) {
    Timing_OneIterComm=0;  /* initialize time */
    tstart=clock();
    /* check for zero beta */
    if (sqrt(cAbs2(beta))<EPS_QMR_CS)
      LogError(EC_ERROR,ONE,POSIT,"QMR_CS fails: zero beta (%.2g%+.2g).",beta[re],beta[im]);
    /* A.v_k; alpha_k=v_k(*).(A.v_k) */
    MatVec(v,Avecbuffer,NULL,false);
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
    s_old[im]=-s_old[im];             /* s_old only used here, hence can be changed */
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
      zeta[re]=zetaabs;
      zeta[im]=0.0;
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
      SWAP_PCMPLX(p_old,p_new);
    }
    /* tau_k=c_k*tau~_k */
    cMultReal(c_new,tautilda,tau);
    /* tau~_k+1=-s_k*tau~_k */
    cMult(s_new,tautilda,temp1);
    cInvSign2(temp1,tautilda);
    /* x_k=x_k-1+tau_k*p_k */
    nIncrem01_cmplx(x,p_new,tau,NULL);
    /* v_k+1=v~_k+1/beta_k+1 */
    cInv(beta,temp1);
    nMultSelf_cmplx(vtilda,temp1);
    SWAP_PCMPLX(v,vtilda);     /* vtilda is as v_k-1 at next iteration */
    /* r_k=r_k-1+(c_k*tau~_k+1/omega_k+1)*v_k+1 */
    cMultReal(c_new/omega_new,tautilda,temp1);
    nIncrem11_d_c(r,v,cAbs2(s_new),temp1,&inprodRplus1);
    /* check progress */
    ProgressReport(inprodRplus1,count,&counter);
    count++;
    Timing_OneIter=clock()-tstart;
  } /* end of the big while loop */

  Timing_OneIterCalc = Timing_OneIter - Timing_OneIterComm;

  /* so, let us check for convergence */
  if (count>maxiter) return -1;  /* no convergence, too many iterations */
  else if (counter>max_count) return -2; /* no convergence, residuals increase too much */
  else return count; /* convergence */
}

/*============================================================*/

int iterative_solver(int method)
    /* choose required iterative method
       do common initialization part */
{
  int i,res,mc;
  double temp;
  char tmp_str[50];

  extern doublecomplex cc_sqrt[MAXNMAT][3];

  /* instead of solving system (I+D.C).x=b , C - diagonal matrix with couple constants
   *                                         D - symmetric interaction matrix of Green's tensor
   * we solve system (I+S.D.S).(S.x)=(S.b), S=sqrt(C), them
   * total interaction matrix is symmetric and Jacobi-preconditioned for any discribution of m
   */
  /* p=b=(S.Einc) is right part of the linear system; used only here, in iteration methods themselves
         p is completely different vector */
  nMult_mat(p,Einc,cc_sqrt);

  temp=nNorm2(p); /* |r_0|^2 when x_0=0 */
  resid_scale=1/temp;
  /* calculate A.(x_0=b), r_0=b-A.(x_0=b) and |r_0|^2 */
  MatVec(p,Avecbuffer,NULL,false);
  nSubtr(r,p,Avecbuffer,&inprodR);
  /* check which x_0 is better */
  if (temp<inprodR) {  /* use x_0=0 */
    nInit(x);
    /* r=p, but faster than copy, p is not used afterwards */
    SWAP_PCMPLX(r,p);
    inprodR=temp;
    strcpy(tmp_str,"x_0 = 0\n");
  }
  else {              /* use x_0=Einc */
    /* x=p, but faster than copy, p is not used afterwards */
    SWAP_PCMPLX(x,p);
    strcpy(tmp_str,"x_0 = E_inc\n");
  }
  epsB=eps*eps*inprodR;
  /* print start values */
  if (ringid==ROOT) {
    prev_err=sqrt(resid_scale*inprodR);
    sprintf(tmp_str+strlen(tmp_str),"RE_000 = %1.10e\n",prev_err);
    if (!orient_avg) {
      fprintf(logfile,"%s",tmp_str);
      fflush(logfile);
    }
    printf("%s",tmp_str);
    fflush(stdout);
  }
  /* call appropriate iterative method */
  if (method==IT_CGNR) res=CGNR(mc=MAXCOUNT_CGNR);
  else if (method==IT_BICGSTAB) res=BiCGStab(mc=MAXCOUNT_BICGSTAB);
  else if (method==IT_BICG_CS) res=BiCG_CS(mc=MAXCOUNT_BICG_CS);
  else if (method==IT_QMR_CS) res=QMR_CS(mc=MAXCOUNT_QMR_CS);
  /* error output */
  if (res==-2) LogError(EC_ERROR,ONE,POSIT,
    "Iterations haven't converged in maximum allowed number of iterations (%d)",maxiter);
  else if (res==-1) LogError(EC_ERROR,ONE,POSIT,
    "Residual norm haven't decreased for maximum allowed number of iterations (%d)",mc);
  /* postprocessing */
  /* x is a solution of a modified system, not exactly internal field - should be used further
     except fot adaptive technique as starting vector for next system */
  nMult_mat(p,x,cc_sqrt);   /* p is now vector of polarizations -  */
                            /* faster to calculate ,e.g. scattered field */
  return res;
}
