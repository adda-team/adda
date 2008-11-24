/* FILE : crosssec.c
 * AUTH : Maxim Yurkin
 * DESCR: All the functions to calculate scattering qunatities (except Mueller matrix).
 *        Functions to read different parameters from files.
 *        Initialization of orientation of the particle.
 *
 *        Previous versions by Martijn Frijlink
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdlib.h>
#include <time.h>
#include "vars.h"
#include "cmplx.h"
#include "const.h"
#include "Romberg.h"
#include "crosssec.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"
#include "io.h"

/* SEMI-GLOBAL VARIABLES */

/* defined and initialized in calculator.c */
extern double *E2_alldir,*E2_alldir_buffer;
extern const doublecomplex cc[][3];
/* defined and initialized in param.c */
extern const double prop_0[3],incPolX_0[3],incPolY_0[3];
extern const int ScatRelation;
/* defined and initialized in timing.c */
extern clock_t Timing_EField_ad, Timing_calc_EField_ad, Timing_comm_EField_ad,
               Timing_EField_sg, Timing_calc_EField_sg, Timing_comm_EField_sg;

/* used in CalculateE.c */
Parms_1D phi_sg;
/* used in calculator.c */
Parms_1D parms_alpha;   /* parameters of integration over alpha */
Parms_1D parms[2];      /* parameters for integration over theta,phi or beta,gamma */
angle_set beta_int,gamma_int,theta_int,phi_int; /* sets of angles */
/* used in param.c */
char avg_string[MAX_PARAGRAPH]; /* string for output of function that reads averaging parameters */

/*=====================================================================*/

INLINE int AlldirIndex(const int theta,const int phi)
    /* Convert the (theta,phi) couple into a linear array index */
{
  return (theta*phi_int.N + phi);
}

/*=====================================================================*/

void InitRotation (void)
   /* initialize matrices used for reference frame transformation */
{
  double ca,sa,cb,sb,cg,sg;
  double beta_matr[3][3];
  double alph,bet,gam;     /* in radians */

  /* initialization of angle values in radians */
  alph=Deg2Rad(alph_deg);
  bet=Deg2Rad(bet_deg);
  gam=Deg2Rad(gam_deg);
  /* calculation of rotation matrix */
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
  /* rotation of incident field */
  MatrVec(beta_matr,prop_0,prop);
  MatrVec(beta_matr,incPolY_0,incPolY);
  MatrVec(beta_matr,incPolX_0,incPolX);
}

/*=====================================================================*/
            /* currently not used */
static int ReadLine(FILE *,const char *,char *,const int) ATT_UNUSED;

static int ReadLine(FILE *file,const char *fname,  /* opened file and filename */
                    char *buf,const int buf_size)  /* buffer for line and its size */
    /* reads the first uncommented line; returns 1 if EOF reached */
{
  while (!feof(file)) {
    fgets(buf,buf_size,file);
    if (*buf!='#') {  /* if uncommented */
      if (strstr(buf,"\n")==NULL)  LogError(EC_ERROR,ONE_POS,
            "Buffer overflow while reading '%s' (size of uncommented line > %d)",
            fname,buf_size-1);
      else return 0;   /* complete line is read */
    }                         /* finish reading the commented line */
    else while (strstr(buf,"\n")==NULL  && !feof(file)) fgets(buf,buf_size,file);
  }
  return 1;
}

/*=====================================================================*/

static void ReadLineStart(FILE *file,const char *fname,  /* opened file and filename */
                          char *buf,const int buf_size,  /* buffer for line and its size */
                          const char *start)             /* beginning of the line to search */
    /* reads the first line that starts with 'start' */
{
  while (!feof(file)) {
    fgets(buf,buf_size,file);
    if (strstr(buf,start)==buf) { /* if correct beginning */
      if (strstr(buf,"\n")==NULL) LogError(EC_ERROR,ONE_POS,
            "Buffer overflow while reading '%s' (size of essential line > %d)",
            fname,buf_size-1);
      else return;  /* line found and fits into buffer */
    }                               /* finish reading unmatched line */
    else while (strstr(buf,"\n")==NULL && !feof(file)) fgets(buf,buf_size,file);
  }
  LogError(EC_ERROR,ONE_POS,
           "String '%s' is not found (in correct place) in file '%s'",start,fname);
}

/*=====================================================================*/

static int ScanIntegrParms(
     FILE *file,const char *fname,  /* opened file and filename */
     angle_set *a,                  /* pointer to angle set */
     Parms_1D *b,                   /* pointer to parameters of integration */
     const int ifcos,               /* if space angles equally in cos */
     char *buf,char* temp,          /* 2 buffers */
     const int buf_size)            /* and their size */
   /* scan integration parameters for angles from file */
{
  int i;
  double unit;

  /* scan file */
  ReadLineStart(file,fname,buf,buf_size,"min=");
  if (sscanf(buf,"min=%lf",&(a->min))!=1) return 1;
  ReadLineStart(file,fname,buf,buf_size,"max=");
  if (sscanf(buf,"max=%lf",&(a->max))!=1) return 1;
  ReadLineStart(file,fname,buf,buf_size,"Jmax=");
  if (sscanf(buf,"Jmax=%d",&(b->Jmax))!=1) return 1;
  ReadLineStart(file,fname,buf,buf_size,"K=");
  if (sscanf(buf,"K=%d",&(b->K))!=1) return 1;
  ReadLineStart(file,fname,buf,buf_size,"eps=");
  if (sscanf(buf,"eps=%lf",&(b->eps))!=1) return 1;

  ReadLineStart(file,fname,buf,buf_size,"equiv=");
  if (sscanf(buf,"equiv=%s",temp)!=1) return 1;
  if (strcmp(temp,"true")==0) b->equival=TRUE;
  else if (strcmp(temp,"false")==0) b->equival=FALSE;
  else LogError(EC_ERROR,ONE_POS,"Wrong argument of 'equiv' option in file %s",fname);

  /* fill all parameters */
  if (a->min==a->max) {
    a->N=b->Grid_size=1;
    b->Jmax=1;
  }
  else {
    /* consistency check */
    if (a->min>a->max) LogError(EC_ERROR,ONE_POS,
      "Wrong range (min=%g, max=%g) in file %s (max must be >= min)",a->min,a->max,fname);
    if (b->Jmax<ROMB_KMIN) LogError(EC_ERROR,ONE_POS,
      "Wrong Jmax (%d) in file %s (must be >=%d)",b->Jmax,fname,ROMB_KMIN);
    if (b->K<1) LogError(EC_ERROR,ONE_POS,
      "Wrong K (%d) in file %s (must be >=1)",b->K,fname);
    if (b->eps<0) LogError(EC_ERROR,ONE_POS,
      "Wrong eps (%g) in file %s (must be >=0)",b->eps,fname);

    a->N=b->Grid_size=(1 << (b->Jmax-1)) + 1;
    if (b->equival && a->N>1) (a->N)--;
  }
  /* initialize points of integration */
  if ((a->val=(double *) malloc(a->N*sizeof(double)))==NULL)
    LogError(EC_ERROR,ONE_POS,"Could not malloc integration array");
  memory += a->N*sizeof(double);

  if (ifcos) {                          /* make equal intervals in cos(angle) */
    /* consistency check */
    if (a->min<0) LogError(EC_ERROR,ONE_POS,
      "Wrong min (%g) in file %s (must be >=0 for this angle)",a->min,fname);
    if (a->max>180) LogError(EC_ERROR,ONE_POS,
      "Wrong max (%g) in file %s (must be <=180 for this angle)",a->max,fname);
    b->min=cos(Deg2Rad(a->max));
    if (fabs(b->min)<ROUND_ERR) b->min=0; /* just for convenience */
    b->max=cos(Deg2Rad(a->min));
    if (b->Grid_size==1) a->val[0]=a->min;
    else {
      unit = (b->max - b->min)/(b->Grid_size-1);
      for (i=0;i<a->N;i++) a->val[i] = Rad2Deg(acos(b->min+unit*i));
    }
  }
  else {			/* make equal intervals in angle */
    b->min=Deg2Rad(a->min);
    b->max=Deg2Rad(a->max);
    if (b->Grid_size==1) a->val[0]=a->min;
    else {
      unit = (a->max - a->min)/(b->Grid_size-1);
      for (i=0;i<a->N;i++) a->val[i] = a->min + unit*i;
    }
  }

  return 0;
}

/*=====================================================================*/

static int ScanAngleSet(FILE *file,const char *fname, /* opened file and filename */
                        angle_set *a,                 /* pointers to angle set */
                        char *buf,char *temp,         /* 2 buffers */
                        const int buf_size)           /* and their size */
   /* scan range or set of angles (theta or phi) from file (used for scat_grid) */
{
  int i;
  double unit;

  ReadLineStart(file,fname,buf,buf_size,"type=");
  if (sscanf(buf,"type=%s",temp)!=1) return -1;
  ReadLineStart(file,fname,buf,buf_size,"N=");
  if (sscanf(buf,"N=%d",&(a->N))!=1) return -1;
  /* initialize angle array */
  if ((a->val=(double *) malloc(a->N*sizeof(double)))==NULL)
    LogError(EC_ERROR,ONE_POS,"Could not malloc angle array");
  memory += a->N*sizeof(double);

  if (strcmp(temp,"range")==0) {
    ReadLineStart(file,fname,buf,buf_size,"min=");
    if (sscanf(buf,"min=%lf",&(a->min))!=1) return -1;
    ReadLineStart(file,fname,buf,buf_size,"max=");
    if (sscanf(buf,"max=%lf",&(a->max))!=1) return -1;
    if (a->N==1) a->val[0]=(a->max + a->min)/2;
    else {
      unit = (a->max - a->min)/(a->N - 1);
      for (i=0;i<a->N;i++) a->val[i] = a->min + unit*i;
    }
    return SG_RANGE;
  }
  else if (strcmp(temp,"values")==0) {
    ReadLineStart(file,fname,buf,buf_size,"values=");
    for (i=0;i<a->N;i++) {
      fgets(buf,buf_size,file);
      if (strstr(buf,"\n")==NULL) LogError(EC_ERROR,ONE_POS,
        "Buffer overflow while scanning lines in file '%s' (line size > %d)",fname,buf_size-1);
      if (sscanf(buf,"%lf\n",a->val+i)!=1) LogError(EC_ERROR,ONE_POS,
        "Failed scanning values from line '%s' in file '%s'",buf,fname);
    }
    return SG_VALUES;
  }
  else LogError(EC_ERROR,ONE_POS,"Unknown type '%s' in file '%s'",temp,fname);
  /* not actually reached */
  return -1;
}

/*=====================================================================*/

/* for convenience common error in functions
   ReadAvgParms, ReadAlldirParms, and Read ScatGridParms */
#define READ_ERROR LogError(EC_ERROR,ONE_POS,"Wrong format of file '%s'",fname);

void ReadAvgParms(const char *fname)
  /* read parameters of orientation averaging from a file */
{
  FILE *input;
  char buf[BUF_LINE],temp[BUF_LINE];

  /* open file */
  input=FOpenErr(fname,"r",ALL_POS);
  /*scan file */
  ReadLineStart(input,fname,buf,BUF_LINE,"alpha:");
  if (ScanIntegrParms(input,fname,&alpha_int,&parms_alpha,FALSE,buf,temp,BUF_LINE)) READ_ERROR;
  ReadLineStart(input,fname,buf,BUF_LINE,"beta:");
  if (ScanIntegrParms(input,fname,&beta_int,&parms[THETA],TRUE,buf,temp,BUF_LINE)) READ_ERROR;
  ReadLineStart(input,fname,buf,BUF_LINE,"gamma:");
  if (ScanIntegrParms(input,fname,&gamma_int,&parms[PHI],FALSE,buf,temp,BUF_LINE)) READ_ERROR;
  /* close file */
  FCloseErr(input,fname,ALL_POS);
  /* print info to string */
  SPRINTZ(avg_string,
    "alpha: from %g to %g in %d steps\n"\
    "beta: from %g to %g in (up to) %d steps (equally spaced in cosine values)\n"\
    "gamma: from %g to %g in (up to) %d steps\n"\
    "see file 'log_orient_avg' for details\n",
    alpha_int.min,alpha_int.max,alpha_int.N,beta_int.min,beta_int.max,beta_int.N,
    gamma_int.min,gamma_int.max,gamma_int.N);

  D("ReadAvgParms complete");
}
/*=====================================================================*/

void ReadAlldirParms(const char *fname)
   /* read integration parameters for asymmetry-paramter & C_sca
      should not be used together with orientation averaging because
      they use the same storage space - parms */
{
  FILE *input;
  char buf[BUF_LINE],temp[BUF_LINE];

  /* open file */
  input=FOpenErr(fname,"r",ALL_POS);
  /*scan file */
  ReadLineStart(input,fname,buf,BUF_LINE,"theta:");
  if (ScanIntegrParms(input,fname,&theta_int,&parms[THETA],TRUE,buf,temp,BUF_LINE)) READ_ERROR;
  ReadLineStart(input,fname,buf,BUF_LINE,"phi:");
  if (ScanIntegrParms(input,fname,&phi_int,&parms[PHI],FALSE,buf,temp,BUF_LINE)) READ_ERROR;
  /* close file */
  FCloseErr(input,fname,ALL_POS);
  /* print info */
  FPRINTZ(logfile,
    "\nScattered field is calculated for all directions (for integrated scattering quantities)\n"\
    "theta: from %g to %g in (up to) %d steps (equally spaced in cosine values)\n"\
    "phi: from %g to %g in (up to) %d steps\n"\
    "see files 'log_int_***' for details\n\n",
    theta_int.min,theta_int.max,theta_int.N,phi_int.min,phi_int.max,phi_int.N);

  D("ReadAlldirParms complete");
}

/*=====================================================================*/

void ReadScatGridParms(const char *fname)
   /* read parameters of the grid on which to calculate scattered field */
{
  FILE *input;
  char buf[BUF_LINE],temp[BUF_LINE];
  int theta_type,phi_type,i;

  /* open file */
  input=FOpenErr(fname,"r",ALL_POS);
  /* scan file */
  ReadLineStart(input,fname,buf,BUF_LINE,"global_type=");
  if (sscanf(buf,"global_type=%s",temp)!=1) READ_ERROR;
  if (strcmp(temp,"grid")==0) {
    angles.type = SG_GRID;
    ReadLineStart(input,fname,buf,BUF_LINE,"theta:");
    if ((theta_type=ScanAngleSet(input,fname,&(angles.theta),buf,temp,BUF_LINE))<0) READ_ERROR;
    if (phi_integr) {
      ReadLineStart(input,fname,buf,BUF_LINE,"phi_integr:");
      if (ScanIntegrParms(input,fname,&(angles.phi),&phi_sg,FALSE,buf,temp,BUF_LINE)) READ_ERROR;
      phi_type = SG_RANGE;
    }
    else {
      ReadLineStart(input,fname,buf,BUF_LINE,"phi:");
      if ((phi_type=ScanAngleSet(input,fname,&(angles.phi),buf,temp,BUF_LINE))<0) READ_ERROR;
    }
    angles.N=angles.theta.N*angles.phi.N;
  }
  else if (strcmp(temp,"pairs")==0) {
    if (phi_integr)
      LogError(EC_ERROR,ONE_POS,"Integration over phi can't be done with 'global_type=pairs'");
    angles.type = SG_PAIRS;
    ReadLineStart(input,fname,buf,BUF_LINE,"N=");
    if (sscanf(buf,"N=%d",&(angles.N))!=1) READ_ERROR;
    angles.theta.N=angles.phi.N=angles.N;
    /* malloc angle arrays */
    if ((angles.theta.val=(double *) malloc(angles.N*sizeof(double)))==NULL)
      LogError(EC_ERROR,ONE_POS,"Could not malloc angles.theta.val");
    if ((angles.phi.val=(double *) malloc(angles.N*sizeof(double)))==NULL)
      LogError(EC_ERROR,ONE_POS,"Could not malloc angles.phi.val");
    memory += 2*angles.N*sizeof(double);

    ReadLineStart(input,fname,buf,BUF_LINE,"pairs=");
    for (i=0;i<angles.N;i++) {
      fgets(buf,BUF_LINE,input);
      if (strstr(buf,"\n")==NULL) LogError(EC_ERROR,ONE_POS,
        "Buffer overflow while scanning lines in file '%s' (line size > %d)",fname,BUF_LINE-1);
      if (sscanf(buf,"%lf %lf\n",angles.theta.val+i,angles.phi.val+i)!=2)
        LogError(EC_ERROR,ONE_POS,"Failed scanning values from line '%s' in file '%s'",buf,fname);
    }
  }
  else LogError(EC_ERROR,ONE_POS,"Unknown global_type '%s' in file '%s'",temp,fname);
  /* close file */
  FCloseErr(input,fname,ALL_POS);
  /* print info */
  if (ringid==ROOT) {
    fprintf(logfile,"\nScattered field is calculated for multiple directions\n");
    if (angles.type==SG_GRID) {
      if (theta_type==SG_RANGE)
        fprintf(logfile,"theta: from %g to %g in %d steps\n",
                angles.theta.min,angles.theta.max,angles.theta.N);
      else if (theta_type==SG_VALUES)
        fprintf(logfile,"theta: %d given values\n",angles.theta.N);
      if (phi_type==SG_RANGE) {
        fprintf(logfile,"phi: from %g to %g in %d steps\n",
                angles.phi.min,angles.phi.max,angles.phi.N);
        if (phi_integr) fprintf(logfile,"(Mueller matrix is integrated over phi)\n");
      }
      else if (phi_type==SG_VALUES)
        fprintf(logfile,"phi: %d given values\n",angles.phi.N);
    }
    else if (angles.type==SG_PAIRS)
      fprintf(logfile,"Total %d given (theta,phi) pairs\n",angles.N);
    fprintf(logfile,"\n");
  }
  D("ReadScatGridParms complete");
#undef READ_ERROR
}

/*=====================================================================*/

void CalcField (doublecomplex *ebuff,  /* where to write calculated scattering amplitude */
	        const double *n)       /* scattering direction */
  /*  Near-optimal routine to compute the scattered fields at one specific
      angle (more exactly - scattering amplitude) */
{
  double kr,kkk;
  doublecomplex a,m2,dpr;
  doublecomplex sum[3],tbuff[3];
  int i,j,jjj;
  double temp, na;
  doublecomplex mult_mat[MAX_NMAT];
  const int scat_avg=TRUE;

  if (ScatRelation==SQ_SO) {
    /* calculate correction coefficient */
    if (scat_avg) na=0;
    else na=DotProd(n,prop);
    temp=kd*kd/24;
    for(i=0;i<Nmat;i++) {
      cSquare(ref_index[i],m2);
      /* mult_mat=1-(kd^2/24)(m^2-2(n.a)m+1) */
      mult_mat[i][RE]=1-temp*(m2[RE]-2*na*ref_index[i][RE]+1);
      mult_mat[i][IM]=temp*(2*na*ref_index[i][IM]-m2[IM]);
    }
  }
  for(i=0;i<3;i++) sum[i][RE]=sum[i][IM]=0.0;

  for (j=0;j<local_nvoid_Ndip;++j) {
    jjj=3*j;
    /* kr=k*r.n */
    kr=WaveNum*DotProd(DipoleCoord+3*j,n);
    /* a=exp(-ikr.n) */
    cExp(-kr,a);
                          /* multiply by a correction coefficient */
    if (ScatRelation==SQ_SO) cMultSelf(a,mult_mat[material[j]]);
    /* sum(P*exp(-ik*r.n)) */
    for(i=0;i<3;i++) {
      sum[i][RE]+=pvec[jjj+i][RE]*a[RE]-pvec[jjj+i][IM]*a[IM];
      sum[i][IM]+=pvec[jjj+i][RE]*a[IM]+pvec[jjj+i][IM]*a[RE];
    }
  } /* end for j */
  /* ebuff=(I-nxn).sum=sum-n*(n.sum) */
  crDotProd(sum,n,dpr);
  cScalMultRVec(n,dpr,tbuff);
  cvSubtr(sum,tbuff,ebuff);

  /* multiply it by (-i*k^3) */
  kkk=WaveNum*WaveNum*WaveNum;
  for(i=0;i<3;i++) {
    temp=ebuff[i][RE];
    ebuff[i][RE]=ebuff[i][IM]*kkk;
    ebuff[i][IM]=-temp*kkk;
  }
}

/*=====================================================================*/

double ExtCross(const double *incPol)
   /* Calculate the Extinction cross-section */
{
  doublecomplex ebuff[3],tmp;
  double sum;

  CalcField (ebuff,prop);
  crDotProd(ebuff,incPol,tmp);    /* incPol is real, so no conjugate is needed */
  sum=tmp[RE];
  MyInnerProduct(&sum,double_type,1);
  return 4*PI*sum/(WaveNum*WaveNum);
}

/*=====================================================================*/

double AbsCross(void)
  /* Calculate the Absorption cross-section for process 0 */
{
  int dip,index,i,j;
  unsigned char mat;
  double sum, dummy, temp1,temp2;
  doublecomplex m2;
  double *m; /* not doublecomplex=double[2] to allow assignment to it */
  double cc_inv_im[MAX_NMAT][3];   /* -Im(1/cc)=Im(cc)/|cc|^2 */
  double mult_mat[MAX_NMAT];

  if (ScatRelation==SQ_DRAINE) {
    /* calculate constant and cc_inv_im */
    dummy = 2*WaveNum*WaveNum*WaveNum/3;
    for (i=0;i<Nmat;i++) for (j=0;j<3;j++) cc_inv_im[i][j]=cc[i][j][IM]/cAbs2(cc[i][j]);
    /* main cycle */
    for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip) {
      mat=material[dip];
      index=3*dip;
      /* Im(P.Eexc(*))-(2/3)k^3*|P|^2=|P|^2*(-Im(1/cc)-(2/3)k^3) */
      for(i=0;i<3;i++) sum+=(cc_inv_im[mat][i] - dummy)*cAbs2(pvec[index+i]);
    }
  }
  else if (ScatRelation==SQ_SO) {
    /* calculate constants */
    temp1=kd*kd/6;
    temp2=4*PI/(gridspace*gridspace*gridspace);
    for (i=0;i<Nmat;i++) {
      m=ref_index[i];
      cSquare(m,m2);
      m2[RE]-=1;
        /* mult_mat=-Im(1/hi)*(1+(kd*Im(m))^2)/d^3;  hi=(m^2-1)/(4*PI)  */
      mult_mat[i]=temp2*m2[IM]*(1+temp1*m[IM]*m[IM])/cAbs2(m2);
    }
    /* main cycle */
    for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip)
      sum+=mult_mat[material[dip]]*cvNorm2(pvec+3*dip);
  }
  MyInnerProduct(&sum,double_type,1);
  return 4*PI*WaveNum*sum;
}

/*=====================================================================*/
                 
void CalcAlldir(void)
   /* calculate scattered field in many directions */
{
  int index,npoints,point,i,j;
  clock_t tstart,tstart2;
  double robserver[3],incPolpar[3],incPolper[3],cthet,sthet,cphi,sphi,th,ph;
  doublecomplex ebuff[3];

  /* Calculate field */
  tstart = clock();
  npoints = theta_int.N*phi_int.N;
  PRINTZ("Calculating scattered field for the whole solid angle:\n");
  for (i=0,point=0;i<theta_int.N;++i) {
    th=Deg2Rad(theta_int.val[i]);
    cthet=cos(th);
    sthet=sin(th);
    for (j=0;j<phi_int.N;++j) {
      ph=Deg2Rad(phi_int.val[j]);
      cphi=cos(ph);
      sphi=sin(ph);
      /* robserver = cos(theta)*prop + sin(theta)*[cos(phi)*incPolX + sin(phi)*incPolY]; */
      LinComb(incPolX,incPolY,cphi,sphi,robserver);
      LinComb(prop,robserver,cthet,sthet,robserver);
      /* calculate scattered field - main bottleneck */
      CalcField(ebuff,robserver);
      /* set Epar and Eper - use E2_alldir array to store them
         this is done to decrease communications in 1.5 times */

      /* incPolper = sin(phi)*incPolX - cos(phi)*incPolY; */
      LinComb(incPolX,incPolY,sphi,-cphi,incPolper);
      /* incPolpar = -sin(theta)*prop + cos(theta)*[cos(phi)*incPolX + sin(phi)*incPolY]; */
      LinComb(incPolX,incPolY,cphi,sphi,incPolpar);
      LinComb(prop,incPolpar,-sthet,cthet,incPolpar);

      index=2*point;
      crDotProd(ebuff,incPolper,((doublecomplex*)E2_alldir)[index]);
      crDotProd(ebuff,incPolpar,((doublecomplex*)E2_alldir)[index+1]);

      point++;
      if (((10*point)%npoints)<10) {
	PRINTZ(" %d%%",100*point/npoints);
        FFLUSHZ(stdout);
      }
    }
  }
  /* accumulate fields */
  tstart2 = clock();
  Accumulate(E2_alldir,4*npoints,E2_alldir_buffer);
  Timing_comm_EField_ad = clock() - tstart2;
  /* calculate square of the field */
  for (point=0;point<npoints;point++)
    E2_alldir[point] = cAbs2(((doublecomplex*)E2_alldir)[2*point]) +
                       cAbs2(((doublecomplex*)E2_alldir)[2*point+1]);
  PRINTZ("  done\n");
  FFLUSHZ(stdout);
  /* timing */
  Timing_EField_ad = clock() - tstart;
  Timing_calc_EField_ad = Timing_EField_ad - Timing_comm_EField_ad;
  Timing_EField += Timing_EField_ad;
}

/*=====================================================================*/

void CalcScatGrid(const char which)
   /* calculate scattered field in many directions */
{
  int index,point,i,j,n;
  clock_t tstart;
  double robserver[3],incPolpar[3],incPolper[3],cthet,sthet,cphi,sphi,th,ph;
  doublecomplex ebuff[3];
  doublecomplex *Egrid; /* either EgridX or EgridY */

  /* Calculate field */
  tstart = clock();
  /* choose which array to fill */
  if (which=='X') Egrid=EgridX;
  else if (which=='Y') Egrid=EgridY;
  /* set type of cycling through angles */
  if (angles.type==SG_GRID) n=angles.phi.N;
  else if (angles.type==SG_PAIRS) n=1;
  PRINTZ("Calculating grid of scattered field:\n");
  /* main cycle */
  for (i=0,point=0;i<angles.theta.N;++i) {
    th=Deg2Rad(angles.theta.val[i]);
    cthet=cos(th);
    sthet=sin(th);
    for (j=0;j<n;++j) {
      if (angles.type==SG_GRID) ph=Deg2Rad(angles.phi.val[j]);
      else if (angles.type==SG_PAIRS) ph=Deg2Rad(angles.phi.val[i]);
      cphi=cos(ph);
      sphi=sin(ph);
      /* robserver = cos(theta)*prop + sin(theta)*[cos(phi)*incPolX + sin(phi)*incPolY]; */
      LinComb(incPolX,incPolY,cphi,sphi,robserver);
      LinComb(prop,robserver,cthet,sthet,robserver);
      /* calculate scattered field - main bottleneck */
      CalcField(ebuff,robserver);
      /* set Epar and Eper - use Egrid array to store them
         this is done to decrease communications in 1.5 times */

      /* incPolper = sin(phi)*incPolX - cos(phi)*incPolY; */
      LinComb(incPolX,incPolY,sphi,-cphi,incPolper);
      /* incPolpar = -sin(theta)*prop + cos(theta)*[cos(phi)*incPolX + sin(phi)*incPolY]; */
      LinComb(incPolX,incPolY,cphi,sphi,incPolpar);
      LinComb(prop,incPolpar,-sthet,cthet,incPolpar);

      index=2*point;
      crDotProd(ebuff,incPolper,Egrid[index]);
      crDotProd(ebuff,incPolpar,Egrid[index+1]);

      point++;
      if (((10*point)%angles.N)<10) {
	PRINTZ(" %d%%",100*point/angles.N);
        FFLUSHZ(stdout);
      }
    }
  }
  /* accumulate fields; timing */
  Timing_calc_EField_sg = clock() - tstart;
  tstart = clock();
  Accumulate((double *)Egrid,4*angles.N,Egrid_buffer);
  PRINTZ("  done\n");
  FFLUSHZ(stdout);
  Timing_comm_EField_sg = clock() - tstart;
  Timing_EField_sg = Timing_calc_EField_sg + Timing_comm_EField_sg;
  Timing_EField += Timing_EField_sg;
}

/*=====================================================================*/

static void CscaIntegrand(const int theta,const int phi,double *res)
  /* function that is transferred to integration module
     when calculating Csca */
{
  res[0]=E2_alldir[AlldirIndex(theta,phi)];
}

/*=====================================================================*/

double ScaCross(void)
  /* Calculate the scattering cross section 
   * from the integral */
{
  clock_t tstart;
  char fname[MAX_FNAME];
  double res;

  sprintf(fname,"%s/" F_LOG_INT_CSCA,directory);

  tstart = clock();
  Romberg2D(parms,CscaIntegrand,1,&res,fname);
  res*=4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
  return res;
}

/*=====================================================================*/

static void gIntegrand(const int theta,const int phi,double *res)
  /* function that is transferred to integration module
     when calculating g */
{
  double E_square,th,ph;
  th=Deg2Rad(theta_int.val[theta]);
  ph=Deg2Rad(phi_int.val[phi]);

  E_square=E2_alldir[AlldirIndex(theta,phi)];
  res[0] = E_square*sin(th)*cos(ph);
  res[1] = E_square*sin(th)*sin(ph);
  res[2] = E_square*cos(th);
}
 
/*=====================================================================*/

void AsymParm(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */
{
  int comp;
  clock_t tstart;
  char log_int[MAX_FNAME];

  sprintf(log_int,"%s/" F_LOG_INT_ASYM,directory);

  tstart = clock();
  Romberg2D(parms,gIntegrand,3,vec,log_int);
  for (comp=0;comp<3;++comp) vec[comp]*=4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

static void gxIntegrand(const int theta,const int phi,double *res)
  /* function that is transferred to integration module
     when calculating g_x */
{
  res[0]=E2_alldir[AlldirIndex(theta,phi)]*sin(Deg2Rad(theta_int.val[theta]))
                                           *cos(Deg2Rad(phi_int.val[phi]));
}

/*=====================================================================*/

void AsymParm_x(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */
{
  clock_t tstart;
  char log_int[MAX_FNAME];

  sprintf(log_int,"%s/" F_LOG_INT_ASYM F_LOG_X,directory);
 
  tstart = clock();
  Romberg2D(parms,gxIntegrand,1,vec,log_int);
  vec[0] *= 4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

static void gyIntegrand(const int theta,const int phi,double *res)
  /* function that is transferred to integration module
     when calculating g_y */
{
  res[0]=E2_alldir[AlldirIndex(theta,phi)]*sin(Deg2Rad(theta_int.val[theta]))
                                           *sin(Deg2Rad(phi_int.val[phi]));
}

/*=====================================================================*/

void AsymParm_y(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */
{
  clock_t tstart;
  char log_int[MAX_FNAME];
 
  sprintf(log_int,"%s/" F_LOG_INT_ASYM F_LOG_Y,directory);
 
  tstart = clock();
  Romberg2D(parms,gyIntegrand,1,vec,log_int);
  vec[0] *= 4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

static void gzIntegrand(const int theta,const int phi,double *res)
  /* function that is transferred to integration module
     when calculating g_z */
{
  res[0]=E2_alldir[AlldirIndex(theta,phi)]*cos(Deg2Rad(theta_int.val[theta]));
}

/*=====================================================================*/

void AsymParm_z(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */

{
  clock_t tstart;
  char log_int[MAX_FNAME];
 
  sprintf(log_int,"%s/" F_LOG_INT_ASYM F_LOG_Z ,directory);
 
  tstart = clock();
  Romberg2D(parms,gzIntegrand,1,vec,log_int);
  vec[0] *= 4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

void Frp_mat(double Fsca_tot[3],double *Fsca,
	     double Finc_tot[3],double *Finc,
	     double Frp_tot[3],double *Frp)
   /* Calculate the Radiation Pressure by direct calculation
    * of the scattering force. Per dipole the force of 
    * the incoming photons, the scattering force and the 
    * radiation pressure are calculated as intermediate results */
{
  int j,l,i,comp,index;
  int local_nvoid_d0, local_nvoid_d1;
  int *nvoid_array;
  unsigned char *materialT;
  double *rdipT;
  doublecomplex *pT;
  doublecomplex temp;
  doublecomplex dummy,_E_inc;
  int lll;
  double r,r2;      /* (squared) absolute distance */
  doublecomplex
    n[3],          /* unit vector in the direction of r_{jl}
			* complex will always be zero */
    a,ab1,ab2,     /* see chapter ... */
    c1[3],c2[3],   /* idem */
    x_cg[3],       /* complex conjungate P*_j */
    Pn_j,  /* n_jl.P_l */
    Pn_l,  /* P*_j.n_jl */
    inp;   /* P*_j.P_l */


  for (comp=0;comp<3;++comp) Fsca_tot[comp]=Finc_tot[comp]=Frp_tot[comp]=0.0;
  /* Convert internal fields to dipole moments;
     Calculate incoming force per dipole */
  for (j=0;j<local_nvoid_Ndip;++j) {
    dummy[RE]=dummy[IM]=0.0;
    for (comp=0;comp<3;++comp) {
      index = 3*j+comp;
      /* Im(P.E*inc) */
      _E_inc[RE] = Einc[index][RE];
      _E_inc[IM] = -Einc[index][IM];
      cMult(pvec[index],_E_inc,temp);
      cAdd(dummy,temp,dummy);
    }
    Finc[3*j+2] = WaveNum*dummy[IM]/2;
    Finc_tot[2] += Finc[3*j+2];
  }

  /* Because of the parallelisation by row-block decomposition
     the distributed arrays involved need to be gathered on each node
     a) material -> materialT
     b) DipoleCoord -> rdipT
     c) pvec -> pT
     */
  /* initialize local_nvoid_d0 and local_nvoid_d1 */
  nvoid_array=iVector(0,nprocs-1);
  nvoid_array[ringid]=local_nvoid_Ndip;
  AllGather(nvoid_array+ringid,nvoid_array,int_type,nprocs);
  local_nvoid_d0=0;
  for (i=0;i<ringid;i++) local_nvoid_d0+=nvoid_array[i];
  local_nvoid_d1=local_nvoid_d0+local_nvoid_Ndip;
  free(nvoid_array);
  /* requires a lot of additional memory */
  if ((materialT = (unsigned char *) malloc(nvoid_Ndip*sizeof(char)))==NULL)
    LogError(EC_ERROR,ALL_POS,"Could not malloc materialT");
  if ((rdipT = (double *) malloc(3*nvoid_Ndip*sizeof(double)))==NULL)
    LogError(EC_ERROR,ALL_POS,"Could not malloc rdipT");
  if ((pT = (doublecomplex *) malloc(3*nvoid_Ndip*sizeof(doublecomplex)))==NULL)
    LogError(EC_ERROR,ALL_POS,"Could not malloc pT");

  memcpy(materialT+local_nvoid_d0,material,local_nvoid_Ndip*sizeof(char));
  memcpy(pT+3*local_nvoid_d0,pvec,3*local_nvoid_Ndip*sizeof(doublecomplex));
  memcpy(rdipT+3*local_nvoid_d0,DipoleCoord,3*local_nvoid_Ndip*sizeof(double));

  AllGather(materialT+local_nvoid_d0,materialT,char_type,local_nvoid_Ndip);
  AllGather(pT+3*local_nvoid_d0,pT,cmplx_type,3*local_nvoid_Ndip);
  AllGather(rdipT+3*local_nvoid_d0,rdipT,double_type,3*local_nvoid_Ndip);

  /* Calculate scattering force per dipole */
  for (j=local_nvoid_d0;j<local_nvoid_d1;++j) {
    int jjj = 3*j;

    for (l=0;l<nvoid_Ndip;++l) if (j!=l) {
      lll = 3*l;
      r2 = 0;
      Pn_j[RE]=Pn_j[IM]=Pn_l[RE]=Pn_l[IM]=inp[RE]=inp[IM]=0.0;

      /* Set distance related variables */
      for (comp=0;comp<3;++comp) {
        n[comp][IM] = 0;
        n[comp][RE] = rdipT[jjj+comp] - rdipT[lll+comp];
	r2 += n[comp][RE]*n[comp][RE];
      }
      r = sqrt(r2);
      n[0][RE]/=r; n[1][RE]/=r; n[2][RE]/=r;

      /* Set the scalar products a.b1 and a.b2 */
      a[RE] = cos(WaveNum*r);
      a[IM] = sin(WaveNum*r);
      ab1[RE] = 3/(r2*r2) - WaveNum*WaveNum/r2;
      ab2[RE] = -WaveNum*WaveNum/r2;
      ab1[IM] = -3*WaveNum/(r*r2);
      ab2[IM] = WaveNum*WaveNum*WaveNum/r;
      cMultSelf(ab1,a);
      cMultSelf(ab2,a);

      /* Prepare c1 and c2 */
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
	/* Set c1 */
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

	/* Set c2 */
	cMult(Pn_j,Pn_l,temp);
	cMult(n[comp],temp,c2[comp]);
        c2[comp][RE] *= -1;
	c2[comp][IM] *= -1;

	cMult(inp,n[comp],temp);
	cAdd(c2[comp],temp,c2[comp]);

	/* Fsca_{jl} = ... */
	cMultSelf(c1[comp],ab1);
	cMultSelf(c2[comp],ab2);
	Fsca[jjj-3*local_d0+comp] += (c1[comp][RE] + c2[comp][RE])/2;
      }
    } /* end l-loop */

    /* Concluding */
    for (comp=0;comp<3;++comp) {
      Fsca_tot[comp] += Fsca[jjj-3*local_d0+comp];

      Frp[jjj-3*local_d0+comp] = Finc[jjj-3*local_d0+comp] + Fsca[jjj-3*local_d0+comp];

      Frp_tot[comp] += Frp[jjj-3*local_d0+comp];
    }
  } /* end j-loop */

  /* Accumulate the total forces on all nodes */
  MyInnerProduct(Finc_tot+2,double_type,1);
  MyInnerProduct(Fsca_tot,double_type,3);
  MyInnerProduct(Frp_tot,double_type,3);

  free(materialT);
  free(rdipT);
  free(pT);
}
