/* FILE : crosssec.c
 * AUTH : Martijn Frijlink
 * DESCR: 1) The functions Ext_cross, Abs_cross and Frp_mat are self-explanatory.
 *        2) The function set_parms reads the integration-parameters from file.
 *        3) The function calc_alldir calculates the scattered field in a
 *           predefined set of directions (predefined by the integration-
 *           parameters). The function fill_tab calculates the trigonometric
 *           function-values belonging to these directions.
 *        4) The functions Sca_cross calls the Romberg-routine for integrating
 *           the scattered intensity resulting in the scattering cross section.
 *        5) The functions Asym_parm_x, Asym_parm_y and Asym_parm_z each
 *           calculate the corresponding component of the asymmetry-parameter,
 *           not yet divided by the scattering cross section. It can therefore
 *           also be interpreted as the scattering force in cross section units.
 *
 *        Currently is developed by Maxim Yurkin
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "cmplx.h"
#include "const.h"
#include "Romberg.h"
#include "crosssec.h"
#include "comm.h"
#include "debug.h"
#include "memory.h"

#define X 0
#define Y 1
#define Z 2
#define R3 3

/* Convert the (Theta,Phi) couple into a linear array index */
#define grid_index(T,P) (R3*((T)*parms[PHI].Grid_size + (P)))

double *si[Narg],*co[Narg];
Parms_1D parms[Narg];
double *n;

extern doublecomplex *Egrid;
extern double *Egrid_buffer; /* defined in main */
extern doublecomplex *x;     /* solution of the internal field */
extern doublecomplex *p;     /* polarization */
extern doublecomplex *Einc;          /* defined in main */
extern char *material;           /* defined calculator */
extern int Nmat;                /* defined in main */
extern double *DipoleCoord;
extern doublecomplex cc[MAXNMAT][3]; /* defined in main */
extern FILE *logfile;           /* defined in main */
extern char directory[200];     /* defined in main */
extern int store_all_dir;       /* defined in main */
extern clock_t            /* defined in main */
  Timing_EField,
  Timing_EField_ad,
  Timing_calc_EField_ad,
  Timing_comm_EField_ad,
  Timing_Integration,
  Timing_FileIO;
extern double WaveNum;
extern char Romb_parms[200];
extern doublecomplex ref_index[MAXNMAT];
extern double kd;
extern double gridspace;
extern double prop[3], incPolX[3], incPolY[3];
extern int ScatRelation;

double beta_matr[3][3];

extern double alph_deg,bet_deg,gam_deg;
double alph,bet,gam;     /* in radians */

integr_parms alpha_int, beta_int, gamma_int;
extern Parms_1D parms_alpha;

extern int equivalent;

/*=====================================================================*/

void init_rotation (void)
   /* initializa matrices used for reference frame transformation */
{
  double ca,sa,cb,sb,cg,sg;
  extern double prop_0[3], incPolY_0[3], incPolX_0[3];
  /* initialization of angle values in radians */
  alph=deg2rad(alph_deg);
  bet=deg2rad(bet_deg);
  gam=deg2rad(gam_deg);
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

void calc_field (doublecomplex *ebuff,  /* where to write calculated scattering amplitude */
	         double *n)             /* scattering direction */
{
  /*  Near-optimal routine to compute the scattered fields at one specific
   *  angle. (more exactly - scattering amplitude)
   */

  double kr,kkk;
  doublecomplex a,m2,dpr;
  doublecomplex sum[3],tbuff[3];
  int i,j,jjj;
  double temp, na;
  double f[3][3];
  doublecomplex mult_mat[MAXNMAT];

  if (ScatRelation==SOrd) {
    /* calculate correction coefficient */
    temp=kd*kd/24;
    for(i=0;i<Nmat;i++) {
      na=DotProd(n,prop);
      cSquare(ref_index[i],m2);
      temp=kd*kd/24;
      mult_mat[i][re]=1-temp*(m2[re]-2*na*ref_index[i][re]+1);
      mult_mat[i][im]=temp*(2*na*ref_index[i][im]-m2[im]);       /* mult_mat=1-(kd^2/24)(m^2-2(n.a)m+1) */
    }
  }
  for(i=0;i<3;i++) sum[i][re]=sum[i][im]=0.0;

  for (j=0;j<local_nvoid_Ndip;++j) {
    jjj=3*j;
    /* kr=k*r.n */
    kr=WaveNum*DotProd(DipoleCoord+3*j,n);
    /* a=exp(-ikr.n) */
    cExp(-kr,a);
                          /* multiply by a correction coefficient */
    if (ScatRelation==SOrd) cMultSelf(a,mult_mat[material[j]]);
    /* sum(P*exp(-ik*r.n)) */
    for(i=0;i<3;i++) {
      sum[i][re]+=p[jjj+i][re]*a[re]-p[jjj+i][im]*a[im];
      sum[i][im]+=p[jjj+i][re]*a[im]+p[jjj+i][im]*a[re];
    }
  } /* end for j */
  /* ebuff=(I-nxn).sum=sum-n*(n.sum) */
  crDotProd(sum,n,dpr);
  cScalMultRVec(n,dpr,tbuff);
  cvSubtr(sum,tbuff,ebuff);

  /* multiply it by (-i*k^3) */
  kkk=WaveNum*WaveNum*WaveNum;
  for(i=0;i<3;i++) {
    temp=ebuff[i][re];
    ebuff[i][re]=ebuff[i][im]*kkk;
    ebuff[i][im]=-temp*kkk;
  }
}

/*=====================================================================*/

double Ext_cross(double *incPol)
   /* Calculate the Extinction cross-section */
{
  doublecomplex ebuff[3],tmp;
  double sum;
  int i;

  calc_field (ebuff,prop);
  crDotProd(ebuff,incPol,tmp);                      /* incPol is real, so no conjugate is needed */
  sum=tmp[re];
  my_inner_product(&sum,double_type,1);
  return 4*PI*sum/(WaveNum*WaveNum);
}

/*=====================================================================*/

double Abs_cross(void)
  /* Calculate the Absorption cross-section for process 0 */
{
  int dip,index,i,j;
  char mat;
  double sum, dummy, tmp, temp1,temp2;
  doublecomplex m2;
  double *m; /* not doublecomplex=double[2] to allow assignment to it */
  double cc_inv_im[MAXNMAT][3];   /* -Im(1/cc)=Im(cc)/|cc|^2 */
  double mult_mat[MAXNMAT];

  if (ScatRelation==DRAINE) {
    /* calculate constant and cc_inv_im */
    dummy = 2*WaveNum*WaveNum*WaveNum/3;
    for (i=0;i<Nmat;i++) for (j=0;j<3;j++) cc_inv_im[i][j]=cc[i][j][im]/cAbs2(cc[i][j]);
    /* main cycle */
    for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip) {
      mat=material[dip];
      index=3*dip;
      /* Im(P.Eexc(*))-(2/3)k^3*|P|^2=|P|^2*(-Im(1/cc)-(2/3)k^3) */
      for(i=0;i<3;i++) sum+=(cc_inv_im[mat][i] - dummy)*cAbs2(p[index+i]);
    }
  }
  else if (ScatRelation==SOrd) {
    /* calculate constants */
    temp1=kd*kd/6;
    temp2=4*PI/(gridspace*gridspace*gridspace);
    for (i=0;i<Nmat;i++) {
      m=ref_index[i];
      cSquare(m,m2);
      m2[re]-=1;
        /* mult_mat=-Im(1/hi)*(1+(kd*Im(m))^2)/d^3;  hi=(m^2-1)/(4*PI)  */
      mult_mat[i]=temp2*m2[im]*(1+temp1*m[im]*m[im])/cAbs2(m2);
    }
    /* main cycle */
    for (dip=0,sum=0;dip<local_nvoid_Ndip;++dip)
      sum+=mult_mat[material[dip]]*cvNorm2(p+3*dip);
  }
  my_inner_product(&sum,double_type,1);
  return 4*PI*WaveNum*sum;
}

/*=====================================================================*/

void set_Parms(void)
   /* Set integration parameters for asymmetry-paramter & C_sca */
{
  int arg;
  FILE *input;
  double dummy[3];

  if (ringid==ROOT) {
    /* default Robm_parms = "Romb5.dat", as set in ADDAmain.c */
    /* use -Romb_input 'name' to change name of input file */
    if ((input=fopen(Romb_parms,"r"))==NULL)
      LogError(EC_ERROR,ONE,POSIT,"Failed to open '%s'",Romb_parms);

    for (arg=0;arg<Narg;++arg) {
      if (arg == THETA) 
        fscanf(input,
              "theta\tmin\t%lg\n\tmax\t%lg\n\tJMAX\t%d\n"\
              "\tepsilon\t%lg\n\tK\t%d\n",
              &dummy[0],&dummy[1],&(parms[arg].JMAX),
              &dummy[2],&(parms[arg].K));
      else 
        fscanf(input,
              "phi\tmin\t%lg\n\tmax\t%lg\n\tJMAX\t%d\n"\
              "\tepsilon\t%lg\n\tK\t%d\n",
              &dummy[0],&dummy[1],&(parms[arg].JMAX),
              &dummy[2],&(parms[arg].K));

      parms[arg].min = (double) dummy[0];
      parms[arg].max = (double) dummy[1];
      parms[arg].INT_EPS = (double) dummy[2];
      parms[arg].Grid_size = (1 << (parms[arg].JMAX-1)) + 1;
    }
    fclose(input);
  }
  
  Bcast_parms(parms);
  fprintz(logfile,
          "Integration parameters\n"\
          "\t\t\t\tPHI\tTHETA\n"\
          "\tEPS\t\t\t%.8lg\t%.8lg\n"\
          "\tMaximal number of\t%d\t%d\n\trefinement-stages\n"\
          "\tNumber of evaluations\t%d\t%d\n\tfor an extrapolation\n"\
          "\tlower boundary\t\t%.8lg\t%.8lg\n"\
          "\tupper boundary\t\t%.8lg\t%.8lg\n",
          (double) parms[PHI].INT_EPS,(double) parms[THETA].INT_EPS,
          parms[PHI].JMAX,parms[THETA].JMAX,
          parms[PHI].K,parms[THETA].K,
          (double) parms[PHI].min,(double) parms[THETA].min,
          (double) parms[PHI].max,(double) parms[THETA].max);
}

/*=====================================================================*/

void fill_tab(void)
   /* calculate in advance the sines and cosines
    * for calculating E field for all_dir
    * and integrating the asymmetry-parameter
    *     and the scattering cross-section */
{
  int i,arg;
  double d;
 
  for (arg=0;arg<Narg;++arg) {
    /* Memory-allocation */
    if ((si[arg] = (double *) malloc(parms[arg].Grid_size*sizeof(double)))==NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc si[%d]",arg);
    if ((co[arg] = (double *) malloc(parms[arg].Grid_size*sizeof(double)))==NULL)
      LogError(EC_ERROR,ALL,POSIT,"Could not malloc si[%d]",arg);

    d = (parms[arg].max - parms[arg].min)/((double)parms[arg].Grid_size-1);
    for (i=0;i<parms[arg].Grid_size;++i) {
      if (arg==THETA) {
        co[arg][i] = parms[arg].min+d*i;
        si[arg][i] = sqrt(1-co[arg][i]*co[arg][i]);
      }
      else {
        si[arg][i] = sin(parms[arg].min+d*i);
        co[arg][i] = cos(parms[arg].min+d*i);
      }
    }
  }
}

/*=====================================================================*/

void finish_int(void)
{
  free(si[PHI]);
  free(co[PHI]);
  free(si[THETA]);
  free(co[THETA]);
}

/*=====================================================================*/
                 
void calc_alldir(char which)
   /* calculate scattered field in many directions */  
{
  int index,npoints,point,theta,phi;
  clock_t tstart;
  double robserver[R3];
  char name[200];
  FILE *scatfile;
  double d_cos_theta, dphi;

  npoints = parms[THETA].Grid_size*parms[PHI].Grid_size;
  fill_tab();
  
  /* Calculate field */
  tstart = clock();
  for (theta=0,point=0;theta<parms[THETA].Grid_size;++theta)
    for (phi=0;phi<parms[PHI].Grid_size;++phi,++point)
      {
	robserver[X] = si[THETA][theta]*co[PHI][phi];
	robserver[Y] = si[THETA][theta]*si[PHI][phi];
	robserver[Z] = co[THETA][theta];
	calc_field(Egrid+grid_index(theta,phi),robserver);
	if ((point % (npoints/10)) == 0)
          {
	    printz("Calculated %d percent of the scattered field\n",
		   100*point/npoints);
            fflushz(stdout);
          }
      }
  Timing_calc_EField_ad = clock() - tstart;
  tstart = clock();
  accumulate((double *)Egrid,2*R3*npoints,Egrid_buffer);
  Timing_comm_EField_ad = clock() - tstart;
  Timing_EField_ad = Timing_calc_EField_ad + Timing_comm_EField_ad;
  Timing_EField += Timing_EField_ad;

  /* Store in file */
  if (ringid==ROOT && store_all_dir)
    {
      tstart = clock();

      strcpy(name,directory);

      if (which == 'X') strcat(name,"/EgridX");
      else if (which == 'Y') strcat(name,"/EgridY");

      if ((scatfile=fopen(name,"w"))==NULL)      
        LogError(EC_ERROR,ONE,POSIT,"Failed to open '%s'",name); 


      d_cos_theta = (parms[THETA].max - parms[THETA].min)/((double)parms[THETA].Grid_size-1);
      dphi = rad2deg((parms[PHI].max - parms[PHI].min)/((double)parms[PHI].Grid_size-1));

      for (theta=0;theta<parms[THETA].Grid_size;++theta)
	for (phi=0;phi<parms[PHI].Grid_size;++phi)
	  {
	    index = grid_index(theta,phi);
	    fprintf(scatfile,"%.4f\t %.4f\t %.7e\t %.7e\t %.7e\t %.7e\t %.7e\t %.7e\n",
		    rad2deg(acos(theta*d_cos_theta+parms[THETA].min)),phi*dphi,
                    (double) Egrid[index+X][re],(double) Egrid[index+X][im],
                    (double) Egrid[index+Y][re],(double) Egrid[index+Y][im],
                    (double) Egrid[index+Z][re],(double) Egrid[index+Z][im]);
	  }
      fclose(scatfile);
	  
      Timing_FileIO += clock() - tstart;
    }
}

/*=====================================================================*/

void C_sca_integrand(int theta,int phi, double *res)
{
  res[0]=cvNorm2(&Egrid[grid_index(theta,phi)]);    /* this norm should be calculated once for all */
  						    /* functions that use it   */
}

/*=====================================================================*/

double Sca_cross(void)
  /* Calculate the scattering cross section 
   * from the integral */
{
  clock_t tstart;
  char fname[200];
  double res;

  sprintf(fname,"%s/log_Csca_int.dat",directory);

  tstart = clock();
  Romberg2D(parms,C_sca_integrand,1,&res,fname);
  res*=4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
  return res;
}

/*=====================================================================*/

void g_integrand(int theta,int phi,double *res)
{
  double E_square;
  
  E_square=cvNorm2(&Egrid[grid_index(theta,phi)]);
  res[X] = E_square*si[THETA][theta]*co[PHI][phi];
  res[Y] = E_square*si[THETA][theta]*si[PHI][phi];
  res[Z] = E_square*co[THETA][theta];
}
 
/*=====================================================================*/

void Asym_parm(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */
{
  int comp;
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_int.dat",directory);

  tstart = clock();
  Romberg2D(parms,g_integrand,3,vec,log_int);
  for (comp=0;comp<3;++comp) vec[comp]*=4*PI/(WaveNum*WaveNum); 
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

void g_x_integrand(int theta,int phi,double *res)
{
  res[0]=cvNorm2(&Egrid[grid_index(theta,phi)])*si[THETA][theta]*co[PHI][phi];
}

/*=====================================================================*/

void Asym_parm_x(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */
{
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_x_int.dat",directory);
 
  tstart = clock();
  Romberg2D(parms,g_x_integrand,1,vec,log_int);
  vec[0] *= 4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

void g_y_integrand(int theta,int phi,double *res)
{
  res[0]=cvNorm2(&Egrid[grid_index(theta,phi)])*si[THETA][theta]*si[PHI][phi];
}

/*=====================================================================*/

void Asym_parm_y(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */
{
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_y_int.dat",directory);
 
  tstart = clock();
  Romberg2D(parms,g_y_integrand,1,vec,log_int);
  vec[0] *= 4*PI/(WaveNum*WaveNum);
  Timing_Integration += clock() - tstart;
}
/*=====================================================================*/

void g_z_integrand(int theta,int phi,double *res)
{
  res[0]=cvNorm2(&Egrid[grid_index(theta,phi)])*si[THETA][theta]*co[PHI][phi];
}

/*=====================================================================*/

void Asym_parm_z(double *vec)
  /* Calculate the unnormalized asymmetry parameter,
   * i.e. not yet normalized by Csca */

{
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_z_int.dat",directory);
 
  tstart = clock();
  Romberg2D(parms,g_z_integrand,1,vec,log_int);
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
  int j,l,i,comp,index,jjj;
  int local_nvoid_d0, local_nvoid_d1;
  int *nvoid_array;
  char *materialT;
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
    dummy[re]=dummy[im]=0.0;
    for (comp=0;comp<3;++comp) {
      index = 3*j+comp;
      /* Im(P.E*inc) */
      _E_inc[re] = Einc[index][re];
      _E_inc[im] = -Einc[index][im];
      cMult(p[index],_E_inc,temp);
      cAdd(dummy,temp,dummy);
    }
    Finc[3*j+Z] = WaveNum*dummy[im]/2;
    Finc_tot[Z] += Finc[3*j+Z];
  }

  /* Because of the parallelisation by row-block decomposition
     the distributed arrays involved need to be gathered on each node
     a) material -> materialT
     b) DipoleCoord -> rdipT
     c) p -> pT
     */
  /* initialize local_nvoid_d0 and local_nvoid_d1 */
  nvoid_array=ivector(0,nprocs-1);
  nvoid_array[ringid]=local_nvoid_Ndip;
  all_gather(nvoid_array+ringid,nvoid_array,int_type,nprocs);
  local_nvoid_d0=0;
  for (i=0;i<ringid;i++) local_nvoid_d0+=nvoid_array[i];
  local_nvoid_d1=local_nvoid_d0+local_nvoid_Ndip;
  free(nvoid_array);
  /* requires a lot of additional memory */
  if ((materialT = (char *) malloc(nvoid_Ndip*sizeof(char)))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc materialT");
  if ((rdipT = (double *) malloc(3*nvoid_Ndip*sizeof(double)))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc rdipT");
  if ((pT = (doublecomplex *) malloc(3*nvoid_Ndip*sizeof(doublecomplex)))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc pT");

  memcpy(materialT+local_nvoid_d0,material,local_nvoid_Ndip*sizeof(char));
  memcpy(pT+3*local_nvoid_d0,p,3*local_nvoid_Ndip*sizeof(doublecomplex));
  memcpy(rdipT+3*local_nvoid_d0,DipoleCoord,3*local_nvoid_Ndip*sizeof(double));

  all_gather(materialT+local_nvoid_d0,materialT,int_type,local_nvoid_Ndip);
  all_gather(pT+3*local_nvoid_d0,pT,cmplx_type,3*local_nvoid_Ndip);
  all_gather(rdipT+3*local_nvoid_d0,rdipT,double_type,3*local_nvoid_Ndip);

  /* Calculate scattering force per dipole */
  for (j=local_nvoid_d0;j<local_nvoid_d1;++j) {
    int jjj = 3*j;

    for (l=0;l<nvoid_Ndip;++l) if (j!=l) {
      lll = 3*l;
      r2 = 0;
      Pn_j[re]=Pn_j[im]=Pn_l[re]=Pn_l[im]=inp[re]=inp[im]=0.0;

      /* Set distance related variables */
      for (comp=0;comp<3;++comp) {
        n[comp][im] = 0;
        n[comp][re] = rdipT[jjj+comp] - rdipT[lll+comp];
	r2 += n[comp][re]*n[comp][re];
      }
      r = sqrt(r2);
      n[X][re]/=r; n[Y][re]/=r; n[Z][re]/=r;

      /* Set the scalar products a.b1 and a.b2 */
      a[re] = cos(WaveNum*r);
      a[im] = sin(WaveNum*r);
      ab1[re] = 3/(r2*r2) - WaveNum*WaveNum/r2;
      ab2[re] = -WaveNum*WaveNum/r2;
      ab1[im] = -3*WaveNum/(r*r2);
      ab2[im] = WaveNum*WaveNum*WaveNum/r;
      cMultSelf(ab1,a);
      cMultSelf(ab2,a);

      /* Prepare c1 and c2 */
      for (comp=0;comp<3;++comp) {
	x_cg[comp][re] = pT[jjj+comp][re];
	x_cg[comp][im] = -pT[jjj+comp][im];
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
	c1[comp][re] *= -5;
	c1[comp][im] *= -5;

	cMult(inp,n[comp],temp);
	cAdd(c1[comp],temp,c1[comp]);
	cMult(Pn_j,pT[lll+comp],temp);
	cAdd(c1[comp],temp,c1[comp]);
	cMult(x_cg[comp],Pn_l,temp);
	cAdd(c1[comp],temp,c1[comp]);

	/* Set c2 */
	cMult(Pn_j,Pn_l,temp);
	cMult(n[comp],temp,c2[comp]);
        c2[comp][re] *= -1;
	c2[comp][im] *= -1;

	cMult(inp,n[comp],temp);
	cAdd(c2[comp],temp,c2[comp]);

	/* Fsca_{jl} = ... */
	cMultSelf(c1[comp],ab1);
	cMultSelf(c2[comp],ab2);
	Fsca[jjj-3*local_d0+comp] += (c1[comp][re] + c2[comp][re])/2;
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
  my_inner_product(Finc_tot+Z,double_type,1);
  my_inner_product(Fsca_tot,double_type,3);
  my_inner_product(Frp_tot,double_type,3);

  free(materialT);
  free(rdipT);
  free(pT);
}

/*=====================================================================*/
char *scan_parms(char *buf, integr_parms *a, Parms_1D *b, int ifcos)

{
  int i;
  double unit;
  
  if ((buf=strstr(buf,"min="))==NULL) return NULL;
  if (sscanf(buf,"min=%lf",&(a->min))!=1) return NULL;
  if ((buf=strstr(buf,"max="))==NULL) return NULL;
  if (sscanf(buf,"max=%lf",&(a->max))!=1) return NULL;
  if ((buf=strstr(buf,"JMAX="))==NULL) return NULL;
  if (sscanf(buf,"JMAX=%d",&(b->JMAX))!=1) return NULL;
  if ((buf=strstr(buf,"K="))==NULL) return NULL;
  if (sscanf(buf,"K=%d",&(b->K))!=1) return NULL;
  if ((buf=strstr(buf,"eps="))==NULL) return NULL;
  if (sscanf(buf,"eps=%lf",&(b->INT_EPS))!=1) return NULL;

  if (a->min==a->max) {
    a->N=b->Grid_size=1;
    b->JMAX=1;
  }
  else a->N=b->Grid_size=(1 << (b->JMAX-1)) + 1;
  /* initialize points of integration */
  if ((a->val=(double *) malloc(a->N*sizeof(double)))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"Could not malloc integration array"); 
  if (ifcos) {                          /* make equal intervals in cos(angle) */
    b->min=cos(deg2rad(a->max));
    if (fabs(b->min)<1e-15) b->min=0; /* just for convenience */
    b->max=cos(deg2rad(a->min));
    if (a->N==1) a->val[0]=a->min;
    else {
      unit = (b->max - b->min)/(a->N-1);
      for (i=0;i<a->N;i++) a->val[i] = rad2deg(acos(b->min+unit*i));
    }
  }
  else {			/* make equal intervals in angle */
    b->min=deg2rad(a->min);
    b->max=deg2rad(a->max);
    if (a->N==1) a->val[0]=a->min;
    else {
      unit = (a->max - a->min)/(a->N-1);
      for (i=0;i<a->N;i++) a->val[i] = a->min + unit*i;
    }
  }

  return buf;
}

/*=====================================================================*/

void read_avg_parms(char *fname) 
  /* read parameters of orientation averaging from a file */
{
  FILE *input;
  char *buf,*temp,*begin;
  unsigned int buf_size=500;
  extern char avg_string[];

  if ((begin=buf=(char *)malloc(buf_size*sizeof(char))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc buf");
  if ((temp=(char *)malloc(buf_size*sizeof(char))) == NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc temp");

  if ((input=fopen(fname,"r"))==NULL) 
    LogError(EC_ERROR,ONE,POSIT,"Failed to open file '%s'",fname);

  strcpy(buf,"");
  while(!feof(input)) {
    fgets(temp,buf_size,input);
    if(*temp!='#') {
      if (buf_size < strlen(buf) + strlen(temp))
        LogError(EC_ERROR,ONE,POSIT,"Buffer overflow while reading '%s'",fname);
      strcat(buf,temp);
    }
  }
  if((buf=strstr(buf,"alpha"))==NULL) LogError(EC_ERROR,ONE,POSIT,"Wrong format of file '%s'",fname);
  if((buf=scan_parms(buf,&alpha_int,&parms_alpha,false))==NULL) 
    LogError(EC_ERROR,ONE,POSIT,"Wrong format of file '%s'",fname);
  
  if((buf=strstr(buf,"beta"))==NULL) LogError(EC_ERROR,ONE,POSIT,"Wrong format of file '%s'",fname);
  if((buf=scan_parms(buf,&beta_int,&parms[THETA],true))==NULL) 
    LogError(EC_ERROR,ONE,POSIT,"Wrong format of file '%s'",fname);
  
  if((buf=strstr(buf,"gamma"))==NULL) LogError(EC_ERROR,ONE,POSIT,"Wrong format of file '%s'",fname);
  if((buf=scan_parms(buf,&gamma_int,&parms[PHI],false))==NULL)
    LogError(EC_ERROR,ONE,POSIT,"Wrong format of file '%s'",fname);
  
  if (equivalent && alpha_int.N>1) alpha_int.N--;
  if (equivalent && gamma_int.N>1) gamma_int.N--;

  free(begin);
  free(temp);
 
  sprintz(avg_string,
    "alpha: from %g to %g in %d steps\n"\
    "beta: from %g to %g in (up to) %d steps (equally spaced in cosine values)\n"\
    "gamma: from %g to %g in (up to) %d steps\n"\
    "see file 'romb_log' for details\n",
    alpha_int.min,alpha_int.max,alpha_int.N,beta_int.min,beta_int.max,beta_int.N,
    gamma_int.min,gamma_int.max,gamma_int.N);
  D("read_avg_parms complete");
}

/*=====================================================================*/

void free_parms(void)
  /* free values of angles for integration */
{
  free(alpha_int.val);
  free(beta_int.val);
  free(gamma_int.val);
}
