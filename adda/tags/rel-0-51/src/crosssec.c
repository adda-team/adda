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
#include <time.h>
#include <math.h>
#include "cmplx.h"
#include "types.h"
#include "const.h"
#include "Romberg.h"
#include "crosssec.h"
#include "comm.h"
#include "debug.h"

#define X 0
#define Y 1
#define Z 2
#define R3 3

/* Convert the (Theta,Phi) tuple into a linear array index */
/*#define grid_index(T,P) (R3*((T)*parms[THETA].Grid_size + (P)))original code Martijn, wrong!!!*/
#define grid_index(T,P) (R3*((T)*parms[PHI].Grid_size + (P)))

doublecomplex **dCmatrix (int nrl, int nrh, int ncl, int nch);

double *si[Narg],*co[Narg];
Parms_1D parms[Narg];
double *n;

extern doublecomplex *Egrid;         /* defined in main */
extern doublecomplex *x;		/* defined in main */
extern doublecomplex *Einc;          /* defined in main */
extern int *material;           /* defined calculator */
extern int Nmat;                /* defined in main */
extern doublecomplex cc[10][3]; /* defined in main */
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
extern char Romb_parms[200];
extern doublecomplex ref_index[10];
extern double kd;
extern double prop[3], incPolX[3], incPolY[3];
extern int ScatRelation;

double beta_matr[3][3];

extern double alph_deg,bet_deg,gam_deg;
double alph,bet,gam;     /* in radians */

integr_parms alpha_int, beta_int, gamma_int;
Parms_1D parms_alpha;

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

void calc_field (doublecomplex *x,   /* internal fields */
 	    doublecomplex *ebuff,    /* where to write calculated scattering amplitude */
            double **rdip,           /* location of dipoles */
	    double *n,               /* scattering direction */
            double kk,               /* k, wave number */
	    int nlocalDip            /* number of dipoles on current processor */
	    )
     
{
  /*  Near-optimal routine to compute the scattered fields at one specific
   *  angle. (more exactly - scattering amplitude)
   *
   *  Michel Grimminck 1995
   *  Changed by Maxim Yurkin 2005
   */
  
  double kr,rt00, rt01, rt02, rt11, rt12, rt22, kkk;
  doublecomplex a,m2;
  doublecomplex sum[3];
  int k,l,m,j,jjj;
  double temp, na;
  double f[3][3];
  doublecomplex tbuff[3];

  
  /* calculate projection matrix (I-nxn) */
  rt00=n[1]*n[1]+n[2]*n[2];
  rt11=n[0]*n[0]+n[2]*n[2];
  rt22=n[0]*n[0]+n[1]*n[1];
  rt01=-n[1]*n[0];
  rt02=-n[2]*n[0];
  rt12=-n[2]*n[1];
  
  for(m=0;m<Nmat;m++) { /* for each refraction index */
    for(k=0;k<3;k++) sum[k].r=sum[k].i=0.0;
    
    for (j = 0; j < nlocalDip; ++j) 
      if (material[j]==m) { /* for each dipole */
	/* calculate field in cartesian coordinates */
	jjj = 3 * j;
	
	/* calculate phase of field */
	
	kr=kk*DotProd(rdip[j],n);
	/* kr=k*r.n */
      
	a.r = cos (kr);         /* a=exp(-ikr.n) */
	a.i = - sin (kr);
	
	for(k=0;k<3;k++) {
	  sum[k].r+=x[jjj+k].r*a.r-x[jjj+k].i*a.i;      /* sum(E*exp(-ik*r.n)) */
	  sum[k].i+=x[jjj+k].r*a.i+x[jjj+k].i*a.r;
	}
      } /* end for j */

    /* multiply with coupleconstant */
    for (k = 0; k < 3; ++k) cMult(sum[k],cc[m][k],&(sum[k]));

    f[0][0] = rt00;
    f[0][1] = f[1][0] = rt01;
    f[0][2] = f[2][0] = rt02;   /* Construct full matrix */
    f[1][1] = rt11;             /* f=I-nxn               */
    f[1][2] = f[2][1] = rt12;  
    f[2][2] = rt22;
    
    for (k = 0; k < 3; ++k) {
      tbuff[k].r=tbuff[k].i=0.0;
      for (l = 0; l < 3; ++l) {
	tbuff[k].r += f[k][l] * sum[l].r;      /* f.sum */
	tbuff[k].i += f[k][l] * sum[l].i;
      }
    }
    kkk=kk*kk*kk;
    if (ScatRelation==SOrd) {   /* multiply by a correction coefficient */
      na=DotProd(n,prop);
      cSquare(ref_index[m],&m2);
      temp=kd*kd/24;
      a.r=1-temp*(m2.r-2*na*ref_index[m].r+1);
      a.i=temp*(2*na*ref_index[m].i-m2.i);       /* a=1-(kd^2/24)(m^2-2(n.a)m+1) */
      for(k = 0; k < 3; ++k) cMult(tbuff[k],a,&(tbuff[k]));
    }
    /* multiply it by (-i*k^3); add it to the E field */
    for(k = 0; k < 3; ++k) {
      kkk=kk*kk*kk;
      ebuff[k].r+=tbuff[k].i*kkk;
      ebuff[k].i-=tbuff[k].r*kkk;
    }
  } /* end m */
}

/*=====================================================================*/

int n_non_void_sites(/* gather the number of non-void sites,
                        i.e. actual dipoles, in the particle */
                     int mat_count[],
                     int Nmat)
{       
  double ndips;

  ndips = (double) local_Ndip - mat_count[Nmat-1];
  my_inner_product(&ndips);

  return (int) ndips;
}

/*=====================================================================*/

double Ext_cross(/* Calculate the Extinction cross-section */
		 doublecomplex *x,   /* dipole-moments */
		 double *incPol,    /* incident polarization */
		  double **rdip,    
		 double k)           /* wave-number incoming wave */
{
  doublecomplex ebuff[3],tmp;
  double sum;

  ebuff[0].r = ebuff[0].i = 0.0;
  ebuff[1].r = ebuff[1].i = 0.0;
  ebuff[2].r = ebuff[2].i = 0.0;
  calc_field (x,ebuff,rdip,prop,k,local_Ndip);
  crDotProd(ebuff,incPol,&tmp);                      /* incPol is real, so no conjugate is needed */
  sum=tmp.r;
  my_inner_product(&sum);
  return 4*PI*sum/(k*k);
}

/*=====================================================================*/

double Abs_cross(/* Calculate the Absorption cross-section for proces 0 */ 
		 doublecomplex *x,     /* dipole-moment */
		 double k)        /* wave-number incoming wave */
{
  int dip,index,mat,i;
  double sum, dummy, tmp, temp;
  doublecomplex m, m2;

  dummy = 2*k*k*k/3;
  temp = kd*kd/6;
  for (dip=0,sum=0;dip<local_Ndip;++dip)
    {
      mat=material[dip];
      if (mat<Nmat-1)  {
        index=3*dip;
        if (ScatRelation==DRAINE) 
          for(i=0; i<3; i++) sum+=(cc[mat][i].i - dummy*cAbs2(cc[mat][i]))*cAbs2(x[index+i]);
        else if (ScatRelation==SOrd) {
          m=ref_index[mat];
          cSquare(m,&m2);
          m2.r-=1;
          for(i=0,tmp=0; i<3; i++) tmp+=cAbs2(cc[mat][i])*cAbs2(x[index+i]);
          tmp*=4*PI*(-m2.i)*(1+temp*m.i*m.i)/cAbs2(m2);
          sum+=tmp;
        }
      }
   }
  my_inner_product(&sum);
  return 4*PI*k*sum;
}

/*=====================================================================*/

void set_Parms(/* Set integration parameters for asymmetry-paramter & C_sca */
	       void)
{
  int 
    arg;
  clock_t
    tstart;
  char
    *name[Narg],
    *buffer;
  FILE 
    *input;

  /*name[THETA] = "theta";
  name[PHI] = "phi";*/

  if (ringid == 0)
    {
      tstart = clock();
      /* default Robm_parms = "Romb5.dat", as set in CD2main.c */
      /* use -Romb_input 'name' to change name of input file */
      if ((input=fopen(Romb_parms,"r"))==NULL)
        LogError(EC_ERROR,ONE,POSIT,"Failed to open '%s'",Romb_parms); 

      for (arg=0;arg<Narg;++arg)
        {
	  double
	    dummy[3];
	
        if (arg == THETA)
          {
            fscanf(input,
                  "theta\tmin\t%lg\n\tmax\t%lg\n\tJMAX\t%d\n"\
                  "\tepsilon\t%lg\n\tK\t%d\n",
                  &dummy[0],&dummy[1],&(parms[arg].JMAX),
                  &dummy[2],&(parms[arg].K));
          }
        else
          {
            fscanf(input,
                  "phi\tmin\t%lg\n\tmax\t%lg\n\tJMAX\t%d\n"\
                  "\tepsilon\t%lg\n\tK\t%d\n",
                  &dummy[0],&dummy[1],&(parms[arg].JMAX),
                  &dummy[2],&(parms[arg].K));
          }
 
	
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

void fill_tab(/* calculate in advance the sines and cosines
	       * for integrating the asymmetry-parameter
	       * and the scattering cross-section */
	      void)
{
  int
    i,arg;
  double
    d;
 
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
                 
void calc_alldir(doublecomplex *x,   /* dipole moments */
		 double **rdip,   /* dipole coordinates */
		 double k,         /* wave number */
		 char which)     
{
  int
    comp,
    index,
    npoints,point,
    theta,phi;
  clock_t
    tstart;
  double robserver[R3];
  char 
    name[200],
    command[200];
  FILE *scatfile;
  double dtheta, dphi;

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
	Egrid[grid_index(theta,phi)].r=Egrid[grid_index(theta,phi)].i=
	  Egrid[grid_index(theta,phi)+1].r=Egrid[grid_index(theta,phi)+1].i=
	  Egrid[grid_index(theta,phi)+2].r=Egrid[grid_index(theta,phi)+2].i=0;
	calc_field(x,&Egrid[grid_index(theta,phi)],rdip, 
		   robserver,k,local_Ndip);
	if ((point % (npoints/10)) == 0)
          {
	    printz("Calculated %d percent of the scattered field\n",
		   100*point/npoints);fflush(stdout);
          }
      }
  Timing_calc_EField_ad = clock() - tstart;
  tstart = clock();
  accumulate((double *)Egrid,2*R3*npoints);
  Timing_comm_EField_ad = clock() - tstart;
  Timing_EField_ad = Timing_calc_EField_ad + Timing_comm_EField_ad;
  Timing_EField += Timing_EField_ad;

  /* Store in file */
  if ((ringid == 0) && store_all_dir)
    {
      tstart = clock();

      strcpy(name,directory);

      if (which == 'X') strcat(name,"/EgridX");
      else if (which == 'Y') strcat(name,"/EgridY");

      if ((scatfile=fopen(name,"w"))==NULL)      
        LogError(EC_ERROR,ONE,POSIT,"Failed to open '%s'",name); 
 
      dtheta = (180.0/PI)*(parms[THETA].max - parms[THETA].min)/((double)parms[THETA].Grid_size-1);
      dphi = (180.0/PI)*(parms[PHI].max - parms[PHI].min)/((double)parms[PHI].Grid_size-1);
      /* fprintf(scatfile,"theta\tphi\tEx.r\tEx.i\tEy.r\tEy.i\tEz.r\tEz.i\t\n\n"); */
      for (theta=0;theta<parms[THETA].Grid_size;++theta)
	for (phi=0;phi<parms[PHI].Grid_size;++phi)
	  {
	    index = grid_index(theta,phi);
	    fprintf(scatfile,"%.4f\t %.4f\t %.7e\t %.7e\t %.7e\t %.7e\t %.7e\t %.7e\n",
		    theta*dtheta,phi*dphi,
                    (double) Egrid[index+X].r,(double) Egrid[index+X].i,
                    (double) Egrid[index+Y].r,(double) Egrid[index+Y].i,
                    (double) Egrid[index+Z].r,(double) Egrid[index+Z].i);
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

void Sca_cross(/* Calculate the scattering cross section 
		* from the integral */
	       double k,        /* wave-number incoming wave */
	       double *res)     /* result */
{
  clock_t tstart;
  char fname[200];

  sprintf(fname,"%s/log_Csca_int.dat",directory);

  tstart = clock();
  Romberg2D(parms,&C_sca_integrand,1,res,fname);
  res[0]*=4*PI/(k*k);
  Timing_Integration += clock() - tstart;
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

void Asym_parm(/* Calculate the unnormalized asymmetry parameter,
                * i.e. not yet normalized by Csca */
               double k,
               double vec[3])
{
  int comp;
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_int.dat",directory);
 
  tstart = clock();
  Romberg2D(parms,&g_integrand,3,vec,log_int);
  for (comp=0;comp<3;++comp) vec[comp]*=4*PI/(k*k); 
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

void g_x_integrand(int theta,int phi,double *res)
{
  res[0]=cvNorm2(&Egrid[grid_index(theta,phi)])*si[THETA][theta]*co[PHI][phi];
}

/*=====================================================================*/

void Asym_parm_x(/* Calculate the unnormalized asymmetry parameter,
                  * i.e. not yet normalized by Csca */
                 double k,
                 double *vec)
{
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_x_int.dat",directory);
 
  tstart = clock();
  Romberg2D(parms,&g_x_integrand,1,vec,log_int);
  vec[0] *= 4*PI/(k*k);
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

void g_y_integrand(int theta,int phi,double *res)
{
  res[0]=cvNorm2(&Egrid[grid_index(theta,phi)])*si[THETA][theta]*si[PHI][phi];
}

/*=====================================================================*/

void Asym_parm_y(/* Calculate the unnormalized asymmetry parameter,
                  * i.e. not yet normalized by Csca */
                 double k,
                 double *vec)
{
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_y_int.dat",directory);
 
  tstart = clock();
  Romberg2D(parms,&g_y_integrand,1,vec,log_int);
  vec[0] *= 4*PI/(k*k);
  Timing_Integration += clock() - tstart;
}
/*=====================================================================*/

void g_z_integrand(int theta,int phi,double *res)
{
  res[0]=cvNorm2(&Egrid[grid_index(theta,phi)])*si[THETA][theta]*co[PHI][phi];
}

/*=====================================================================*/

void Asym_parm_z(/* Calculate the unnormalized asymmetry parameter,
                  * i.e. not yet normalized by Csca */
                 double k,
                 double *vec)
{
  clock_t tstart;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_z_int.dat",directory);
 
  tstart = clock();
  Romberg2D(parms,&g_z_integrand,1,vec,log_int);
  vec[0] *= 4*PI/(k*k);
  Timing_Integration += clock() - tstart;
}

/*=====================================================================*/

void Frp_mat(/* Calculate the Radiation Pressure by direct calculation
	        of the scattering force. Per dipole the force of 
		the incoming photons, the scattering force and the 
		radiation pressure are calculated as intermediate results */
	     double Fsca_tot[3],double *Fsca,
	     double Finc_tot[3],double *Finc,
	     double Frp_tot[3],double *Frp,
	     doublecomplex *x,double **rdip,double k)
{
  int 
    j,l,i,
    comp;
  int *materialT;
  double *rdipT;
  doublecomplex *xT;
  doublecomplex temp; 
 
  for (comp=0;comp<3;++comp)
    {
      Fsca_tot[comp] = Finc_tot[comp] = Frp_tot[comp] = 0;
    }

  /* Convert internal fields to dipole moments;
     Calculate incoming force per dipole */
  for (j=0;j<local_Ndip;++j)
    if (material[j]<Nmat-1)
      {
	doublecomplex dummy = {0,0};
	
	for (comp=0;comp<3;++comp)
	  {
	    int index = 3*j+comp;
	    doublecomplex cc_d,_E_inc;
	    
	    /* Conversion */
	    cMult(x[index],cc[material[j]][comp],&(x[index]));

	    /* Im(P.E*inc) */
	    _E_inc.r = Einc[index].r;
	    _E_inc.i = -Einc[index].i;
	    cMult(x[index],_E_inc,&temp);
	    cAdd(dummy,temp,&(dummy));
	  }

	Finc[3*j+Z] = k*dummy.i/2; 
	Finc_tot[Z] += Finc[3*j+Z];
      }

  /* Because of the parallelisation by row-block decomposition 
     the distributed arrays involved need to be gathered on each node 
     a) material -> materialT
     b) rdip -> rdipT
     c) x -> xT
     */
  if ((materialT = (int *) malloc(Ndip*sizeof(int)))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc materialT");
  if ((rdipT = (double *) malloc(3*Ndip*sizeof(double)))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc rdipT");
  if ((xT = (doublecomplex *) malloc(3*Ndip*sizeof(doublecomplex)))==NULL)
    LogError(EC_ERROR,ALL,POSIT,"Could not malloc xT");

  memcpy(&materialT[local_d0],material,local_Ndip*sizeof(int));
  memcpy(&xT[3*local_d0],x,3*local_Ndip*sizeof(doublecomplex));
  for (j=0;j<local_Ndip;++j)
    for (comp=0;comp<3;++comp)
      rdipT[3*(j+local_d0)+comp] = rdip[j][comp];

  all_gather(&materialT[local_d0],materialT,"int",local_Ndip);
  all_gather(&xT[3*local_d0],xT,"doublecomplex",3*local_Ndip);
  all_gather(&rdipT[3*local_d0],rdipT,"double",3*local_Ndip);

  /* Calculate scattering force per dipole */
  for (j=local_d0;j<local_d1;++j)
    if (materialT[j]<Nmat-1)
      {
	int jjj = 3*j;

	for (l=0;l<Ndip;++l)
	  if (materialT[l]<Nmat-1 && j!=l)
	    {
	      int 
		lll = 3*l,
		comp;
	      double
		r,r2 = 0;      /* (squared) absolute distance */
	      doublecomplex
		n[3],          /* unit vector in the direction of r_{jl}
				* complex will always be zero */ 
		a,ab1,ab2,     /* see chapter ... */
		c1[3],c2[3],   /* idem */
		x_cg[3],       /* complex conjungate P*_j */
		Pn_j = {0,0},  /* n_jl.P_l */
		Pn_l = {0,0},  /* P*_j.n_jl */
		inp = {0,0};   /* P*_j.P_l */

	      /* Set distance related variables */
	      for (comp=0;comp<3;++comp)
		{
		  n[comp].i = 0;
		  n[comp].r = rdipT[jjj+comp] - rdipT[lll+comp];
		  r2 += n[comp].r*n[comp].r;
		}
	      r = sqrt(r2);
	      n[X].r/=r; n[Y].r/=r; n[Z].r/=r;

	      /* Set the scalair products a.b1 and a.b2 */
	      a.r = cos(k*r);
	      a.i = sin(k*r);
	      ab1.r = 3/(r2*r2) - k*k/r2;  ab2.r = -k*k/r2;
	      ab1.i = -3*k/(r*r2);         ab2.i = k*k*k/r;
	      cMult(a,ab1,&ab1);           cMult(a,ab2,&ab2);

	      /* Prepare c1 and c2 */
	      for (comp=0;comp<3;++comp)
		{
		  x_cg[comp].r = xT[jjj+comp].r;
		  x_cg[comp].i = -xT[jjj+comp].i;
		  cMult(x_cg[comp],n[comp],&temp);
		  cAdd(Pn_j,temp,&Pn_j);
		  cMult(n[comp],xT[lll+comp],&temp);
		  cAdd(Pn_l,temp,&Pn_l);
		  cMult(x_cg[comp],xT[lll+comp],&temp);
		  cAdd(inp,temp,&inp);
		}

	      for (comp=0;comp<3;++comp)
		{
		  /* Set c1 */
		  cMult(Pn_j,Pn_l,&temp);
		  cMult(n[comp],temp,&(c1[comp]));
		  c1[comp].r *= -5;
		  c1[comp].i *= -5;

		  cMult(inp,n[comp],&temp);
		  cAdd(c1[comp],temp,&(c1[comp]));
		  cMult(Pn_j,xT[lll+comp],&temp);
		  cAdd(c1[comp],temp,&(c1[comp]));
		  cMult(x_cg[comp],Pn_l,&temp);
		  cAdd(c1[comp],temp,&(c1[comp]));

		  /* Set c2 */
	  	  cMult(Pn_j,Pn_l,&temp);
		  cMult(n[comp],temp,&(c2[comp]));
                  c2[comp].r *= -1;
		  c2[comp].i *= -1;		  

		  cMult(inp,n[comp],&temp);
		  cAdd(c2[comp],temp,&(c2[comp]));

		  /* Fsca_{jl} = ... */
		  cMult(ab1,c1[comp],&(c1[comp]));
		  cMult(ab2,c2[comp],&(c2[comp]));
		  Fsca[jjj-3*local_d0+comp] += (c1[comp].r + c2[comp].r)/2;
		}
	    } /* end j-loop */
	
	/* Concluding */
	for (comp=0;comp<3;++comp)
	  {
	    Fsca_tot[comp] += Fsca[jjj-3*local_d0+comp];

	    Frp[jjj-3*local_d0+comp] = Finc[jjj-3*local_d0+comp] + Fsca[jjj-3*local_d0+comp];

	    Frp_tot[comp] += Frp[jjj-3*local_d0+comp];
	  }
      } /* end l-loop */

  /* Accumulate the total forces on all nodes */
  my_inner_product(&Finc_tot[Z]);
  for (comp=0;comp<3;++comp) 
    { 
      my_inner_product(&Fsca_tot[comp]);
      my_inner_product(&Frp_tot[comp]);
    }

  free(materialT);
  free(rdipT);
  free(xT);
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
  int buf_size=500;
  extern char avg_string[];
  
  begin=buf=(char *)malloc(buf_size*sizeof(char));
  temp=(char *)malloc(buf_size*sizeof(char));

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

