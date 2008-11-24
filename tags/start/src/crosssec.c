/* FILE : crosssec.c
 * Author : Martijn Frijlink
 * Descr : 
 *   1) The functions Ext_cross, Abs_cross and Frp_mat are self-explanatory.
 *   2) The function set_parms reads the integration-parameters from file.
 *   3) The function calc_alldir calculates the scattered field in a 
 *      predefined set of directions (predefined by the integration-
 *      parameters). The function fill_tab calculates the trigonometric 
 *      function-values belonging to these directions.
 *   4) The functions Sca_cross calls the Romberg-routine for integrating
 *      the scattered intensity resulting in the scattering cross section.
 *   5) The functions Asym_parm_x, Asym_parm_y and Asym_parm_z each 
 *      calculate the corresponding component of the asymmetry-parameter,
 *      not yet divided by the scattering cross section. It can therefore
 *      also be interpreted as the scattering force in cross section units.
 */

#include <stdio.h>
#include <math.h>
#include "cmplx.h"
#include "types.h"
#include "const.h"
#include "Romberg.h"
#include "crosssec.h"
#include "comm.h"

#define X 0
#define Y 1
#define Z 2
#define R3 3

/* Convert the (Theta,Phi) tuple into a linear array index */
/*#define grid_index(T,P) (R3*((T)*parms[THETA].Grid_size + (P)))original code Martijn, wrong!!!*/
#define grid_index(T,P) (R3*((T)*parms[PHI].Grid_size + (P)))

dcomplex **dCmatrix (int nrl, int nrh, int ncl, int nch);

REAL *si[Narg],*co[Narg];
Parms_1D parms[Narg];
Rvector n;

extern dcomplex *Egrid;         /* defined in main */
extern dcomplex *x;		/* defined in main */
extern dcomplex *Einc;          /* defined in main */
extern int *material;           /* defined calculator */
extern int Nmat;                /* defined in main */
extern doublecomplex cc[10];    /* defined in main */
extern FILE *logfile;           /* defined in main */
extern char directory[200];     /* defined in main */
extern int store_all_dir;       /* defined in main */
extern unsigned long            /* defined in main */ 
  Timing_EField,
  Timing_calc_EField,
  Timing_comm_EField,
  Timing_Integration,
  Timing_FileIO;
extern char Romb_parms[200];

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

double Ext_cross(/* Calculate the Extinction cross-section */
		 dcomplex *x,     /* dipole-moments */
		 dcomplex *Einc,  /* incoming wave */
		 double k)        /* wave-number incoming wave */
{
  int
    dip,index;
  double
    sum;
  doublecomplex
    inpr;

  for (dip=0,sum=0;dip<local_Ndip;++dip)
    {
      if (material[dip]<Nmat-1)
	{
	  index=3*dip+X;
	  inpr.r = Einc[index].r*x[index].r + Einc[index].i*x[index].i;
	  inpr.i = Einc[index].i*x[index].r - Einc[index].r*x[index].i;

	  index=3*dip+Y;
	  inpr.r += Einc[index].r*x[index].r + Einc[index].i*x[index].i;
	  inpr.i += Einc[index].i*x[index].r - Einc[index].r*x[index].i;

	  index=3*dip+Z;
	  inpr.r += Einc[index].r*x[index].r + Einc[index].i*x[index].i;
	  inpr.i += Einc[index].i*x[index].r - Einc[index].r*x[index].i;

	  sum += cc[material[dip]].i*inpr.r - cc[material[dip]].r*inpr.i;
	}
    }
  my_inner_product(&sum);

  return 4*PI*k*sum;
}

double Sca_diff(/* Evaluate the formal expression for 
		   the difference between the integrated and 
		   direct scattering cross section (see scriptie) */
		dcomplex *x,     /* dipole-moments */
		dcomplex *r,  /* incoming wave */
		double k)        /* wave-number incoming wave */
{
  int
    dip,index;
  double
    sum;
  doublecomplex
    inpr;

  for (dip=0,sum=0;dip<local_Ndip;++dip)
    {
      if (material[dip]<Nmat-1)
	{
	  index=3*dip+X;
	  inpr.r = x[index].r*r[index].r + x[index].i*r[index].i;
	  inpr.i = x[index].i*r[index].r - x[index].r*r[index].i;

	  index=3*dip+Y;
	  inpr.r += x[index].r*r[index].r + x[index].i*r[index].i;
	  inpr.i += x[index].i*r[index].r - x[index].r*r[index].i;

	  index=3*dip+Z;
	  inpr.r += x[index].r*r[index].r + x[index].i*r[index].i;
	  inpr.i += x[index].i*r[index].r - x[index].r*r[index].i;

	  sum += cc[material[dip]].i*inpr.r + cc[material[dip]].r*inpr.i;
	}
    }
  my_inner_product(&sum);

  return -4*PI*k*sum;
}

double Abs_cross(/* Calculate the Absorption cross-section for proces 0 */ 
		 dcomplex *x,     /* dipole-moment */
		 double k)        /* wave-number incoming wave */
{
  int
    dip,
    index;
  double
    inpr,
    sum;
  double
    dummy;
    
  dummy = 2*k*k*k/3;

  for (dip=0,sum=0;dip<local_Ndip;++dip)
    if (material[dip]<Nmat-1)
      {
	index = 3*dip+X;
	inpr = x[index].r*x[index].r + x[index].i*x[index].i;
	index = 3*dip+Y;
	inpr += x[index].r*x[index].r + x[index].i*x[index].i;
	index = 3*dip+Z;
	inpr += x[index].r*x[index].r + x[index].i*x[index].i;

	index = material[dip];
	sum += inpr * (cc[index].i - 
	  dummy*(cc[index].r*cc[index].r + cc[index].i*cc[index].i));
      }

  my_inner_product(&sum);

  return 4*PI*k*sum;
}




void set_Parms(/* Set integration parameters for asymmetry-paramter & C_sca */
	       void)
{
  int 
    arg;
  unsigned long
    tstart,tstop;
  char
    *name[Narg],
    *buffer;
  FILE 
    *input;

  /*name[THETA] = "theta";
  name[PHI] = "phi";*/

  if (ringid == 0)
    {
      tstart = extime();
      /* default Robm_parms = "Romb5.dat", as set in CD2main.c */
      /* use -Romb_input 'name' to change name of input file */
      input = fopen(Romb_parms,"r");

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
 
	
	parms[arg].min = (REAL) dummy[0];
	parms[arg].max = (REAL) dummy[1];
	parms[arg].INT_EPS = (REAL) dummy[2];
        parms[arg].Grid_size = (1 << parms[arg].JMAX) + 1;
      }
 
      fclose(input);
      tstop = extime();
      Timing_FileIO += tstop - tstart;
    }
  
  Bcast_parms(parms);
  fprintz(logfile,
          "Integration parameters\n"\
          "\t\t\t\tPHI\tTHETA\n"\
          "\tEPS\t\t\t%lg\t%lg\n"\
          "\tMaximal number of\t%d\t%d\n\trefinement-stages\n"\
          "\tNumber of evaluations\t%d\t%d\n\tfor an extrapolation\n"\
          "\tlower boundary\t\t%lg\t%lg\n"\
          "\tupper boundary\t\t%lg\t%lg\n",
          (double) parms[PHI].INT_EPS,(double) parms[THETA].INT_EPS,
          parms[PHI].JMAX,parms[THETA].JMAX,
          parms[PHI].K,parms[THETA].K,
          (double) parms[PHI].min,(double) parms[THETA].min,
          (double) parms[PHI].max,(double) parms[THETA].max);
}

void fill_tab(/* calculate in advance the sines and cosines
	       * for integrating the asymmetry-parameter
	       * and the scattering cross-section */
	      void)
{
  int
    i,arg;
  REAL
    d;
 
  for (arg=0;arg<Narg;++arg)
    {
      /* Memory-allocation */
      si[arg] = (REAL *) malloc(parms[arg].Grid_size*sizeof(REAL));
      co[arg] = (REAL *) malloc(parms[arg].Grid_size*sizeof(REAL));

      d = (parms[arg].max - parms[arg].min)/((REAL)parms[arg].Grid_size-1);
      for (i=0;i<parms[arg].Grid_size;++i)
	{
          si[arg][i] = sin(d*i);
          co[arg][i] = cos(d*i);
        }
    }
}

void finish_int(void)
{
  free(si[PHI]);
  free(co[PHI]);
  free(si[THETA]);
  free(co[THETA]);
}

void calc_alldir(dcomplex *x,   /* dipole moments */
		 REAL **rdip,   /* dipole coordinates */
		 double k, char which)      /* wave number */
{
  int
    comp,
    index,
    npoints,point,
    theta,phi;
  unsigned long
    tstart,tstop;
  double robserver[R3];
  dcomplex **fint;
  char 
    name[200],
    command[200];
  FILE *scatfile;
  REAL dtheta, dphi;
  
  if ((fint = dCmatrix (0, 2, 0, 2)) == NULL) {
    LogError (EC_ERROR, "crosssec.c",
              "ringID %d, could not malloc fint",
	      ringid);
    AbortServer (1);
  }

  npoints = parms[THETA].Grid_size*parms[PHI].Grid_size;
  fill_tab();
  
  /* Calculate field */
  tstart = extime();
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
		   robserver,k,fint,local_Ndip);
	if ((point % (npoints/10)) == 0)
          {
	    printz("Calculated %d percent of the scattered field\n",
		   100*point/npoints);fflush(stdout);
          }
      }
  tstop = extime();
  Timing_calc_EField = tstop - tstart;
  tstart = extime();
  accumulate((REAL *)Egrid,2*R3*npoints);
  tstop = extime();
  Timing_comm_EField = tstop - tstart;
  Timing_EField = Timing_calc_EField + Timing_comm_EField;

  /* Store in file */
  if ((ringid == 0) && store_all_dir)
    {
      tstart = extime();

      strcpy(name,directory);

      if(which == 'X')
        strcat(name,"/EgridX");
      else
        strcat(name,"/EgridY");

      scatfile = fopen(name,"w");      
 
      dtheta = (180.0/PI)*(parms[THETA].max - parms[THETA].min)/((REAL)parms[THETA].Grid_size-1);
      dphi = (180.0/PI)*(parms[PHI].max - parms[PHI].min)/((REAL)parms[PHI].Grid_size-1);
      /* fprintf(scatfile,"theta\tphi\tEx.r\tEx.i\tEy.r\tEy.i\tEz.r\tEz.i\t\n\n"); */
      for (theta=0;theta<parms[THETA].Grid_size;++theta)
	for (phi=0;phi<parms[PHI].Grid_size;++phi)
	  {
	    index = grid_index(theta,phi);
	    fprintf(scatfile,"%.4f\t %.4f\t %e\t %e\t %e\t %e\t %e\t %e\n",
		    theta*dtheta,phi*dphi,
                    (double) Egrid[index+X].r,(double) Egrid[index+X].i,
                    (double) Egrid[index+Y].r,(double) Egrid[index+Y].i,
                    (double) Egrid[index+Z].r,(double) Egrid[index+Z].i);
	  }
      fclose(scatfile);
	  
      tstop = extime();
      Timing_FileIO += tstop - tstart;
    }
}



Rvector C_sca_integrand(int theta,int phi)
{
  int
    index,
    comp;

  index = grid_index(theta,phi);
  for (comp=0,n.x[0] = 0;comp<R3;++comp)
    {
      n.x[0] += 
	Egrid[index+comp].r*Egrid[index+comp].r + 
	Egrid[index+comp].i*Egrid[index+comp].i;
    }
  n.x[0] *= si[THETA][theta];

  return n;
}

void Sca_cross(/* Calculate the scattering cross section 
		* from the integral */
	       double k,        /* wave-number incoming wave */
	       double *res)     /* result */
{
  Rvector
    dummy;
  unsigned long
    tstart,tstop;
  char logfile[200];

  sprintf(logfile,"%s/log_Csca_int.dat",directory);
  n.x = (REAL *)malloc(sizeof(REAL));
  dummy.x = (REAL *) malloc(sizeof(REAL));

  tstart = extime();
  dummy = outer_qromb(parms,C_sca_integrand,1,logfile);
  tstop = extime();
  Timing_Integration += tstop - tstart;

  *res = dummy.x[0];

  free(n.x);
  free(dummy.x);
}







Rvector g_integrand(int theta,int phi)
{
  int
    index,
    comp;
  REAL E_square;
 
  index = grid_index(theta,phi);
  for (comp=0,E_square=0;comp<R3;++comp)
    {
      E_square +=
        Egrid[index+comp].r*Egrid[index+comp].r +
        Egrid[index+comp].i*Egrid[index+comp].i;
    }
  n.x[X] = E_square*si[THETA][theta]*si[THETA][theta]*co[PHI][phi];
  n.x[Y] = E_square*si[THETA][theta]*si[THETA][theta]*si[PHI][phi];
  n.x[Z] = E_square*si[THETA][theta]*co[THETA][theta];
 
  return n;
}
 
void Asym_parm(/* Calculate the unnormalized asymmetry parameter,
                * i.e. not yet normalized by Csca */
               double k,
               double vec[3])
{
  int
    comp;
  unsigned long
    tstart,tstop;
  Rvector
    dummy;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_int.dat",directory);
  n.x = (REAL *) malloc(3*sizeof(REAL));
  dummy.x = (REAL *) malloc(3*sizeof(REAL));
 
  tstart = extime();
  dummy = outer_qromb(parms,g_integrand,3,log_int);
  tstop = extime();
  Timing_Integration += tstop - tstart;
 
  for (comp=0;comp<3;++comp)
    vec[comp] = dummy.x[comp];
 
  free(n.x);
  free(dummy.x);
}




Rvector g_x_integrand(int theta,int phi)
{
  int
    index,
    comp;
  REAL E_square;
 
  index = grid_index(theta,phi);
  for (comp=0,E_square=0;comp<R3;++comp)
    {
      E_square +=
        Egrid[index+comp].r*Egrid[index+comp].r +
        Egrid[index+comp].i*Egrid[index+comp].i;
    }
  n.x[0] = E_square*si[THETA][theta]*si[THETA][theta]*co[PHI][phi];
 
  return n;
}

void Asym_parm_x(/* Calculate the unnormalized asymmetry parameter,
                  * i.e. not yet normalized by Csca */
                 double k,
                 double *vec)
{
  int
    comp;
  unsigned long
    tstart,tstop;
  Rvector
    dummy;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_x_int.dat",directory);
  n.x = (REAL *) malloc(sizeof(REAL));
  dummy.x = (REAL *) malloc(sizeof(REAL));
 
  tstart = extime();
  dummy = outer_qromb(parms,g_x_integrand,1,log_int);
  tstop = extime();
  Timing_Integration += tstop - tstart;
 
  *vec = dummy.x[0];
 
  free(n.x);
  free(dummy.x);
}

Rvector g_y_integrand(int theta,int phi)
{
  int
    index,
    comp;
  REAL E_square;
 
  index = grid_index(theta,phi);
  for (comp=0,E_square=0;comp<R3;++comp)
    {
      E_square +=
        Egrid[index+comp].r*Egrid[index+comp].r +
        Egrid[index+comp].i*Egrid[index+comp].i;
    }
  n.x[0] = E_square*si[THETA][theta]*si[THETA][theta]*si[PHI][phi];
 
  return n;
}

void Asym_parm_y(/* Calculate the unnormalized asymmetry parameter,
                  * i.e. not yet normalized by Csca */
                 double k,
                 double *vec)
{
  int
    comp;
  unsigned long
    tstart,tstop;
  Rvector
    dummy;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_y_int.dat",directory);
  n.x = (REAL *) malloc(sizeof(REAL));
  dummy.x = (REAL *) malloc(sizeof(REAL));
 
  tstart = extime();
  dummy = outer_qromb(parms,g_y_integrand,1,log_int);
  tstop = extime();
  Timing_Integration += tstop - tstart;
 
  *vec = dummy.x[0];
 
  free(n.x);
  free(dummy.x);
}

Rvector g_z_integrand(int theta,int phi)
{
  int
    index,
    comp;
  REAL E_square;
 
  index = grid_index(theta,phi);
  for (comp=0,E_square=0;comp<R3;++comp)
    {
      E_square +=
        Egrid[index+comp].r*Egrid[index+comp].r +
        Egrid[index+comp].i*Egrid[index+comp].i;
    }
  n.x[0] = E_square*si[THETA][theta]*co[THETA][theta];
 
  return n;
}

void Asym_parm_z(/* Calculate the unnormalized asymmetry parameter,
                  * i.e. not yet normalized by Csca */
                 double k,
                 double *vec)
{
  int
    comp;
  unsigned long
    tstart,tstop;
  Rvector
    dummy;
  char log_int[200];
 
  sprintf(log_int,"%s/log_asym_z_int.dat",directory);
  n.x = (REAL *) malloc(sizeof(REAL));
  dummy.x = (REAL *) malloc(sizeof(REAL));
 
  tstart = extime();
  dummy = outer_qromb(parms,g_z_integrand,1,log_int);
  tstop = extime();
  Timing_Integration += tstop - tstart;
 
  *vec = dummy.x[0];
 
  free(n.x);
  free(dummy.x);
}

dcomplex cmultiply(dcomplex a,dcomplex b)
{
  dcomplex c;

  c.r = a.r*b.r - a.i*b.i;
  c.i = a.r*b.i + a.i*b.r;

  return c;
}

dcomplex cadd(dcomplex a,dcomplex b)
{
  dcomplex c;

  c.r = a.r + b.r;
  c.i = a.i + b.i;

  return c;
}

void Frp_mat(/* Calculate the Radiation Pressure by direct calculation
	        of the scattering force. Per dipole the force of 
		the incoming photons, the scattering force and the 
		radiation pressure are calculated as intermediate results */
	     double Fsca_tot[3],REAL *Fsca,
	     double Finc_tot[3],REAL *Finc,
	     double Frp_tot[3],REAL *Frp,
	     dcomplex *x,REAL **rdip,double k)
{
  int 
    j,l,i,
    comp;
  int
    *materialT;
  REAL
    *rdipT;
  dcomplex
    *xT;
  char *Fname = __FILE__;
 
  for (comp=0;comp<3;++comp)
    {
      Fsca_tot[comp] = Finc_tot[comp] = Frp_tot[comp] = 0;
    }

  /* Convert internal fields to dipole moments;
     Calculate incoming force per dipole */
  for (j=0;j<local_Ndip;++j)
    if (material[j]<Nmat-1)
      {
	dcomplex dummy = {0,0};
	
	for (comp=0;comp<3;++comp)
	  {
	    int index = 3*j+comp;
	    dcomplex cc_d,_E_inc;
	    
	    /* Conversion */
	    cc_d.r = (REAL) cc[material[j]].r;
	    cc_d.i = (REAL) cc[material[j]].i;
	    x[index] = cmultiply(x[index],cc_d);

	    /* Im(P.E*inc) */
	    _E_inc.r = Einc[index].r;
	    _E_inc.i = -Einc[index].i;
	    dummy = cadd(dummy,cmultiply(x[index],_E_inc));
	  }

	Finc[3*j+Z] = k*dummy.i/2; 
	Finc_tot[Z] += Finc[3*j+Z];
      }

  /* Because of the parallelisation by row-block decomposition 
     the distributed arrays involved need to be gathered on each node 
     a) material -> material
     b) rdip -> rdipT
     c) x -> xT
     */
  if ((materialT = (int *) malloc(Ndip*sizeof(int)))==NULL)
    {
      LogError (EC_ERROR, Fname,
                "processor %d, ringID %d, could not malloc materialT",
                MyProcId, RingId);
      AbortServer (1);
    }
  if ((rdipT = (REAL *) malloc(3*Ndip*sizeof(REAL)))==NULL)
    {
      LogError (EC_ERROR, Fname,
                "processor %d, ringID %d, could not malloc rdipT",
                MyProcId, RingId);
      AbortServer (1);
    }
  if ((xT = (dcomplex *) malloc(3*Ndip*sizeof(dcomplex)))==NULL)
    {
      LogError (EC_ERROR, Fname,
                "processor %d, ringID %d, could not malloc xT",
                MyProcId, RingId);
      AbortServer (1);
    }

  memcpy(&materialT[local_d0],material,local_Ndip*sizeof(int));
  memcpy(&xT[3*local_d0],x,3*local_Ndip*sizeof(dcomplex));
  for (j=0;j<local_Ndip;++j)
    for (comp=0;comp<3;++comp)
      rdipT[3*(j+local_d0)+comp] = rdip[j][comp];

  all_gather(&materialT[local_d0],materialT,"int",local_Ndip);
  all_gather(&xT[3*local_d0],xT,"dcomplex",3*local_Ndip);
  all_gather(&rdipT[3*local_d0],rdipT,"REAL",3*local_Ndip);

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
	      REAL
		r,r2 = 0;      /* (squared) absolute distance */
	      dcomplex
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
	      ab1 = cmultiply(a,ab1);      ab2 = cmultiply(a,ab2);

	      /* Prepare c1 and c2 */
	      for (comp=0;comp<3;++comp)
		{
		  x_cg[comp].r = xT[jjj+comp].r;
		  x_cg[comp].i = -xT[jjj+comp].i;

		  Pn_j = cadd(Pn_j,cmultiply(x_cg[comp],n[comp]));
		  Pn_l = cadd(Pn_l,cmultiply(n[comp],   xT[lll+comp]));
		  inp  = cadd(inp, cmultiply(x_cg[comp],xT[lll+comp]));
		}

	      for (comp=0;comp<3;++comp)
		{
		  /* Set c1 */
		  c1[comp] = cmultiply(n[comp],cmultiply(Pn_j,Pn_l));
		  c1[comp].r *= -5;
		  c1[comp].i *= -5;

		  c1[comp] = cadd(c1[comp],cmultiply(inp,       n[comp]));
		  c1[comp] = cadd(c1[comp],cmultiply(Pn_j,      xT[lll+comp]));
		  c1[comp] = cadd(c1[comp],cmultiply(x_cg[comp],Pn_l));

		  /* Set c2 */
		  c2[comp] = cmultiply(n[comp],cmultiply(Pn_j,Pn_l));
		  c2[comp].r *= -1;
		  c2[comp].i *= -1;		  

		  c2[comp] = cadd(c2[comp],cmultiply(inp,n[comp]));

		  /* Fsca_{jl} = ... */
		  c1[comp] = cmultiply(ab1,c1[comp]);
		  c2[comp] = cmultiply(ab2,c2[comp]);
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
