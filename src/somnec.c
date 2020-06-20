/* last change:  pgm   8 nov 2000    1:04 pm */
/* program somnec(input,output,tape21) */

/* Original code is NEC2C (in public domain) obtained from http://www.si-list.net/swindex.html,
 * http://sharon.esrac.ele.tue.nl/users/pe1rxq/nec2c.rxq/nec2c.rxq-0.2.tar.gz
 * This code is incorporated into ADDA with several changes, which are explicitly indicated by comments. Those include:
 * - no generation of interpolation grid, only single run
 * - numerical precision was changed to double
 * - conjugation (that was in place to couple with other parts of nec2 code) was removed
 * - a few cosmetic changes to remove compiler warnings (with -pedantic flag)
 */

/* TODO: Systematic accuracy study of this code is required. At least 7 digits of precision are desired (for test runs)
 * However, the urgent thing is that for some reason the accuracy of this code is very bad (few percent errors for
 * rho=0, exactly). (issue 176)
 */

// ADDA: the cycling to generate interpolation grid was removed, now it is a single-run routine
/* program to generate nec interpolation grids for fields due to */
/* ground.  field components are computed by numerical evaluation */
/* of modified sommerfeld integrals. */
// ADDA: double precision is now used throughout
/* somnec2d is a double precision version of somnec for use with */
/* nec2d.  an alternate version (somnec2sd) is also provided in which */
/* computation is in single precision but the output file is written */
/* in double precision for use with nec2d.  somnec2sd runs about twic */
/* as fast as the full double precision somnec2d.  the difference */
/* between nec2d results using a for021 file from this code rather */
/* than from somnec2sd was insignficant in the cases tested. */

/* changes made by j bergervoet, 31-5-95: */
/* parameter 0. --> 0.d0 in calling of routine test */
/* status of output files set to 'unknown' */

// ADDA: The following definitions are extracted from the file nec2c.h (only those that are used here)
#include <complex.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#ifndef	TRUE
#define	TRUE	1
#endif
#ifndef	FALSE
#define	FALSE	0
#endif

/* common constants */
#define	PI	3.141592654
#define	TP	6.283185308
#define	PTP	.6283185308
#define PI10	31.41592654
#define	GAMMA	.5772156649
#define C1	-.02457850915
#define C2	.3674669052
#define C3	.7978845608
#define P10	.0703125
#define P20	.1121520996
#define Q10	.125
#define Q20	.0732421875
#define P11	.1171875
#define P21	.1441955566
#define Q11	.375
#define Q21	.1025390625
#define POF	.7853981635
/* MAXH=100 seems sufficient to obtain convergence in all cases; larger values doesn't seem to improve the accuracy
 * even when smaller CRIT is used. Probably other approximations are employed somewhere.
 */
#define MAXH	100
#define CRIT	1.0E-4
#define NM	131072
#define NTS	4

#define cmplx(r, i) ((r)+(i)*I)

void 	som_init(complex double epscf);
static void 	bessel(complex double z, complex double *j0, complex double *j0p);
void 	evlua(double zphIn,double rhoIn,complex double *erv, complex double *ezv,
	complex double *erh, complex double *eph);
static void 	gshank(complex double start, complex double dela, complex double *sum,
	int nans, complex double *seed, int ibk, complex double bk, complex double delb);
static void 	hankel(complex double z, complex double *h0, complex double *h0p);
static void 	lambda(double t, complex double *xlam, complex double *dxlam);
static void 	rom1(int n, complex double *sum, int nx);
static void 	saoa( double t, complex double *ans);

static void 	test(double f1r, double f2r, double *tr, double f1i,
	double f2i, double *ti, double dmin);
static void 	abort_on_error(int why);

/* common /evlcom/ */
static int jh;
static double ck2, ck2sq, tkmag, tsmag, ck1r, zph, rho;
static complex double ct1, ct2, ct3, ck1, ck1sq, cksm;

/* common /cntour/ */
static complex double a, b;

/*-----------------------------------------------------------------------*/

void som_init(complex double epscf)
{
  complex double erv, ezv;

  ck2=TP;
  ck2sq=ck2*ck2;

  /* sommerfeld integral evaluation uses exp(-jwt) */
  // ADDA: removed conjugation of epscf and in evaluation of integrals below

  ck1sq=ck2sq*epscf;
  ck1=csqrt(ck1sq);
  ck1r=creal(ck1);
  tkmag=100.*cabs(ck1);
  tsmag=100.*ck1*conj(ck1);
  cksm=ck2sq/(ck1sq+ck2sq);
  ct1=.5*(ck1sq-ck2sq);
  erv=ck1sq*ck1sq;
  ezv=ck2sq*ck2sq;
  ct2=.125*(erv-ezv);
  erv *= ck1sq;
  ezv *= ck2sq;
  ct3=.0625*(erv-ezv);

  return;
}

/*-----------------------------------------------------------------------*/

/* bessel evaluates the zero-order bessel function */
/* and its derivative for complex argument z. */
static void bessel(complex double z, complex double *j0, complex double *j0p )
{
  int k, i, ib, iz, miz;
  static int m[101], init = FALSE;
  static double a1[25], a2[25];
  double tst, zms;
  complex double p0z, p1z, q0z, q1z, zi, zi2, zk, cz, sz, j0x, j0px;

  /* initialization of constants */
  if( ! init )
  {
    for( k = 1; k <= 25; k++ )
    {
      i = k-1;
      a1[i]=-.25/(k*k);
      a2[i]=1.0/(k+1.0);
    }

    for( i = 1; i <= 101; i++ )
    {
      tst=1.0;
      for( k = 0; k < 24; k++ )
      {
	init = k;
	tst *= -i*a1[k];
	if( tst < 1.0e-6 )
	  break;
      }

      m[i-1] = init+1;
    } /* for( i = 1; i<= 101; i++ ) */

    init = TRUE;
  } /* if(init == 0) */

  zms=z*conj(z);
  if(zms <= 1.e-12)
  {
    *j0=1;
    *j0p=-.5*z;
    return;
  }

  ib=0;
  if(zms <= 37.21)
  {
    if(zms > 36.)
      ib=1;

    /* series expansion */
    iz=zms;
    miz=m[iz];
    *j0=1;
    *j0p=*j0;
    zk=*j0;
    zi=z*z;

    for( k = 0; k < miz; k++ )
    {
      zk *= a1[k]*zi;
      *j0 += zk;
      *j0p += a2[k]*zk;
    }
    *j0p *= -.5*z;

    if(ib == 0)
      return;

    j0x=*j0;
    j0px=*j0p;
  }

  /* asymptotic expansion */
  zi=1./z;
  zi2=zi*zi;
  p0z=1.+(P20*zi2-P10)*zi2;
  p1z=1.+(P11-P21*zi2)*zi2;
  q0z=(Q20*zi2-Q10)*zi;
  q1z=(Q11-Q21*zi2)*zi;
  zk=cexp(I*(z-POF));
  zi2=1./zk;
  cz=.5*(zk+zi2);
  sz=I*.5*(zi2-zk);
  zk=C3*csqrt(zi);
  *j0=zk*(p0z*cz-q0z*sz);
  *j0p=-zk*(p1z*sz+q1z*cz);

  if(ib == 0)
    return;

  zms=cos((sqrt(zms)-6.)*PI10);
  *j0=.5*(j0x*(1.+zms)+ *j0*(1.-zms));
  *j0p=.5*(j0px*(1.+zms)+ *j0p*(1.-zms));

  return;
}

/*-----------------------------------------------------------------------*/

/* evlua controls the integration contour in the complex */
/* lambda plane for evaluation of the sommerfeld integrals */
void evlua(double zphIn,double rhoIn, complex double *erv, complex double *ezv,
    complex double *erh, complex double *eph )
{
  int i, jump;
  static double del, slope, rmis;
  static complex double cp1, cp2, cp3, bk, delta, delta2, sum[6], ans[6];

  zph=zphIn;
  rho=rhoIn;

  del=zph;
  if( rho > del )
    del=rho;

  if(zph >= 2.*rho)
  {
    /* bessel function form of sommerfeld integrals */
    jh=0;
    a=0;
    del=1./del;

    if( del > tkmag)
    {
      b=cmplx(.1*tkmag,-.1*tkmag);
      rom1(6,sum,2);
      a=b;
      b=cmplx(del,-del);
      rom1 (6,ans,2);
      for( i = 0; i < 6; i++ )
	sum[i] += ans[i];
    }
    else
    {
      b=cmplx(del,-del);
      rom1(6,sum,2);
    }

    delta=PTP*del;
    gshank(b,delta,ans,6,sum,0,b,b);
    ans[5] *= ck1;

    /* ADDA: conjugate was removed */
    *erv=ck1sq*ans[2];
    *ezv=ck1sq*(ans[1]+ck2sq*ans[4]);
    *erh=ck2sq*(ans[0]+ans[5]);
    *eph=-ck2sq*(ans[3]+ans[5]);

    return;

  } /* if(zph >= 2.*rho) */

  /* hankel function form of sommerfeld integrals */
  jh=1;
  cp1=cmplx(0.0,.4*ck2);
  cp2=cmplx(.6*ck2,-.2*ck2);
  cp3=cmplx(1.02*ck2,-.2*ck2);
  a=cp1;
  b=cp2;
  rom1(6,sum,2);
  a=cp2;
  b=cp3;
  rom1(6,ans,2);

  for( i = 0; i < 6; i++ )
    sum[i]=-(sum[i]+ans[i]);

  /* path from imaginary axis to -infinity */
  if(zph > .001*rho)
    slope=rho/zph;
  else
    slope=1000.;

  del=PTP/del;
  delta=cmplx(-1.0,slope)*del/sqrt(1.+slope*slope);
  delta2=-conj(delta);
  gshank(cp1,delta,ans,6,sum,0,bk,bk);
  rmis=rho*(creal(ck1)-ck2);

  jump = FALSE;
  if( (rmis >= 2.*ck2) && (rho >= 1.e-10) )
  {
    if(zph >= 1.e-10)
    {
      bk=cmplx(-zph,rho)*(ck1-cp3);
      rmis=-creal(bk)/fabs(cimag(bk));
      if(rmis > 4.*rho/zph)
	jump = TRUE;
    }

    if( ! jump )
    {
      /* integrate up between branch cuts, then to + infinity */
      cp1=ck1-cmplx(.1,.2);
      cp2=cp1+.2;
      bk=cmplx(0.,del);
      gshank(cp1,bk,sum,6,ans,0,bk,bk);
      a=cp1;
      b=cp2;
      rom1(6,ans,1);
      for( i = 0; i < 6; i++ )
	ans[i] -= sum[i];

      gshank(cp3,bk,sum,6,ans,0,bk,bk);
      gshank(cp2,delta2,ans,6,sum,0,bk,bk);
    }

    jump = TRUE;

  } /* if( (rmis >= 2.*ck2) || (rho >= 1.e-10) ) */
  else
    jump = FALSE;

  if( ! jump )
  {
    /* integrate below branch points, then to + infinity */
    for( i = 0; i < 6; i++ )
      sum[i]=-ans[i];

    rmis=creal(ck1)*1.01;
    if( (ck2+1.) > rmis )
      rmis=ck2+1.;

    bk=cmplx(rmis,.99*cimag(ck1));
    delta=bk-cp3;
    delta *= del/cabs(delta);
    gshank(cp3,delta,ans,6,sum,1,bk,delta2);

  } /* if( ! jump ) */

  ans[5] *= ck1;

  /* ADDA: conjugate was removed */
  *erv=ck1sq*ans[2];
  *ezv=ck1sq*(ans[1]+ck2sq*ans[4]);
  *erh=ck2sq*(ans[0]+ans[5]);
  *eph=-ck2sq*(ans[3]+ans[5]);

  return;
}

/*-----------------------------------------------------------------------*/

/* gshank integrates the 6 sommerfeld integrals from start to */
/* infinity (until convergence) in lambda.  at the break point, bk, */
/* the step increment may be changed from dela to delb.  shank's */
/* algorithm to accelerate convergence of a slowly converging series */
/* is used */
static void gshank( complex double start, complex double dela, complex double *sum,
    int nans, complex double *seed, int ibk, complex double bk, complex double delb )
{
  int ibx, j, i, jm, intx, inx, brk=0;
  static double rbk, amg, den, denm;
  static complex double a1, a2, as1, as2, del, aa;
  static complex double q1[6][MAXH], q2[6][MAXH], ans1[6], ans2[6];

  rbk=creal(bk);
  del=dela;
  if(ibk == 0)
    ibx=1;
  else
    ibx=0;

  for( i = 0; i < nans; i++ )
    ans2[i]=seed[i];

  b=start;
  for( intx = 1; intx <= MAXH; intx++ )
  {
    inx=intx-1;
    a=b;
    b += del;

    if( (ibx == 0) && (creal(b) >= rbk) )
    {
      /* hit break point.  reset seed and start over. */
      ibx=1;
      b=bk;
      del=delb;
      rom1(nans,sum,2);
      if( ibx != 2 )
      {
	for( i = 0; i < nans; i++ )
	  ans2[i] += sum[i];
	intx = 0;
	continue;
      }

      for( i = 0; i < nans; i++ )
	ans2[i]=ans1[i]+sum[i];
      intx = 0;
      continue;

    } /* if( (ibx == 0) && (creal(b) >= rbk) ) */

    rom1(nans,sum,2);
    for( i = 0; i < nans; i++ )
      ans1[i] = ans2[i]+sum[i];
    a=b;
    b += del;

    if( (ibx == 0) && (creal(b) >= rbk) )
    {
      /* hit break point.  reset seed and start over. */
      ibx=2;
      b=bk;
      del=delb;
      rom1(nans,sum,2);
      if( ibx != 2 )
      {
	for( i = 0; i < nans; i++ )
	  ans2[i] += sum[i];
	intx = 0;
	continue;
      }

      for( i = 0; i < nans; i++ )
	ans2[i] = ans1[i]+sum[i];
      intx = 0;
      continue;

    } /* if( (ibx == 0) && (creal(b) >= rbk) ) */

    rom1(nans,sum,2);
    for( i = 0; i < nans; i++ )
      ans2[i]=ans1[i]+sum[i];

    den=0.;
    for( i = 0; i < nans; i++ )
    {
      as1=ans1[i];
      as2=ans2[i];

      if(intx >= 2)
      {
	for( j = 1; j < intx; j++ )
	{
	  jm=j-1;
	  aa=q2[i][jm];
	  a1=q1[i][jm]+as1-2.*aa;

	  if( (creal(a1) != 0.) || (cimag(a1) != 0.) )
	  {
	    a2=aa-q1[i][jm];
	    a1=q1[i][jm]-a2*a2/a1;
	  }
	  else
	    a1=q1[i][jm];

	  a2=aa+as2-2.*as1;
	  if( (creal(a2) != 0.) || (cimag(a2) != 0.) )
	    a2=aa-(as1-aa)*(as1-aa)/a2;
	  else
	    a2=aa;

	  q1[i][jm]=as1;
	  q2[i][jm]=as2;
	  as1=a1;
	  as2=a2;

	} /* for( j = 1; i < intx; i++ ) */

      } /* if(intx >= 2) */

      q1[i][intx-1]=as1;
      q2[i][intx-1]=as2;
      amg=fabs(creal(as2))+fabs(cimag(as2));
      if(amg > den)
	den=amg;

    } /* for( i = 0; i < nans; i++ ) */

    denm=1.e-3*den*CRIT;
    jm=intx-3;
    if(jm < 1)
      jm=1;

    for( j = jm-1; j < intx; j++ )
    {
      brk = FALSE;
      for( i = 0; i < nans; i++ )
      {
	a1=q2[i][j];
	den=(fabs(creal(a1))+fabs(cimag(a1)))*CRIT;
	if(den < denm)
	  den=denm;
	a1=q1[i][j]-a1;
	amg=fabs(creal(a1)+fabs(cimag(a1)));
	if(amg > den)
	{
	  brk = TRUE;
	  break;
	}

      } /* for( i = 0; i < nans; i++ ) */

      if( brk ) break;

    } /* for( j = jm-1; j < intx; j++ ) */

    if( ! brk )
    {
      for( i = 0; i < nans; i++ )
	sum[i]=.5*(q1[i][inx]+q2[i][inx]);
      return;
    }

  } /* for( intx = 1; intx <= maxh; intx++ ) */

  /* No convergence */
  printf("z=%g, rho=%g\n",zph,rho);
  abort_on_error(-6);
}

/*-----------------------------------------------------------------------*/

/* hankel evaluates hankel function of the first kind,   */
/* order zero, and its derivative for complex argument z */
static void hankel( complex double z, complex double *h0, complex double *h0p )
{
  int i, k, ib, iz, miz;
  static int m[101], init = FALSE;
  static double a1[25], a2[25], a3[25], a4[25], psi, tst, zms;
  complex double clogz, j0, j0p, p0z, p1z, q0z, q1z, y0, y0p, zi, zi2, zk;

  /* initialization of constants */
  if( ! init )
  {
    psi=-GAMMA;
    for( k = 1; k <= 25; k++ )
    {
      i = k-1;
      a1[i]=-.25/(k*k);
      a2[i]=1.0/(k+1.0);
      psi += 1.0/k;
      a3[i]=psi+psi;
      a4[i]=(psi+psi+1.0/(k+1.0))/(k+1.0);
    }

    for( i = 1; i <= 101; i++ )
    {
      tst=1.0;
      for( k = 0; k < 24; k++ )
      {
	init = k;
	tst *= -i*a1[k];
	if(tst*a3[k] < 1.e-6)
	  break;
      }
      m[i-1]=init+1;
    }

    init = TRUE;

  } /* if( ! init ) */

  zms=z*conj(z);
  if(zms == 0.)
    abort_on_error(-7);

  ib=0;
  if(zms <= 16.81)
  {
    if(zms > 16.)
      ib=1;

    /* series expansion */
    iz=zms;
    miz=m[iz];
    j0=1;
    j0p=j0;
    y0=0;
    y0p=y0;
    zk=j0;
    zi=z*z;

    for( k = 0; k < miz; k++ )
    {
      zk *= a1[k]*zi;
      j0 += zk;
      j0p += a2[k]*zk;
      y0 += a3[k]*zk;
      y0p += a4[k]*zk;
    }

    j0p *= -.5*z;
    clogz=clog(.5*z);
    y0=(2.*j0*clogz-y0)/PI+C2;
    y0p=(2./z+2.*j0p*clogz+.5*y0p*z)/PI+C1*z;
    *h0=j0+I*y0;
    *h0p=j0p+I*y0p;

    if(ib == 0)
      return;

    y0=*h0;
    y0p=*h0p;

  } /* if(zms <= 16.81) */

  /* asymptotic expansion */
  zi=1./z;
  zi2=zi*zi;
  p0z=1.+(P20*zi2-P10)*zi2;
  p1z=1.+(P11-P21*zi2)*zi2;
  q0z=(Q20*zi2-Q10)*zi;
  q1z=(Q11-Q21*zi2)*zi;
  zk=cexp(I*(z-POF))*csqrt(zi)*C3;
  *h0=zk*(p0z+I*q0z);
  *h0p=I*zk*(p1z+I*q1z);

  if(ib == 0)
    return;

  zms=cos((sqrt(zms)-4.)*31.41592654);
  *h0=.5*(y0*(1.+zms)+ *h0*(1.-zms));
  *h0p=.5*(y0p*(1.+zms)+ *h0p*(1.-zms));

  return;
}

/*-----------------------------------------------------------------------*/

/* compute integration parameter xlam=lambda from parameter t. */
static void lambda( double t, complex double *xlam, complex double *dxlam )
{
  *dxlam=b-a;
  *xlam=a+*dxlam*t;
  return;
}

/*-----------------------------------------------------------------------*/

/* rom1 integrates the 6 sommerfeld integrals from a to b in lambda. */
/* the method of variable interval width romberg integration is used. */
static void rom1( int n, complex double *sum, int nx )
{
  int jump, lstep, nogo, i, ns, nt;
  static double z, ze, s, ep, zend, dz=0., dzot=0., tr, ti;
  static complex double t00, t11, t02;
  static complex double g1[6], g2[6], g3[6], g4[6], g5[6], t01[6], t10[6], t20[6];

  lstep=0;
  z=0.;
  ze=1.;
  s=1.;
  ep=s/(1.e4*NM);
  zend=ze-ep;
  for( i = 0; i < n; i++ )
    sum[i]=0;
  ns=nx;
  nt=0;
  saoa(z,g1);

  jump = FALSE;
  while( TRUE )
  {
    if( ! jump )
    {
      dz=s/ns;
      if( (z+dz) > ze )
      {
	dz=ze-z;
	if( dz <= ep )
	  return;
      }

      dzot=dz*.5;
      saoa(z+dzot,g3);
      saoa(z+dz,g5);

    } /* if( ! jump ) */

    nogo=FALSE;
    for( i = 0; i < n; i++ )
    {
      t00=(g1[i]+g5[i])*dzot;
      t01[i]=(t00+dz*g3[i])*.5;
      t10[i]=(4.*t01[i]-t00)/3.;

      /* test convergence of 3 point romberg result */
      test( creal(t01[i]), creal(t10[i]), &tr, cimag(t01[i]), cimag(t10[i]), &ti, 0. );
      if( (tr > CRIT) || (ti > CRIT) )
	nogo = TRUE;
    }

    if( ! nogo )
    {
      for( i = 0; i < n; i++ )
	sum[i] += t10[i];

      nt += 2;
      z += dz;
      if(z > zend)
	return;

      for( i = 0; i < n; i++ )
	g1[i]=g5[i];

      if( (nt >= NTS) && (ns > nx) )
      {
	ns=ns/2;
	nt=1;
      }

      jump = FALSE;
      continue;

    } /* if( ! nogo ) */

    saoa(z+dz*.25,g2);
    saoa(z+dz*.75,g4);
    nogo=FALSE;
    for( i = 0; i < n; i++ )
    {
      t02=(t01[i]+dzot*(g2[i]+g4[i]))*.5;
      t11=(4.*t02-t01[i])/3.;
      t20[i]=(16.*t11-t10[i])/15.;

      /* test convergence of 5 point romberg result */
      test( creal(t11), creal(t20[i]), &tr, cimag(t11), cimag(t20[i]), &ti, 0. );
      if( (tr > CRIT) || (ti > CRIT) )
	nogo = TRUE;
    }

    if( ! nogo )
    {
      for( i = 0; i < n; i++ )
	sum[i] += t20[i];

      nt++;
      z += dz;
      if(z > zend)
	return;

      for( i = 0; i < n; i++ )
	g1[i]=g5[i];

      if( (nt >= NTS) && (ns > nx) )
      {
	ns=ns/2;
	nt=1;
      }

      jump = FALSE;
      continue;

    } /* if( ! nogo ) */

    nt=0;
    if(ns < NM)
    {
      ns *= 2;
      dz=s/ns;
      dzot=dz*.5;

      for( i = 0; i < n; i++ )
      {
	g5[i]=g3[i];
	g3[i]=g2[i];
      }

      jump = TRUE;
      continue;

    } /* if(ns < nm) */

    if( ! lstep )
    {
      lstep = TRUE;
      lambda( z, &t00, &t11 );
    }

    for( i = 0; i < n; i++ )
      sum[i] += t20[i];

    nt++;
    z += dz;
    if(z > zend)
      return;

    for( i = 0; i < n; i++ )
      g1[i]=g5[i];

    if( (nt >= NTS) && (ns > nx) )
    {
      ns /= 2;
      nt=1;
    }

    jump = FALSE;

  } /* while( TRUE ) */

}

/*-----------------------------------------------------------------------*/

/* saoa computes the integrand for each of the 6 sommerfeld */
/* integrals for source and observer above ground */
static void saoa( double t, complex double *ans)
{
  double xlr, sign;
  static complex double xl, dxl, cgam1, cgam2, b0, b0p, com, dgam, den1, den2;

  lambda(t, &xl, &dxl);
  if( jh == 0 )
  {
    /* bessel function form */
    bessel(xl*rho, &b0, &b0p);
    b0  *=2.;
    b0p *=2.;
    cgam1=csqrt(xl*xl-ck1sq);
    cgam2=csqrt(xl*xl-ck2sq);
    if(creal(cgam1) == 0.)
      cgam1=cmplx(0.,-fabs(cimag(cgam1)));
    if(creal(cgam2) == 0.)
      cgam2=cmplx(0.,-fabs(cimag(cgam2)));
  }
  else
  {
    /* hankel function form */
    hankel(xl*rho, &b0, &b0p);
    com=xl-ck1;
    cgam1=csqrt(xl+ck1)*csqrt(com);
    if(creal(com) < 0. && cimag(com) >= 0.)
      cgam1=-cgam1;
    com=xl-ck2;
    cgam2=csqrt(xl+ck2)*csqrt(com);
    if(creal(com) < 0. && cimag(com) >= 0.)
      cgam2=-cgam2;
  }

  xlr=xl*conj(xl);
  if(xlr >= tsmag)
  {
    if(cimag(xl) >= 0.)
    {
      xlr=creal(xl);
      if(xlr >= ck2)
      {
	if(xlr <= ck1r)
	  dgam=cgam2-cgam1;
	else
	{
	  sign=1.;
	  dgam=1./(xl*xl);
	  dgam=sign*((ct3*dgam+ct2)*dgam+ct1)/xl;
	}
      }
      else
      {
	sign=-1.;
	dgam=1./(xl*xl);
	dgam=sign*((ct3*dgam+ct2)*dgam+ct1)/xl;
      } /* if(xlr >= ck2) */

    } /* if(cimag(xl) >= 0.) */
    else
    {
      sign=1.;
      dgam=1./(xl*xl);
      dgam=sign*((ct3*dgam+ct2)*dgam+ct1)/xl;
    }

  } /* if(xlr < tsmag) */
  else
    dgam=cgam2-cgam1;

  den2=cksm*dgam/(cgam2*(ck1sq*cgam2+ck2sq*cgam1));
  den1=1./(cgam1+cgam2)-cksm/cgam2;
  com=dxl*xl*cexp(-cgam2*zph);
  ans[5]=com*b0*den1/ck1;
  com *= den2;

  if(rho != 0.)
  {
    b0p=b0p/rho;
    ans[0]=-com*xl*(b0p+b0*xl);
    ans[3]=com*xl*b0p;
  }
  else
  {
    ans[0]=-com*xl*xl*.5;
    ans[3]=ans[0];
  }

  ans[1]=com*cgam2*cgam2*b0;
  ans[2]=-ans[3]*cgam2*rho;
  ans[4]=com*b0;

  return;
}

/*-----------------------------------------------------------------------*/

/* test for convergence in numerical integration */
// ADDA: copied from nec2c.c
static void test( double f1r, double f2r, double *tr,
    double f1i, double f2i, double *ti, double dmin )
{
  double den;

  den= fabs( f2r);
  *tr= fabs( f2i);

  if( den < *tr)
    den= *tr;
  if( den < dmin)
    den= dmin;

  if( den < 1.0e-37)
  {
    *tr=0.;
    *ti=0.;
    return;
  }

  *tr= fabs(( f1r- f2r)/ den);
  *ti= fabs(( f1i- f2i)/ den);

  return;
}

/*------------------------------------------------------------------------*/

/*  abort_on_error()
 *
 *  prints an error message and exits
 */

static void abort_on_error( int why )
{
  switch( why )
  {
    case -1 : /* abort if input file name too long */
      fprintf( stderr, "%s\n",
	  "nec2c: input file name too long - aborting" );
      break;

    case -2 : /* abort if output file name too long */
      fprintf( stderr, "%s\n",
	  "nec2c: output file name too long - aborting" );
      break;

    case -3 : /* abort on input file read error */
      fprintf( stderr, "%s\n",
	  "nec2c: error reading input file - aborting" );
      break;

    case -4 : /* Abort on malloc failure */
      fprintf( stderr, "%s\n",
	  "nec2c: A memory allocation request has failed - aborting" );
      break;

    case -5 : /* Abort if a GF card is read */
      fprintf( stderr, "%s\n",
	  "nec2c: NGF solution option not supported - aborting" );
      break;

    case -6: /* No convergence in gshank() */
            fprintf( stderr, "%s\n",
	  "nec2c: No convergence in gshank() - aborting" );
      break;

    case -7: /* Error in hankel() */
            fprintf( stderr, "%s\n",
	  "nec2c: hankel not valid for z=0. - aborting" );
      break;
  }  /* switch( why ) */

  /* clean up and quit */
  exit( why );

} /* end of abort_on_error() */
