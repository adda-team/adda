/* Original code (somnec.c) is part of NEC2C v.1.3 (in public domain) obtained
 * from https://www.qsl.net/5b4az/pkg/nec2/nec2c/nec2c-1.3.tar.bz2
 * This code is incorporated into ADDA with several changes, which are explicitly indicated by comments. Those include:
 * - no generation of interpolation grid, only single run
 * - conjugation (that was in place to couple with other parts of nec2 code) was removed
 * - a few cosmetic changes to remove compiler warnings (with -pedantic flag)
 * - replacement of 20 to MAXH in couple of places
 *
 * Further updates
 */

/* TODO: Systematic accuracy study of this code is required. At least 7 digits of precision are desired (for test runs)
 * However, the urgent thing is that for some reason the accuracy of this code is very bad (few percent errors for
 * rho=0, exactly). (issue 176)
 */

// ADDA: the cycling to generate interpolation grid was removed, now it is a single-run routine
/* program to generate nec interpolation grids for fields due to
 ground.  field components are computed by numerical evaluation
 of modified sommerfeld integrals.

 somnec2d is a double precision version of somnec for use with
 nec2d.  an alternate version (somnec2sd) is also provided in which
 computation is in single precision but the output file is written
 in double precision for use with nec2d.  somnec2sd runs about twice
 as fast as the full double precision somnec2d.  the difference
 between nec2d results using a for021 file from this code rather
 than from somnec2sd was insignificant in the cases tested.

 changes made by j bergervoet, 31-5-95:
 parameter 0. --> 0.d0 in calling of routine test
 status of output files set to 'unknown' */

#include "somnec.h" // corresponding header

// ADDA: The following definitions are extracted from the file nec2c.h (only those that are used here)
#include <math.h>
#include <stdbool.h> // for bool
#include <stdio.h>
#include <stdlib.h>

#ifndef	TRUE
#define	TRUE	1
#endif
#ifndef	FALSE
#define	FALSE	0
#endif

/* common constants */
#define PI    3.141592654
#define TP    6.283185308
#define PTP   0.6283185308
#define PI10  31.41592654
#define GAMMA 0.5772156649
#define C1    -0.02457850915
#define C2    0.3674669052
#define C3    0.7978845608
#define P10   0.0703125
#define P20   0.1121520996
#define Q10   0.125
#define Q20   0.0732421875
#define P11   0.1171875
#define P21   0.1441955566
#define Q11   0.375
#define Q21   0.1025390625
#define POF   0.7853981635
/* MAXH=100 seems sufficient to obtain convergence in all cases; larger values doesn't seem to improve the accuracy
 * even when smaller CRIT is used. Probably other approximations are employed somewhere.
 * However, larger MAXH is required for distances larger than several wavelengths to reach specified CRIT. In practice,
 * such CRIT is not needed for large distances, since Sommerfeld integrals are only a correction to the analytical
 * reflection term. Thus, we increased MAXH to 500 for now, but final solution is left for
 * https://github.com/adda-team/adda/issues/176
 * If you get convergence errors, increase MAXH further
 * Also, tsmag below may affect the accuracy. Currently, the relative errors due to tsmag is <= 3e-7
 * TODO: increase accuracy of constants above to double
 */
#define MAXH 500
#define CRIT 1.0E-4
#define NM   131072
#define NTS  4

#define EPS_TAN 0.1 // tangent of minimum angle between straight segment and direction to the pole

#define cmplx(r, i) ((r)+(i)*I)

static void bessel(complex double z, complex double *j0, complex double *j0p);
static void gshank(complex double start, complex double dela, complex double *sum,
	int nans, complex double *seed, int ibk, complex double bk, complex double delb);
static void hankel(complex double z, complex double *h0, complex double *h0p);
static void lambda(double t, complex double *xlam, complex double *dxlam);
static void rom1(int n, complex double *sum, int nx);
static void saoa( double t, complex double *ans);

static void test(double f1r, double f2r, double *tr, double f1i,
	double f2i, double *ti, double dmin);
static void abort_on_error(int why);

/* common /evlcom/ */
static int jh;
static double ck2, ck2sq, tkmag, tsmag, ck1r, zph, rho;
static complex double ct1, ct2, ct3, ck1, ck1sq, cksm;
static complex xl0; // location of surface-plasmon pole
static bool pole;   // whether the surface-plasmon pole is real (lies on the principal Riemann sheet)
static bool denser; // whether the substrate material is denser than the upper medium (vacuum)

/* common /cntour/ */
static complex double a, b;

/*-----------------------------------------------------------------------*/

// a single function from cmplx.h not to include the whole file
static inline double cAbs2(const complex double z)
// square of absolute value of complex number; |z|^2
{
	return creal(z)*creal(z) + cimag(z)*cimag(z);
}

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
	// TODO: the largest of |ck1|, |ck2| should be used to define tsmag
	tsmag=100.*cAbs2(ck1);
	cksm=ck2sq/(ck1sq+ck2sq);
	ct1=.5*(ck1sq-ck2sq);
	erv=ck1sq*ck1sq;
	ezv=ck2sq*ck2sq;
	ct2=.125*(erv-ezv);
	erv *= ck1sq;
	ezv *= ck2sq;
	ct3=.0625*(erv-ezv);
	// determine pole (of D2) position and whether it is real
	xl0=ck1*ck2/csqrt(ck1sq+ck2sq);
	pole=(creal(xl0)>ck2);
	denser=(creal(ck1)>ck2);

	return;
}

//======================================================================================================================

static void residual(complex double *ans)
// increments answer by residual of integrals, involving Hankel functions, over the pole xl0
{
	complex double cgam2,h0,h0p,com,den2;

	hankel(xl0*rho,&h0,&h0p); // assumes that rho!=0
	cgam2=I*xl0*ck2/ck1; // here we employ known value of xl0, and that Re(xl0)>Re(k2)

	den2=-TP*ck1*ck2/(ck1sq*ck1sq-ck2sq*ck2sq); // residual of D2/2
	com=den2*xl0*cexp(-cgam2*zph);
	// ans[5] is not affected

	ans[0]-=com*xl0*((h0p/rho)+h0*xl0); // h0pp = -(h0 + h0p/(rho*xl))
	ans[1]+=com*cgam2*cgam2*h0;
	ans[2]-=com*cgam2*xl0*h0p;
	ans[3]+=com*xl0*h0p/rho;
	ans[4]+=com*h0;

	return;
}

/*-----------------------------------------------------------------------*/

/* bessel evaluates the zero-order bessel function */
/* and its derivative for complex argument z. */
static void bessel( complex double z, complex double *j0, complex double *j0p )
{
	int k, ib;
	static int m[101], init = FALSE;
	static double a1[25], a2[25];
	double zms;
	complex double p0z, p1z, q0z, q1z, zi, zi2, zk, cz, sz, j0x=0, j0px=0;

	/* initialization of constants */
	if( ! init )
	{
		int i;
		for( k = 1; k <= 25; k++ )
		{
			i = k-1;
			a1[i]=-.25/(k*k);
			a2[i]=1.0/(k+1.0);
		}

		for( i = 1; i <= 101; i++ )
		{
			double tst = 1.0;
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

	zms=cAbs2(z);
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
		int iz=(int)zms;
		int miz=m[iz];
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
/* the logic for choosing integration contours follows Section 4.2 "Numerical Evaluation of the Sommerfeld Integrals" of
 * http://users.tpg.com.au/micksw012/files/gnec-theory.pdf (gNEC Theory v. 0.9.15, 02.02.2018) with the only change of
 * conjugating all complex variables due to a different time-harmonic convention. In the following "on the figure"
 * refers to contour descriptions in this section.
 * Improvements upon this logic are noted explicitly.
 */
void evlua(double zphIn,double rhoIn, complex double *erv, complex double *ezv,
	complex double *erh, complex double *eph, int mode)
/* mode determines the choice of contour: automatic (0), Bessel (1), shorter Hankel (2), longer Hankel (3)
 * mode=2 is effectively used is some other value is specified (TODO: add a meaningful test for that)
 */
{
	int i, jump;
	static double del, slope, rmis;
	static complex double cp1, cp2, cp3, bk, delta, delta2, sum[6], ans[6];

	zph=zphIn;
	rho=rhoIn;

	del=zph;
	if( rho > del )
		del=rho;

	if (mode==1 || (mode==0 && zph >= 2.*rho))
	{
		/* Bessel-function form of Sommerfeld integrals.
		 * Sufficiently large Z guarantees rapid exponential decay through the contour
		 * TODO: for large rho, it makes sense to move the inflection point closer to the real axis (due to exponential
		 *       increase of Bessel function with |Im(argument)|
		 */
		jh=0;
		a=0;
		del=1./del;

		if( del > tkmag)
		{
			/* when both rho and Z are small, the contour from zero to inflection point is divided into two parts,
			 * probably to better control the accuracy under the large variation of the integrand over the contour
			 * (at least, of the J0'). This is not explained in the docs, and it is not clear, why the k1 (not k2) is
			 * used for the boundary value tkmag. Still, any changes will require extensive testing to justify.
			 */
			b=cmplx(.1*tkmag,-.1*tkmag);
			rom1(6,sum,2); // from zero to intermediate point
			a=b;
			b=cmplx(del,-del);
			rom1 (6,ans,2); // from intermediate point to the inflection one on the figure
			for( i = 0; i < 6; i++ )
				sum[i] += ans[i];
		}
		else
		{
			b=cmplx(del,-del);
			rom1(6,sum,2); // from zero to inflection point on the figure
		}

		delta=PTP*del;
		gshank(b,delta,ans,6,sum,0,b,b); // horizontal line to infinity on the figure
		ans[5] *= ck1;

		/* ADDA: conjugate was removed */
		*erv=ck1sq*ans[2];
		*ezv=ck1sq*(ans[1]+ck2sq*ans[4]);
		*erh=ck2sq*(ans[0]+ans[5]);
		*eph=-ck2sq*(ans[3]+ans[5]);

		return;

	} /* if(zph >= 2.*rho) */

	/* Hankel-function form of Sommerfeld integrals, based on H0^(1)(rho*xl) */
	jh=1;
	// the following ensures that the contour never crosses the branch cut from ck1
	cp1=(creal(ck1)+cimag(ck1)>0.41*ck2) ? I*0.4*ck2 : 0.99*ck1;
	// bottom position of the contour (imaginary part); should be small for large rho to avoid loss of precision
	double imB=(ck2*rho < 20) ? -.2*ck2 : -4/rho;
	cp2=cmplx(.4*ck2-imB,imB);
	/* for poles close to ck1 (when denser is true), we shift the point c slightly to the right of the pole - this works
	 * fine for all other contours. Otherwise (if denser is false), we need to consider the pole more carefully below
	 */
	cp3=((pole && denser) ? creal(xl0) + 0.02*ck2 : 1.02*ck2)  + I*imB;
	a=cp1;
	b=cp2;
	rom1(6,sum,2); // from a to b on the figure
	a=cp2;
	b=cp3;
	rom1(6,ans,2); // from b to c on the figure

	/* minus signs are used somewhat complicatedly in the following, motivated by the fact that some integrals are
	 * computed from infinity to a point, i.e., in reverse direction to that assumed by gshank
	 */
	for( i = 0; i < 6; i++ )
		sum[i]=-(sum[i]+ans[i]);

	/* path from imaginary axis to -infinity */
	/* Overall, the paths to infinity are chosen so that the exp(-gamma2*Z)*H0(1)(rho*xl) asymptotically decays as
	 * exp(-a*|xl|) with real a>0 (for large |xl|). Below, the difference between delta and delta2 is due to their usage
	 * in the second and first quadrant, respectively, where gamma2 ~ -z and z.
	 */
	if(zph > .001*rho)
		slope=rho/zph;
	// this avoids very large numbers below, but convergence is guaranteed to be excellent along this direction
	else
		slope=1000.;

	del=PTP/del;
	delta=cmplx(-1.0,slope)*del/sqrt(1.+slope*slope);
	delta2=-conj(delta);
	//TODO: replace bk with 0 in last two arguments to gshank, whenever they are not used (i.e., when ibk=0)
	gshank(cp1,delta,ans,6,sum,0,bk,bk); // from infinity to a on the figure
	// by this point ans hold minus integral from infinity to c through a and b
	/* TODO: make the use of point e explicit here, and also test whether a straight line will be sufficient (as is done
	 * below for simpler contour). However, this test may be already implicit in the following slope comparisons.
	 */
	rmis=rho*creal(ck1-0.1-cp3);

	// comparing with ck2 have different units on two sides, TODO: replace with ck2->2pi
	if (mode==3 || (mode==0 && rmis >= 2.*ck2 && rho >= 1.e-10))
	{
		/* this case will lead to nonsense results for rmis<0, in particular when Re(k1)<k2
		 * (can only happen when overridden by mode)
		 */
		jump = TRUE;
		if(zph >= 1.e-10)
		{
			/* Here we choose between two alternatives: (1) integrate directly between two branch points (c to d on
			 * figure) or (2) from c to infinity and then back to d along two vertical lines near the corresponding
			 * branch cuts. The choice of one over another depends on the |arg(a)| in asymptotic expression exp(-a*|xl|)
			 * for the integrand. We want to minimize this argument.
			 * In these cases, |tg[arg(a)]| is approximately tg[arg((Z-i*rho)*(d-c))] and Z/rho, respectively. Inverse
			 * of these quantities are compared in the following.
			 * Factor 4 seems to be an empirical one to balance the two options (either one finite integral or two
			 * infinite ones). However, this specific values is not discussed in any docs, moreover, based on the
			 * balance arguments, it is more logical for it to be 1/4 or 1/2. The only possible argument in favor of 4
			 * is that a finite integral does not necessarily reaches the asymptotic behavior, and thus is more
			 * complicated for adaptive routines.
			 * Anyway, changing this factor can only be done after careful testing of computational times to reach the
			 * prescribed accuracy.
			 */
			bk=cmplx(-zph,rho)*(ck1-cp3);
			rmis=-creal(bk)/fabs(cimag(bk));
			// TODO: better to use inverse comparison, since Z can be 0
			if(rmis > 4.*rho/zph && mode!=3)
				jump = FALSE;
		}

		if( jump )
		{
			/* integrate up to infinity and back between branch cuts, then to + infinity on the right side */
			// TODO: the shifts from the poke should be proportional to either ck1 or ck2
			cp1=ck1-cmplx(.1,.2);
			cp2=cp1+.2;
			bk=cmplx(0.,del);
			gshank(cp1,bk,sum,6,ans,0,bk,bk); // from e to infinity on the figure
			a=cp1;
			b=cp2;
			rom1(6,ans,1);
			// invert sign from previous parts of contour
			for( i = 0; i < 6; i++ )
				ans[i] -= sum[i];

			gshank(cp3,bk,sum,6,ans,0,bk,bk); // from c to infinity on the figure
			gshank(cp2,delta2,ans,6,sum,0,bk,bk); // from f to infinity on the figure
		}

	} /* if( (rmis >= 2.*ck2) && (rho >= 1.e-10) ) */
	else
		jump = FALSE;

	if( ! jump ) // this is always enabled when mode==2 (or any other unrecognized value)
	{
		/* integrate below branch points, then to + infinity */
		// invert sign from previous parts of contour
		for( i = 0; i < 6; i++ )
			sum[i]=-ans[i];

		if (denser) { // here the SP pole is already accounted for by shift of point c
			complex double cp4,d34;
			/* TODO: the shift from the zero may be inadequate for very different real and imaginary parts
			 * but any other empirics must ensure that the slope of the path is always positive
			 */
			cp4=1.01*creal(ck1)+0.99*cimag(ck1)*I;
			d34=cp4-cp3;
			/* test if a direct line from cp3 will hit the branch cut from k1; if yes, integrate from from c to d and
			 * then to infinity, otherwise - directly from ñ to infinity
			 */
			if (cimag(d34)<slope*creal(d34)) gshank(cp3,del*d34/cabs(d34),ans,6,sum,1,cp4,delta2);
			else gshank(cp3,delta2,ans,6,sum,0,0,0);
		}
		else { // here we do not need to worry about the branch cut from k1, but need to account for the SP pole
			if (pole) { // this can be optimized by singularity extraction
				complex double rot=(1+I*EPS_TAN)/(1+EPS_TAN*EPS_TAN); // rotation multiplier for EPS_TAN
				complex double dP3=xl0-cp3;
				// tangent of angle from direction to the pole to the slope (rho/Z)
				double tanA=(creal(dP3)*slope-cimag(dP3))/(creal(dP3)+cimag(dP3)*slope);
				if (tanA>EPS_TAN) { // leaves pole to the right at sufficient distance
					gshank(cp3,delta2,ans,6,sum,0,0,0);
					residual(ans);
				}
				// leaves pole to the left at sufficient distance
				else if (tanA<-EPS_TAN) gshank(cp3,delta2,ans,6,sum,0,0,0);
				/* The following chooses one of the paths around the pole, shifted by angle EPS from the exact direction
				 * to the pole. Generally, we choose the closest one, but ensure that the slope of the first segment is
				 * strictly between 0 and infinity
				 */
				else if ( ( tanA>0 && cimag(dP3)*EPS_TAN<creal(dP3) ) || creal(dP3)*EPS_TAN>=cimag(dP3) ) {
					gshank(cp3,del*rot*dP3/cabs(dP3),ans,6,sum,1,cp3+dP3*rot,delta2);
					residual(ans);
				}
				else gshank(cp3,del*conj(rot)*dP3/cabs(dP3),ans,6,sum,1,cp3+dP3*conj(rot),delta2);
			}
			else gshank(cp3,delta2,ans,6,sum,0,0,0); // from ñ directly to infinity on the figure
		}
	} /* if( ! jump ) */

	ans[5] *= ck1;

	/* ADDA: conjugate was removed */
	// not clear what is the benefit of using 6 integrands instead of 4
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
static void gshank( complex double start, complex double dela,
	complex double *sum, int nans, complex double *seed,
	int ibk, complex double bk, complex double delb )
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
			for( i = 0; i < nans; i++ )
				ans2[i] += sum[i];
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
	int k, ib;
	static int m[101], init = FALSE;
	static double a1[25], a2[25], a3[25], a4[25], psi, tst, zms;
	complex double clogz, j0, j0p, p0z, p1z, q0z, q1z, y0=0, y0p=0, zi, zi2, zk;

	/* initialization of constants */
	if( ! init )
	{
		int i;
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

	zms=cAbs2(z);
	if(zms == 0.)
		abort_on_error(-7);

	ib=0;
	if(zms <= 16.81)
	{
		if(zms > 16.)
			ib=1;

		/* series expansion */
		int iz=(int)zms;
		int miz=m[iz];
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

// ADDA: the following was present in NEC2C v.0.2 but disappeared in v.1.3. Not clear, what is the implication
//		if( (nt >= NTS) && (ns > nx) )
//		{
//			ns /= 2;
//			nt=1;
//		}

		jump = FALSE;

	} /* while( TRUE ) */

}

/*-----------------------------------------------------------------------*/

/* saoa computes the integrand for each of the 6 sommerfeld */
/* integrals for source and observer above ground */
static void saoa( double t, complex double *ans)
{
	double xlr;
	static complex double xl, dxl, cgam1, cgam2, b0, b0p, com, dgam, den1, den2;

	lambda(t, &xl, &dxl);
	// evaluate gamma1 and gamma2 (cgam1, cgam2)
	if( jh == 0 )
	{
		/* bessel function form */
		/* Assuming k1,2 to belong to the first quadrant of the complex plane, the following expressions have branch
		 * cuts from k to i*inf in an arc that stays within the first quadrant with Im(z)>=Im(k) and symmetrically from
		 * -k to -i*inf.
		 * For Re(k)=0, this is exactly the vertical line from k to k+i*inf (as assumed in the theory).
		 * For Im(k)=0, it is a corner from k to 0 and further to i*inf.
		 * This is sufficient to avoid integration contour, used for the Bessel form (concentrated fully in the fourth
		 * quadrant). The exact zero falls on the branch  cut for Im(k)=0, but the function value is zero in this point
		 * anyway. The last part of the code keeps the functions continuous at the boundary of the fourth quadrant, when
		 * Im(k)=0, although that does not seem strictly necessary for the employed contours.
		 */
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
		/* Assuming k1,2 to belong to the first quadrant of the complex plane, the following expressions have branch
		 * cuts from k vertically to k+i*inf and from -k horizontally to -k-inf. The first one is the same as in the
		 * manuals of NEC code (where the integration contours are discussed), while the second one is different. But
		 * it still stays in the third quadrant, which is sufficient for the used integration contours.
		 */
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

	/* evaluate dgam=gamma2-gamma1 avoiding loss of precision for large |xl|
	 * It uses the expansion of gamma=sqrt(z^2-k^2) with error of O(k^8/z^7), which is valid (with sign=1) for |z|>>|k|
	 * and Re(z)>Re(k)sign(Im(z)), assuming k is in the first quadrant.
	 * If sign=-1, then the expansion works for |z|>>|k| and Re(z)<Re(k)sign(Im(z))
	 */
	xlr=cAbs2(xl);
	if(xlr >= tsmag)
	{
		double sign;
		if(cimag(xl) >= 0.)
		{
			xlr=creal(xl);
			if(xlr >= ck2)
			{
				/* the original expression is fine between the two branch cuts, since then gamma1 is close to -gamma2
				 * instead of gamma2 (for large |z|), i.e., they do not cancel each other
				 */
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
			/* this covers the whole second quadrant (which seems sufficient for all used contours),
			 * but the values for 0<Re(xl)<Re(k2) are only correct only if Re(k1)>=Re(k2). Otherwise,
			 * additionally Re(xl)<Re(k1) should be tested (leading to dgam=cgam2-cgam1, as above - TODO).
			 */
			{
				sign=-1.;
				dgam=1./(xl*xl);
				dgam=sign*((ct3*dgam+ct2)*dgam+ct1)/xl;
			} /* if(xlr >= ck2) */

		} /* if(cimag(xl) >= 0.) */
		else
		// this simple condition is sufficient for Im(xl)<0 if we assume that third quadrant is never reached
		{
			sign=1.;
			dgam=1./(xl*xl);
			dgam=sign*((ct3*dgam+ct2)*dgam+ct1)/xl;
		}

	} /* if(xlr < tsmag) */
	else
		dgam=cgam2-cgam1;

	// calculate integrand multipliers D1 and D2 (up to a factor of 2). D2 is rewritten to avoid loss of precision
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
		ans[0]=-com*xl*xl*.5*b0; // here b0=2 due to the logic for Bessel contour above
		ans[3]=ans[0];
	}

	ans[1]=com*cgam2*cgam2*b0;
	ans[2]=-ans[3]*cgam2*rho;
	ans[4]=com*b0;

	return;
}

/*-----------------------------------------------------------------------*/

/* test for convergence in numerical integration */
// ADDA: copied from calculation.c
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
// ADDA: copied from misc.c, first 5 cases are redundant
static void abort_on_error( int why )
{
	switch( why )
	{
		case -1 : /* abort if input file name too long */
			fprintf(stderr, "%s\n",
				"nec2c: input file name too long - aborting");
			break;

		case -2 : /* abort if output file name too long */
			fprintf(stderr, "%s\n",
				"nec2c: output file name too long - aborting");
			break;

		case -3 : /* abort on input file read error */
			fprintf(stderr, "%s\n",
				"nec2c: error reading input file - aborting");
			break;

		case -4 : /* Abort on malloc failure */
			fprintf(stderr, "%s\n",
				"nec2c: A memory allocation request has failed - aborting");
			break;

		case -5 : /* Abort if a GF card is read */
			fprintf(stderr, "%s\n",
				"nec2c: NGF solution option not supported - aborting");
			break;

		// ADDA: updated comment to be specific to ADDA
		case -6: /* No convergence in gshank() */
			fprintf(stderr,"No convergence in gshank() - aborting. Try to increase MAXH in somnec.c and recompile\n");
			break;

		case -7: /* Error in hankel() */
			fprintf(stderr, "%s\n",
				"nec2c: hankel not valid for z=0. - aborting");
			break;
	}  /* switch( why ) */

	/* clean up and quit */
	exit( why );

} /* end of abort_on_error() */
