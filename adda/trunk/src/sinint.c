/* FILE: sinint.c
 * AUTH: Maxim Yurkin
 * DESCR: Function for calculating sine and cosine integrals. Double versions of
 *        routine given in "Numerical Recipes in C" 2nd ed.
 *
 * Copyright (C) 2007 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <math.h>
#include <float.h>        /* for DBL_MIN and DBL_MAX */
#include "const.h"
#include "cmplx.h"
#include "io.h"

#define EPS DBL_EPSILON     /* Relative error, or absolute error near a zero of Ci(x) */
#define MAXIT 150           /* Maximum number of iterations allowed */
#define TMIN 2.0            /* Dividing line between using the series and the continued fraction */

/*============================================================*/

void cisi(const double x,double *ci,double *si)
  /* Computes the cosine and sine integrals Ci(x) and Si(x). Ci(0) is returned
     as a large negative number and no error message is generated. For x<0
     routine returns Ci(-x), while actually Si(-x)=-Si(x), Ci(-x)=Ci(x)-i*pi. */

{
  int i,k,odd;
  double a,err,fact,sign,sum,sumc,sums,t,term;
  doublecomplex h,b,c,d,del,tmp;

  t=fabs(x);
  /* special case */
  if (x==0) {
    *si=0;
    *ci=-DBL_MAX;
    return;
  }
  /* Evaluate continued fraction by modified Lentz's method */
  if (t>TMIN) {
    b[RE]=1;
    b[IM]=t;
    c[RE]=c[IM]=0;
    cInv(b,d);
    cEqual(d,h);
    for (i=2;i<=MAXIT;i++) {
      a=-(i-1)*(i-1);
      b[RE]+=2;
      /* d=1/(a*d+b) */
      cMultReal(a,d,d);
      cAdd(b,d,d);
      cInv(d,d);
      /* c=b+a/c */
      if (i==2) cEqual(b,c);  /* special case, then initially c=0 */
      else {
        cInv(c,c);
        cMultReal(a,c,c);
        cAdd(b,c,c);
      }
      cMult(c,d,del);
      cMultSelf(h,del);
      if (fabs(del[RE]-1)+fabs(del[IM])<EPS) break;
    }
    if (i > MAXIT) LogError(EC_ERROR,ALL_POS,
                            "Failed to converge during calculation of sine integral of %g",x);
    imExp(-t,tmp);
    cMultSelf(h,tmp);
    *ci=-h[RE];
    *si=PI_OVER_TWO+h[IM];
  }
  else {  /* Evaluate both series simultaneously */
    /* Special case: avoid failure of convergence test because of underflow */
    if (t<sqrt(DBL_MIN)) {
      sumc=0;
      sums=t;
    }
    else {
      sum=sums=sumc=0;
      sign=fact=1;
      odd=TRUE;
      for (k=1;k<=MAXIT;k++) {
        fact*=t/k;
        term=fact/k;
        sum+=sign*term;
        err=term/fabs(sum);
        if (odd) {
          sign=-sign;
          sums=sum;
          sum=sumc;
        }
        else {
          sumc=sum;
          sum=sums;
        }
        if (err<EPS) break;
        odd=!odd;
      }
      if (k > MAXIT) LogError(EC_ERROR,ALL_POS,
                            "Failed to converge during calculation of sine integral of %g",x);
    }
    *si=sums;
    *ci=sumc+log(t)+EULER;
  }
  if (x<0) *si=-(*si);
}
