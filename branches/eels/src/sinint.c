/* File: sinint.c
 * $Date::                            $
 * Descr: function for calculating sine and cosine integrals
 *
 *        It was originaly based on routine given in "Numerical Recipes in C" 2nd ed. and then was slightly corrected
 *        according to the 3rd ed. of the same book.
 *
 * Copyright (C) 2007-2010,2012-2013 ADDA contributors
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
#include "io.h"
// system headers
#include <float.h> // for DBL_MIN and DBL_MAX
#include <math.h>

#define EPS DBL_EPSILON // Relative error, or absolute error near a zero of Ci(x)
#define MAXIT 100       // Maximum number of iterations allowed
#define TMIN 2.0        // Dividing line between using the series and the continued fraction
#define BIG DBL_MAX*EPS // A number near machine overflow limit
#define FPMIN DBL_MIN*4 // Number close to the smallest representable number

//======================================================================================================================

void cisi(const double x,double *ci,double *si)
/* Computes the cosine and sine integrals Ci(x) and Si(x). Ci(0) is returned as a large negative number and no error
 * message is generated. For x<0 routine returns -Si(-x) [correct] and Ci(-x), while actually Ci(x)=Ci(-x)-i*pi.
 */

{
	int i,k;
	bool odd;
	double a,err,fact,sign,sum,sumc,sums,t,term;
	doublecomplex h,b,c,d,del;

	t=fabs(x);
	// special case
	if (x==0) {
		*si=0;
		*ci=-BIG;
		return;
	}
	// Evaluate continued fraction by modified Lentz's method
	if (t>TMIN) {
		b = 1 + I*t;
		c=BIG;
		h=d=1/b;
		for (i=1;i<MAXIT;i++) {
			a=-i*i;
			b+=2;
			d=1/(a*d+b);
			c=b+a/c; // for i=1 c=+inf, so careful division should be performed to avoid overflows
			del=c*d;
			h*=del;
			if (fabs(creal(del)-1)+fabs(cimag(del))<=EPS) break;
		}
		if (i>=MAXIT) LogError(ALL_POS,"Failed to converge during calculation of sine integral of "GFORMDEF,x);
		h*=imExp(-t);
		*ci=-creal(h);
		*si=PI_OVER_TWO+cimag(h);
	}
	else { // Evaluate both series simultaneously
		// Special case: avoid failure of convergence test because of underflow
		if (t<sqrt(FPMIN)) {
			sumc=0;
			sums=t;
		}
		else {
			sum=sums=sumc=0;
			sign=fact=1;
			odd=true;
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
			if (k>MAXIT) LogError(ALL_POS,"Failed to converge during calculation of sine integral of "GFORMDEF,x);
		}
		*si=sums;
		*ci=sumc+log(t)+EULER;
	}
	if (x<0) *si=-(*si);
}
