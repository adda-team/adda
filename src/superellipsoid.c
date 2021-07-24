/* Routines for determining volume ratio of superellipsoids.
 *
 * Copyright (C) ADDA contributors
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
#include "io.h"
// system headers
#include <math.h>

// LOCAL VARIABLES
// Coefficients for Lanczos approximation for Gamma function taken from
// Press, Teukolsky, Vertterling, and Flannery, Numerical Recipes, 3rd Ed.,
// (2007), Ch. 6.1
static const double coeffs[] = {
  0.999999999999997092,
  57.1562356658629235,
  -59.5979603554754912,
  14.1360979747417471,
  -0.491913816097620199,
  .339946499848118887e-4,
  .465236289270485756e-4,
  -.983744753048795646e-4,
  .158088703224912494e-3,
  -.210264441724104883e-3,
  .217439618115212643e-3,
  -.164318106536763890e-3,
  .844182239838527433e-4,
  -.261908384015814087e-4,
  .368991826595316234e-5
} ;

static const double g = 4.7421875000000;

//======================================================================================================================

static double Gamma(double x)
// Calculate the gamma function.
// Argument must be > 0, but no check performed here.
// Algorithm used here works on the complex plane, but only need it on the real line.
{
  int i;
  const int N = 14; // Number of terms in Lanczos sum
  double t, result;

  if (fmod(x, 1) == 0) {
    // Calculate the factorial for an integer input.
    // Recall that Gamma(z + 1) = z! for integer z.
    result = 1.;
    for (i = 1; i < x; i++) { // Terminate early!
      result *= i;
    }
    return result;
  }
  else {
    /* For non-integer inputs, use the Lanczos approximation
    as described in Numerical Recipes (2007), ch. 6.1.
    Note that in their code g = 671/128 - 1/2. */

    // Partial sum
    result = coeffs[0];

    for (i = 1; i <= N; i++){
      result += coeffs[i]/(x-1+i);
    }

    t = x+g-0.5;

    return sqrt(TWO_PI) * pow(t,x-0.5) * exp(-t) * result;
  }
}

//======================================================================================================================

static double Beta(double x,double y)
// Calculate the beta function
{
  return Gamma(x)*Gamma(y)/Gamma(x+y);
}

//======================================================================================================================

double SuperellipsoidVolumeRatio(double aspY,double aspZ,double n,double e)
// Calculate volume ratio for a superellipsoid.
// Volume analytically given by Wriedt (2002) eq. 4, derived in Ref. 40.
// Note that the relevant box volume is 8a^3. So, the ratio is
// ratio = (1/4) (b/a) (c/a) n B(n/2 +1, n) e B(e/2, e/2).
{
  double nterm,eterm;

  // lim_{n -> 0} n B(n/2 + 1, n) = 1, and lim_{e -> 0} e B(e/2, e/2) = 4.
  // But direct calculation doesn't work if n and e are too small.
  // Below a threshold, set them to limiting values.
  // The threshold is lower than what is probably practical for DDA calculations
  // since we use finite-size dipoles.

  const double veps = 1e-5;
  if (n<veps) nterm = 1;
  else nterm = n*Beta(n/2+1,n);

  if (e<veps) eterm = 4;
  else eterm = e*Beta(e/2,e/2);

  return 0.25 * aspY * aspZ * nterm * eterm;
}
