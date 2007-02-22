/* FILE: Romberg.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of Romberg routines
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __Romberg_h
#define __Romberg_h

#include "types.h"  /* needed fot Parms_1D */
/* indexes of array */
#define THETA 0
#define PHI   1
/* minimum order to start extrapolation from */
#define ROMB_KMIN  3

double Romberg1D(Parms_1D param,int size,double *data,double *ss);

void Romberg2D(Parms_1D parms_input[2],void (*func_input)(int theta,int phi,double *res),
	       int dim_input, double *ss, char *fname);

#endif /*__Romberg_h*/
