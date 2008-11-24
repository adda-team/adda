/* FILE: Romberg.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of Romberg routines
 *
 * Copyright (C) 2006,2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __Romberg_h
#define __Romberg_h

#include "types.h" // needed for Parms_1D
// indexes of array
#define THETA 0
#define PHI   1

double Romberg1D(Parms_1D param,int size,const double *data,double *ss);

void Romberg2D(const Parms_1D parms_input[2],double (*func_input)(int theta,int phi,double *res),
               int dim_input, double *res, const char *fname);

#endif // __Romberg_h
