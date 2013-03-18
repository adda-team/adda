/* File: Romberg.h
 * $Date::                            $
 * Descr: definitions of Romberg routines
 *
 * Copyright (C) 2006,2008,2010,2012-2013 ADDA contributors
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
#ifndef __Romberg_h
#define __Romberg_h

// project headers
#include "types.h" // needed for Parms_1D
// indexes of array
#define THETA 0
#define PHI   1

double Romberg1D(Parms_1D param,int size,const double * restrict data,double * restrict ss);

void Romberg2D(const Parms_1D parms_input[2],double (*func_input)(int theta,int phi,double * restrict res),
	int dim_input, double * restrict res, const char * restrict fname);

#endif // __Romberg_h
