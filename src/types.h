/* File: types.h
 * $Date::                            $
 * Descr: definitions of various structures
 *
 * Copyright (C) 2006-2010,2012-2013 ADDA contributors
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
#ifndef __types_h
#define __types_h

// project headers
#include "const.h"   // for enum types
// system headers
#include <complex.h>
#include <stdbool.h> // for bool
#include <stddef.h> // for size_t

/* complex numbers; they are defined here so that headers that refer to doublecomplex can include only this small file
 * instead of large cmplx.h
 */
typedef double complex doublecomplex;

typedef struct	      // integration parameters
{
	double eps;       // convergence criterion
	int Jmin;         // minimal number of refinements
	int Jmax;         // maximal number of refinements
	double min;       // minimum
	double max;       // maximum
	size_t Grid_size; // number of grid points
	bool equival;     // whether max and min points are equivalent
	bool periodic;    // whether integrated function is periodic
} Parms_1D;

typedef struct	 // values of angles
{	             // !!! All angles are in degrees
	double min;  // minimum; for convenience (not really needed)
	double max;  // maximum; for convenience (not really needed)
	size_t N;    // number of points
	double * restrict val; // values of points; restrict should be minded in the code !!!
} angle_set;

typedef struct	        // integration parameters
{	                    // !!! All angles are in degrees
	enum scatgrid type; // if pairs are used or grid
	size_t N;           // total number of pairs (grid points)
	angle_set theta;    // values of theta
	angle_set phi;      // values of phi
} scat_grid_angles;

#endif // __types_h
