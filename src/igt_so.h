/* Functions implementing approximate averaging of the Green's tensor over the voxel volume. 
 * Designed to be conveniently used separately from ADDA as well.
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
// system headers
#include <complex.h>

typedef double complex doublecomplex_t;

/* the main function, which implements averaging (integral divided by volume) of the scaled Green's tensor (dyadic):
 * 4pi*k^2*G(r,r') = (k^2*I + nablaXnabla)[exp(ikR)/R]=exp(ikR)*R^-3*[(kR)^2*(I-Rproj)+(ikR-1)(I-3Rproj)]
 * with r' varied over the rectangular paralelipiped (voxel or dipole) with sizes ds_x, ds_y, ds_z.
 * Here I is unit tensor, R=r-r' (vector), and Rproj is tensor with components R_mu*R_nu/R^2 (projector on R),
 * r is given as rvec in input, and k - as wave_num. Result is a symmetric matrix - its 6 nontrivial components are
 * returned (Gxx,Gxy,Gxz,Gyy,Gyz,Gzz).
 *
 * When r is inside the voxel, the singular integral is considered as principal value (with spherical exclusion volume)
 * but the following scaled L-term for this exclusion volume is additionally added: -(4pi/3)I/Vd,
 * where Vd = ds_x*ds_y*ds_z (voxel volume). Thus the result is discontinuous when r crosses the voxel boundary, and 
 * loss of precision is expected near the boundary. Otherwise, the function should be robust for any r (TODO!!!).
 *
 * Analytic (relatively fast) formulaes are used with relative accuracy of O[(kd)^4] for any r,
 * where d is the largest voxel dimension
 */
void CalcIGTso(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6]);


/* the following three functions have the same API, but are used internally for different ranges of r
 * they should only be used for testing or benchmarking
 */
/* TODO: make another argument - mode (AUTO, NEAR, MEDIUM, FAR) defined by macros here (like in const.h)
 * then only one function will be needed (with AUTO recommended for general use) and in the main .c code we can
 * avoid repetition of common parts (like InitIGTvars) and wasting of computation of r2
 */
void CalcIGTso_near(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6]);
void CalcIGTso_medium(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6]);
void CalcIGTso_far(const double rvec[static restrict 3],const double wave_num,const double ds_x,const double ds_y,
	const double ds_z,doublecomplex_t result[static restrict 6]);

