/* Inline complex functions, functions on length-3 real and complex vectors, and several auxiliary functions
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
#include <complex.h>

typedef double complex doublecomplex_t;

void calc_igt_so(double qvec[static 3], const double wave_num, const double ds_x, const double ds_y, const double ds_z, doublecomplex_t result[static restrict 6]);
void calc_igt_so_near(double qvec[static 3], const double wave_num, const double ds_x, const double ds_y, const double ds_z, doublecomplex_t result[static restrict 6]);
void calc_igt_so_middle(double qvec[static 3], const double wave_num, const double ds_x, const double ds_y, const double ds_z, doublecomplex_t result[static restrict 6]);
void calc_igt_so_far(double qvec[static 3], const double wave_num, const double ds_x, const double ds_y, const double ds_z, doublecomplex_t result[static restrict 6]);

