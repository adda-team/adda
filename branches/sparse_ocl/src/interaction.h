/* FILE : interaction.h
 * Descr: the functions used to calculate the interaction term
 *
 * Copyright (C) 2011-2013 ADDA contributors
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
#ifndef __interaction_h
#define __interaction_h

#include "cmplx.h"

/* Calculates interaction term between two dipoles; given integer distance vector {i,j,k} (in units of d). The acting
 * dipole is placed in the origin, while the field is calculated at position given as argument. All six components of
 * the symmetric matrix are computed at once. The elements in result are: [G11, G12, G13, G22, G23, G33]
 */
void (*InterTerm_int)(const int i,const int j,const int k,doublecomplex result[static restrict 6]);
// same as above, but distance is passed as a double vector (in um)
void (*InterTerm_real)(const double qvec[static restrict 3],doublecomplex result[static restrict 6]);

/* Calculates reflection term between two dipoles; given integer distance vector {i,j,k} (in units of d). k is the _sum_
 * of dipole indices along z with respect to the center of bottom dipoles of the particle. Bottom is considered for the
 * current processor (position) and the whole particle (position_full) in FFT and SPARSE modes respectively. The latter
 * behavior is determined by ZsumShift.
 * The acting dipole is placed in the origin, while the field is calculated at position given as argument.
 * Six components of the matrix are computed at once: [GR11, GR12, GR13, GR22, GR23, GR33].
 * The matrix is not symmetric, but satisfies: GR21=GR12, GR31=-GR13, GR32=-GR23. However the large matrix GR, which
 * acts on the total vector of dipole polarizations is still complex-symmetric, since GR[i,j]=GR^T[j,i] (interchange of
 * i and j changes only the sign of x and y components, but not z, which leads to sign change of 13,31,23,32 components)
 */
void (*ReflTerm_int)(const int i,const int j,const int k,doublecomplex result[static restrict 6]);
/* same as above, but distance is passed as a double vector (in um) and its z-component is the sum of heights of
 * source and probe points above the surface. So qvec_in is the actual distance between probe point and the image of the
 * source point.
 */
void (*ReflTerm_real)(const double qvec[static restrict 3],doublecomplex result[static restrict 6]);

void InitInteraction(void);
void FreeInteraction(void);

#endif //__interaction_h
