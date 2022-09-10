/* Header for calculation of Sommerfeld integrals
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
#ifndef __somnec_h
#define __somnec_h

// system headers
/* ADDA uses doublecomplex macro instead of standard 'complex double', but we use the latter in somnec.c/h to
 * keep them standalone
 */
#include <complex.h>

void som_init(complex double epscf);
void evlua(double zphIn,double rhoIn,complex double *erv,complex double *ezv,complex double *erh,complex double *eph,
	int mode);

#endif // __somnec_h
