/* FILE : interaction.h
 * Descr: the functions used to calculate the interaction term
 *
 * Copyright (C) 2006-2012 ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef __interaction_h
#define __interaction_h

#include "cmplx.h"

extern double gridspace, WaveNum;

#ifdef USE_SSE3
extern __m128d c1, c2, c3, zo, inv_2pi, p360, prad_to_deg;
extern __m128d exptbl[361];
#endif //USE_SSE3

extern void (*CalcInterTerm)(const int i,const int j,const int k,doublecomplex * restrict result);
extern void InitInteraction(void);

//============================================================

#endif //__interaction_h
