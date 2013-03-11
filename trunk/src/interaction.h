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

void (*CalcInterTerm)(const int i,const int j,const int k,doublecomplex result[static restrict 6]);
void InitInteraction(void);
void FreeInteraction(void);

#endif //__interaction_h
