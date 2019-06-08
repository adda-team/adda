/* File: bess.h
 * $Date::                            $
 * Descr:
 *
 * Copyright (C) 2006,2008-2013
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double BESSJ0 (double X);
double Sign(double X, double Y);
double BESSJ1 (double X);

void BESSJCS (int N2, double X, double ARRJ[5]);
void BESSJLP (int N, double X, double ARRJ[2]);
