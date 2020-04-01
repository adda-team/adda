/* File: bess.h
 * $Date::      17/07/2019                      $
 * Descr:		Two functions "bessjn2" and "bessjn5" for computing 2 and 5 orders of Bessel functions respectively
 * Reference:	C++ Release 1.0 By J-P Moreau, Paris.(www.jpmoreau.fr) 
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

double bessj0 (double X);
double Sign(double X, double Y);
double bessj1 (double X);

//void bessjn5 (int N2, double X, double ARRJ[5]);
void bessjn2 (int N, double X, double ARRJ[2]);
