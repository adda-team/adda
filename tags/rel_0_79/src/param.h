/* File: param.h
 * $Author$
 * $Date::                            $
 * Descr: inline routines for testing of input parameters
 *
 * Copyright (C) 2006,2008 University of Amsterdam
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
#ifndef __param_h
#define __param_h

#include "function.h" // needed for INLINE and function attributes

typedef struct {
	int l1; // first level index
	int l2; // second level index
} opt_index;

extern opt_index opt; // defined in param.c

void PrintErrorHelp(const char *fmt, ... ) ATT_PRINTF(1,2) ATT_NORETURN;

//============================================================

INLINE void TestPositive(const double val,const char *name)
// check if val is positive, otherwise produces error message
{
	if (val<=0) PrintErrorHelp("Illegal %s (%g), must be positive",name,val);
}

//============================================================

INLINE void TestNonNegative(const double val,const char *name)
// check if val is positive, otherwise produces error message
{
	if (val<0) PrintErrorHelp("Illegal %s (%g), must be nonnegative",name,val);
}

//============================================================

INLINE void TestPositive_i(const int val,const char *name)
// check if val (int) is positive, otherwise produces error message
{
	if (val<=0) PrintErrorHelp("Illegal %s (%d), must be positive",name,val);
}

//============================================================
/* In following 4 functions, one of two letters means either Including or Not-including (left and
 * right point of the interval respectively)
 */

INLINE void TestRangeII(const double val,const char *name,const double min,const double max)
// check if val is in interval [min,max], otherwise produces error message
{
	if (val<min || val>max)
		PrintErrorHelp("Illegal %s (%g), must belong to the interval [%g,%g]",name,val,min,max);
}

//============================================================

INLINE void TestRangeNI(const double val,const char *name,const double min,const double max)
// checks if val is in interval (min,max], otherwise produces error message
{
	if (val<=min || val>max)
		PrintErrorHelp("Illegal %s (%g), must belong to the interval (%g,%g]",name,val,min,max);
}

//============================================================

INLINE void TestRangeIN(const double val,const char *name,const double min,const double max)
// checks if val is in interval [min,max), otherwise produces error message
{
	if (val<min || val>=max)
		PrintErrorHelp("Illegal %s (%g), must belong to the interval [%g,%g)",name,val,min,max);
}

//============================================================

INLINE void TestRangeNN(const double val,const char *name,const double min,const double max)
// checks if val is in interval (min,max), otherwise produces error message
{
	if (val<=min || val>=max)
		PrintErrorHelp("Illegal %s (%g), must belong to the interval (%g,%g)",name,val,min,max);
}

//============================================================

INLINE void TestRange_i(const int val,const char *name,const int min,const int max)
// checks if val (int) is in interval [min,max], otherwise produces error message
{
	if (val<min || val>max)
		PrintErrorHelp("Illegal %s (%d), must belong to the interval [%d,%d]",name,val,min,max);
}

#endif // __param_h
