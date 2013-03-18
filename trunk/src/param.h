/* File: param.h
 * $Date::                            $
 * Descr: inline routines for testing of input parameters
 *
 * Copyright (C) 2006,2008,2010-2013 ADDA contributors
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
#ifndef __param_h
#define __param_h

// project headers
#include "function.h" // needed for function attributes
// system headers
#include <limits.h> // for INT_MIN and INT_MAX
#include <math.h>

typedef struct {
	int l1; // first level index
	int l2; // second level index
} opt_index;

extern opt_index opt; // defined in param.c

void PrintErrorHelp(const char * restrict fmt, ... ) ATT_PRINTF(1,2) ATT_NORETURN;

//======================================================================================================================

static inline void TestPositive(const double val,const char * restrict name)
// check if val is positive, otherwise produces error message
{
	if (val<=0) PrintErrorHelp("Illegal %s ("GFORMDEF"), must be positive",name,val);
}

//======================================================================================================================

static inline void TestNonNegative(const double val,const char * restrict name)
// check if val is nonnegative, otherwise produces error message
{
	if (val<0) PrintErrorHelp("Illegal %s ("GFORMDEF"), must be nonnegative",name,val);
}

//======================================================================================================================

static inline void TestPositive_i(const int val,const char * restrict name)
// check if val (int) is positive, otherwise produces error message
{
	if (val<=0) PrintErrorHelp("Illegal %s (%d), must be positive",name,val);
}

//======================================================================================================================

static inline void TestNonNegative_i(const int val,const char * restrict name)
// check if val (int) is nonnegative, otherwise produces error message
{
	if (val<0) PrintErrorHelp("Illegal %s (%d), must be nonnegative",name,val);
}

//======================================================================================================================
/* In following 4 functions, one of two letters means either Including or Not-including (left and right point of the
 * interval respectively)
 */

static inline void TestRangeII(const double val,const char * restrict name,const double min,const double max)
// check if val is in interval [min,max], otherwise produces error message
{
	if (val<min || val>max)
		PrintErrorHelp("Illegal %s ("GFORMDEF"), must belong to the interval ["GFORMDEF","GFORMDEF"]",name,val,min,max);
}

//======================================================================================================================

static inline void TestRangeNI(const double val,const char * restrict name,const double min,const double max)
// checks if val is in interval (min,max], otherwise produces error message
{
	if (val<=min || val>max)
		PrintErrorHelp("Illegal %s ("GFORMDEF"), must belong to the interval ("GFORMDEF","GFORMDEF"]",name,val,min,max);
}

//======================================================================================================================

static inline void TestRangeIN(const double val,const char * restrict name,const double min,const double max)
// checks if val is in interval [min,max), otherwise produces error message
{
	if (val<min || val>=max)
		PrintErrorHelp("Illegal %s ("GFORMDEF"), must belong to the interval ["GFORMDEF","GFORMDEF")",name,val,min,max);
}

//======================================================================================================================

static inline void TestRangeNN(const double val,const char * restrict name,const double min,const double max)
// checks if val is in interval (min,max), otherwise produces error message
{
	if (val<=min || val>=max)
		PrintErrorHelp("Illegal %s ("GFORMDEF"), must belong to the interval ("GFORMDEF","GFORMDEF")",name,val,min,max);
}

//======================================================================================================================

static inline void ConvertToInteger(double val,const char * restrict name,int *res)
/* converts val to res, but first checks if val is really an integer and in the bounds, otherwise produces an error
 * message
 */
{
	if (val != floor(val)) PrintErrorHelp("Illegal %s ("GFORMDEF"), must be an integer",name,val);
	if (val<INT_MIN || val>INT_MAX) PrintErrorHelp("Illegal %s ("GFORMDEF"), must be inside integer bounds",name,val);
	*res=(int)val;
}

//======================================================================================================================

static inline void TestRange_i(const int val,const char * restrict name,const int min,const int max)
// checks if val (int) is in interval [min,max], otherwise produces error message
{
	if (val<min || val>max) PrintErrorHelp("Illegal %s (%d), must belong to the interval [%d,%d]",name,val,min,max);
}

//======================================================================================================================

static inline void TestGreaterThan_i(const int val,const char * restrict name,const int min)
// checks if val (int) is greater than min, otherwise produces error message
{
	if (val<=min) PrintErrorHelp("Illegal %s (%d), must be greater than %d",name,val,min);
}

#endif // __param_h
