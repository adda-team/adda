/* File: prec_time.c
 * $Date::                            $
 * Descr: precision timing routines (OS dependent); definitions (including inline) - in prec_time.h
 *
 * Copyright (C) 2006,2008,2010,2012-2014 ADDA contributors
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

/* This file should be compiled only in precise timing mode, hence the following declaration is redundant. However, it
 * helps proper syntax checking in IDE, such as Eclipse.
 */
#ifndef PRECISE_TIMING
#  define PRECISE_TIMING
#endif

#include "prec_time.h" // corresponding header, there all other needed includes are added

// LOCAL VARIABLES

#ifdef WINDOWS
static double inv_freq;
#endif

//======================================================================================================================

void InitTime(SYSTEM_TIME * restrict t)
// set time to zero
{
#ifdef WINDOWS
	t->QuadPart=0;
#elif defined(POSIX)
	t->tv_sec=t->tv_usec=0;
#endif
}

//======================================================================================================================

void SetTimerFreq(void)
// set frequency of windows timer; should be called once before running TimerToSec
{
#ifdef WINDOWS
	LARGE_INTEGER freq;

	QueryPerformanceFrequency(&freq);
	inv_freq=1/(double)freq.QuadPart;
#endif
}

//======================================================================================================================

double TimerToSec(const SYSTEM_TIME * restrict t)
// timer to seconds
{
#ifdef WINDOWS
	return (inv_freq*t->QuadPart);
#elif defined(POSIX)
	return (t->tv_sec+MICRO*t->tv_usec);
#endif
}

