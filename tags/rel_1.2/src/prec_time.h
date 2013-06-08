/* File: prec_time.h
 * $Date::                            $
 * Descr: definitions of inline functions for precise timing
 *
 * Copyright (C) 2006-2008,2010,2012-2013 ADDA contributors
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
#ifndef __prec_time_h
#define __prec_time_h

/* Precise timing gives an accuracy of order micro_sec. It gives extensive information on timing of FFT initialization,
 * D-matrix initialization, and Matrix Vector multiplication. It is optimized to consume as little time as possible by
 * itself. It is used mostly for locating and optimizing the bottlenecks of the code execution. It is not ANSI C,
 * therefore is system dependent, though is expected to work for most.
 */

#ifdef PRECISE_TIMING

// project headers
#include "os.h"       // for OS definitions
#include "function.h" // for static inline and function attributes

#ifdef WINDOWS
#	define SYSTEM_TIME LARGE_INTEGER
#elif defined(POSIX)
#	include <sys/time.h> // for timeval and gettimeofday
#	include <stdio.h>    // needed for definition of NULL
#	define SYSTEM_TIME struct timeval
#else
#	error "Unknown operation system. Precise timing is not supported."
#endif

#define FFORMPT "%.4f" // format for precise timing results

void InitTime(SYSTEM_TIME * restrict t);
void SetTimerFreq(void);
double TimerToSec(const SYSTEM_TIME * restrict t) ATT_PURE;
double DiffSec(const SYSTEM_TIME * restrict t1,const SYSTEM_TIME * restrict t2) ATT_PURE;

//======================================================================================================================

static inline void Elapsed(const SYSTEM_TIME * restrict t1,const SYSTEM_TIME * restrict t2,SYSTEM_TIME * restrict res)
// compute time difference
{
#ifdef WINDOWS
	res->QuadPart=t2->QuadPart-t1->QuadPart;
#elif defined(POSIX)
	res->tv_sec=t2->tv_sec-t1->tv_sec;
	res->tv_usec=t2->tv_usec-t1->tv_usec;
#endif
}

//======================================================================================================================

static inline void ElapsedInc(const SYSTEM_TIME * restrict t1,const SYSTEM_TIME * restrict t2,
	SYSTEM_TIME * restrict res)
// compute time difference, increment result by this value
{
#ifdef WINDOWS
	res->QuadPart+=(t2->QuadPart-t1->QuadPart);
#elif defined(POSIX)
	res->tv_sec+=(t2->tv_sec-t1->tv_sec);
	res->tv_usec+=(t2->tv_usec-t1->tv_usec);
#endif
}

//======================================================================================================================

static inline void GetTime(SYSTEM_TIME * restrict t)
// get current time
{
#ifdef WINDOWS
	QueryPerformanceCounter(t);
#elif defined(POSIX)
	// gettimeofday is described only in POSIX 1003.1-2001, but it should work for many other systems
	gettimeofday(t,NULL);
#endif
}

#endif // PRECISE_TIMING

#endif // __prec_time_h
