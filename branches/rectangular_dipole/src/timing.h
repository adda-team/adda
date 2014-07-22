/* File: timing.h
 * $Date::                            $
 * Descr: definitions for usual timing; should be completely portable
 *
 * Copyright (C) 2006,2008-2009,2012-2014 ADDA contributors
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
#ifndef __timing_h
#define __timing_h

// project headers
#include "os.h"
#include "parbas.h"

#ifdef ADDA_MPI
#	define TIME_TYPE double
#	define GET_TIME() MPI_Wtime()
#else
#	include <time.h>
#	define TIME_TYPE clock_t
#	define GET_TIME() clock()
#endif

#ifdef WINDOWS
#	define SYSTEM_TIME LARGE_INTEGER
#	define GET_SYSTEM_TIME(t) QueryPerformanceCounter(t)
#elif defined(POSIX)
/* It make sense to switch to clock_gettime, especially with types CLOCK_MONOTONIC or CLOCK_MONOTONIC_RAW, to be
 * independent of wall time synchronization, etc. (as is now for Windows functions). However, that would be not so
 * portable and is probably overkill.
 */
#	include <sys/time.h> // for timeval and gettimeofday
#	include <stdio.h>    // needed for definition of NULL
#	define SYSTEM_TIME struct timeval
// gettimeofday is described only in POSIX 1003.1-2001, but it should work for many other systems
#	define GET_SYSTEM_TIME(t) gettimeofday(t,NULL)
#else
#	include <time.h>
#	define SYSTEM_TIME time_t
#	define GET_SYSTEM_TIME(t) time(t)
#endif

void StartTime(void);
void InitTiming(void);
void FinalStatistics(void);
double DiffSystemTime(const SYSTEM_TIME * restrict t1,const SYSTEM_TIME * restrict t2);

#endif // __timing_h
