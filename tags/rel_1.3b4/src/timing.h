/* File: timing.h
 * $Date::                            $
 * Descr: definitions for usual timing; should be completely portable
 *
 * Copyright (C) 2006,2008-2009,2012-2013 ADDA contributors
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
#include "parbas.h"

#ifdef ADDA_MPI
#	define TIME_TYPE double
#	define GET_TIME MPI_Wtime
#else
#	include <time.h>
#	define TIME_TYPE clock_t
#	define GET_TIME clock
#endif

void StartTime(void);
void InitTiming(void);
void FinalStatistics(void);

#endif // __timing_h
