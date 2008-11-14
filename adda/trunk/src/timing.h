/* FILE: timing.h
 * AUTH: Maxim Yurkin
 * DESCR: Definitions for usual timing; should be completely portable
 *
 * Copyright (C) 2006,2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __timing_h
#define __timing_h

#include "parbas.h"

#ifdef MPI
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
