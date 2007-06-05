/* FILE: timing.h
 * AUTH: Maxim Yurkin
 * DESCR: Definitions for usual timing; should be completely portable
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __timing_h
#define __timing_h

#include "const.h"  /* for MPI */

#ifdef MPI
#include <mpi.h>
#define TIME_TYPE double
#define GET_TIME MPI_Wtime
#else
#include <time.h>
#define TIME_TYPE clock_t
#define GET_TIME clock
#endif

void StartTime(void);
void InitTiming(void);
void FinalStatistics(void);

#endif /* __timing_h */
