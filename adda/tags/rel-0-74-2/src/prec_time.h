/* FILE: prec_time.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of inline functions for
 *        precise timing
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __prec_time_h
#define __prec_time_h

/* Precise timing gives an accuracy of order micro_sec. It gives extensive information
   on timing of FFT init, Dmatrix init, and Matrix Vector multiplication.
   It is optimized not to take by itself as little time as possible.
   It is used mostly for locating and optimizing the bottlenecks of the code execution.
   It is not ANSI C, therefore is system dependent, though is expected to work for most. */

/*#define PRECISE_TIMING  /* uncomment to conduct precise timing */

#ifdef PRECISE_TIMING

#include "const.h"  /* needed for INLINE */

#ifdef _WIN32  /* Windows */
# include <windows.h>
# define SYSTEM_TIME LARGE_INTEGER
#else    /* UNIX */
/* this works for many different systems, if gettimeofday is present */
# include <sys/time.h>
# include <stdio.h>    /* needed for definition of NULL */
# define SYSTEM_TIME struct timeval
#endif

void InitTime(SYSTEM_TIME *t);
void SetTimerFreq(void);
double TimerToSec(const SYSTEM_TIME *t) ATT_PURE;
double DiffSec(const SYSTEM_TIME *t1,const SYSTEM_TIME *t2) ATT_PURE;

/*============================================================*/

INLINE void elapsed(const SYSTEM_TIME *t1,const SYSTEM_TIME *t2,SYSTEM_TIME *res)
     /* compute time difference */
{
#ifdef _WIN32  /* Windows */
  res->QuadPart=t2->QuadPart-t1->QuadPart;
#else          /* UNIX */
  res->tv_sec=t2->tv_sec-t1->tv_sec;
  res->tv_usec=t2->tv_usec-t1->tv_usec;
#endif
}

/*============================================================*/

INLINE void ElapsedInc(const SYSTEM_TIME *t1,const SYSTEM_TIME *t2,SYSTEM_TIME *res)
     /* compute time difference, increment result by this value */
{
#ifdef _WIN32  /* Windows */
  res->QuadPart+=(t2->QuadPart-t1->QuadPart);
#else          /* UNIX */
  res->tv_sec+=(t2->tv_sec-t1->tv_sec);
  res->tv_usec+=(t2->tv_usec-t1->tv_usec);
#endif
}

/*============================================================*/

INLINE void GetTime(SYSTEM_TIME *t)
     /* compute time difference, increment result by this value */
{
#ifdef _WIN32  /* Windows */
  QueryPerformanceCounter(t);
#else          /* UNIX */
  gettimeofday(t,NULL);
#endif
}

#endif /* PRECISE_TIMING */

#endif /*__prec_time_h*/
