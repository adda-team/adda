/* FILE: timing.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of inline functions for
 *        precise timing
 */
#ifndef __timing_h
#define __timing_h

/*#define PRECISE_TIMING    /* uncomment to conduct precise timing */

#ifdef PRECISE_TIMING     /* precise timing definitions and functions
                             expected accuracy of order micro_sec */
#include "const.h"  /* needed for INLINE */

#ifdef _WIN32  /* Windows */
# include <windows.h>
# define SYSTEM_TIME LARGE_INTEGER
#else    /* UNIX */
# include <sys/time.h>
# include <stdio.h>    /* needed for definition of NULL */
# define SYSTEM_TIME struct timeval
#endif

void InitTime(SYSTEM_TIME *t);
void SetTimerFreq(void);
double TimerToSec(SYSTEM_TIME *t);
double DiffSec(SYSTEM_TIME *t1,SYSTEM_TIME *t2);

/*============================================================*/

INLINE void elapsed(SYSTEM_TIME *t1,SYSTEM_TIME *t2,SYSTEM_TIME *res)
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

INLINE void ElapsedInc(SYSTEM_TIME *t1,SYSTEM_TIME *t2,SYSTEM_TIME *res)
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

#endif /*__timing_h*/
