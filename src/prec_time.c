/* FILE: prec_time.c
 * AUTH: Maxim Yurkin
 * DESCR: precision timing routines (OS dependent)
 *        definitions and inlines are in timing.h
 */

#include "prec_time.h"

#ifdef PRECISE_TIMING /* following is only for precise timing */

#ifdef _WIN32     /* Windows */
double inv_freq;
#else             /* UNIX */
# define MICRO 1E-6
#endif

/*============================================================*/

void InitTime(SYSTEM_TIME *t)
     /* set time to zero */
{
#ifdef _WIN32  /* Windows */
  t->QuadPart=0;
#else          /* UNIX */
  t->tv_sec=t->tv_usec=0;
#endif
}

/*============================================================*/

void SetTimerFreq(void)
   /* set frequency of windows timer;
   should be called once before running TimerToSec or DiffSec*/
{
#ifdef _WIN32  /* Windows */
  LARGE_INTEGER freq;

  QueryPerformanceFrequency(&freq);
  inv_freq=1/(double)freq.QuadPart;
#endif
}

/*============================================================*/

double TimerToSec(SYSTEM_TIME *t)
   /* timer to seconds */
{
#ifdef _WIN32   /* Windows */
  return (inv_freq*t->QuadPart);
#else           /* UNIX */
  return (t->tv_sec+MICRO*t->tv_usec);
#endif
}

/*============================================================*/

double DiffSec(SYSTEM_TIME *t1,SYSTEM_TIME *t2)
   /* difference between two times in seconds */
{
  SYSTEM_TIME res;

  elapsed(t1,t2,&res);
  return TimerToSec(&res);
}

#endif /* PRECISE_TIMING */

