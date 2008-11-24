/* FILE: io.h
 * AUTH: Maxim Yurkin
 * DESCR: io routines
 */
#ifndef __io_h
#define __io_h

#include <stdio.h>

/* file locking is made quite robust, however it is a complex operation that
   can cause unexpected behaviour (permanent locks) especially when
   program is terminated externally (e.g. because of MPI failure).
   Moreover, it is not ANSI C, hence may have problems on some particular systems. */

/*#define NOT_USE_LOCK     /* uncomment to disable file locking */

#ifndef NOT_USE_LOCK
# define USE_LOCK
#endif

void LogError(int ErrCode,int who,char *FName,int line,char *Format,...);
void PrintBoth(FILE *file,char *fmt, ... );

void InitVariables(void);
void ParseParameters(int argc,char **argv);
void VariablesInterconnect(void);
void DirectoryLog(int argc,char **argv);
void PrintInfo(void);

#endif /* __io_h */




