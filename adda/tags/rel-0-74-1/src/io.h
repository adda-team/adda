/* FILE: io.h
 * AUTH: Maxim Yurkin
 * DESCR: io and error handling routines
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __io_h
#define __io_h

#include <stdio.h>

/* file locking is made quite robust, however it is a complex operation that
   can cause unexpected behaviour (permanent locks) especially when
   program is terminated externally (e.g. because of MPI failure).
   Moreover, it is not ANSI C, hence may have problems on some particular systems.
   Currently file locking functions are only in param.c */

/*#define NOT_USE_LOCK     /* uncomment to disable file locking */

#ifndef NOT_USE_LOCK
# define USE_LOCK
#endif

void LogError(int ErrCode,int who,const char *fname,int line,const char *fmt,...);
void PrintError(const char *fmt, ... );
void LogPending(void);
void PrintBoth(FILE *file,const char *fmt, ... );

FILE *FOpenErr(const char *fname,const char *mode,int who,const char *err_fname,int lineN);
void FCloseErr(FILE *file,const char *fname,int who,const char *err_fname,int lineN);
void RemoveErr(const char *fname,int who,const char *err_fname,int lineN);
void MkDirErr(const char *dirname,int who,const char *err_fname,int lineN);

#endif /* __io_h */




