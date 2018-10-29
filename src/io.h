/* File: io.h
 * $Date::                            $
 * Descr: i/o and error handling routines
 *
 * Copyright (C) 2006,2008-2013
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
#ifndef __io_h
#define __io_h

// project headers
#include "const.h"    // for enum types
#include "function.h" // for function attributes
// system headers
#include <stdio.h>    // for file
#include <stdarg.h>   // for va_list

/* File locking is made quite robust, however it is a complex operation that can cause unexpected behavior (permanent
 * locks) especially when program is terminated externally (e.g. because of MPI failure). Moreover, it is not ANSI C,
 * hence may have problems on some particular systems. Currently file locking functions are only in param.c
 */

//#define NOT_USE_LOCK  // uncomment to disable file locking
//#define ONLY_LOCKFILE // uncomment to use only lock file, without file locking over NFS

#ifndef NOT_USE_LOCK
#	define USE_LOCK
#	ifndef ONLY_LOCKFILE
#		define LOCK_FOR_NFS // currently this works only for POSIX
#	endif
#endif

// Common parts of function declaration and calls; they are passed to ProcessError and DebugPrintf
#define ERR_LOC_DECL const enum enwho who,const char * restrict srcfile,const int srcline
#define ERR_LOC_CALL who,srcfile,srcline

// A way of calling snprintf and vsnprintf resistant to buffer overflows, but without errors
#define SNPRINTF_SHIFT_ROBUST(shift,tmp,str,size,...) { \
	tmp=snprintf(str+shift,size-shift,__VA_ARGS__); \
	if (tmp>0) { shift+=tmp; if (shift>=size) shift=size-1; } \
}
#define VSNPRINTF_SHIFT_ROBUST(shift,tmp,str,size,...) { \
	tmp=vsnprintf(str+shift,size-shift,__VA_ARGS__); \
	if (tmp>0) { shift+=tmp; if (shift>=size) shift=size-1; } \
}

char *dyn_sprintf(const char *format, ...) ATT_PRINTF(1,2) ATT_MALLOC;
char *rea_sprintf(char *str,const char *format, ...) ATT_PRINTF(2,3) ATT_MALLOC;
void WrapLines(char * restrict str);
char *WrapLinesCopy(const char * restrict str);
void LogError(ERR_LOC_DECL,const char * restrict fmt,...) ATT_PRINTF(4,5) ATT_NORETURN;
void LogWarning(enum ec code,ERR_LOC_DECL,const char * restrict fmt,...) ATT_PRINTF(5,6);
size_t SnprintfErr(ERR_LOC_DECL,char * restrict str,const size_t size,const char * restrict fmt,...) ATT_PRINTF(6,7);
size_t SnprintfShiftErr(ERR_LOC_DECL,const size_t shift,char * restrict str,const size_t size,const char * restrict fmt,
	...) ATT_PRINTF(7,8);
size_t VsnprintfErr(ERR_LOC_DECL,char * restrict str,const size_t size,const char * restrict fmt,va_list args);
void PrintError(const char * restrict fmt, ... ) ATT_PRINTF(1,2) ATT_NORETURN;
void LogPending(void);
void PrintBoth(FILE * restrict file,const char * restrict fmt, ... ) ATT_PRINTF(2,3);

FILE *FOpenErr(const char * restrict fname,const char * restrict mode,ERR_LOC_DECL);
void FCloseErr(FILE * restrict file,const char * restrict fname,ERR_LOC_DECL);
void RemoveErr(const char * restrict fname,ERR_LOC_DECL);
void MkDirErr(const char * restrict dirname,ERR_LOC_DECL);

char *FGetsError(FILE * restrict file,const char * restrict fname,size_t *line,char * restrict buf,const int buf_size,
	ERR_LOC_DECL);
size_t SkipNLines(FILE * restrict file,const size_t n);
size_t SkipComments(FILE * restrict file);

#endif // __io_h
