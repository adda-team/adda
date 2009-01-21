/* File: io.h
 * $Author$
 * $Date::                            $
 * Descr: i/o and error handling routines
 *
 * Copyright (C) 2006,2008 University of Amsterdam
 * Copyright (C) 2009 Institute of Chemical Kinetics and Combustion & University of Amsterdam
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU
 * General Public License as published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
 * Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
#ifndef __io_h
#define __io_h

#include <stdio.h>    // for file
#include "function.h" // for function attributes
#include "const.h"    // for enum types

/* File locking is made quite robust, however it is a complex operation that can cause unexpected
 * behavior (permanent locks) especially when program is terminated externally (e.g. because of MPI
 * failure). Moreover, it is not ANSI C, hence may have problems on some particular systems.
 * Currently file locking functions are only in param.c
 */

//#define NOT_USE_LOCK  // uncomment to disable file locking
//#define ONLY_LOCKFILE // uncomment to use only lock file, without file locking over NFS

#ifndef NOT_USE_LOCK
#	define USE_LOCK
#	ifndef ONLY_LOCKFILE
#		define LOCK_FOR_NFS // currently this works only for POSIX
#	endif
#endif

void WrapLines(char *str);
char *WrapLinesCopy(const char *str);
void LogError(enum ec ErrCode,enum enwho who,const char *fname,int line,
	const char *fmt,...) ATT_PRINTF(5,6);
void PrintError(const char *fmt, ... ) ATT_PRINTF(1,2) ATT_NORETURN;
void LogPending(void);
void PrintBoth(FILE *file,const char *fmt, ... ) ATT_PRINTF(2,3);

FILE *FOpenErr(const char *fname,const char *mode,enum enwho who,const char *err_fname,
               int lineN) ATT_MALLOC;
void FCloseErr(FILE *file,const char *fname,enum enwho who,const char *err_fname,int lineN);
void RemoveErr(const char *fname,enum enwho who,const char *err_fname,int lineN);
void MkDirErr(const char *dirname,enum enwho who,const char *err_fname,int lineN);

#endif // __io_h
