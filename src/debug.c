/* File: debug.c
 * $Date::                            $
 * Descr: functions for printing debugging information when compiling with option -DDEBUGFULL
 *
 * Copyright (C) 2006-2008,2010,2012-2013 ADDA contributors
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

/* This file should be compiled only in full debug mode, hence the following declaration is redundant. However, it helps
 * proper syntax checking in IDE, such as Eclipse.
 */
#ifndef DEBUGFULL
#  define DEBUGFULL
#endif

#include "const.h" // keep this first
#include "debug.h" // corresponding header
// project headers
#include "cmplx.h"
#include "comm.h"
#include "io.h"
#include "types.h"
#include "vars.h"
// system headers
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

//======================================================================================================================

void DebugPrintf(ERR_LOC_DECL,const char * restrict fmt, ... )
/* Prints anything to stdout with additional debug information
 *
 * Not thread-safe! Should not be called in parallel from multiple threads (e.g. OpenMP)
 */
{
	va_list args;
	static char msg[MAX_PARAGRAPH]; // not to allocate it at every call

	if (who==ALL || IFROOT) { // controls whether output should be produced
		va_start(args,fmt);
		VsnprintfErr(ERR_LOC_CALL,msg,MAX_PARAGRAPH,fmt,args);
#ifdef PARALLEL
		if (who==ALL) printf("(ringID=%i) DEBUG: %s:%d: %s \n",ringid,srcfile,srcline,msg);
		else
#endif
		printf("DEBUG: %s:%d: %s \n",srcfile,srcline,msg);
		fflush(stdout);
		va_end(args);
	}
}

//======================================================================================================================

void FieldPrint (doublecomplex * restrict x)
/* print current field at certain dipole -- not used; left for deep debug; NOT ROBUST, since DipoleCoord is not always
 * available
 */
{
	int i=9810;

	i*=3;
	fprintf(logfile,"Dipole coordinates = "GFORM3V"\n",COMP3V(DipoleCoord+i));
	fprintf(logfile,"E = "CFORM3V,REIM3V(x+i));
}
