/* File: debug.c
 * $Author$
 * $Date::                            $
 * Descr: functions for printing debugging information when compiling with option -DDEBUG
 *
 *        Previous versions by "vesseur"
 *
 * Copyright (C) 2006-2008 University of Amsterdam
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
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "types.h"
#include "debug.h"
#include "const.h"

#ifdef DEBUG
//============================================================

void DebugPrintf(const char *fname,const int line,const char *fmt, ... )
{
	va_list args;
	char pr_line[MAX_PARAGRAPH];
	extern int ringid;

	va_start(args, fmt);
	sprintf(pr_line,"DEBUG: %s:%d: ",fname,line);
	vsprintf(pr_line+strlen(pr_line),fmt,args);
	strcat(pr_line, "\n");
	printf("(ringID=%i) %s",ringid,pr_line);
	fflush(stdout);
	va_end(args);
}

//=======================================================

void FieldPrint (doublecomplex *x)
/* print current field at certain dipole -- not used; left for deep debug; NOT ROBUST, since
 * DipoleCoord is not always available
 */
{
	extern FILE *logfile;
	extern double *DipoleCoord;
	int i=9810;

	i*=3;
	fprintf(logfile,"Dipole coordinates = "GFORM3V"\n",
		DipoleCoord[i],DipoleCoord[i+1],DipoleCoord[i+2]);
	fprintf(logfile,"E = "CFORM3V,x[i][RE],x[i][IM],x[i+1][RE],x[i+1][IM],x[i+2][RE],x[i+2][IM]);
}

#endif // DEBUG
