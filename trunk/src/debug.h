/* File: debug.h
 * $Author$
 * $Date::                            $
 * Descr: definitions for debug functions
 *
 *        Previous versions by "vesseur"
 *
 * Copyright (C) 2006,2008 University of Amsterdam
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
#ifndef __debug_h
#define __debug_h

/* Debugging implies turning on additional information messages during the code execution. A simple
 * and convenient tool to generate such messages is used.
 */

//#define DEBUG // uncomment to degug

#ifdef DEBUG

#	include "function.h" // for function attributes
#	define D(p) DebugPrintf(__FILE__,__LINE__,p)
#	define D2(p,a) DebugPrintf(__FILE__,__LINE__,p,a)
#	define D2z(p,a) if (ringid==ROOT) DebugPrintf(__FILE__,__LINE__,p,a)
void DebugPrintf(const char *fname,int line,const char *fmt,...) ATT_PRINTF(3,4);
void FieldPrint(doublecomplex *x) ATT_UNUSED;
#else
#	define D(p)
#	define D2(p,a)
#	define D2z(p,a)
#endif

#endif // __debug_h
