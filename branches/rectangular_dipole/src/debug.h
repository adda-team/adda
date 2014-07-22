/* File: debug.h
 * $Date::                            $
 * Descr: definitions for debug functions
 *
 * Copyright (C) 2006,2008-2010,2012-2013 ADDA contributors
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
#ifndef __debug_h
#define __debug_h

/* Debugging implies turning on additional information messages during the code execution. A simple and convenient tool
 * to generate such messages is used.
 */

// this are used not always, but only for manual debug (hence can be always defined)
#define PRINT_REAL(x) printf(#x"=%g\n",x)
#define PRINT_COMPLEX(x) printf(#x"=%g%+gi\n",REIM(x))

#ifdef DEBUGFULL

#	include "const.h"    // for POS
#	include "function.h" // for function attributes
#	include "io.h"       // for ERR_LOC_DECL
#	include "types.h"    // for doublecomplex
#	define D(...) DebugPrintf(ALL_POS,__VA_ARGS__)
#	define Dz(...) DebugPrintf(ONE_POS,__VA_ARGS__)
void DebugPrintf(ERR_LOC_DECL,const char * restrict fmt,...) ATT_PRINTF(4,5);
void FieldPrint(doublecomplex * restrict x) ATT_UNUSED;
#else
#	define D(...)
#	define Dz(...)
#endif

#endif // __debug_h
