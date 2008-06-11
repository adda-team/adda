/* FILE: debug.h
 * AUTH: Maxim Yurkin
 * DESCR: Definitions for debug functions
 *
 *        Previous versions by "vesseur"
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __debug_h
#define __debug_h

/* Debugging implies turning on additional information messages
   during the code execution. A simple and convenient tool to
   generate such messages is used. */

/*#define DEBUG     /* uncomment to degug */

#ifdef DEBUG

#include "function.h"  /* for function attributes */

# define D(p)        DebugPrintf(__FILE__, __LINE__, p)
# define D2(p,a)     DebugPrintf(__FILE__, __LINE__, p,a)
# define D2z(p,a)    if (ringid==ROOT) DebugPrintf(__FILE__, __LINE__, p,a)
void DebugPrintf(const char *fname,int line,const char *fmt, ...) ATT_PRINTF(3,4);
void FieldPrint(doublecomplex *x) ATT_UNUSED;
#else
# define D(p)
# define D2(p,a)
# define D2z(p,a)
#endif

#endif /*__debug_h*/
