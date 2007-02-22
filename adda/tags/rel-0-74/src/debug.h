/* FILE: debug.h
 * AUTH: Maxim Yurkin
 * DESCR: Definitions for debug functions
 *
 *        Previous versions by "vesseur"
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __debug_h
#define __debug_h

/* Debugging implies turning on additional information messages
   during the code execution. A simple and convenient tool to
   generate such messages is used. */

/*#define DEBUG     /* uncomment to degug */

#ifdef DEBUG
# define D(p)        DebugPrintf(__FILE__, __LINE__, p)
# define D2(p,a)     DebugPrintf(__FILE__, __LINE__, p,a)
# define D3(p,a,b)   DebugPrintf(__FILE__, __LINE__, p,a,b)
# define D4(p,a,b,c) DebugPrintf(__FILE__, __LINE__, p,a,b,c)
void DebugPrintf(const char *fname,int line,const char *fmt, ...);
void FieldPrint(doublecomplex *x);
#else
# define D(p)
# define D2(p,a)
# define D3(p,a,b)
# define D4(p,a,b,c)
#endif

#endif /*__debug_h*/
