/* FILE: debug.h
 * AUTH: vesseur
 * DESCR: Definitions for debug functions
 *
 *        Currently is developed by Maxim Yurkin
 */
#ifndef __debug_h
#define __debug_h

/* Debugging implies turning on additional information messages
   during the code execution. A simple and convenient tool to
   generate such messages is used. */

/*#define DEBUG     /* uncomment to degug */

#ifdef DEBUG
#  define D(p)        DebugPrintf(__FILE__, __LINE__, p)
#  define D2(p,a)     DebugPrintf(__FILE__, __LINE__, p,a)
#  define D3(p,a,b)   DebugPrintf(__FILE__, __LINE__, p,a,b)
#  define D4(p,a,b,c) DebugPrintf(__FILE__, __LINE__, p,a,b,c)

void DebugPrintf(char *fname, int line, char *fmt, ...);

#else
#  define D(p)
#  define D2(p,a)
#  define D3(p,a,b)
#  define D4(p,a,b,c)
#endif

#endif /*__debug_h*/
