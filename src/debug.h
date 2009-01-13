/* FILE: debug.h
 * AUTH: vesseur
 * DESCR: Definitions for debug functions
 *
 *        Currently is developed by Maxim Yurkin
 */
#ifndef __debug_h
#define __debug_h

/*#define DEBUG     /* uncomment to degug */

#ifndef DEBUG
#  ifdef DEBUG_ALL
#    define DEBUG
#  endif
#endif

#ifdef DEBUG
#  define D(p)        debug_printf(__FILE__, __LINE__, p)
#  define D2(p,a)     debug_printf(__FILE__, __LINE__, p,a)
#  define D3(p,a,b)   debug_printf(__FILE__, __LINE__, p,a,b)
#  define D4(p,a,b,c) debug_printf(__FILE__, __LINE__, p,a,b,c)

void debug_printf(char *fname, int line, char *fmt, ...);

#else
#  define D(p)
#  define D2(p,a)
#  define D3(p,a,b)
#  define D4(p,a,b,c)
#endif

#endif /*__debug_h*/
