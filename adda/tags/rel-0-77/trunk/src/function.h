/* FILE: function.h
 * AUTH: Maxim Yurkin
 * DESCR: INLINE definition and function attributes
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __function_h
#define __function_h

/* specify to inline some functions */
/* if problems with compiler change to "static" */
#define INLINE static __inline

/* some optimization shortcuts; can be easily removed if not supported by compiler */
#ifdef __GNUC__
# define ATT_NORETURN      __attribute__ ((__noreturn__))
# define ATT_PRINTF(a,b)   __attribute__ ((__format__(__printf__,a,b)))
# define ATT_UNUSED        __attribute__ ((__unused__))
# define ATT_PURE          __attribute__ ((__pure__))
# define ATT_MALLOC        __attribute__ ((__malloc__))
#else
# define ATT_NORETURN
# define ATT_PRINTF(a,b)
# define ATT_UNUSED
# define ATT_PURE
# define ATT_MALLOC
#endif

#endif /*__function_h*/
