/* File: function.h
 * $Date::                            $
 * Descr: INLINE definition and function attributes
 *
 * Copyright (C) 2006,2008,2010-2011 ADDA contributors
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
#ifndef __function_h
#define __function_h

// specify to inline some functions; if there are problems with compiler change to "static"
#define INLINE static __inline

// attribute options for GCC compilers (Intel compiler may also recognize them)
#ifdef __GNUC__
	// sets a macro for testing GCC version (copied from _mingw.h)
#	ifdef __GNUC_MINOR__
#		define GCC_PREREQ(major,minor) \
			(__GNUC__ > (major) || (__GNUC__ == (major) && __GNUC_MINOR__ >= (minor)))
#	else
#		define GCC_PREREQ(major, minor)  0
#	endif
	// The following chooses between __printf__ and __gnu_printf__ attributes
#	if GCC_PREREQ(4,4)
		/* for newer gcc. In particular, it removes warnings about %z in printf-type functions on
		 * Windows. Native Windows libraries do not support this format specifier, however MinGW
		 * contains replacement functions to take care of it. So when using MinGW, system
		 * limitations are not that relevant.
		 * This causes certain problems with icc 12.0. It claims to support gcc 4.5 but does not
		 * recognize this format. However, this is dealt with disabling corresponding warning in
		 * Makefile.
		 */
#		define ATT_PRINTF(a,b) __attribute__ ((__format__(__gnu_printf__,a,b)))
#	else // for older gcc
#		define ATT_PRINTF(a,b) __attribute__ ((__format__(__printf__,a,b)))
#	endif
#	if GCC_PREREQ(3,0) // pure and malloc require gcc 3.0
#		define ATT_PURE        __attribute__ ((__pure__))
#		define ATT_MALLOC      __attribute__ ((__malloc__))
#	else
#		define ATT_PURE
#		define ATT_MALLOC
#	endif
#	define ATT_NORETURN    __attribute__ ((__noreturn__))
#	define ATT_UNUSED      __attribute__ ((__unused__))
#else
#	define ATT_PRINTF(a,b)
#	define ATT_PURE
#	define ATT_MALLOC
#	define ATT_NORETURN
#	define ATT_UNUSED
#endif

#endif // __function_h
