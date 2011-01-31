/* File: function.h
 * $Date::                            $
 * Descr: INLINE definition and function attributes
 *
 * Copyright (C) 2006,2008,2010 ADDA contributors
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

// some optimization shortcuts; can be easily removed if not supported by compiler
#ifdef __GNUC__ // this should also work with Intel
#	define ATT_NORETURN    __attribute__ ((__noreturn__))
#	define ATT_PRINTF(a,b) __attribute__ ((__format__(__printf__,a,b)))
#	define ATT_UNUSED      __attribute__ ((__unused__))
#	define ATT_PURE        __attribute__ ((__pure__))
#	define ATT_MALLOC      __attribute__ ((__malloc__))
#else
#	define ATT_NORETURN
#	define ATT_PRINTF(a,b)
#	define ATT_UNUSED
#	define ATT_PURE
#	define ATT_MALLOC
#endif

#endif // __function_h
