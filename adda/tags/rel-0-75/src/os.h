/* FILE: portable.h
 * AUTH: Maxim Yurkin
 * DESCR: determines which operation system is used
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __os_h
#define __os_h

/* If not WINDOWS or POSIX is found, some parts of the program, such as Precise Timing and
   File Locking, will fail to compile */
#ifdef _WIN32
# define WINDOWS
# include <windows.h>
/* this list is not exhaustive. gcc always defines __POSIX__ on POSIX-compliant systems,
   however other compilers do not necessarily do the same. You may define it manually */
#elif defined(__POSIX__) || defined(unix) || defined (__unix) || defined (__unix__)
# define POSIX
#endif

#endif /*__os_h*/
