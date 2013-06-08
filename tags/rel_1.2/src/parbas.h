/* File: parbas.h
 * $Date::                            $
 * Descr: parallel basics; includes necessary headers and checks version of the standard.
 *
 * Copyright (C) 2007-2010,2013 ADDA contributors
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
#ifndef __parbas_h
#define __parbas_h

#ifdef ADDA_MPI
#	include <mpi.h>
// define required version of MPI
#	define MPI_VER_REQ 2
#	define MPI_SUBVER_REQ 0
// check MPI version for conformity during compilation
#	if !defined(MPI_VERSION) || !defined(MPI_SUBVERSION)
#		error "Can not determine MPI version, hence MPI is too old."
#	elif (MPI_VERSION<MPI_VER_REQ) || ((MPI_VERSION==MPI_VER_REQ) && (MPI_SUBVERSION<MPI_SUBVER_REQ))
#		error "MPI version is too old."
#	endif

/* Hopefully MPI_SIZE_T will be defined in the future MPI versions. As of version 2.2 there is only MPI_AINT, which is
 * by design similar to size_t. However, there are no guarantees of the precise correspondence. The code below is lame,
 * but should work almost always. Another alternative is to use overloaded functions, but we want to stick with C.
 */
#ifndef MPI_SIZE_T
#	include <limits.h>
#	include <stdint.h>
#	if (SIZE_MAX == UINT_MAX)
#		define MPI_SIZE_T MPI_UNSIGNED
#	elif (SIZE_MAX == ULONG_MAX)
#		define MPI_SIZE_T MPI_UNSIGNED_LONG
#	elif (SIZE_MAX == ULLONG_MAX)
#		define MPI_SIZE_T MPI_UNSIGNED_LONG_LONG
#	else
#		error "No MPI alternative for size_t. Create an issue at http://code.google.com/p/a-dda/issues/"
#	endif
#endif

#endif // ADDA_MPI

#endif // __parbas_h
