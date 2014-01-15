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
#	include "os.h" // for awareness of WINDOWS
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

// We require only version 2.0, but use the some extensions from 2.2, if available

/* While MPI 2.2 fully supports bool complex datatypes (including reduce operations), there is lack of the support of
 * reduction on Windows. The most advanced implementation (for which binaries are available) is MPICH 1.4.1p1 - it
 * doesn't support reduce on complex numbers, although the same Unix version does.
 * http://trac.mpich.org/projects/mpich/ticket/1525
 * Moreover, MPICH2 1.5 and further can't even be compiled on Windows. http://trac.mpich.org/projects/mpich/ticket/1557
 *
 * OpenMPI may also be affected by similar inconsistencies - see http://svn.boost.org/trac/ompi/ticket/3127
 */
/* whether there is support of reduction for advanced datatypes (like bool and complex):
 * At least 2.2 and check that MPICH2 version is not deficient (better than 1.4.1p1 on Windows)
 * TODO: test if that is enough for OpenMPI, or whether a version check should be added
 */
#define EXT_MPI_REDUCE ( ((MPI_VERSION>2) || ((MPI_VERSION==2) && (MPI_SUBVERSION>=2))) && \
		!( defined(MPICH2) && defined(WINDOWS) && (MPICH2_NUMVERSION<=10401301) ) )

#if defined(MPI_C_BOOL) && EXT_MPI_REDUCE
#	define EXT_MPI_22
#	define SUPPORT_MPI_BOOL
#	define mpi_bool MPI_C_BOOL
#else
/* this is not perfectly portable, but should work on most hardware. These datatypes do not need to be fully compatible
 * but should only have the same size and map 0 to 0 and !0 to !0.
 */
#	define mpi_bool MPI_SIGNED_CHAR
#endif

#ifdef MPI_C_DOUBLE_COMPLEX
#	define EXT_MPI_22
#	define SUPPORT_MPI_COMPLEX
#	if EXT_MPI_REDUCE
#		define SUPPORT_MPI_COMPLEX_REDUCE
#	endif
#endif

/* If any extensions are used at compile time, the runtime requirements are incremented to 2.2 or to the MPI version
 * used for compilation (whichever is smaller). So, if we are using the same MPI for compilation and runtime, this will
 * introduce no limitations. However, if runtime MPI is different from the one used for compilation, this will prevent
 * at least some of problems.
 */
#ifdef MPI_EXT_22
#	if (MPI_VERSION==2) && (MPI_SUBVERSION<2)
#		define RUN_MPI_SUBVER_REQ MPI_SUBVERSION
#	else
#		define RUN_MPI_SUBVER_REQ 2
#	endif
#else // if no extensions, the requirement is the same as during compilation
#	define RUN_MPI_SUBVER_REQ MPI_SUBVER_REQ
#endif

/* Hopefully MPI_SIZE_T will be defined in the future MPI versions. As of version 3.0 there is only MPI_AINT, which is
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
