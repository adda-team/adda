/* File: parbas.h
 * $Date::                            $
 * Descr: parallel basics; includes necessary headers and checks version of the standard.
 *
 * Copyright (C) 2007-2010,2013-2014 ADDA contributors
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

#include "const.h" // for GREATER_EQ2
#include "os.h"    // for awareness of WINDOWS
#include <mpi.h>
/* This minimum requirement (2.1) is based on functions MPI_Allgather(v), which are used for radiation forces and sparse
 * mode. Should not be a problem, since is supported by OpenMPI since 1.3, and MPICH2 since 1.1. If those functions are
 * removed, then even MPI 2.0 will do.
 */
#define MPI_VER_REQ 2
#define MPI_SUBVER_REQ 1
// check MPI version for conformity during compilation
#if !defined(MPI_VERSION) || !defined(MPI_SUBVERSION)
#	error "Can not determine MPI version, hence MPI is too old."
#else
#	define MPI_PREREQ(major,minor) GREATER_EQ2(MPI_VERSION,MPI_SUBVERSION,major,minor)
#	if !MPI_PREREQ(MPI_VER_REQ,MPI_SUBVER_REQ)
#		error "MPI version is too old."
#	endif
#endif

/* We use some extensions from 2.2, if available. Namely MPI_C_BOOL and MPI_C_DOUBLE_COMPLEX. Those types may be
 * available in earlier implementations, but we actually use them only when implementation declares itself conforming to
 * MPI 2.2 (to avoid some subtle problems). For OpenMPI this is since version 1.7.3. If this version is found, it is
 * further required during runtime.
 *
 * While MPI 2.2 fully supports bool & complex datatypes (including reduce operations), there is lack of the support of
 * reduction on Windows. The most advanced implementation (for which binaries are available) is MPICH2 1.4.1p1 - it
 * doesn't support reduce on complex numbers, although the same Unix version does.
 * http://trac.mpich.org/projects/mpich/ticket/1525 . Moreover, MPICH2 1.5 and further can't even be compiled on
 * Windows. http://trac.mpich.org/projects/mpich/ticket/1557
 *
 * Thus, we use a special test of MPICH2 deficiency further on.
 */
#define DEFICIENT_MPICH2 ( defined(MPICH2) && defined(WINDOWS) && (MPICH2_NUMVERSION<=10401301) )

// Runtime (library) requirements depend on the MPI used for compilation
#define RUN_MPI_VER_REQ MPI_VER_REQ

#if MPI_PREREQ(2,2)
#	define RUN_MPI_SUBVER_REQ 2
//	Complex is used either partly or fully (with reduce), bool is used only when fully supported.
#	define SUPPORT_MPI_COMPLEX
#	if !DEFICIENT_MPICH2
#		define SUPPORT_MPI_COMPLEX_REDUCE
#		define SUPPORT_MPI_BOOL
#	endif
#else
#	define RUN_MPI_SUBVER_REQ MPI_SUBVER_REQ
#endif

#ifdef SUPPORT_MPI_BOOL
#	define mpi_bool MPI_C_BOOL
#else
/* this is not perfectly portable, but should work on most hardware. These datatypes do not need to be fully compatible
 * but should only have the same size and map 0 to 0 and !0 to !0.
 */
#	define mpi_bool MPI_SIGNED_CHAR
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
