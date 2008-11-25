/* File: parbas.h
 * $Author: yurkin $
 * $Date:: 2008-11-25 10:26:04 +0600 #$
 * Descr: parallel basics; includes necessary headers and checks version of the standard.
 *
 * Copyright (C) 2007,2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __parbas_h
#define __parbas_h

#ifdef MPI
#	include <mpi.h>
// define required version of MPI
#	define MPI_VER_REQ 1
#	define MPI_SUBVER_REQ 2
// check MPI version for conformity during compilation
#	if !defined(MPI_VERSION) || !defined(MPI_SUBVERSION)
#		error *** Can not determine MPI version, hence MPI is too old. ***
#	elif (MPI_VERSION<MPI_VER_REQ) || ((MPI_VERSION==MPI_VER_REQ) && (MPI_SUBVERSION<MPI_SUBVER_REQ))
#		error *** MPI version is too old. ***
#	endif
#endif

#endif // __parbas_h
