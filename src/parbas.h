/* File: parbas.h
 * $Date::                            $
 * Descr: parallel basics; includes necessary headers and checks version of the standard.
 *
 * Copyright (C) 2007-2009 ADDA contributors
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
#ifndef __parbas_h
#define __parbas_h

#ifdef ADDA_MPI
#	include <mpi.h>
// define required version of MPI
#	define MPI_VER_REQ 2
#	define MPI_SUBVER_REQ 0
// check MPI version for conformity during compilation
#	if !defined(MPI_VERSION) || !defined(MPI_SUBVERSION)
#		error *** Can not determine MPI version, hence MPI is too old. ***
#	elif (MPI_VERSION<MPI_VER_REQ) || ((MPI_VERSION==MPI_VER_REQ) && (MPI_SUBVERSION<MPI_SUBVER_REQ))
#		error *** MPI version is too old. ***
#	endif
#endif

#endif // __parbas_h
