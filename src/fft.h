/* File: fft.h
 * $Date::                            $
 * Descr: definitions of FFT parameters and routines
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
#ifndef __fft_h
#define __fft_h

/* Temperton FFT is a simple one, its source code is supplied together with ADDA. The only
 * inconvenience is that it is in Fortran (not easily incorporated into a project under Windows
 * (using any C/C++ developing tool) and should be compiled separately.
 * FFTW3 requires separate installation of package from http://www.fftw.org, however it is highly
 * optimized to the particular hardware and is generally significantly faster. Therefore, it is the
 * default.
 */

//#define FFT_TEMPERTON // uncomment to use Temperton FFT

#ifndef FFT_TEMPERTON
#	define FFTW3 // FFTW3 is default
#elif defined(NO_FORTRAN)
// this is checked in Makefile, but additional check is always good
#	error Tempertron FFT is implemented in Fortran, hence is incompatible with NO_FORTRAN option
#endif

// direction of FFT and transpose; complies with FFTW3 definition
#define FFT_FORWARD -1
#define FFT_BACKWARD 1

void fftX(int isign);
void fftY(int isign);
void fftZ(int isign);
void TransposeYZ(int direction);
void InitDmatrix(void);
void Free_FFT_Dmat(void);
int fftFit(int size, int _div);
void CheckNprocs(void);

#ifdef OPENCL
void clfftX(const clFFT_Direction);
void clfftY(const clFFT_Direction);
void clfftZ(const clFFT_Direction);
#endif

#endif // __fft_h
