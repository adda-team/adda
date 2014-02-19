/* File: fft.h
 * $Date::                            $
 * Descr: definitions of FFT parameters and routines; void in sparse mode
 *
 * Copyright (C) 2006,2008,2010-2013 ADDA contributors
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
#ifndef SPARSE

#ifndef __fft_h
#define __fft_h

#ifndef FFT_TEMPERTON
#	define FFTW3 // FFTW3 is default
#endif
#ifdef OPENCL
#	ifndef CLFFT_APPLE
#		define CLFFT_AMD // CLFFT_AMD is default for OPENCL
#	endif
#endif

// direction of FFT and transpose; complies with definitions of FFTW3, TempertonFFT, Apple and AMD OpenCL FFTs
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

#endif // __fft_h

#endif // !SPARSE
