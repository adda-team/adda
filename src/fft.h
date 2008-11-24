/* FILE: fft.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of fft parameters and routines
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __fft_h
#define __fft_h

/* Temperton FFT is a simple one, its source code is supplied together with ADDA.
   The only inconvenience is that it is in fortran (not easily incorporated into a
   project under Windows (using any C/C++ developing tool) and should be
   compiled separately.
   FFTW3 requires separate installation of package from http://www.fftw.org, however
   it is hoghly optimized to the particular hardware and is generally significantly faster.
   Therfore, it is the default. */  

/*#define FFT_TEMPERTON        /* uncomment to use Temperton FFT      */

#ifndef FFT_TEMPERTON
#  define FFTW3          /* FFTW3 is default */
#endif

/* direction of FFT and transpose */
/* complies with FFTW3 definition */
#define FFT_FORWARD -1
#define FFT_BACKWARD 1

void fftX(int isign);
void fftY(int isign);
void fftZ(int isign);
void TransposeYZ(int direction);
extern void InitDmatrix(void);
extern void Free_FFT_Dmat(void);
int fftFit(int size, int div);
void CheckNprocs(void);

#endif /*__fft_h*/
