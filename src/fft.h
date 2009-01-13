/* FILE: fft.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of fft parameters and routines
 */
#ifndef __fft_h
#define __fft_h

/*#define TEMPERTON        /* uncomment to use Temperton FFT      */

#ifndef TEMPERTON
#  define FFTW3          /* FFTW3 is default */
#endif

#ifdef FFTW3
# include <fftw3.h>
/* define level of planning for usual and Dmatrix (DM) FFT */
/* FFTW_ESTIMATE (heuristics), FFTW_MEASURE (def), FTW_PATIENT, or FFTW_EXHAUSTIVE */
# define PLAN_FFTW FFTW_MEASURE
# define PLAN_FFTW_DM FFTW_ESTIMATE
#endif

/* direction of FFT and transpose */
/* complies with FFTW3 definition */
#define FORWARD -1
#define BACKWARD 1

/* for transpose YZ */
#define BLOCK 64

void fftX(int isign);
void fftY(int isign);
void fftZ(int isign);
void transposeYZ(int direction);

#endif /*__fft_h*/
