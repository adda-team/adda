/* FILE: MatVec_nim.c
 * AUTH: Alfons Hoekstra
 * DATE: december 1990
 * HIST: some revisions done in september 1992
 * DESC: this file contains the code for (MatVec_nim, MatVecAndInp_nim)
 *       and (MatVecHer_nim, MatVecAndInpHer_nim).
 */

/* calculate local matrix vector product of decomposed interaction
 * matrix with rk or pk. The matrix elements are calculated when they
 * are needed.
 */

#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include "debug.h"
#include "cmplx.h"

#define PI 3.14159265359

extern void couple_matrix ( );
extern void free_dCmatrix ();
extern dcomplex **dCmatrix ( );

static char *Fname = "MatVec_nim";


/********* USER INTERFACE TO ABOVE MATVEC FUNCTIONS *********/

void MatVec_nim (argvec, resultvec, nldip, ndip, procid, nproc,
		 wavenum)
     dcomplex *argvec, *resultvec;
     double   wavenum;
     int      nldip, ndip, procid, nproc;
{
  both_MatVec(argvec, resultvec, NULL, nldip, ndip,
	      procid, nproc, wavenum,0);
}


void MatVecAndInp_nim (argvec, resultvec, inprod, nldip, ndip,
		       procid, nproc, wavenum)   
     dcomplex *argvec, *resultvec;
     double   *inprod, wavenum;
     int      nldip, ndip, procid, nproc;
     
{
  both_MatVec(argvec, resultvec, inprod, nldip, ndip,
	      procid, nproc, wavenum,0);
}

void MatVecHer_nim (argvec, resultvec, nldip, ndip, procid, nproc,
		    wavenum)
     dcomplex *argvec, *resultvec;
     double   wavenum;
     int      nldip, ndip, procid, nproc;
{
  both_MatVec(argvec, resultvec, NULL, nldip, ndip,
	      procid, nproc, wavenum,1);
}


void MatVecAndInpHer_nim (argvec, resultvec, inprod, nldip, ndip,
			  procid, nproc, wavenum)   
     dcomplex *argvec, *resultvec;
     double   *inprod, wavenum;
     int      nldip, ndip, procid, nproc;
     
{
  both_MatVec(argvec, resultvec, inprod, nldip, ndip,
	      procid, nproc, wavenum,1);
}



