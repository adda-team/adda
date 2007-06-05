/* FILE: MatVec_nim.c
 * AUTH: Alfons Hoekstra
 * DESCR: this file contains the code for (MatVec_nim, MatVecAndInp_nim)
 *        and (MatVecHer_nim, MatVecAndInpHer_nim). Just interface for 
 *        calling function both_MatVec
 *
 *        Currently is developed by Maxim Yurkin
 */

#include <stdio.h>
#include "cmplx.h"

extern void couple_matrix ( );
extern void free_dCmatrix ();
extern doublecomplex **dCmatrix ( );

/*==========================================================*/
/********* USER INTERFACE TO ABOVE MATVEC FUNCTIONS *********/

void MatVec_nim (argvec, resultvec, nldip, ndip, procid, nproc,
		 wavenum)
     doublecomplex *argvec, *resultvec;
     double   wavenum;
     int      nldip, ndip, procid, nproc;
{
  both_MatVec(argvec, resultvec, NULL, nldip, ndip,
	      procid, nproc, wavenum,0);
}

/*==========================================================*/

void MatVecAndInp_nim (argvec, resultvec, inprod, nldip, ndip,
		       procid, nproc, wavenum)   
     doublecomplex *argvec, *resultvec;
     double   *inprod, wavenum;
     int      nldip, ndip, procid, nproc;
     
{
  both_MatVec(argvec, resultvec, inprod, nldip, ndip,
	      procid, nproc, wavenum,0);
}

/*==========================================================*/

void MatVecHer_nim (argvec, resultvec, nldip, ndip, procid, nproc,
		    wavenum)
     doublecomplex *argvec, *resultvec;
     double   wavenum;
     int      nldip, ndip, procid, nproc;
{
  both_MatVec(argvec, resultvec, NULL, nldip, ndip,
	      procid, nproc, wavenum,1);
}

/*==========================================================*/

void MatVecAndInpHer_nim (argvec, resultvec, inprod, nldip, ndip,
			  procid, nproc, wavenum)   
     doublecomplex *argvec, *resultvec;
     double   *inprod, wavenum;
     int      nldip, ndip, procid, nproc;
     
{
  both_MatVec(argvec, resultvec, inprod, nldip, ndip,
	      procid, nproc, wavenum,1);
}
