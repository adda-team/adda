/* File: linalg.h
 * $Author$
 * $Date::                            $
 * Descr: definitions for linear algebra operations on large vectors; see source (linalg.c) for
 *        details
 *
 * Copyright (C) 2006,2008 University of Amsterdam
 * Copyright (C) 2010 Institute of Chemical Kinetics and Combustion & University of Amsterdam
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
#ifndef __linalg_h
#define __linalg_h

#include "types.h"    // for doublecomplex
#include "function.h" // for function attributes
#include "timing.h"   // for TIME_TYPE

void nInit(doublecomplex *a);
void nCopy(doublecomplex *a,doublecomplex *b);
double nNorm2(doublecomplex *a,TIME_TYPE *comm_timing) ATT_PURE;
void nDotProd(doublecomplex *a,doublecomplex *b,doublecomplex c,TIME_TYPE *comm_timing);
void nDotProd_conj(doublecomplex *a,doublecomplex *b,doublecomplex c,TIME_TYPE *comm_timing);
void nDotProdSelf_conj(doublecomplex *a,doublecomplex c,TIME_TYPE *comm_timing);
void nDotProdSelf_conj_Norm2(doublecomplex *a,doublecomplex c,double *norm,TIME_TYPE *comm_timing);
void nIncrem110_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2);
void nIncrem011_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2);
void nIncrem111_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2,const doublecomplex c3);
void nIncrem01(doublecomplex *a,doublecomplex *b,const double c,double *inprod,
	TIME_TYPE *comm_timing);
void nIncrem10(doublecomplex *a,doublecomplex *b,const double c,double *inprod,
	TIME_TYPE *comm_timing);
void nIncrem11_d_c(doublecomplex *a,doublecomplex *b,const double c1,const doublecomplex c2,
	double *inprod,TIME_TYPE *comm_timing);
void nIncrem01_cmplx(doublecomplex *a,doublecomplex *b,const doublecomplex c,double *inprod,
	TIME_TYPE *comm_timing);
void nIncrem10_cmplx(doublecomplex *a,doublecomplex *b,const doublecomplex c,double *inprod,
	TIME_TYPE *comm_timing);
void nLinComb_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2,double *inprod,TIME_TYPE *comm_timing);
void nLinComb1_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	double *inprod,TIME_TYPE *comm_timing);
void nSubtr(doublecomplex *a,doublecomplex *b,doublecomplex *c,double *inprod,
	TIME_TYPE *comm_timing);
void nMult_cmplx(doublecomplex *a,doublecomplex *b,const doublecomplex c);
void nMultSelf_cmplx(doublecomplex *a,const doublecomplex c);
void nMult_mat(doublecomplex *a,doublecomplex *b,doublecomplex c[][3]);
void nMultSelf_mat(doublecomplex *a,doublecomplex c[][3]);
void nConj(doublecomplex *a);

#endif // __linalg_h
