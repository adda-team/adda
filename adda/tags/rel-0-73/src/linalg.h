/* FILE: linalg.h
 * AUTH: Maxim Yurkin
 * DESCR: Definitions for linear algebra operations on large vectors
 *        see source (linalg.c) for description
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#ifndef __linalg_h
#define __linalg_h

void nInit(doublecomplex *a);
void nCopy(doublecomplex *a,doublecomplex *b);
double nNorm2(doublecomplex *a);
void nDotProd(doublecomplex *a,doublecomplex *b,doublecomplex c);
void nDotProd_conj(doublecomplex *a,doublecomplex *b,doublecomplex c);
void nDotProdSelf_conj(doublecomplex *a,doublecomplex c);
void nDotProdSelf_conj_Norm2(doublecomplex *a,doublecomplex c,double *norm);
void nIncrem110_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2);
void nIncrem011_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2);
void nIncrem111_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2,doublecomplex c3);
void nIncrem01(doublecomplex *a,doublecomplex *b,double c,double *inprod);
void nIncrem10(doublecomplex *a,doublecomplex *b,double c,double *inprod);
void nIncrem11_d_c(doublecomplex *a,doublecomplex *b,double c1,doublecomplex c2,double *inprod);
void nIncrem01_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c,double *inprod);
void nIncrem10_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c,double *inprod);
void nLinComb_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                     doublecomplex c1, doublecomplex c2,double *inprod);
void nLinComb1_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                     doublecomplex c1,double *inprod);
void nSubtr(doublecomplex *a, doublecomplex*b,doublecomplex *c,double *inprod);
void nMult_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c);
void nMultSelf_cmplx(doublecomplex *a,doublecomplex c);
void nMult_mat(doublecomplex *a,doublecomplex *b,doublecomplex c[][3]);
void nMultSelf_mat(doublecomplex *a,doublecomplex c[][3]);
void nConj(doublecomplex *a);

#endif /* __linalg_h */
