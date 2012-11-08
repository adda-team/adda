/* File: linalg.c
 * $Date::                            $
 * Descr: different linear algebra operations for use with iterative solvers; highly specialized
 *
 *        'const' can be used for many more function variables, however it doesn't work in
 *        combination with 'doublecomplex *' or more nested lists. That seems to be a principal
 *        limitation of C standard (some compilers may work, some produce warnings).
 *        Common feature of many functions is accepting timing argument. If it is not NULL, it is
 *        incremented by the time used for communication.
 *
 * Copyright (C) 2006-2008,2010-2012 ADDA contributors
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
#include "linalg.h" // corresponding header
// project headers
#include "comm.h"
#include "types.h"
#include "vars.h"
// system headers
#include <string.h>

/* There are several optimization ideas used in this file:
 * 1) 'restrict' keyword tells the compiler that used doublecomplex array do not alias, however it
 * is unclear whether the compiler will also assume that they do not alias with input pointers to
 * double or single doublecomplex (also passed as a pointer to double). Moreover, these single
 * variables are the ones that are most often accessed inside the loop. So we define them as local
 * register variables, and connect them to the input values outside of the main loop.
 * 2) pragmas, indicating that the loop is supposed to be very large, are used. However, they are
 * probably understood only by the Intel compiler.
 * 3) If usage of some function has coinciding arguments, than a special function for such case is
 * created. In particular, this allows consistent usage of 'restrict' keyword almost for all
 * function arguments.
 * 4) Deeper optimizations, such as loop unrolling, are left to the compiler.
 */
//============================================================

void nInit(doublecomplex * restrict a)
// initialize vector a with null values
{
	register size_t i;
	register const size_t n=local_nRows;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) a[i][RE]=a[i][IM]=0;
}

//============================================================

void nCopy(doublecomplex * restrict a,doublecomplex * restrict b)
// copy vector b to a (a=b); !!! they must not alias !!!
{
	memcpy(a,b,local_nRows*sizeof(doublecomplex));
}

//============================================================

double nNorm2(doublecomplex * restrict a,TIME_TYPE *comm_timing)
// squared norm of a large vector a
{
	register size_t i;
	register const size_t n=local_nRows;
	register double sum=0;
	double inprod; // needed for MyInnerProduct, since register variable can't be dereferenced
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];
	// this function is not called inside the main iteration loop
	inprod=sum;
	MyInnerProduct(&inprod,double_type,1,comm_timing);
	return inprod;
}

//============================================================

void nDotProd(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex c,
	TIME_TYPE *comm_timing)
/* dot product of two large vectors; c=a.b; here the dot implies conjugation
 * !!! a and b must not alias !!! (to enforce use of function nNorm2 for coinciding arguments)
 */
{
	register size_t i;
	register const size_t n=local_nRows;
	register double cre=0,cim=0;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		cre += a[i][RE]*b[i][RE] + a[i][IM]*b[i][IM];
		cim += a[i][IM]*b[i][RE] - a[i][RE]*b[i][IM];
	}
	c[RE]=cre;
	c[IM]=cim;
	MyInnerProduct(c,cmplx_type,1,comm_timing);
}

//============================================================

void nDotProd_conj(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex c,
	TIME_TYPE *comm_timing)
/* conjugate dot product of two large vectors; c=a.b*=b.a*; here the dot implies conjugation
 * !!! a and b must not alias !!! (to enforce use of function nDotProdSelf_conj for coinciding
 * arguments)
 */
{
	register size_t i;
	register const size_t n=local_nRows;
	register double cre=0,cim=0;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		cre += a[i][RE]*b[i][RE] - a[i][IM]*b[i][IM];
		cim += a[i][IM]*b[i][RE] + a[i][RE]*b[i][IM];
	}
	c[RE]=cre;
	c[IM]=cim;
	MyInnerProduct(c,cmplx_type,1,comm_timing);
}

//============================================================

void nDotProdSelf_conj(doublecomplex * restrict a,doublecomplex c,TIME_TYPE *comm_timing)
// conjugate dot product of vector on itself; c=a.a*; here the dot implies conjugation
{
	register size_t i;
	register const size_t n=local_nRows;
	register double cre=0,cim=0;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		cre+=a[i][RE]*a[i][RE]-a[i][IM]*a[i][IM];
		cim+=a[i][RE]*a[i][IM];
	}
	c[RE]=cre;
	c[IM]=cim;
	MyInnerProduct(c,cmplx_type,1,comm_timing);
	c[IM]*=2;
}

//============================================================

void nDotProdSelf_conj_Norm2(doublecomplex * restrict a,doublecomplex c,double * restrict norm,
	TIME_TYPE *comm_timing)
/* Computes both conjugate dot product of vector on itself (c=a.a*)
 * and its Hermitian squared norm=||a||^2
 */
{
	register size_t i;
	register const size_t n=local_nRows;
	register double s0=0,s1=0,s2=0;
	double buf[3];

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		s0+=a[i][RE]*a[i][RE];
		s1+=a[i][IM]*a[i][IM];
		s2+=a[i][RE]*a[i][IM];
	}
	buf[0]=s0;
	buf[1]=s1;
	buf[2]=s2;
	MyInnerProduct(buf,double_type,3,comm_timing);
	*norm=buf[0]+buf[1];
	c[RE]=buf[0]-buf[1];
	c[IM]=2*buf[2];
}

//============================================================

void nIncrem110_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex * restrict c,const doublecomplex c1,const doublecomplex c2)
// a=c1*a+c2*b+c; !!! a,b,c must not alias !!!
{
	register size_t i;
	register const size_t n=local_nRows;
	register double tmp;
	register const double c1re=c1[RE],c1im=c1[IM],c2re=c2[RE],c2im=c2[IM];

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		tmp=a[i][RE];
		a[i][RE] = c1re*a[i][RE] - c1im*a[i][IM] + c2re*b[i][RE] - c2im*b[i][IM] + c[i][RE];
		a[i][IM] = c1re*a[i][IM] + c1im*tmp + c2re*b[i][IM] + c2im*b[i][RE] + c[i][IM];
	}
}

//============================================================

void nIncrem011_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex * restrict c,const doublecomplex c1,const doublecomplex c2)
// a+=c1*b+c2*c; !!! a,b,c must not alias !!!
{
	register size_t i;
	register const size_t n=local_nRows;
	register const double c1re=c1[RE],c1im=c1[IM],c2re=c2[RE],c2im=c2[IM];

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		a[i][RE] += c1re*b[i][RE] - c1im*b[i][IM] + c2re*c[i][RE] - c2im*c[i][IM];
		a[i][IM] += c1re*b[i][IM] + c1im*b[i][RE] + c2re*c[i][IM] + c2im*c[i][RE];
	}
}

//============================================================

void nIncrem110_d_c_conj(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex * restrict c,const double c1,const doublecomplex c2,double * restrict inprod,
	TIME_TYPE *comm_timing)
/* a=c1*a(*)+c2*b(*)+c; one constant is real, another - complex,
 * vectors a and c are conjugated during the evaluation; !!! a,b,c must not alias !!!
 */
{
	register const size_t n=local_nRows;
	register size_t i;
	// extra variable for c1 is probably redundant, but won't harm
	register const double cd=c1,cre=c2[RE],cim=c2[IM];
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = cd*a[i][RE] + cre*b[i][RE] + cim*b[i][IM] + c[i][RE]; // a=c1*a(*)+c2*b(*)+c
			a[i][IM] = -cd*a[i][IM] + cim*b[i][RE] - cre*b[i][IM] + c[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = cd*a[i][RE] + cre*b[i][RE] + cim*b[i][IM] + c[i][RE]; // a=c1*a(*)+c2*b(*)+c
			a[i][IM] = -cd*a[i][IM] + cim*b[i][RE] - cre*b[i][IM] + c[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nIncrem111_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex * restrict c,const doublecomplex c1,const doublecomplex c2,const doublecomplex c3)
// a=c1*a+c2*b+c3*c; !!! a,b,c must not alias !!!
{
	register size_t i;
	register const size_t n=local_nRows;
	register double tmp;
	register const double c1re=c1[RE],c1im=c1[IM],c2re=c2[RE],c2im=c2[IM],c3re=c3[RE],c3im=c3[IM];

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		tmp=a[i][RE];
		a[i][RE] = c1re*a[i][RE] - c1im*a[i][IM] + c2re*b[i][RE] - c2im*b[i][IM]
		         + c3re*c[i][RE] - c3im*c[i][IM];
		a[i][IM] = c1re*a[i][IM] + c1im*tmp + c2re*b[i][IM] + c2im*b[i][RE]
		         + c3re*c[i][IM] + c3im*c[i][RE];
	}
}

//============================================================

void nIncrem(doublecomplex * restrict a,doublecomplex * restrict b,double * restrict inprod,
	TIME_TYPE *comm_timing)
// a+=b, inprod=|a|^2; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] += b[i][RE]; // a+=b
			a[i][IM] += b[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] += b[i][RE]; // a+=b
			a[i][IM] += b[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nDecrem(doublecomplex * restrict a,doublecomplex * restrict b,double * restrict inprod,
	TIME_TYPE *comm_timing)
// a-=b, inprod=|a|^2; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] -= b[i][RE]; // a-=b
			a[i][IM] -= b[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] -= b[i][RE]; // a-=b
			a[i][IM] -= b[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nIncrem01(doublecomplex * restrict a,doublecomplex * restrict b,const double c,
	double * restrict inprod,TIME_TYPE *comm_timing)
// a=a+c*b, inprod=|a|^2; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cd=c; // extra variable for c is probably redundant, but won't harm
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] += cd*b[i][RE]; // a+=cd*b
			a[i][IM] += cd*b[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] += cd*b[i][RE]; // a+=cd*b
			a[i][IM] += cd*b[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nIncrem10(doublecomplex * restrict a,doublecomplex * restrict b,const double c,
	double * restrict inprod,TIME_TYPE *comm_timing)
// a=c*a+b, inprod=|a|^2; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cd=c; // extra variable for c is probably redundant, but won't harm
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = cd*a[i][RE] + b[i][RE]; // a=cd*a+b
			a[i][IM] = cd*a[i][IM] + b[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = cd*a[i][RE] + b[i][RE]; // a=cd*a+b
			a[i][IM] = cd*a[i][IM] + b[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nIncrem11_d_c(doublecomplex * restrict a,doublecomplex * restrict b,const double c1,
	const doublecomplex c2,double * restrict inprod,TIME_TYPE *comm_timing)
/* a=c1*a+c2*b, inprod=|a|^2 , one constant is double, another - complex;
 * !!! a and b must not alias !!!
 */
{
	register const size_t n=local_nRows;
	register size_t i;
	// extra variable for c1 is probably redundant, but won't harm
	register const double c1d=c1,c2re=c2[RE],c2im=c2[IM];
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = c1d*a[i][RE] + c2re*b[i][RE] - c2im*b[i][IM]; // a=c1d*a+c2*b
			a[i][IM] = c1d*a[i][IM] + c2re*b[i][IM] + c2im*b[i][RE];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = c1d*a[i][RE] + c2re*b[i][RE] - c2im*b[i][IM]; // a=c1d*a+c2*b
			a[i][IM] = c1d*a[i][IM] + c2re*b[i][IM] + c2im*b[i][RE];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nIncrem01_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,const doublecomplex c,
	double * restrict inprod,TIME_TYPE *comm_timing)
// a=a+c*b, inprod=|a|^2; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cre=c[RE],cim=c[IM];
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] += cre*b[i][RE] - cim*b[i][IM]; // a+=c*b
			a[i][IM] += cre*b[i][IM] + cim*b[i][RE];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] += cre*b[i][RE] - cim*b[i][IM]; // a+=c*b
			a[i][IM] += cre*b[i][IM] + cim*b[i][RE];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nIncrem10_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,const doublecomplex c,
	double * restrict inprod,TIME_TYPE *comm_timing)
// a=c*a+b, inprod=|a|^2; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cre=c[RE],cim=c[IM];
	register double sum=0,tmp;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			tmp=a[i][RE]; // a=c*a+b
			a[i][RE] = cre*a[i][RE] - cim*a[i][IM] + b[i][RE];
			a[i][IM] = cre*a[i][IM] + cim*tmp + b[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			tmp=a[i][RE]; // a=c*a+b
			a[i][RE] = cre*a[i][RE] - cim*a[i][IM] + b[i][RE];
			a[i][IM] = cre*a[i][IM] + cim*tmp + b[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nLinComb_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex * restrict c,const doublecomplex c1,const doublecomplex c2,
	double * restrict inprod,TIME_TYPE *comm_timing)
// a=c1*b+c2*c, inprod=|a|^2; !!! a,b,c must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double c1re=c1[RE],c1im=c1[IM],c2re=c2[RE],c2im=c2[IM];
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			// a=c1*b+c2*c
			a[i][RE] = c1re*b[i][RE] - c1im*b[i][IM] + c2re*c[i][RE] - c2im*c[i][IM];
			a[i][IM] = c1re*b[i][IM] + c1im*b[i][RE] + c2re*c[i][IM] + c2im*c[i][RE];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			// a=c1*b+c2*c
			a[i][RE] = c1re*b[i][RE] - c1im*b[i][IM] + c2re*c[i][RE] - c2im*c[i][IM];
			a[i][IM] = c1re*b[i][IM] + c1im*b[i][RE] + c2re*c[i][IM] + c2im*c[i][RE];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nLinComb1_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex * restrict c,const doublecomplex c1,double * restrict inprod,
	TIME_TYPE *comm_timing)
// a=c1*b+c, inprod=|a|^2; !!! a,b,c must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double c1re=c1[RE],c1im=c1[IM];
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			// a=c1*b+c
			a[i][RE] = c1re*b[i][RE] - c1im*b[i][IM] + c[i][RE];
			a[i][IM] = c1re*b[i][IM] + c1im*b[i][RE] + c[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			// a=c1*b+c
			a[i][RE] = c1re*b[i][RE] - c1im*b[i][IM] + c[i][RE];
			a[i][IM] = c1re*b[i][IM] + c1im*b[i][RE] + c[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nLinComb1_cmplx_conj(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex * restrict c,const doublecomplex c1,double * restrict inprod,
	TIME_TYPE *comm_timing)
// a=c1*b(*)+c, inprod=|a|^2; !!! a,b,c must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double c1re=c1[RE],c1im=c1[IM];
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = c1re*b[i][RE] + c1im*b[i][IM] + c[i][RE]; // a=c1*b(*)+c
			a[i][IM] = c1im*b[i][RE] - c1re*b[i][IM] + c[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = c1re*b[i][RE] + c1im*b[i][IM] + c[i][RE]; // a=c1*b(*)+c
			a[i][IM] = c1im*b[i][RE] - c1re*b[i][IM] + c[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}


//============================================================

void nSubtr(doublecomplex * restrict a,doublecomplex * restrict b,doublecomplex * restrict c,
	double *inprod,TIME_TYPE *comm_timing)
// a=b-c, inprod=|a|^2; !!! a,b,c must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register double sum=0;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = b[i][RE] - c[i][RE]; // a=b-c
			a[i][IM] = b[i][IM] - c[i][IM];
		}
	}
	else {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<n;i++) {
			a[i][RE] = b[i][RE] - c[i][RE]; // a=b-c
			a[i][IM] = b[i][IM] - c[i][IM];
			sum += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // sum=|a|^2
		}
		(*inprod)=sum;
		MyInnerProduct(inprod,double_type,1,comm_timing);
	}
}

//============================================================

void nMult(doublecomplex * restrict a,doublecomplex * restrict b,const double c)
// multiply vector by a real constant; a=c*b; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cd=c; // extra variable for c is probably redundant, but won't harm

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		a[i][RE] = cd*b[i][RE]; // a=c*b
		a[i][IM] = cd*b[i][IM];
	}
}

//============================================================

void nMult_cmplx(doublecomplex * restrict a,doublecomplex * restrict b,const doublecomplex c)
// multiply vector by a complex constant; a=c*b; !!! a and b must not alias !!!
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cre=c[RE],cim=c[IM];

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		a[i][RE] = cre*b[i][RE] - cim*b[i][IM]; // a=c*b
		a[i][IM] = cre*b[i][IM] + cim*b[i][RE];
	}
}

//============================================================

void nMultSelf(doublecomplex * restrict a,const double c)
// multiply vector by a real constant; a*=c
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cd=c; // extra variable for c is probably redundant, but won't harm

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		a[i][RE] *= cd; // a*=c
		a[i][IM] *= cd;
	}
}
//============================================================

void nMultSelf_conj(doublecomplex * restrict a,const double c)
// conjugate vector and multiply it by a real constant; a=c*a(*)
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cd=c; // extra variable for c is probably redundant, but won't harm

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		a[i][RE] *= cd; // a=c*a(*)
		a[i][IM] *= -cd;
	}
}

//============================================================

void nMultSelf_cmplx(doublecomplex * restrict a,const doublecomplex c)
// multiply vector by a complex constant; a*=c
{
	register const size_t n=local_nRows;
	register size_t i;
	register const double cre=c[RE],cim=c[IM];
	register double tmp;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) {
		tmp=a[i][RE];
		a[i][RE] = cre*a[i][RE] - cim*a[i][IM]; // a[i]*=c
		a[i][IM] = cre*a[i][IM] + cim*tmp;
	}
}

//============================================================

void nMult_mat(doublecomplex * restrict a,doublecomplex * restrict b,
	doublecomplex (* restrict c)[3])
/* multiply by a function of material of a dipole and component; a[3*i+j]=c[mat[i]][j]*b[3*i+j]
 * !!! a,b,c must not alias !!!
 */
{
	register const size_t nd=local_nvoid_Ndip; // name 'nd' to distinguish with 'n' used elsewhere
	register size_t i,k;
	register int j;
	/* Hopefully, the following declaration is enough to allow efficient loop unrolling. So the
	 * compiler should understand that none of the used vectors alias. Otherwise, deeper
	 * optimization should be used.
	 */
	doublecomplex * restrict val;

#pragma loop count (100000)
#pragma ivdep
	for (i=0,k=0;i<nd;i++) {
		val=c[material[i]];
		for (j=0;j<3;j++) { // we assume that compiler will unroll this loop
			a[k][RE] = val[j][RE]*b[k][RE] - val[j][IM]*b[k][IM];
			a[k][IM] = val[j][RE]*b[k][IM] + val[j][IM]*b[k][RE];
			k++;
		}
	}
}

//============================================================

void nMultSelf_mat(doublecomplex * restrict a,doublecomplex (* restrict c)[3])
/* multiply by a function of material of a dipole and component; a[3*i+j]*=c[mat[i]][j]
 * !!! a and c must not alias !!!
 */
{
	register const size_t nd=local_nvoid_Ndip; // name 'nd' to distinguish with 'n' used elsewhere
	register size_t i,k;
	register int j;
	register double tmp;
	/* Hopefully, the following declaration is enough to allow efficient loop unrolling. So the
	 * compiler should understand that none of the used vectors alias. Otherwise, deeper
	 * optimization should be used.
	 */
	doublecomplex * restrict val;

#pragma loop count (100000)
#pragma ivdep
	for (i=0,k=0;i<nd;i++) {
		val=c[material[i]];
		for (j=0;j<3;j++) { // we assume that compiler will unroll this loop
			tmp=a[k][RE];
			a[k][RE] = val[j][RE]*a[k][RE] - val[j][IM]*a[k][IM];
			a[k][IM] = val[j][RE]*a[k][IM] + val[j][IM]*tmp;
			k++;
		}
	}
}

//============================================================

void nConj(doublecomplex * restrict a)
// complex conjugate of the vector
{
	register const size_t n=local_nRows;
	register size_t i;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<n;i++) a[i][IM]=-a[i][IM];
}
