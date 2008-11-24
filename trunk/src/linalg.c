/* FILE: linalg.c
 * AUTH: Maxim Yurkin
 * DESCR: Different linear algebra operations for use with iterative solvers
 *        Highly specialized for DDA
 *
 *        'const' can be used for many more function variables, however it doesn't
 *        work in combination with 'doublecomplex *' or more nested lists. That seems
 *        to be a principal limitation of C standard (some compilers may work, some produce
 *        warnings)
 *
 * Copyright (C) 2006-2008 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <string.h>
#include "vars.h"
#include "types.h"
#include "comm.h"
#include "linalg.h"

//============================================================

void nInit(doublecomplex *a)
// initialize vector a with null values
{
	size_t i;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;i++) a[i][RE]=a[i][IM]=0.0;
}

//============================================================

void nCopy(doublecomplex *a,doublecomplex *b)
// copy vector b to a
{
	memcpy(a,b,nlocalRows*sizeof(doublecomplex));
}

//============================================================

double nNorm2(doublecomplex *a,TIME_TYPE *timing)
// squared norm of a large vector a
{
	size_t i;
	double inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) inprod += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];
	// this function is not called inside the main iteration loop
	MyInnerProduct(&inprod,double_type,1,timing);
	return inprod;
}

//============================================================

void nDotProd(doublecomplex *a,doublecomplex *b,doublecomplex c,TIME_TYPE *timing)
// dot product of two large vectors; c=a.b
{
	size_t i;

	c[RE]=c[IM]=0.0;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		c[RE] += a[i][RE]*b[i][RE] + a[i][IM]*b[i][IM];
		c[IM] += a[i][IM]*b[i][RE] - a[i][RE]*b[i][IM];
	}
	MyInnerProduct(c,cmplx_type,1,timing);
}

//============================================================

void nDotProd_conj(doublecomplex *a,doublecomplex *b,doublecomplex c,TIME_TYPE *timing)
// conjugate dot product of two large vectors; c=a.b*=b.a*
{
	size_t i;

	c[RE]=c[IM]=0.0;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		c[RE] += a[i][RE]*b[i][RE] - a[i][IM]*b[i][IM];
		c[IM] += a[i][IM]*b[i][RE] + a[i][RE]*b[i][IM];
	}
	MyInnerProduct(c,cmplx_type,1,timing);
}

//============================================================

void nDotProdSelf_conj(doublecomplex *a,doublecomplex c,TIME_TYPE *timing)
// conjugate dot product of vector on itself; c=a.a*
{
	size_t i;

	c[RE]=c[IM]=0.0;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		c[RE]+=a[i][RE]*a[i][RE]-a[i][IM]*a[i][IM];
		c[IM]+=a[i][RE]*a[i][IM];
	}
	MyInnerProduct(c,cmplx_type,1,timing);
	c[IM]*=2;
}

//============================================================

void nDotProdSelf_conj_Norm2(doublecomplex *a,doublecomplex c,double *norm,TIME_TYPE *timing)
/* Computes both conjugate dot product of vector on itself (c=a.a*)
 * and its Hermitian squared norm=||a||^2
 */
{
	size_t i;
	double buf[3];

	buf[0]=buf[1]=buf[2]=0.0;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		buf[0]+=a[i][RE]*a[i][RE];
		buf[1]+=a[i][IM]*a[i][IM];
		buf[2]+=a[i][RE]*a[i][IM];
	}
	MyInnerProduct(buf,double_type,3,timing);
	*norm=buf[0]+buf[1];
	c[RE]=buf[0]-buf[1];
	c[IM]=2*buf[2];
}

//============================================================

void nIncrem110_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2)
// a=c1*a+c2*b+c
{
	size_t i;
	double tmp;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		tmp=a[i][RE];
		a[i][RE] = c1[RE]*a[i][RE] - c1[IM]*a[i][IM] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM] + c[i][RE];
		a[i][IM] = c1[RE]*a[i][IM] + c1[IM]*tmp + c2[RE]*b[i][IM] + c2[IM]*b[i][RE] + c[i][IM];
	}
}

//============================================================

void nIncrem011_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2)
// a+=c1*b+c2*c
{
	size_t i;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		a[i][RE] += c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c2[RE]*c[i][RE] - c2[IM]*c[i][IM];
		a[i][IM] += c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c2[RE]*c[i][IM] + c2[IM]*c[i][RE];
	}
}

//============================================================

void nIncrem111_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2,const doublecomplex c3)
// a=c1*a+c2*b+c3*c
{
	size_t i;
	double tmp;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		tmp=a[i][RE];
		a[i][RE] = c1[RE]*a[i][RE] - c1[IM]*a[i][IM] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM]
		         + c3[RE]*c[i][RE] - c3[IM]*c[i][IM];
		a[i][IM] = c1[RE]*a[i][IM] + c1[IM]*tmp + c2[RE]*b[i][IM] + c2[IM]*b[i][RE]
		         + c3[RE]*c[i][IM] + c3[IM]*c[i][RE];
	}
}

//============================================================

void nIncrem01(doublecomplex *a,doublecomplex *b,const double c,double *inprod,TIME_TYPE *timing)
// a=a+c*b, inprod=|a|^2
{
	size_t i;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] += c*b[i][RE]; // a+=c*b
			a[i][IM] += c*b[i][IM];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] += c*b[i][RE]; // a+=c*b
			a[i][IM] += c*b[i][IM];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nIncrem10(doublecomplex *a,doublecomplex *b,const double c,double *inprod,TIME_TYPE *timing)
// a=c*a+b, inprod=|a|^2
{
	size_t i;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] = c*a[i][RE] + b[i][RE]; // a=c*a+b
			a[i][IM] = c*a[i][IM] + b[i][IM];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] = c*a[i][RE] + b[i][RE]; // a=c*a+b
			a[i][IM] = c*a[i][IM] + b[i][IM];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nIncrem11_d_c(doublecomplex *a,doublecomplex *b,const double c1,const doublecomplex c2,
	double *inprod,TIME_TYPE *timing)
// a=c1*a+c2*b, inprod=|a|^2 , one constant is double, another - complex
{
	size_t i;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] = c1*a[i][RE] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM]; // a=c1*a+c2*b
			a[i][IM] = c1*a[i][IM] + c2[RE]*b[i][IM] + c2[IM]*b[i][RE];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] = c1*a[i][RE] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM]; // a=c1*a+c2*b
			a[i][IM] = c1*a[i][IM] + c2[RE]*b[i][IM] + c2[IM]*b[i][RE];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nIncrem01_cmplx(doublecomplex *a,doublecomplex *b,const doublecomplex c,double *inprod,
	TIME_TYPE *timing)
// a=a+c*b, inprod=|a|^2
{
	size_t i;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] += c[RE]*b[i][RE] - c[IM]*b[i][IM]; // a+=c*b
			a[i][IM] += c[RE]*b[i][IM] + c[IM]*b[i][RE];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] += c[RE]*b[i][RE] - c[IM]*b[i][IM]; // a+=c*b
			a[i][IM] += c[RE]*b[i][IM] + c[IM]*b[i][RE];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nIncrem10_cmplx(doublecomplex *a,doublecomplex *b,const doublecomplex c,double *inprod,
	TIME_TYPE *timing)
// a=c*a+b, inprod=|a|^2
{
	size_t i;
	double tmp;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			tmp=a[i][RE]; // a=c*a+b
			a[i][RE] = c[RE]*a[i][RE] - c[IM]*a[i][IM] + b[i][RE];
			a[i][IM] = c[RE]*a[i][IM] + c[IM]*tmp + b[i][IM];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			tmp=a[i][RE]; // a=c*a+b
			a[i][RE] = c[RE]*a[i][RE] - c[IM]*a[i][IM] + b[i][RE];
			a[i][IM] = c[RE]*a[i][IM] + c[IM]*tmp + b[i][IM];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nLinComb_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	const doublecomplex c2,double *inprod,TIME_TYPE *timing)
// a=c1*b+c2*c, inprod=|a|^2
{
	size_t i;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			// a=c1*b+c2*c
			a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c2[RE]*c[i][RE] - c2[IM]*c[i][IM];
			a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c2[RE]*c[i][IM] + c2[IM]*c[i][RE];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			// a=c1*b+c2*c
			a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c2[RE]*c[i][RE] - c2[IM]*c[i][IM];
			a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c2[RE]*c[i][IM] + c2[IM]*c[i][RE];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nLinComb1_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,const doublecomplex c1,
	double *inprod,TIME_TYPE *timing)
// a=c1*b+c, inprod=|a|^2
{
	size_t i;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			// a=c1*b+c
			a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c[i][RE];
			a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c[i][IM];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			// a=c1*b+c
			a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c[i][RE];
			a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c[i][IM];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nSubtr(doublecomplex *a,doublecomplex *b,doublecomplex *c,double *inprod,TIME_TYPE *timing)
// a=b-c, inprod=|a|^2
{
	size_t i;

	if (inprod==NULL) {
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] = b[i][RE] - c[i][RE]; // a=b-c
			a[i][IM] = b[i][IM] - c[i][IM];
		}
	}
	else {
		*inprod=0.0;
#pragma loop count (100000)
#pragma ivdep
		for (i=0;i<nlocalRows;++i) {
			a[i][RE] = b[i][RE] - c[i][RE]; // a=b-c
			a[i][IM] = b[i][IM] - c[i][IM];
			(*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM]; // *inprod=|a|^2
		}
		MyInnerProduct(inprod,double_type,1,timing);
	}
}

//============================================================

void nMult_cmplx(doublecomplex *a,doublecomplex *b,const doublecomplex c)
// multiply vector by a complex constant; a=c*b
{
	size_t i;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		a[i][RE] = c[RE]*b[i][RE] - c[IM]*b[i][IM]; // a[i]=c*b[i]
		a[i][IM] = c[RE]*b[i][IM] + c[IM]*b[i][RE];
	}
}

//============================================================

void nMultSelf_cmplx(doublecomplex *a,const doublecomplex c)
// multiply vector by a complex constant; a*=c
{
	size_t i;
	double tmp;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) {
		tmp=a[i][RE];
		a[i][RE] = c[RE]*a[i][RE] - c[IM]*a[i][IM]; // a[i]*=c
		a[i][IM] = c[RE]*a[i][IM] + c[IM]*tmp;
	}
}

//============================================================

void nMult_mat(doublecomplex *a,doublecomplex *b,doublecomplex c[][3])
// multiply by a function of material of a dipole and component; a[3*i+j]=c[mat[i]][j]*b[3*i+j]
{
	size_t i,k;
	int j;
	doublecomplex *val;

	k=0;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<local_nvoid_Ndip;++i) {
		val=c[material[i]];
		for (j=0;j<3;j++) {
			a[k][RE] = val[j][RE]*b[k][RE] - val[j][IM]*b[k][IM];
			a[k][IM] = val[j][RE]*b[k][IM] + val[j][IM]*b[k][RE];
			k++;
		}
	}
}

//============================================================

void nMultSelf_mat(doublecomplex *a,doublecomplex c[][3])
// multiply by a function of material of a dipole and component; a[3*i+j]*=c[mat[i]][j]
{
	size_t i,k;
	int j;
	double tmp;
	doublecomplex *val;

	k=0;
#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<local_nvoid_Ndip;++i) {
		val=c[material[i]];
		for (j=0;j<3;j++) {
			tmp=a[k][RE];
			a[k][RE] = val[j][RE]*a[k][RE] - val[j][IM]*a[k][IM];
			a[k][IM] = val[j][RE]*a[k][IM] + val[j][IM]*tmp;
			k++;
		}
	}
}

//============================================================

void nConj(doublecomplex *a)
// complex conjugate of the vector
{
	size_t i;

#pragma loop count (100000)
#pragma ivdep
	for (i=0;i<nlocalRows;++i) a[i][IM]=-a[i][IM];
}
