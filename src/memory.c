/* File: memory.c
 * $Date::                            $
 * Descr: allocation and freeing of different vectors and matrices; checks for 'out of memory';
 *        resistant to 0 sizes in allocation and NULL in freeing
 *
 * Copyright (C) 2006-2008,2010,2012-2013 ADDA contributors
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
#include "const.h" // keep this first
#include "memory.h" // corresponding header
// project headers
#include "fft.h"
#include "io.h"
#include "types.h"
// system headers
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef FFTW3
#	include <fftw3.h> // for fftw_malloc; types.h should be defined before (to match C99 complex type)
#endif

// common error check
#define MALLOC_ERROR LogError(ERR_LOC_CALL,"Could not malloc %s",name)
#define CHECK_NULL(size,v) if ((size)!=0 && (v)==NULL) MALLOC_ERROR
#define CHECK_SIZE(size,type) if ((SIZE_MAX/sizeof(type))<(size)) MALLOC_ERROR
#define IF_FREE(v) if((v)!=NULL) free(v)
#define OVERFLOW LogError(ERR_LOC_CALL,"Integer overflow in '%s'",name);

//======================================================================================================================

void CheckOverflow(const double size,OTHER_ARGUMENTS)
// checks if size can fit into size_t type, otherwise overflow will happen before memory allocation
{
	if (size>SIZE_MAX) OVERFLOW;
}

//======================================================================================================================

size_t MultOverflow(const size_t a,const size_t b,OTHER_ARGUMENTS)
// multiplies two integers and checks for overflow
{
	if ((SIZE_MAX/a)<b) OVERFLOW;
	return(a*b);
}

//======================================================================================================================
doublecomplex *complexVector(const size_t size,OTHER_ARGUMENTS)
// allocates complex vector
{
	doublecomplex * restrict v;

	CHECK_SIZE(size,doublecomplex);
#ifdef FFTW3
	v=(doublecomplex *)fftw_malloc(size*sizeof(doublecomplex));
#else
	v=(doublecomplex *)malloc(size*sizeof(doublecomplex));
#endif
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

double **doubleMatrix(const size_t rows,const size_t cols,OTHER_ARGUMENTS)
// allocates double matrix (rows x cols)
{
	register size_t i;
	double ** restrict m;

	CHECK_SIZE(rows,double *);
	CHECK_SIZE(cols,double);
	m=(double **)malloc(rows*sizeof(double *));
	CHECK_NULL(rows,m);
	for (i=0;i<rows;i++) {
		m[i]=(double *)malloc(cols*sizeof(double));
		CHECK_NULL(cols,m[i]);
	}
	return m;
}

//======================================================================================================================

double *doubleVector(const size_t size,OTHER_ARGUMENTS)
// allocates double vector
{
	double * restrict v;

	CHECK_SIZE(size,double);
	v=(double *)malloc(size*sizeof(double));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

double *doubleRealloc(double *ptr,const size_t size,OTHER_ARGUMENTS)
// reallocates double vector ptr to a larger size
{
	double *v;

	CHECK_SIZE(size,double);
	v=(double *)realloc(ptr,size*sizeof(double));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

double *doubleVector2(const size_t nl,const size_t nh,OTHER_ARGUMENTS)
// allocates double vector with indices from nl to nh; all arguments must be non-negative and nh>=nl
{
	double * restrict v;
	size_t size;

	if (nh<nl || nh-nl==SIZE_MAX) MALLOC_ERROR;
	else size=nh-nl+1;
	CHECK_SIZE(size,double);
	v=(double *)malloc(size*sizeof(double));
	CHECK_NULL(size,v);
	v-=nl;
	return v;
}

//======================================================================================================================

int **intMatrix(const size_t nrl,const size_t nrh,const size_t ncl,const size_t nch,OTHER_ARGUMENTS)
// allocates integer matrix with indices [nrl,nrh]x[ncl,nch]; all arguments must be non-negative; nrh>=nrl; nch>=ncl
{
	register size_t i;
	size_t rows,cols;
	int ** restrict m;

	if (nrh<nrl || nrh-nrl==SIZE_MAX) MALLOC_ERROR;
	else rows=nrh-nrl+1;
	if (nch<ncl || nch-ncl==SIZE_MAX) MALLOC_ERROR;
	else cols=nch-ncl+1;
	CHECK_SIZE(rows,int *);
	CHECK_SIZE(cols,int);
	m=(int **)malloc(rows*sizeof(int *));
	CHECK_NULL(rows,m);
	m-=nrl;
	for (i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc(cols*sizeof(int));
		CHECK_NULL(cols,m[i]);
		m[i]-=ncl;
	}
	return m;
}

//======================================================================================================================

int *intVector(const size_t size,OTHER_ARGUMENTS)
// allocates integer vector
{
	int * restrict v;

	CHECK_SIZE(size,int);
	v=(int *)malloc(size*sizeof(int));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

unsigned short *ushortVector(const size_t size,OTHER_ARGUMENTS)
// allocates unsigned short vector
{
	unsigned short * restrict v;

	CHECK_SIZE(size,short);
	v=(unsigned short *)malloc(size*sizeof(short));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

char *charVector(const size_t size,OTHER_ARGUMENTS)
// allocates unsigned char vector
{
	char * restrict v;

	CHECK_SIZE(size,char);
	v=(char *)malloc(size*sizeof(char));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

char *charRealloc(char *ptr,const size_t size,OTHER_ARGUMENTS)
// reallocates char vector ptr to a larger size
{
	char *v;

	CHECK_SIZE(size,char);
	v=(char *)realloc(ptr,size*sizeof(char));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

unsigned char *ucharVector(const size_t size,OTHER_ARGUMENTS)
// allocates unsigned char vector
{
	unsigned char * restrict v;

	CHECK_SIZE(size,char);
	v=(unsigned char *)malloc(size*sizeof(char));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

bool *boolVector(const size_t size,OTHER_ARGUMENTS)
// allocates bool vector
{
	bool * restrict v;

	CHECK_SIZE(size,bool);
	v=(bool *)malloc(size*sizeof(bool));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

size_t *sizetVector(const size_t size,OTHER_ARGUMENTS)
// allocates bool vector
{
	size_t * restrict v;

	CHECK_SIZE(size,size_t);
	v=(size_t *)malloc(size*sizeof(size_t));
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

void *voidVector(const size_t size,OTHER_ARGUMENTS)
// allocates void vector
{
	void * restrict v;

	v=malloc(size);
	CHECK_NULL(size,v);
	return v;
}

//======================================================================================================================

void Free_cVector (doublecomplex * restrict v)
// frees complex vector
{
#ifdef FFTW3
	if (v!=NULL) fftw_free(v);
#else
	IF_FREE(v);
#endif
}

//======================================================================================================================

void Free_dMatrix(double ** restrict m,const size_t rows)
// frees double matrix (rows x cols)
{
	register size_t i;

	for (i=0;i<rows;i++) IF_FREE(m[i]);
	IF_FREE(m);
}

//======================================================================================================================

void Free_dVector2(double * restrict v,const size_t nl)
// frees double vector with indices from nl; all arguments must be non-negative
{
	IF_FREE(v+nl);
}

//======================================================================================================================

void Free_iMatrix(int ** restrict m,const size_t nrl,const size_t nrh,const size_t ncl)
// frees integer matrix with indices [nrl,nrh]x[ncl,...]; all arguments must be non-negative
{
	register size_t i;

	for (i=nrh;i>=nrl;i--) IF_FREE(m[i]+ncl);
	IF_FREE(m+nrl);
}

//======================================================================================================================

void Free_general(void * restrict v)
// frees general vector; kept in a special function for future development
// !!! Must not be used for complex vectors - use Free_cVector instead !!!
{
	IF_FREE(v);
}
