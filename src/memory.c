/* FILE: memory.c
 * AUTH: Maxim Yurkin
 * DESCR: allocation and freeing of different vectors and matrices, checks for
 *        out of memory performed;
 *        Resistant to 0 sizes in allocation and NULL in freeing
 *
 *        Previous versions by Alfons Hoekstra
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <stdlib.h>
#include "types.h"
#include "memory.h"
#include "fft.h"
#include "io.h"
#include "const.h"

#ifdef FFTW3
# include <fftw3.h>   /* for fftw_malloc */
#endif

/* common error check */
#define CHECK_NULL(size,v) if ((size)!=0 && (v)==NULL) LogError(EC_ERROR,who,fname,line, \
                                                                "Could not malloc %s",name)
#define IF_FREE(v) if((v)!=NULL) free(v)

/*============================================================*/

doublecomplex *complexVector(const int size,OTHER_ARGUMENTS)
  /* allocates complex vector */
{
  doublecomplex *v;

#ifdef FFTW3
  v=(doublecomplex *)fftw_malloc(size*sizeof(doublecomplex));
#else
  v=(doublecomplex *)malloc(size*sizeof(doublecomplex));
#endif
  CHECK_NULL(size,v);
  return v;
}

/*============================================================*/

double **doubleMatrix(const int rows,const int cols,OTHER_ARGUMENTS)
  /* allocates double matrix (rows x cols) */
{
  register int i;
  double **m;

  m=(double **)malloc(rows*sizeof(double));
  CHECK_NULL(rows,m);
  for (i=0;i<rows;i++) {
    m[i]=(double *)malloc(cols*sizeof(double));
    CHECK_NULL(cols,m[i]);
  }
  return m;
}

/*============================================================*/

double *doubleVector(const int size,OTHER_ARGUMENTS)
   /* allocates double vector */
{
  double *v;

  v=(double *)malloc(size*sizeof(double));
  CHECK_NULL(size,v);
  return v;
}

/*============================================================*/

double *doubleVector2(const int nl,const int nh,OTHER_ARGUMENTS)
  /* allocates double vector with indices from nl to nh */
{
  double *v;

  v=(double *)malloc((nh-nl+1)*sizeof(double));
  CHECK_NULL(nh-nl+1,v);
  v-=nl;
  return v;
}

/*============================================================*/

int **intMatrix(const int nrl,const int nrh,const int ncl,const int nch,OTHER_ARGUMENTS)
  /* allocates integer matrix with indices [nrl,nrh]x[ncl,nch] */
{
  register int i;
  int **m;

  m=(int **)malloc((nrh-nrl+1)*sizeof(int));
  CHECK_NULL(nrh-nrl+1,m);
  m-=nrl;
  for (i=nrl;i<=nrh;i++) {
    m[i]=(int *)malloc((nch-ncl+1)*sizeof(int));
    CHECK_NULL(nch-ncl+1,m[i]);
    m[i]-=ncl;
  }
  return m;
}

/*============================================================*/

int *intVector(const int size,OTHER_ARGUMENTS)
  /* allocates integer vector */
{
  int *v;

  v=(int *)malloc(size*sizeof(int));
  CHECK_NULL(size,v);
  return v;
}

/*============================================================*/

unsigned short *ushortVector(const int size,OTHER_ARGUMENTS)
  /* allocates unsigned short vector */
{
  unsigned short *v;

  v=(unsigned short *)malloc(size*sizeof(short));
  CHECK_NULL(size,v);
  return v;
}

/*============================================================*/

char *charVector(const int size,OTHER_ARGUMENTS)
  /* allocates unsigned char vector */
{
  char *v;

  v=(char *)malloc(size*sizeof(char));
  CHECK_NULL(size,v);
  return v;
}

/*============================================================*/

unsigned char *ucharVector(const int size,OTHER_ARGUMENTS)
  /* allocates unsigned char vector */
{
  unsigned char *v;

  v=(unsigned char *)malloc(size*sizeof(char));
  CHECK_NULL(size,v);
  return v;
}

/*============================================================*/

void *voidVector(const int size,OTHER_ARGUMENTS)
  /* allocates void vector */
{
  void *v;

  v=malloc(size);
  CHECK_NULL(size,v);
  return v;
}

/*============================================================*/

void Free_cVector (doublecomplex *v)
   /* frees complex vector */
{
#ifdef FFTW3
  if (v!=NULL) fftw_free(v);
#else
  IF_FREE(v);
#endif
}

/*============================================================*/

void Free_dMatrix(double **m,const int rows)
  /* frees double matrix (rows x cols) */
{
  register int i;

  for (i=0;i<rows;i++) IF_FREE(m[i]);
  IF_FREE(m);
}

/*============================================================*/

void Free_dVector2(double *v,const int nl)
  /* frees double vector with indices from nl */
{
  IF_FREE(v+nl);
}

/*============================================================*/

void Free_iMatrix(int **m,const int nrl,const int nrh,const int ncl)
   /* frees integer matrix with indices [nrl,nrh]x[ncl,...] */
{
  register int i;

  for (i=nrh;i>=nrl;i--) IF_FREE(m[i]+ncl);
  IF_FREE(m+nrl);
}

/*============================================================*/

void Free_general(void *v)
  /* frees general vector; kept in a special function for future development */
{
  IF_FREE(v);
}

