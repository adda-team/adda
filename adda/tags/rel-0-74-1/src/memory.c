/* FILE: memory.c
 * AUTH: Alfons Hoekstra
 * DESCR: allocation and freeing of real, complex, and integer
 *        vectors and matrices
 *
 *        Currently is developed by Maxim Yurkin
 *
 * Copyright (C) 2006 M.A. Yurkin and A.G. Hoekstra
 * This code is covered by the GNU General Public License.
 */
#include <stdio.h>
#include <stdlib.h>
#include "cmplx.h"
#include "memory.h"
#include "fft.h"

#ifdef FFTW3
# include <fftw3.h>   /* for fftw_malloc */
#endif

/*============================================================*/
/******************* Complex matrices and vectors ***************/

/* cMatrix(nrl, nrh, ncl, nch)
 * Allocate a double complex matrix with range
 * [nrl..nrh][ncl..nch]
 */
doublecomplex **cMatrix (const int nrl,const int nrh,const int ncl,const int nch)
{
  register int	i;
  doublecomplex	**m;
  
  m = (doublecomplex **) malloc ((nrh-nrl+1)*sizeof(doublecomplex));
  if (m == NULL)
    return (m);
  m -= nrl;
  
  for (i=nrl; i<=nrh; i++) {
    m[i] = (doublecomplex *) malloc ((nch-ncl+1)*sizeof(doublecomplex));
    if (m[i] == NULL)
      return (NULL);
    m[i] -= ncl;
  }
  
  return (m);
}

/*============================================================*/

/* Free_cMatrix(m, nrl, nrh, ncl)
 * Frees a double complex matrix allocated with cMatrix
 */
void Free_cMatrix (doublecomplex **m,const int nrl,const int nrh,const int ncl)
{ 
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}

/*============================================================*/

/* cVector(size)
 *  return a pointer to space for a complex vector
 *  having indices ranging from 0 to size-1
 */
doublecomplex *cVector(const int size)
{
  doublecomplex *v;

#ifdef FFTW3
  v = (doublecomplex *) fftw_malloc(size*sizeof(doublecomplex));
#else
  v = (doublecomplex *) malloc(size*sizeof(doublecomplex));
#endif  
  return v;
}

/*============================================================*/

/* Free_cVector(m)
 * Frees a double complex vector allocated with cVector
 */
void Free_cVector (doublecomplex *v)
{ 
#ifdef FFTW3
  fftw_free((void *)v);
#else  
  free(v);
#endif
}

/*============================================================*/
/********************* Double matrices and vectors *************/

/* dMatrix(nrl, nrh, ncl, nch)
 * Allocate a double matrix with range [nrl..nrh][ncl..nch]
 */
double **dMatrix (const int nrl,const int nrh,const int ncl,const int nch)
{
  register int	i;
  double **m;
  
  m = (double **) malloc ((nrh-nrl+1)*sizeof(double));
  if (m == NULL)
    return (m);
  m -= nrl;
  
  for (i=nrl; i<=nrh; i++) {
    m[i] = (double *) malloc ((nch-ncl+1)*sizeof(double));
    if (m[i] == NULL)
      return (NULL);
    m[i] -= ncl;
  }
  
  return (m);
}

/*============================================================*/

void Free_dMatrix (double **m,const int nrl,const int nrh,const int ncl)
{
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}

/*============================================================*/

/* dVector(nl, nh)
 *  return a pointer to space for a real vector in double precision
 *  having indices ranging from nl to nh
 */
double *dVector(const int nl,const int nh)
{
  double *v;
 
  v = (double *) malloc( (nh-nl+1) * sizeof(double) );
  if ( v == NULL )
    return( v );
  v -= nl;
  return( v );
}

/*============================================================*/

/* Free_dVector(v, nl)
 * Frees a double vector allocated with dVector
 */
void Free_dVector (double *v,const int nl)
{
  free((char *) (v+nl));
}

/******************** Integer Vector and Matrix ***************/

/* iVector(nl, nh)
 *  return a pointer to space for a integer vector
 *  having indices ranging from nl to nh
 */
int *iVector(const int nl,const  int nh)
{ 
  int *v;
  
  v = (int *) malloc( (nh-nl+1) * sizeof(int) );
  if ( v == NULL )
    return( v );
  v -= nl;
  return( v );
}

/*============================================================*/

/* Free_iVector(m, nl)
 * Frees an integer vector allocated with iVector
 */
void Free_iVector (int *v,const int nl)
{ 
  free((char *) (v+nl));
}

/*============================================================*/

/* iMatrix(nrl, nrh, ncl, nch)
 * Allocate an int matrix with range [nrl..nrh][ncl..nch]
 */
int **iMatrix (const int nrl,const int nrh,const int ncl,const int nch)
{
  register int	i;
  int **m;
  
  m = (int **) malloc ((nrh-nrl+1)*sizeof(int));
  if (m == NULL)
    return (m);
  m -= nrl;
  
  for (i=nrl; i<=nrh; i++) {
    m[i] = (int *) malloc ((nch-ncl+1)*sizeof(int));
    if (m[i] == NULL)
      return (NULL);
    m[i] -= ncl;
  }
  
  return (m);
}

/*============================================================*/

void Free_iMatrix (int **m,const int nrl,const int nrh,const int ncl)
{ 
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}
