/* FILE: memory.c
 * AUTH: Alfons Hoekstra
 * DESCR: allocation and freeing of real and complex vectors and matrices
 *
 *        Currently is developed by Maxim Yurkin
 */

#include <stdio.h>
#include <stdlib.h>
#include "cmplx.h"
#include "memory.h"
#include "fft.h"

/*============================================================*/
/******************* Complex matrices and vectors ***************/

/* dCmatrix(nrl, nrh, ncl, nch)
 * Allocate a double complex matrix with range
 * [nrl..nrh][ncl..nch]
 */
doublecomplex **dCmatrix (int nrl, int nrh, int ncl, int nch)
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

/* free_dCmatrix(m, nrl, nrh, ncl)
 * Frees a double complex matrix allocated with dCmatrix
 */
void free_dCmatrix (doublecomplex **m, int nrl, int nrh, int ncl)
{ 
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}

/*============================================================*/

/* dCvector(size)
 *  return a pointer to space for a complex vector
 *  having indices ranging from 0 to size-1
 */
doublecomplex *dCvector(int size)
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

/* free_dCvector(m)
 * Frees a double complex vector allocated with dCvector
 */
void free_dCvector (doublecomplex *v)
{ 
#ifdef FFTW3
  fftw_free((void *)v);
#else  
  free(v);
#endif
}

/*============================================================*/
/********************* Double matrices and vectors *************/

/* dmatrix(nrl, nrh, ncl, nch)
 * Allocate a double matrix with range [nrl..nrh][ncl..nch]
 */
double **dmatrix (int nrl, int nrh, int ncl, int nch)
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

void free_dmatrix (double **m, int nrl, int nrh, int ncl)
{ 
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}

/*============================================================*/

/* dvector(nl, nh)
 *  return a pointer to space for a real vector in double precision
 *  having indices ranging from nl to nh
 */
double *dvector(int nl,int nh)
{
  double *v;
 
  v = (double *) malloc( (nh-nl+1) * sizeof(double) );
  if ( v == NULL )
    return( v );
  v -= nl;
  return( v );
}

/*============================================================*/

/* free_dvector(v, nl)
 * Frees a double vector allocated with dvector
 */
void free_dvector (double *v, int nl)
{
  free((char *) (v+nl));
}

/******************** Integer Vector and Matrix ***************/

/* ivector(nl, nh)
 *  return a pointer to space for a integer vector
 *  having indices ranging from nl to nh
 */
int *ivector(int nl, int nh)
{ 
  int *v;
  
  v = (int *) malloc( (nh-nl+1) * sizeof(int) );
  if ( v == NULL )
    return( v );
  v -= nl;
  return( v );
}

/*============================================================*/

/* free_ivector(m, nl)
 * Frees an integer vector allocated with ivector
 */
void free_ivector (int *v, int nl)
{ 
  free((char *) (v+nl));
}

/*============================================================*/

/* imatrix(nrl, nrh, ncl, nch)
 * Allocate an int matrix with range [nrl..nrh][ncl..nch]
 */
int **imatrix (int nrl, int nrh, int ncl, int nch)
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

void free_imatrix (int **m, int nrl, int nrh, int ncl)
{ 
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}
