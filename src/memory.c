/* FILE: memory.c
 * AUTH: Alfons Hoekstra
 * DATE: 17-10-1990
 * HIST: see logfile ChangeLog.
 */

#include <stdio.h>
#include "cmplx.h"

static char Fheader[]="$Id$";

/******************* Complex matrices and vectors ***************/

/* dCmatrix(nrl, nrh, ncl, nch)
 * Allocate a double complex matrix with range
 * [nrl..nrh][ncl..nch]
 */
dcomplex **dCmatrix (int nrl, int nrh, int ncl, int nch)
{
  register int	i;
  dcomplex	**m;
  
  m = (dcomplex **) malloc ((nrh-nrl+1)*sizeof(dcomplex));
  if (m == NULL)
    return (m);
  m -= nrl;
  
  for (i=nrl; i<=nrh; i++) {
    m[i] = (dcomplex *) malloc ((nch-ncl+1)*sizeof(dcomplex));
    if (m[i] == NULL)
      return (NULL);
    m[i] -= ncl;
  }
  
  return (m);
}



/* free_dCmatrix(m, nrl, nrh, ncl)
 * Frees a double complex matrix allocated with dCmatrix
 */
void free_dCmatrix (dcomplex **m, int nrl, int nrh, int ncl)
{ 
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}




/* dCvector(nlm nh)
 *  return a pointer to space for a complex vector
 *  having indices ranging from nl to nh
 * - Original written in ANSI C by S.V. Schell, 1990 Sep 14
 */
dcomplex *dCvector(int nl,int nh)
{
  dcomplex *v;
 
  v = (dcomplex *) calloc( (nh-nl+1) * sizeof(dcomplex),1 );
  if ( v == NULL )
    return( v );
  v -= nl;
  return( v );
}




/********************* Double matrices and vectors *************/

/* dmatrix(nrl, nrh, ncl, nch)
 * Allocate a double matrix with range [nrl..nrh][ncl..nch]
 */
REAL **dmatrix (int nrl, int nrh, int ncl, int nch)
{
  register int	i;
  REAL **m;
  
  m = (REAL **) malloc ((nrh-nrl+1)*sizeof(REAL));
  if (m == NULL)
    return (m);
  m -= nrl;
  
  for (i=nrl; i<=nrh; i++) {
    m[i] = (REAL *) malloc ((nch-ncl+1)*sizeof(REAL));
    if (m[i] == NULL)
      return (NULL);
    m[i] -= ncl;
  }
  
  return (m);
}

void free_dmatrix (REAL **m, int nrl, int nrh, int ncl)
{ 
  register int	i;
  for (i=nrh; i>=nrl; i--) free((char *)(m[i]+ncl));
  free((char *)(m+nrl));
}


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

/* free_dvector(v, nl)
 * Frees a double vector allocated with dvector
 */
void free_dvector (double *v, int nl)
{
  free((char *) (v+nl));
}



/******************** Integer Vector **********************/


/* ivector(nl, nh)
 *  return a pointer to space for a real vector in double precision
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

char ram_memory(char *mem,int size)
/* make sure that memory addresses mem to mem+size are loaded from
 * virtual memory to ram. (this cost virtually no extra time if it is
 * already in ram)
 */
{
  int i;
  char bogus;
  
  for(i=0;i<size;i+=2048) bogus+=mem[i];
  
  return(bogus); /* return value to prevent compiler from removing it */
}


