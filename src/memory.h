/* FILE: Romberg.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        memory allocation and freeing
 */
#ifndef __memory_h
#define __memory_h

#define MBYTE 1048576.0

doublecomplex **dCmatrix (int nrl, int nrh, int ncl, int nch);
void free_dCmatrix (doublecomplex **m, int nrl, int nrh, int ncl);
doublecomplex *dCvector(int size);
void free_dCvector (doublecomplex *v);

double **dmatrix (int nrl, int nrh, int ncl, int nch);
void free_dmatrix (double **m, int nrl, int nrh, int ncl);
double *dvector(int nl,int nh);
void free_dvector (double *v, int nl);

int *ivector(int nl, int nh);
void free_ivector (int *v, int nl);
int **imatrix (int nrl, int nrh, int ncl, int nch);
void free_imatrix (int **m, int nrl, int nrh, int ncl);

#endif /*__memory_h*/
