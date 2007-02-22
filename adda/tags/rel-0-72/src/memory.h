/* FILE: memory.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        memory allocation and freeing
 */
#ifndef __memory_h
#define __memory_h

#define MBYTE 1048576.0

doublecomplex **cMatrix (int nrl, int nrh, int ncl, int nch);
void Free_cMatrix (doublecomplex **m, int nrl, int nrh, int ncl);
doublecomplex *cVector(int size);
void Free_cVector (doublecomplex *v);

double **dMatrix (int nrl, int nrh, int ncl, int nch);
void Free_dMatrix (double **m, int nrl, int nrh, int ncl);
double *dVector(int nl,int nh);
void Free_dVector (double *v, int nl);

int *iVector(int nl, int nh);
void Free_iVector (int *v, int nl);
int **iMatrix (int nrl, int nrh, int ncl, int nch);
void Free_iMatrix (int **m, int nrl, int nrh, int ncl);

#endif /*__memory_h*/
