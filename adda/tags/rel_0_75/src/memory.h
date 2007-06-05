/* FILE: memory.h
 * AUTH: Maxim Yurkin
 * DESCR: definitions of functions for
 *        memory allocation and freeing
 *
 * Copyright (C) 2006 University of Amsterdam
 * This code is covered by the GNU General Public License.
 */
#ifndef __memory_h
#define __memory_h

#define MBYTE 1048576.0

doublecomplex **cMatrix (int nrl, int nrh, int ncl, int nch) ATT_MALLOC;
void Free_cMatrix (doublecomplex **m, int nrl, int nrh, int ncl);
doublecomplex *cVector(int size) ATT_MALLOC;
void Free_cVector (doublecomplex *v);

double **dMatrix (int nrl, int nrh, int ncl, int nch) ATT_MALLOC;
void Free_dMatrix (double **m, int nrl, int nrh, int ncl);
double *dVector(int nl,int nh) ATT_MALLOC;
void Free_dVector (double *v, int nl);

int *iVector(int nl, int nh) ATT_MALLOC;
void Free_iVector (int *v, int nl);
int **iMatrix (int nrl, int nrh, int ncl, int nch) ATT_MALLOC;
void Free_iMatrix (int **m, int nrl, int nrh, int ncl);

#endif /*__memory_h*/
