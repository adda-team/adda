/* FILE: linalg.c
 * AUTH: Maxim Yurkin
 * DESCR: Different linear algebra operations for use with iterative solvers
 *        Highly specialized for DDA
 */
#include <time.h>
#include "vars.h"
#include "cmplx.h"
#include "comm.h"
#include "linalg.h"

/*============================================================*/

void nInit(doublecomplex *a)
    /* initialize vector a with null values */
{
  size_t i;

  for (i=0;i<nlocalRows;i++) a[i][RE]=a[i][IM]=0.0;
}

/*============================================================*/

void nCopy(doublecomplex *a,doublecomplex *b)
    /* copy vector b to a */
{
  memcpy(a,b,nlocalRows*sizeof(doublecomplex));
}

/*============================================================*/

double nNorm2(doublecomplex *a)
    /* squared norm of a large vector a */
{
  size_t i;
  double inprod=0.0;

  for (i=0;i<nlocalRows;++i)
    inprod += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];
  MyInnerProduct(&inprod,double_type,1);
  return inprod;
}

/*============================================================*/

void nDotProd(doublecomplex *a,doublecomplex *b,doublecomplex c)
    /* dot product of two large vectors; c=a.b */
{
  size_t i;
  doublecomplex dotprod;
  clock_t tstart;

  c[RE]=c[IM]=0.0;
  for (i=0;i<nlocalRows;++i) {
    c[RE] += a[i][RE]*b[i][RE] + a[i][IM]*b[i][IM];
    c[IM] += a[i][IM]*b[i][RE] - a[i][RE]*b[i][IM];
  }
  tstart=clock();
  MyInnerProduct(c,cmplx_type,1);
  Timing_OneIterComm+=clock()-tstart;
}

/*============================================================*/

void nDotProd_conj(doublecomplex *a,doublecomplex *b,doublecomplex c)
    /* conjugate dot product of two large vectors; c=a.b*=b.a* */
{
  size_t i;
  doublecomplex dotprod;
  clock_t tstart;

  c[RE]=c[IM]=0.0;
  for (i=0;i<nlocalRows;++i) {
    c[RE] += a[i][RE]*b[i][RE] - a[i][IM]*b[i][IM];
    c[IM] += a[i][IM]*b[i][RE] + a[i][RE]*b[i][IM];
  }
  tstart=clock();
  MyInnerProduct(c,cmplx_type,1);
  Timing_OneIterComm+=clock()-tstart;
}

/*============================================================*/

void nDotProdSelf_conj(doublecomplex *a,doublecomplex c)
    /* conjugate dot product of vector on itself; c=a.a* */
{
  size_t i;
  clock_t tstart;

  c[RE]=c[IM]=0.0;
  for (i=0;i<nlocalRows;++i) {
    c[RE]+=a[i][RE]*a[i][RE]-a[i][IM]*a[i][IM];
    c[IM]+=a[i][RE]*a[i][IM];
  }
  tstart=clock();
  MyInnerProduct(c,cmplx_type,1);
  Timing_OneIterComm+=clock()-tstart;
  c[IM]*=2;
}

/*============================================================*/

void nDotProdSelf_conj_Norm2(doublecomplex *a,doublecomplex c,double *norm)
    /* Computes both conjugate dot product of vector on itself (c=a.a*)
         and its hermitian Norm squared norm=||a||^2 */
{
  size_t i;
  double buf[3];
  clock_t tstart;

  buf[0]=buf[1]=buf[2]=0.0;
  for (i=0;i<nlocalRows;++i) {
    buf[0]+=a[i][RE]*a[i][RE];
    buf[1]+=a[i][IM]*a[i][IM];
    buf[2]+=a[i][RE]*a[i][IM];
  }
  tstart=clock();
  MyInnerProduct(buf,double_type,3);
  Timing_OneIterComm+=clock()-tstart;
  *norm=buf[0]+buf[1];
  c[RE]=buf[0]-buf[1];
  c[IM]=2*buf[2];
}

/*============================================================*/

void nIncrem110_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2)
     /* a=c1*a+c2*b+c */
{
  size_t i;
  double tmp;

  for (i=0;i<nlocalRows;++i) {
    tmp=a[i][RE];
    a[i][RE] = c1[RE]*a[i][RE] - c1[IM]*a[i][IM] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM] + c[i][RE];
    a[i][IM] = c1[RE]*a[i][IM] + c1[IM]*tmp + c2[RE]*b[i][IM] + c2[IM]*b[i][RE] + c[i][IM];
  }
}

/*============================================================*/

void nIncrem011_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2)
     /* a+=c1*b+c2*c */
{
  size_t i;

  for (i=0;i<nlocalRows;++i) {
    a[i][RE] += c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c2[RE]*c[i][RE] - c2[IM]*c[i][IM];
    a[i][IM] += c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c2[RE]*c[i][IM] + c2[IM]*c[i][RE];
  }
}

/*============================================================*/

void nIncrem111_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2,doublecomplex c3)
     /* a=c1*a+c2*b+c3*c */
{
  size_t i;
  double tmp;

  for (i=0;i<nlocalRows;++i) {
    tmp=a[i][RE];
    a[i][RE] = c1[RE]*a[i][RE] - c1[IM]*a[i][IM] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM] +
               c3[RE]*c[i][RE] - c3[IM]*c[i][IM];
    a[i][IM] = c1[RE]*a[i][IM] + c1[IM]*tmp + c2[RE]*b[i][IM] + c2[IM]*b[i][RE] +
               c3[RE]*c[i][IM] + c3[IM]*c[i][RE];
  }
}

/*============================================================*/

void nIncrem01(doublecomplex *a,doublecomplex *b,double c,double *inprod)
     /* a=a+c*b, inprod=|a|^2 */
{
  size_t i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] += c*b[i][RE];   /* a+=c*b */
      a[i][IM] += c*b[i][IM];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] += c*b[i][RE];   /* a+=c*b */
      a[i][IM] += c*b[i][IM];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];   /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem10(doublecomplex *a,doublecomplex *b,double c,double *inprod)
     /* a=c*a+b, inprod=|a|^2 */
{
  size_t i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] = c*a[i][RE] + b[i][RE];  /* a=c*a+b */
      a[i][IM] = c*a[i][IM] + b[i][IM];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] = c*a[i][RE] + b[i][RE];  /* a=c*a+b */
      a[i][IM] = c*a[i][IM] + b[i][IM];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem11_d_c(doublecomplex *a,doublecomplex *b,double c1,doublecomplex c2,double *inprod)
     /* a=c1*a+c2*b, inprod=|a|^2 , one constant is double, another - complex */
{
  size_t i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] = c1*a[i][RE] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM];     /* a=c1*a+c2*b */
      a[i][IM] = c1*a[i][IM] + c2[RE]*b[i][IM] + c2[IM]*b[i][RE];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] = c1*a[i][RE] + c2[RE]*b[i][RE] - c2[IM]*b[i][IM];     /* a=c1*a+c2*b */
      a[i][IM] = c1*a[i][IM] + c2[RE]*b[i][IM] + c2[IM]*b[i][RE];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem01_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c,double *inprod)
     /* a=a+c*b, inprod=|a|^2 */
{
  size_t i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] += c[RE]*b[i][RE] - c[IM]*b[i][IM];     /* a+=c*b */
      a[i][IM] += c[RE]*b[i][IM] + c[IM]*b[i][RE];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] += c[RE]*b[i][RE] - c[IM]*b[i][IM];     /* a+=c*b */
      a[i][IM] += c[RE]*b[i][IM] + c[IM]*b[i][RE];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem10_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c,double *inprod)
     /* a=c*a+b, inprod=|a|^2 */
{
  size_t i;
  double tmp;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      tmp=a[i][RE];                                         /* a=c*a+b */
      a[i][RE] = c[RE]*a[i][RE] - c[IM]*a[i][IM] + b[i][RE];
      a[i][IM] = c[RE]*a[i][IM] + c[IM]*tmp + b[i][IM];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      tmp=a[i][RE];                                         /* a=c*a+b */
      a[i][RE] = c[RE]*a[i][RE] - c[IM]*a[i][IM] + b[i][RE];
      a[i][IM] = c[RE]*a[i][IM] + c[IM]*tmp + b[i][IM];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nLinComb_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                     doublecomplex c1, doublecomplex c2,double *inprod)
     /* a=c1*b+c2*c, inprod=|a|^2 */
{
  size_t i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c2*c */
      a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c2[RE]*c[i][RE] - c2[IM]*c[i][IM];
      a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c2[RE]*c[i][IM] + c2[IM]*c[i][RE];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c2*c */
      a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c2[RE]*c[i][RE] - c2[IM]*c[i][IM];
      a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c2[RE]*c[i][IM] + c2[IM]*c[i][RE];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nLinComb1_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                     doublecomplex c1,double *inprod)
     /* a=c1*b+c, inprod=|a|^2 */
{
  size_t i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c */
      a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c[i][RE];
      a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c[i][IM];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c */
      a[i][RE] = c1[RE]*b[i][RE] - c1[IM]*b[i][IM] + c[i][RE];
      a[i][IM] = c1[RE]*b[i][IM] + c1[IM]*b[i][RE] + c[i][IM];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nSubtr(doublecomplex *a,doublecomplex *b,doublecomplex *c,double *inprod)
     /* a=b-c, inprod=|a|^2 */
{
  size_t i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] = b[i][RE] - c[i][RE];     /* a=b-c */
      a[i][IM] = b[i][IM] - c[i][IM];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][RE] = b[i][RE] - c[i][RE];     /* a=b-c */
      a[i][IM] = b[i][IM] - c[i][IM];
      (*inprod) += a[i][RE]*a[i][RE] + a[i][IM]*a[i][IM];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    MyInnerProduct(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nMult_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c)
     /* multiply vector by a complex constant; a=c*b */
{
  size_t i;

  for (i=0;i<nlocalRows;++i) {
    a[i][RE] = c[RE]*b[i][RE] - c[IM]*b[i][IM];     /* a[i]=c*b[i] */
    a[i][IM] = c[RE]*b[i][IM] + c[IM]*b[i][RE];
  }
}

/*============================================================*/

void nMultSelf_cmplx(doublecomplex *a,doublecomplex c)
   /* multiply vector by a complex constant; a*=c */
{
  size_t i;
  double tmp;

  for (i=0;i<nlocalRows;++i) {
    tmp=a[i][RE];
    a[i][RE] = c[RE]*a[i][RE] - c[IM]*a[i][IM];     /* a[i]*=c */
    a[i][IM] = c[RE]*a[i][IM] + c[IM]*tmp;
  }
}

/*============================================================*/

void nMult_mat(doublecomplex *a,doublecomplex *b,doublecomplex c[][3])
   /* multiply by a function of material of a dipole and component;
      a[3*i+j]=c[mat[i]][j]*b[3*i+j] */
{
  int i,j,k;
  doublecomplex *val;

  k=0;
  for (i=0;i<local_nvoid_Ndip;++i) {
    val=c[material[i]];
    for (j=0;j<3;j++) {
      a[k][RE] = val[j][RE]*b[k][RE] - val[j][IM]*b[k][IM];
      a[k][IM] = val[j][RE]*b[k][IM] + val[j][IM]*b[k][RE];
      k++;
    }  
  }
}

/*============================================================*/

void nMultSelf_mat(doublecomplex *a,doublecomplex c[][3])
   /* multiply by a function of material of a dipole and component;
      a[3*i+j]*=c[mat[i]][j] */
{
  int i,j,k;
  double tmp;
  doublecomplex *val;

  k=0;
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

/*============================================================*/

void nConj(doublecomplex *a)
   /* complex conjugate of the vector */
{
  size_t i;

  for (i=0;i<nlocalRows;++i) a[i][IM]=-a[i][IM];
}

