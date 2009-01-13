/* FILE: linalg.c
 * AUTH: Maxim Yurkin
 * DESCR: Different linear algebra operations for use with iterative solvers
 *        Highly specialized for DDA
 */
#include "time.h"
#include "cmplx.h"
#include "comm.h"
#include "linalg.h"

extern clock_t Timing_OneIterComm;
extern int nlocalRows;
extern char *material;

/*============================================================*/

void nInit(doublecomplex *a)
    /* initialize vector a with null values */
{
  int i;
  
  for (i=0;i<nlocalRows;i++) a[i][re]=a[i][im]=0.0;
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
  int i;
  double inprod=0.0;

  for (i=0;i<nlocalRows;++i)
    inprod += a[i][re]*a[i][re] + a[i][im]*a[i][im];
  my_inner_product(&inprod,double_type,1);
  return inprod;
}

/*============================================================*/

void nDotProd(doublecomplex *a,doublecomplex *b,doublecomplex c)
    /* dot product of two large vectors; c=a.b */
{
  int i;
  doublecomplex dotprod;
  clock_t tstart;

  c[re]=c[im]=0.0;
  for (i=0;i<nlocalRows;++i) {
    c[re] += a[i][re]*b[i][re] + a[i][im]*b[i][im];
    c[im] += a[i][im]*b[i][re] - a[i][re]*b[i][im];
  }
  tstart=clock();
  my_inner_product(c,cmplx_type,1);
  Timing_OneIterComm+=clock()-tstart;
}

/*============================================================*/

void nDotProd_conj(doublecomplex *a,doublecomplex *b,doublecomplex c)
    /* conjugate dot product of two large vectors; c=a.b*=b.a* */
{
  int i;
  doublecomplex dotprod;
  clock_t tstart;

  c[re]=c[im]=0.0;
  for (i=0;i<nlocalRows;++i) {
    c[re] += a[i][re]*b[i][re] - a[i][im]*b[i][im];
    c[im] += a[i][im]*b[i][re] + a[i][re]*b[i][im];
  }
  tstart=clock();
  my_inner_product(c,cmplx_type,1);
  Timing_OneIterComm+=clock()-tstart;
}

/*============================================================*/

void nDotProdSelf_conj(doublecomplex *a,doublecomplex c)
    /* conjugate dot product of vector on itself; c=a.a* */
{
  int i;
  clock_t tstart;

  c[re]=c[im]=0.0;
  for (i=0;i<nlocalRows;++i) {
    c[re]+=a[i][re]*a[i][re]-a[i][im]*a[i][im];
    c[im]+=a[i][re]*a[i][im];
  }
  tstart=clock();
  my_inner_product(c,cmplx_type,1);
  Timing_OneIterComm+=clock()-tstart;
  c[im]*=2;
}

/*============================================================*/

void nDotProdSelf_conj_Norm2(doublecomplex *a,doublecomplex c,double *norm)
    /* Computes both conjugate dot product of vector on itself (c=a.a*)
         and its hermitian Norm squared norm=||a||^2 */
{
  int i;
  double buf[3];
  clock_t tstart;

  buf[0]=buf[1]=buf[2]=0.0;
  for (i=0;i<nlocalRows;++i) {
    buf[0]+=a[i][re]*a[i][re];
    buf[1]+=a[i][im]*a[i][im];
    buf[2]+=a[i][re]*a[i][im];
  }
  tstart=clock();
  my_inner_product(buf,double_type,3);
  Timing_OneIterComm+=clock()-tstart;
  *norm=buf[0]+buf[1];
  c[re]=buf[0]-buf[1];
  c[im]=2*buf[2];
}

/*============================================================*/

void nIncrem110_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2)
     /* a=c1*a+c2*b+c */
{
  int i;
  double tmp;

  for (i=0;i<nlocalRows;++i) {
    tmp=a[i][re];
    a[i][re] = c1[re]*a[i][re] - c1[im]*a[i][im] + c2[re]*b[i][re] - c2[im]*b[i][im] + c[i][re];
    a[i][im] = c1[re]*a[i][im] + c1[im]*tmp + c2[re]*b[i][im] + c2[im]*b[i][re] + c[i][im];
  }
}

/*============================================================*/

void nIncrem011_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2)
     /* a+=c1*b+c2*c */
{
  int i;

  for (i=0;i<nlocalRows;++i) {
    a[i][re] += c1[re]*b[i][re] - c1[im]*b[i][im] + c2[re]*c[i][re] - c2[im]*c[i][im];
    a[i][im] += c1[re]*b[i][im] + c1[im]*b[i][re] + c2[re]*c[i][im] + c2[im]*c[i][re];
  }
}

/*============================================================*/

void nIncrem111_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                      doublecomplex c1,doublecomplex c2,doublecomplex c3)
     /* a=c1*a+c2*b+c3*c */
{
  int i;
  double tmp;

  for (i=0;i<nlocalRows;++i) {
    tmp=a[i][re];
    a[i][re] = c1[re]*a[i][re] - c1[im]*a[i][im] + c2[re]*b[i][re] - c2[im]*b[i][im] +
               c3[re]*c[i][re] - c3[im]*c[i][im];
    a[i][im] = c1[re]*a[i][im] + c1[im]*tmp + c2[re]*b[i][im] + c2[im]*b[i][re] +
               c3[re]*c[i][im] + c3[im]*c[i][re];
  }
}

/*============================================================*/

void nIncrem01(doublecomplex *a,doublecomplex *b,double c,double *inprod)
     /* a=a+c*b, inprod=|a|^2 */
{
  int i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][re] += c*b[i][re];   /* a+=c*b */
      a[i][im] += c*b[i][im];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][re] += c*b[i][re];   /* a+=c*b */
      a[i][im] += c*b[i][im];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];   /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem10(doublecomplex *a,doublecomplex *b,double c,double *inprod)
     /* a=c*a+b, inprod=|a|^2 */
{
  int i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][re] = c*a[i][re] + b[i][re];  /* a=c*a+b */
      a[i][im] = c*a[i][im] + b[i][im];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][re] = c*a[i][re] + b[i][re];  /* a=c*a+b */
      a[i][im] = c*a[i][im] + b[i][im];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem11_d_c(doublecomplex *a,doublecomplex *b,double c1,doublecomplex c2,double *inprod)
     /* a=c1*a+c2*b, inprod=|a|^2 , one constant is double, another - complex */
{
  int i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][re] = c1*a[i][re] + c2[re]*b[i][re] - c2[im]*b[i][im];     /* a=c1*a+c2*b */
      a[i][im] = c1*a[i][im] + c2[re]*b[i][im] + c2[im]*b[i][re];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][re] = c1*a[i][re] + c2[re]*b[i][re] - c2[im]*b[i][im];     /* a=c1*a+c2*b */
      a[i][im] = c1*a[i][im] + c2[re]*b[i][im] + c2[im]*b[i][re];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem01_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c,double *inprod)
     /* a=a+c*b, inprod=|a|^2 */
{
  int i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][re] += c[re]*b[i][re] - c[im]*b[i][im];     /* a+=c*b */
      a[i][im] += c[re]*b[i][im] + c[im]*b[i][re];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][re] += c[re]*b[i][re] - c[im]*b[i][im];     /* a+=c*b */
      a[i][im] += c[re]*b[i][im] + c[im]*b[i][re];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nIncrem10_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c,double *inprod)
     /* a=c*a+b, inprod=|a|^2 */
{
  int i;
  double tmp;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      tmp=a[i][re];                                         /* a=c*a+b */
      a[i][re] = c[re]*a[i][re] - c[im]*a[i][im] + b[i][re];
      a[i][im] = c[re]*a[i][im] + c[im]*tmp + b[i][im];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      tmp=a[i][re];                                         /* a=c*a+b */
      a[i][re] = c[re]*a[i][re] - c[im]*a[i][im] + b[i][re];
      a[i][im] = c[re]*a[i][im] + c[im]*tmp + b[i][im];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nLinComb_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                     doublecomplex c1, doublecomplex c2,double *inprod)
     /* a=c1*b+c2*c, inprod=|a|^2 */
{
  int i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c2*c */
      a[i][re] = c1[re]*b[i][re] - c1[im]*b[i][im] + c2[re]*c[i][re] - c2[im]*c[i][im];
      a[i][im] = c1[re]*b[i][im] + c1[im]*b[i][re] + c2[re]*c[i][im] + c2[im]*c[i][re];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c2*c */
      a[i][re] = c1[re]*b[i][re] - c1[im]*b[i][im] + c2[re]*c[i][re] - c2[im]*c[i][im];
      a[i][im] = c1[re]*b[i][im] + c1[im]*b[i][re] + c2[re]*c[i][im] + c2[im]*c[i][re];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nLinComb1_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex *c,
                     doublecomplex c1,double *inprod)
     /* a=c1*b+c, inprod=|a|^2 */
{
  int i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c */
      a[i][re] = c1[re]*b[i][re] - c1[im]*b[i][im] + c[i][re];
      a[i][im] = c1[re]*b[i][im] + c1[im]*b[i][re] + c[i][im];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      /* a=c1*b+c */
      a[i][re] = c1[re]*b[i][re] - c1[im]*b[i][im] + c[i][re];
      a[i][im] = c1[re]*b[i][im] + c1[im]*b[i][re] + c[i][im];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nSubtr(doublecomplex *a,doublecomplex *b,doublecomplex *c,double *inprod)
     /* a=b-c, inprod=|a|^2 */
{
  int i;
  clock_t tstart;

  if (inprod==NULL) {
    for (i=0;i<nlocalRows;++i) {
      a[i][re] = b[i][re] - c[i][re];     /* a=b-c */
      a[i][im] = b[i][im] - c[i][im];
    }
  }
  else {
    *inprod=0.0;
    for (i=0;i<nlocalRows;++i) {
      a[i][re] = b[i][re] - c[i][re];     /* a=b-c */
      a[i][im] = b[i][im] - c[i][im];
      (*inprod) += a[i][re]*a[i][re] + a[i][im]*a[i][im];  /* *inprod=|a|^2   */
    }
    tstart=clock();
    my_inner_product(inprod,double_type,1);
    Timing_OneIterComm+=clock()-tstart;
  }
}

/*============================================================*/

void nMult_cmplx(doublecomplex *a,doublecomplex *b,doublecomplex c)
     /* multiply vector by a complex constant; a=c*b */
{
  int i;

  for (i=0;i<nlocalRows;++i) {
    a[i][re] = c[re]*b[i][re] - c[im]*b[i][im];     /* a[i]=c*b[i] */
    a[i][im] = c[re]*b[i][im] + c[im]*b[i][re];
  }
}

/*============================================================*/

void nMultSelf_cmplx(doublecomplex *a,doublecomplex c)
   /* multiply vector by a complex constant; a*=c */
{
  int i;
  double tmp;

  for (i=0;i<nlocalRows;++i) {
    tmp=a[i][re];
    a[i][re] = c[re]*a[i][re] - c[im]*a[i][im];     /* a[i]*=c */
    a[i][im] = c[re]*a[i][im] + c[im]*tmp;
  }
}

/*============================================================*/

void nMult_mat(doublecomplex *a,doublecomplex *b,doublecomplex c[][3])
   /* multiply by a function of material of a dipole and component; a[3*i+j]=c[mat[i]][j]*b[3*i+j] */
{
  int i,j,k;
  doublecomplex *val;

  k=0;
  for (i=0;i<local_nvoid_Ndip;++i) {
    val=c[material[i]];
    for (j=0;j<3;j++) {
      a[k][re] = val[j][re]*b[k][re] - val[j][im]*b[k][im];
      a[k][im] = val[j][re]*b[k][im] + val[j][im]*b[k][re];
      k++;
    }  
  }
}

/*============================================================*/

void nMultSelf_mat(doublecomplex *a,doublecomplex c[][3])
   /* multiply by a function of material of a dipole and component; a[3*i+j]*=c[mat[i]][j] */
{
  int i,j,k;
  double tmp;
  doublecomplex *val;

  k=0;
  for (i=0;i<local_nvoid_Ndip;++i) {
    val=c[material[i]];
    for (j=0;j<3;j++) {
      tmp=a[k][re];
      a[k][re] = val[j][re]*a[k][re] - val[j][im]*a[k][im];
      a[k][im] = val[j][re]*a[k][im] + val[j][im]*tmp;
      k++;
    }
  }
}

/*============================================================*/

void nConj(doublecomplex *a)
   /* complex conjugate of the vector */
{
  int i;

  for (i=0;i<nlocalRows;++i) a[i][im]=-a[i][im];
}

