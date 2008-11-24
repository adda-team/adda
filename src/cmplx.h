/*
 *	cmplx.h
 *
 */

#if defined(SINGLE)
#define REAL float
#elif defined(DOUBLE)
#define REAL double
#endif

typedef struct {
        REAL   r;              /* real part      */
        REAL   i;              /* imaginary part */
}       dcomplex;                /* real-precision Complex */

typedef struct {
        double r;              /* real part      */
        double i;              /* imaginary part */
}       doublecomplex;                /* double-precision Complex */

