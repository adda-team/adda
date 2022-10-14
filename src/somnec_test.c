/* This is a standalone module for testing the calculation of Sommerfeld integrals (in somnec.c/h). It runs somnec in
 * all possible modes (with respect to integration contours), measures the execution time, and compares the results with
 * each other. While many of the contours can work for a wide range of parameters, some may fail for specific values.
 * However, this is not necessarily a failure of somnec, since such combinations should never occur in automatic regime.
 * The only exception from this test-all rationale is not executing the Bessel and Hankel forms for Z=0 and rho=0,
 * respectively. The latter case would lead to only trivial test that automatic regime leads to Hankel contour.
 *
 *
 * It is not yet included as a target in Makefile, so need to be compiled manually by a command like:
 * gcc -O3 -Wall -o somnec_test somnec*.c
 *
 * Copyright (C) ADDA contributors
 * This file is part of ADDA.
 *
 * ADDA is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 * ADDA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with ADDA. If not, see
 * <http://www.gnu.org/licenses/>.
 */
// project headers
#include "somnec.h"
// system headers
#include <math.h>
#include <stdio.h>
#include <time.h>

#include "timing.h"

// useful macros for printing complex numbers
#define CFORM "%.10g%+.10gi"
#define REIM(a) creal(a),cimag(a)

//======================================================================================================================

double DiffSystemTime(const SYSTEM_TIME * restrict t1,const SYSTEM_TIME * restrict t2)
/* compute difference (in seconds) between two system times; not very fast
 * !!! order of arguments is inverse to that in standard difftime (for historical reasons)
 * This function was copied from timing.c to make these test routine standalone.
 */
{
#ifdef WINDOWS
	LARGE_INTEGER freq;
	QueryPerformanceFrequency(&freq);
	return (double)(t2->QuadPart - t1->QuadPart)/(double)(freq.QuadPart);
#elif defined(POSIX)
	return (double)(t2->tv_sec - t1->tv_sec) + MICRO*(double)(t2->tv_usec - t1->tv_usec);
#else // fallback for 1s-precision timer
	return difftime(*t2,*t1);
#endif
}

//======================================================================================================================

double norm2(const complex double a1,const complex double a2,const complex double a3,const complex double a4)
// computes the squared L2 norm of vector from C^4, i.e. ||a||^2
{
	return a1*conj(a1)+a2*conj(a2)+a3*conj(a3)+a4*conj(a4);
}

//======================================================================================================================

double RelError(const complex double a1,const complex double a2,const complex double a3,const complex double a4,
	const complex double b1,const complex double b2,const complex double b3,const complex double b4)
// computes the relative L2 difference between two vectors from C^4, i.e. sqrt(2*||a-b||^2/(||a||^2+||b||^2))
{
	return sqrt( 2*norm2(a1-b1,a2-b2,a3-b3,a4-b4)/(norm2(a1,a2,a3,a4)+norm2(b1,b2,b3,b4)) );
}
//======================================================================================================================

int main(int argc,char **argv)
{
	double mre,mim;
	double rho,Z;
	complex double eps,erv0,ezv0,erh0,eph0,erv1,ezv1,erh1,eph1,erv2,ezv2,erh2,eph2,erv3,ezv3,erh3,eph3;
	SYSTEM_TIME t1,t2;
	// TODO: add meaningful error messages
	if (argc != 5) return 1;
	if (sscanf(argv[1],"%lf",&mre)!=1) return 2;
	if (sscanf(argv[2],"%lf",&mim)!=1) return 2;
	if (sscanf(argv[3],"%lf",&rho)!=1) return 2;
	if (sscanf(argv[4],"%lf",&Z)!=1) return 2;

	eps=mre+I*mim;
	eps=eps*eps;
	som_init(eps);

	GET_SYSTEM_TIME(&t1);
	evlua(Z,rho,&erv0,&ezv0,&erh0,&eph0,0);
	GET_SYSTEM_TIME(&t2);
	printf("Mode 0 (%.6f s):\n",DiffSystemTime(&t1,&t2));
	printf("Irv="CFORM"\n",REIM(erv0));
	printf("Izv="CFORM"\n",REIM(ezv0));
	printf("Irh="CFORM"\n",REIM(erh0));
	printf("Iph="CFORM"\n",REIM(eph0));

	if (Z==0) printf("Mode 1 (Bessel contour) does not work for Z=0\n");
	else {
		GET_SYSTEM_TIME(&t1);
		evlua(Z,rho,&erv1,&ezv1,&erh1,&eph1,1);
		GET_SYSTEM_TIME(&t2);
		printf("Mode 1 (%.6f s):\n",DiffSystemTime(&t1,&t2));
		printf("Irv="CFORM"\n",REIM(erv1));
		printf("Izv="CFORM"\n",REIM(ezv1));
		printf("Irh="CFORM"\n",REIM(erh1));
		printf("Iph="CFORM"\n",REIM(eph1));
	}

	if (rho==0) printf("Modes 2 and 3 (Hankel contours) does not work for rho=0\n");
	else {
		GET_SYSTEM_TIME(&t1);
		evlua(Z,rho,&erv2,&ezv2,&erh2,&eph2,2);
		GET_SYSTEM_TIME(&t2);
		printf("Mode 2 (%.6f s):\n",DiffSystemTime(&t1,&t2));
		printf("Irv="CFORM"\n",REIM(erv2));
		printf("Izv="CFORM"\n",REIM(ezv2));
		printf("Irh="CFORM"\n",REIM(erh2));
		printf("Iph="CFORM"\n",REIM(eph2));

		GET_SYSTEM_TIME(&t1);
		evlua(Z,rho,&erv3,&ezv3,&erh3,&eph3,3);
		GET_SYSTEM_TIME(&t2);
		printf("Mode 3 (%.6f s):\n",DiffSystemTime(&t1,&t2));
		printf("Irv="CFORM"\n",REIM(erv3));
		printf("Izv="CFORM"\n",REIM(ezv3));
		printf("Irh="CFORM"\n",REIM(erh3));
		printf("Iph="CFORM"\n",REIM(eph3));
	}

	if (rho==0) printf("\nRelative error 0vs1: %g\n",RelError(erv0,ezv0,erh0,eph0,erv1,ezv1,erh1,eph1));
	else if (Z==0) {
		printf("\nRelative error 0vs2: %g\n",RelError(erv0,ezv0,erh0,eph0,erv2,ezv2,erh2,eph2));
		printf("\nRelative error 3vs2: %g\n",RelError(erv3,ezv3,erh3,eph3,erv2,ezv2,erh2,eph2));
	}
	else {
		printf("\nRelative error 0vs1: %g\n",RelError(erv0,ezv0,erh0,eph0,erv1,ezv1,erh1,eph1));
		printf("\nRelative error 2vs1: %g\n",RelError(erv2,ezv2,erh2,eph2,erv1,ezv1,erh1,eph1));
		printf("\nRelative error 3vs1: %g\n",RelError(erv3,ezv3,erh3,eph3,erv1,ezv1,erh1,eph1));
	}

	return 0;
}
