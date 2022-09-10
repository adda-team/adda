/* This is a standalone module for testing the calculation of Sommerfeld integrals (in somnec.c/h).
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
#include <stdio.h>

// useful macros for printing complex numbers
#define CFORM "%.10g%+.10gi"
#define REIM(a) creal(a),cimag(a)

//======================================================================================================================

int main(int argc,char **argv)
{
	double mre,mim;
	double rho,Z;
	complex double eps,erv0,ezv0,erh0,eph0,erv1,ezv1,erh1,eph1,erv2,ezv2,erh2,eph2,erv3,ezv3,erh3,eph3;
	// TODO: add meaningful error messages
	if (argc != 5) return 1;
	if (sscanf(argv[1],"%lf",&mre)!=1) return 2;
	if (sscanf(argv[2],"%lf",&mim)!=1) return 2;
	if (sscanf(argv[3],"%lf",&rho)!=1) return 2;
	if (sscanf(argv[4],"%lf",&Z)!=1) return 2;

	eps=mre+I*mim;
	eps=eps*eps;
	som_init(eps);

	evlua(Z,rho,&erv0,&ezv0,&erh0,&eph0,0);
	printf("Mode 0:\n");
	printf("Irv="CFORM"\n",REIM(erv0));
	printf("Izv="CFORM"\n",REIM(ezv0));
	printf("Irh="CFORM"\n",REIM(erh0));
	printf("Iph="CFORM"\n",REIM(eph0));

	evlua(Z,rho,&erv1,&ezv1,&erh1,&eph1,1);
	printf("\nMode 1:\n");
	printf("Irv="CFORM"\n",REIM(erv1));
	printf("Izv="CFORM"\n",REIM(ezv1));
	printf("Irh="CFORM"\n",REIM(erh1));
	printf("Iph="CFORM"\n",REIM(eph1));

	evlua(Z,rho,&erv2,&ezv2,&erh2,&eph2,2);
	printf("\nMode 2:\n");
	printf("Irv="CFORM"\n",REIM(erv2));
	printf("Izv="CFORM"\n",REIM(ezv2));
	printf("Irh="CFORM"\n",REIM(erh2));
	printf("Iph="CFORM"\n",REIM(eph2));

	evlua(Z,rho,&erv3,&ezv3,&erh3,&eph3,3);
	printf("\nMode 3:\n");
	printf("Irv="CFORM"\n",REIM(erv3));
	printf("Izv="CFORM"\n",REIM(ezv3));
	printf("Irh="CFORM"\n",REIM(erh3));
	printf("Iph="CFORM"\n",REIM(eph3));

	printf("\nSum of absolute errors 0vs1: %g\n",cabs(erv0-erv1)+cabs(ezv0-ezv1)+cabs(erh0-erh1)+cabs(eph0-eph1));
	printf("\nSum of absolute errors 2vs1: %g\n",cabs(erv2-erv1)+cabs(ezv2-ezv1)+cabs(erh2-erh1)+cabs(eph2-eph1));
	printf("\nSum of absolute errors 3vs1: %g\n",cabs(erv3-erv1)+cabs(ezv3-ezv1)+cabs(erh3-erh1)+cabs(eph3-eph1));

	return 0;
}
