/* File: bess.c
 * $Date::      17/07/2019                      $
 * Descr: 		Two functions "bessjn2" and "bessjn5" for computing 2 and 5 orders of Bessel functions respectively
 * Reference:	C++ Release 1.0 By J-P Moreau, Paris.(www.jpmoreau.fr)  
 *
 * Copyright (C) 2006-2013 ADDA contributors
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>



/*	This subroutine calculates the First Kind Bessel Function of
	order 0, for any real number X. The polynomial approximation by
	series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.
	REFERENCES:
	M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	VOL.5, 1962.
*/
double bessj0 (double X) { 
	const double
		P1=1.0, P2=-0.1098628627E-2, P3=0.2734510407E-4,
		P4=-0.2073370639E-5, P5= 0.2093887211E-6,
		Q1=-0.1562499995E-1, Q2= 0.1430488765E-3, Q3=-0.6911147651E-5,
		Q4= 0.7621095161E-6, Q5=-0.9349451520E-7,
		R1= 57568490574.0, R2=-13362590354.0, R3=651619640.7,
		R4=-11214424.18, R5= 77392.33017, R6=-184.9052456,
		S1= 57568490411.0, S2=1029532985.0, S3=9494680.718,
		S4= 59272.64853, S5=267.8532712, S6=1.0;
    double
		AX,FR,FS,Z,FP,FQ,XX,Y, TMP;
	if (X==0.0) return 1.0;
	AX = fabs(X);
	if (AX < 8.0) {
		Y = X*X;
		FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
		FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
		TMP = FR/FS;
	}
	else {
		Z = 8./AX;
		Y = Z*Z;
		XX = AX-0.785398164;
		FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
		FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
		TMP = sqrt(0.636619772/AX)*(FP*cos(XX)-Z*FQ*sin(XX));
	}
return TMP;
}

double Sign(double X, double Y) {
	if (Y<0.0) return (-fabs(X));
	else return (fabs(X));
}


/*	This subroutine calculates the First Kind Bessel Function of
	order 1, for any real number X. The polynomial approximation by
	series of Chebyshev polynomials is used for 0<X<8 and 0<8/X<1.	
	REFERENCES:
	M.ABRAMOWITZ,I.A.STEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, 1965.
	C.W.CLENSHAW, NATIONAL PHYSICAL LABORATORY MATHEMATICAL TABLES,
	VOL.5, 1962.
*/
double bessj1 (double X) {
	const double  
		P1=1.0, P2=0.183105E-2, P3=-0.3516396496E-4, P4=0.2457520174E-5,
		P5=-0.240337019E-6,  P6=0.636619772,
		Q1= 0.04687499995, Q2=-0.2002690873E-3, Q3=0.8449199096E-5,
		Q4=-0.88228987E-6, Q5= 0.105787412E-6,
		R1= 72362614232.0, R2=-7895059235.0, R3=242396853.1,
		R4=-2972611.439,   R5=15704.48260,  R6=-30.16036606,
		S1=144725228442.0, S2=2300535178.0, S3=18583304.74,
		S4=99447.43394,    S5=376.9991397,  S6=1.0;

	double AX,FR,FS,Y,Z,FP,FQ,XX, TMP;

	AX = fabs(X);
	if (AX < 8.0) {
		Y = X*X;
		FR = R1+Y*(R2+Y*(R3+Y*(R4+Y*(R5+Y*R6))));
		FS = S1+Y*(S2+Y*(S3+Y*(S4+Y*(S5+Y*S6))));
		TMP = X*(FR/FS);
	}
	else {
		Z = 8.0/AX;
		Y = Z*Z;
		XX = AX-2.35619491;
		FP = P1+Y*(P2+Y*(P3+Y*(P4+Y*P5)));
		FQ = Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*Q5)));
		TMP = sqrt(P6/AX)*(cos(XX)*FP-Z*sin(XX)*FQ)*Sign(S6,X);
	}
	return TMP;
}
/*
void bessjn5 (int N2, double X, double ARRJ[5]) { // Computing 5 orders of Bessel function
	const int A = 40;
	const double BIGNO = 1e10,  BIGNI = 1e-10;

	double TOX,BJM,BJ,BJP,SUM,TMP,BJM1,BJ1,BJP1,SUM1,TMP1;
	int J, JSUM, JSUM1, M, K, N, sgn, sgn1;

	N = (N2 - 2);
	if (X == 0.0) {
		switch(N2) {
			case -2:
			ARRJ[0] = ARRJ[1] = ARRJ[2] = ARRJ[3] = 0.0;
			ARRJ[4] = 1.0;
			break;
			case -1:
			ARRJ[0] = ARRJ[1] = ARRJ[2] = ARRJ[4] = 0.0;
			ARRJ[3] = 1.0;
			break;
			case 0:
			ARRJ[0] = ARRJ[1] = ARRJ[3] = ARRJ[4] = 0.0;
			ARRJ[2] = 1.0;
			break;
			case 1:
			ARRJ[0] = ARRJ[2] = ARRJ[3] = ARRJ[4] = 0.0;
			ARRJ[1] = 1.0;
			break;
			case 2:
			ARRJ[1] = ARRJ[2] = ARRJ[3] = ARRJ[4] = 0.0;
			ARRJ[0] = 1.0;
			break;
			default:
			ARRJ[0] = ARRJ[1] = ARRJ[2] = ARRJ[3] = ARRJ[4] =0.0;
			break;
		}
	}
	else {
		if ((N == -5)||(N == -4)||(N == -3)||(N == -2)||(N == -1)||(N == 0)||(N == 1)) {
			switch(N) {
				case -5:
				ARRJ[4] = -bessj1(X); 
				ARRJ[3] = -(2.0/X)*ARRJ[4]-bessj0(X);
				ARRJ[2] = -(4.0/X)*ARRJ[3]-ARRJ[4];
				ARRJ[1] = -(6.0/X)*ARRJ[2]-ARRJ[3];
				ARRJ[0] = -(8.0/X)*ARRJ[1]-ARRJ[2];
				break;
				case -4:
				ARRJ[4] = bessj0(X); 
				ARRJ[3] = -bessj1(X);
				ARRJ[2] = -(2.0/X)*ARRJ[3]-ARRJ[4];
				ARRJ[1] = -(4.0/X)*ARRJ[2]-ARRJ[3];
				ARRJ[0] = -(6.0/X)*ARRJ[1]-ARRJ[2];
				break;
				case -3:
				ARRJ[4] = bessj1(X); 
				ARRJ[3] = bessj0(X);
				ARRJ[2] = -ARRJ[4];;
				ARRJ[1] = -(2.0/X)*ARRJ[2]-ARRJ[3];
				ARRJ[0] = -(4.0/X)*ARRJ[1]-ARRJ[2];
				break;
				case -2:
				ARRJ[3] = bessj1(X); 
				ARRJ[2] = bessj0(X);
				ARRJ[1] = -ARRJ[4];;
				ARRJ[1] = -(2.0/X)*ARRJ[2]-ARRJ[3];
				ARRJ[0] = -(4.0/X)*ARRJ[1]-ARRJ[2];
				break;
				case -1:
				ARRJ[0] = -bessj1(X);
				ARRJ[1] = bessj0(X);
				ARRJ[2] = -ARRJ[0];
				ARRJ[3] = (2.0/X)*ARRJ[2]-ARRJ[1];
				ARRJ[4] = (4.0/X)*ARRJ[3]-ARRJ[2];
				break;
				case 0:
				ARRJ[0] = bessj0(X);
				ARRJ[1] = bessj1(X);
				ARRJ[2] = (2.0/X)*ARRJ[1]-ARRJ[0];
				ARRJ[3] = (4.0/X)*ARRJ[2]-ARRJ[1];
				ARRJ[4] = (6.0/X)*ARRJ[3]-ARRJ[2];
				break;
				case 1:
				ARRJ[0] = bessj1(X);
				ARRJ[1] = (2.0/X)*ARRJ[0]-bessj0(X);
				ARRJ[2] = (4.0/X)*ARRJ[1]-ARRJ[0];
				ARRJ[3] = (6.0/X)*ARRJ[2]-ARRJ[1];
				ARRJ[4] = (8.0/X)*ARRJ[3]-ARRJ[2];
				break;
			}
		}
		else {
			TOX = 2.0/X;
			sgn = 1;
			sgn1 = 1;
			N = abs(N2-2);
			if ((N2-2 <0)&((N2-2)%2 != 0)) sgn = -1;
			if ((N2-1 <0)&((N2-1)%2 != 0)) sgn1 = -1;
			if (X > 1.0*N) {
				BJM = bessj0(X);
				BJ  = bessj1(X);
				
				for (J=1; J<N; J++) {  // Upward recurrence
					BJP = J*TOX*BJ-BJM;
					BJM = BJ;
					BJ  = BJP;
				}
				
				BJP = N*TOX*BJ-BJM;
				ARRJ[0] = sgn*BJ; 
				ARRJ[1] = sgn1*BJP;
				ARRJ[2] = (N+2)*(2.0/X)*ARRJ[1]-ARRJ[0];
				ARRJ[3] = (N+3)*(2.0/X)*ARRJ[2]-ARRJ[1];
				ARRJ[4] = (N+4)*(2.0/X)*ARRJ[3]-ARRJ[2];
			}
			else {
				M = (int) (2*((N+floor(sqrt(1.0*(A*N))))/2));
				K = (int) (2*((N+1+floor(sqrt(1.0*(A*(N+1)))))/2));
				TMP = 0.0;
				TMP1 = 0.0;
				JSUM = 0;
				JSUM1 = 0;
				SUM = 0.0;
				SUM1 = 0.0;
				BJP = 0.0;
				BJP1 = 0.0;
				BJ  = 1.0;
				BJ1  = 1.0;
				
				for (J=K; J>0; J--) {  // Downward recurrence ( Miller's algorithm )
					BJM1 = J*TOX*BJ1-BJP1;
					BJP1 = BJ1;
					BJ1  = BJM1;
					if (fabs(BJ1) > BIGNO) {
						BJ1  = BJ1*BIGNI;
						BJP1 = BJP1*BIGNI;
						TMP1 = TMP1*BIGNI;
						SUM1 = SUM1*BIGNI;
					}
					if (JSUM1 != 0)  SUM1 += BJ1;
					JSUM1 = 1-JSUM1;
					if (J == N+1)  TMP1 = BJP1;
					if (J <= M) {
						BJM = J*TOX*BJ-BJP;
						BJP = BJ;
						BJ  = BJM;
						if (fabs(BJ) > BIGNO) {
							BJ  = BJ*BIGNI;
							BJP = BJP*BIGNI;
							TMP = TMP*BIGNI;
							SUM = SUM*BIGNI;
						}
						if (JSUM != 0)  SUM += BJ;
						JSUM = 1-JSUM;
						if (J == N)  TMP = BJP;
					}
				}
				
				ARRJ[0] = sgn*TMP/(2.0*SUM-BJ);
				ARRJ[1] = sgn1*TMP1/(2.0*SUM1-BJ1);
				ARRJ[2] = (N+2)*(2.0/X)*ARRJ[1]-ARRJ[0];
				ARRJ[3] = (N+3)*(2.0/X)*ARRJ[2]-ARRJ[1];
				ARRJ[4] = (N+4)*(2.0/X)*ARRJ[3]-ARRJ[2];

			}
		}
	}
}
*/

void bessjn2 (int N, double X, double ARRJ[2]) {  // Computing 2 orders of Bessel function
	const int A = 40;
	const double BIGNO = 1e10,  BIGNI = 1e-10;

	double TOX,BJM,BJ,BJP,SUM,TMP,BJM1,BJ1,BJP1,SUM1,TMP1;
	int J, JSUM, JSUM1, M, K;
	if (X == 0.0) {
		if (N == 0) { ARRJ[0] = 1.0; ARRJ[1] = 0.0;}
		else { ARRJ[0] = 0.0; ARRJ[1] = 0.0;} 
	}
	else {
		if ((N == 0)||(N == 1)) {
			if (N == 0) { ARRJ[0] = bessj0(X); ARRJ[1] = bessj1(X);}
			if (N == 1) { ARRJ[0] = bessj1(X); ARRJ[1] = (2.0/X)*ARRJ[0]-bessj0(X);}
		}
		else {
			TOX = 2.0/X;
			if (X > 1.0*N) {
				BJM = bessj0(X);
				BJ  = bessj1(X);
				
				for (J=1; J<N; J++) {  // Upward recurrence
					BJP = J*TOX*BJ-BJM;
					BJM = BJ;
					BJ  = BJP;
				}
				
				BJP = N*TOX*BJ-BJM;
				ARRJ[0] = BJ; ARRJ[1] = BJP;
			}
			else {
				M = (int) (2*((N+floor(sqrt(1.0*(A*N))))/2));
				K = (int) (2*((N+1+floor(sqrt(1.0*(A*(N+1)))))/2));
				TMP = 0.0;
				TMP1 = 0.0;
				JSUM = 0;
				JSUM1 = 0;
				SUM = 0.0;
				SUM1 = 0.0;
				BJP = 0.0;
				BJP1 = 0.0;
				BJ  = 1.0;
				BJ1  = 1.0;
				
				for (J=K; J>0; J--) {  // Downward recurrence ( Miller's algorithm )
					BJM1 = J*TOX*BJ1-BJP1;
					BJP1 = BJ1;
					BJ1  = BJM1;
					if (fabs(BJ1) > BIGNO) {
						BJ1  = BJ1*BIGNI;
						BJP1 = BJP1*BIGNI;
						TMP1 = TMP1*BIGNI;
						SUM1 = SUM1*BIGNI;
					}
					if (JSUM1 != 0)  SUM1 += BJ1;
					JSUM1 = 1-JSUM1;
					if (J == N+1)  TMP1 = BJP1;
					if (J <= M) {
						BJM = J*TOX*BJ-BJP;
						BJP = BJ;
						BJ  = BJM;
						if (fabs(BJ) > BIGNO) {
							BJ  = BJ*BIGNI;
							BJP = BJP*BIGNI;
							TMP = TMP*BIGNI;
							SUM = SUM*BIGNI;
						}
						if (JSUM != 0)  SUM += BJ;
						JSUM = 1-JSUM;
						if (J == N)  TMP = BJP;
					}
				}
				
				ARRJ[0] = TMP/(2.0*SUM-BJ);
				ARRJ[1] = TMP1/(2.0*SUM1-BJ1);
			}
		}
	}
}
