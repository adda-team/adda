      SUBROUTINE D07HRE(NDIM,WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D07HRE
C***KEYWORDS basic integration rule, degree 7
C***PURPOSE  To initialize a degree 7 basic rule, and null rules.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-05-31
C***DESCRIPTION  D07HRE initializes a degree 7 integration rule,
C            two degree 5 null rules, one degree 3 null rule and one
C            degree 1 null rule for the hypercube [-1,1]**NDIM.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   WTLENG Integer.
C          The number of weights in each of the rules.
C          WTLENG MUST be set equal to 6.
C
C   ON RETURN
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   G      Real array of dimension (NDIM, WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1, J), ..., G(NDIM, J) are the are the generators for the
C          points associated with the Jth weights.
C   ERRCOF Real array of dimension 6.
C          Heuristic error coefficients that are used in the
C          error estimation in BASRUL.
C   RULPTS Real array of dimension WTLENG.
C          A work array.
C
C***REFERENCES A. Genz and A. Malik,
C             "An Imbedded Family of Fully Symmetric Numerical
C              Integration Rules",
C              SIAM J Numer. Anal. 20 (1983), pp. 580-588.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D07HRE
C
C   Global variables
C
      INTEGER NDIM,WTLENG
      DOUBLE PRECISION W(5,WTLENG),G(NDIM,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local Variables
C
      DOUBLE PRECISION RATIO,LAM0,LAM1,LAM2,LAMP,TWONDM
      INTEGER I,J
C
C***FIRST EXECUTABLE STATEMENT D07HRE
C
C
C     Initialize generators, weights and RULPTS
C
      DO 30 J = 1,WTLENG
          DO 10 I = 1,NDIM
              G(I,J) = 0
10        CONTINUE
          DO 20 I = 1,5
              W(I,J) = 0
20        CONTINUE
          RULPTS(J) = 2*NDIM
30    CONTINUE
      TWONDM = 2**NDIM
      RULPTS(WTLENG) = TWONDM
      RULPTS(WTLENG-1) = 2*NDIM* (NDIM-1)
      RULPTS(1) = 1
C
C     Compute squared generator parameters
C
      LAM0 = 0.4707
      LAMP = 0.5625
      LAM1 = 4/ (15-5/LAM0)
      RATIO = (1-LAM1/LAM0)/27
      LAM2 = (5-7*LAM1-35*RATIO)/ (7-35*LAM1/3-35*RATIO/LAM0)
C
C     Compute degree 7 rule weights
C
      W(1,6) = 1/ (3*LAM0)**3/TWONDM
      W(1,5) = (1-5*LAM0/3)/ (60* (LAM1-LAM0)*LAM1**2)
      W(1,3) = (1-5*LAM2/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM2))/
     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(1,5)
      W(1,2) = (1-5*LAM1/3-5*TWONDM*W(1,6)*LAM0* (LAM0-LAM1))/
     +         (10*LAM2* (LAM2-LAM1))
C
C     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules
C
      W(2,6) = 1/ (36*LAM0**3)/TWONDM
      W(2,5) = (1-9*TWONDM*W(2,6)*LAM0**2)/ (36*LAM1**2)
      W(2,3) = (1-5*LAM2/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM2))/
     +         (10*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(2,5)
      W(2,2) = (1-5*LAM1/3-5*TWONDM*W(2,6)*LAM0* (LAM0-LAM1))/
     +         (10*LAM2* (LAM2-LAM1))
      W(3,6) = 5/ (108*LAM0**3)/TWONDM
      W(3,5) = (1-9*TWONDM*W(3,6)*LAM0**2)/ (36*LAM1**2)
      W(3,3) = (1-5*LAMP/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAMP))/
     +         (10*LAM1* (LAM1-LAMP)) - 2* (NDIM-1)*W(3,5)
      W(3,4) = (1-5*LAM1/3-5*TWONDM*W(3,6)*LAM0* (LAM0-LAM1))/
     +         (10*LAMP* (LAMP-LAM1))
      W(4,6) = 1/ (54*LAM0**3)/TWONDM
      W(4,5) = (1-18*TWONDM*W(4,6)*LAM0**2)/ (72*LAM1**2)
      W(4,3) = (1-10*LAM2/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM2))/
     +         (20*LAM1* (LAM1-LAM2)) - 2* (NDIM-1)*W(4,5)
      W(4,2) = (1-10*LAM1/3-10*TWONDM*W(4,6)*LAM0* (LAM0-LAM1))/
     +         (20*LAM2* (LAM2-LAM1))
C
C     Set generator values
C
      LAM0 = SQRT(LAM0)
      LAM1 = SQRT(LAM1)
      LAM2 = SQRT(LAM2)
      LAMP = SQRT(LAMP)
      DO 40 I = 1,NDIM
          G(I,WTLENG) = LAM0
40    CONTINUE
      G(1,WTLENG-1) = LAM1
      G(2,WTLENG-1) = LAM1
      G(1,WTLENG-4) = LAM2
      G(1,WTLENG-3) = LAM1
      G(1,WTLENG-2) = LAMP
C
C     Compute final weight values.
C     The null rule weights are computed from differences between
C     the degree 7 rule weights and lower degree rule weights.
C
      W(1,1) = TWONDM
      DO 70 J = 2,5
          DO 50 I = 2,WTLENG
              W(J,I) = W(J,I) - W(1,I)
              W(J,1) = W(J,1) - RULPTS(I)*W(J,I)
50        CONTINUE
70    CONTINUE
      DO 80 I = 2,WTLENG
          W(1,I) = TWONDM*W(1,I)
          W(1,1) = W(1,1) - RULPTS(I)*W(1,I)
80    CONTINUE
C
C     Set error coefficients
C
      ERRCOF(1) = 5
      ERRCOF(2) = 5
      ERRCOF(3) = 1
      ERRCOF(4) = 5
      ERRCOF(5) = 0.5
      ERRCOF(6) = 0.25
C
C***END D07HRE
C
      END
