      SUBROUTINE DRLHRE(NDIM,CENTER,HWIDTH,WTLENG,G,W,ERRCOF,NUMFUN,
     +                  FUNSUB,SCALES,NORMS,X,NULL,BASVAL,RGNERR,DIRECT)
C***BEGIN PROLOGUE DRLHRE
C***KEYWORDS basic numerical integration rule
C***PURPOSE  To compute basic integration rule values.
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 90-02-06
C***DESCRIPTION DRLHRE computes basic integration rule values for a
C            vector of integrands over a hyper-rectangular region.
C            These are estimates for the integrals. DRLHRE also computes
C            estimates for the errors and determines the coordinate axis
C            where the fourth difference for the integrands is largest.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   WTLENG Integer.
C          The number of weights in the basic integration rule.
C   G      Real array of dimension (NDIM,WTLENG).
C          The fully symmetric sum generators for the rules.
C          G(1,J), ..., G(NDIM,J) are the are the generators for the
C          points associated with the Jth weights.
C   W      Real array of dimension (5,WTLENG).
C          The weights for the basic and null rules.
C          W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C   ERRCOF Real array of dimension 6.
C          The error coefficients for the rules.
C          It is assumed that the error is computed using:
C           IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3)
C             THEN ERROR = ERRCOF(3)*N1
C             ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3)
C           ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6))
C          where N1-N4 are the null rules, EP is the error
C          for the parent
C          subregion and ES is the error for the sibling subregion.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM,X,NUMFUN,FUNVLS).
C           Input Parameters:
C            X      Real array of dimension NDIM.
C                   Defines the evaluation point.
C            NDIM   Integer.
C                   Number of variables for the integrand.
C            NUMFUN Integer.
C                   Number of components for the vector integrand.
C           Output Parameters:
C            FUNVLS Real array of dimension NUMFUN.
C                   The components of the integrand at the point X.
C   SCALES Real array of dimension (3,WTLENG).
C          Scaling factors used to construct new null rules based
C          on a linear combination of two successive null rules
C          in the sequence of null rules.
C   NORMS  Real array of dimension (3,WTLENG).
C          2**NDIM/(1-norm of the null rule constructed by each of the
C          scaling factors.)
C   X      Real Array of dimension NDIM.
C          A work array.
C   NULL   Real array of dimension (NUMFUN, 8)
C          A work array.
C
C   ON RETURN
C
C   BASVAL Real array of dimension NUMFUN.
C          The values for the basic rule for each component
C          of the integrand.
C   RGNERR Real array of dimension NUMFUN.
C          The error estimates for each component of the integrand.
C   DIRECT Real.
C          The coordinate axis where the fourth difference of the
C          integrand values is largest.
C
C***REFERENCES
C   A.C.Genz and A.A.Malik, An adaptive algorithm for numerical
C   integration over an N-dimensional rectangular region,
C   J.Comp.Appl.Math., 6:295-302, 1980.
C
C   T.O.Espelid, Integration Rules, Null Rules and Error
C   Estimation, Reports in Informatics 33, Dept. of Informatics,
C   Univ. of Bergen, 1988.
C
C***ROUTINES CALLED: DFSHRE, FUNSUB
C
C***END PROLOGUE DRLHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER WTLENG,NUMFUN,NDIM
      DOUBLE PRECISION CENTER(NDIM),X(NDIM),HWIDTH(NDIM),BASVAL(NUMFUN),
     +                 RGNERR(NUMFUN),NULL(NUMFUN,8),W(5,WTLENG),
     +                 G(NDIM,WTLENG),ERRCOF(6),DIRECT,SCALES(3,WTLENG),
     +                 NORMS(3,WTLENG)
C
C   Local variables.
C
      DOUBLE PRECISION RGNVOL,DIFSUM,DIFMAX,FRTHDF
      INTEGER I,J,K,DIVAXN
      DOUBLE PRECISION SEARCH,RATIO
C
C***FIRST EXECUTABLE STATEMENT DRLHRE
C
C
C       Compute volume of subregion, initialize DIVAXN and rule sums;
C       compute fourth differences and new DIVAXN (RGNERR is used
C       for a work array here). The integrand values used for the
C       fourth divided differences are accumulated in rule arrays.
C
      RGNVOL = 1
      DIVAXN = 1
      DO 10 I = 1,NDIM
          RGNVOL = RGNVOL*HWIDTH(I)
          X(I) = CENTER(I)
          IF (HWIDTH(I).GT.HWIDTH(DIVAXN)) DIVAXN = I
10    CONTINUE
      CALL FUNSUB(NDIM,X,NUMFUN,RGNERR)
      DO 30 J = 1,NUMFUN
          BASVAL(J) = W(1,1)*RGNERR(J)
          DO 20 K = 1,4
              NULL(J,K) = W(K+1,1)*RGNERR(J)
20        CONTINUE
30    CONTINUE
      DIFMAX = 0
      RATIO = (G(1,3)/G(1,2))**2
      DO 60 I = 1,NDIM
          X(I) = CENTER(I) - HWIDTH(I)*G(1,2)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,5))
          X(I) = CENTER(I) + HWIDTH(I)*G(1,2)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,6))
          X(I) = CENTER(I) - HWIDTH(I)*G(1,3)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,7))
          X(I) = CENTER(I) + HWIDTH(I)*G(1,3)
          CALL FUNSUB(NDIM,X,NUMFUN,NULL(1,8))
          X(I) = CENTER(I)
          DIFSUM = 0
          DO 50 J = 1,NUMFUN
              FRTHDF = 2* (1-RATIO)*RGNERR(J) - (NULL(J,7)+NULL(J,8)) +
     +                 RATIO* (NULL(J,5)+NULL(J,6))
C
C           Ignore differences below roundoff
C
              IF (RGNERR(J)+FRTHDF/4.NE.RGNERR(J)) DIFSUM = DIFSUM +
     +            ABS(FRTHDF)
              DO 40 K = 1,4
                  NULL(J,K) = NULL(J,K) + W(K+1,2)*
     +                        (NULL(J,5)+NULL(J,6)) +
     +                        W(K+1,3)* (NULL(J,7)+NULL(J,8))
40            CONTINUE
              BASVAL(J) = BASVAL(J) + W(1,2)* (NULL(J,5)+NULL(J,6)) +
     +                    W(1,3)* (NULL(J,7)+NULL(J,8))
50        CONTINUE
          IF (DIFSUM.GT.DIFMAX) THEN
              DIFMAX = DIFSUM
              DIVAXN = I
          END IF
60    CONTINUE
      DIRECT = DIVAXN
C
C    Finish computing the rule values.
C
      DO 90 I = 4,WTLENG
          CALL DFSHRE(NDIM,CENTER,HWIDTH,X,G(1,I),NUMFUN,FUNSUB,RGNERR,
     +                NULL(1,5))
          DO 80 J = 1,NUMFUN
              BASVAL(J) = BASVAL(J) + W(1,I)*RGNERR(J)
              DO 70 K = 1,4
                  NULL(J,K) = NULL(J,K) + W(K+1,I)*RGNERR(J)
70            CONTINUE
80        CONTINUE
90    CONTINUE
C
C    Compute errors.
C
      DO 130 J = 1,NUMFUN
C
C    We search for the null rule, in the linear space spanned by two
C    successive null rules in our sequence, which gives the greatest
C    error estimate among all normalized (1-norm) null rules in this
C    space.
C
          DO 110 I = 1,3
              SEARCH = 0
              DO 100 K = 1,WTLENG
                  SEARCH = MAX(SEARCH,ABS(NULL(J,I+1)+SCALES(I,
     +                     K)*NULL(J,I))*NORMS(I,K))
100           CONTINUE
              NULL(J,I) = SEARCH
110       CONTINUE
          IF (ERRCOF(1)*NULL(J,1).LE.NULL(J,2) .AND.
     +        ERRCOF(2)*NULL(J,2).LE.NULL(J,3)) THEN
              RGNERR(J) = ERRCOF(3)*NULL(J,1)
          ELSE
              RGNERR(J) = ERRCOF(4)*MAX(NULL(J,1),NULL(J,2),NULL(J,3))
          END IF
          RGNERR(J) = RGNVOL*RGNERR(J)
          BASVAL(J) = RGNVOL*BASVAL(J)
130   CONTINUE
C
C***END DRLHRE
C
      END
