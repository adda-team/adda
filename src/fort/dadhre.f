      SUBROUTINE DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,
     +                  EPSABS,EPSREL,KEY,RESTAR,NUM,LENW,WTLENG,
     +                  RESULT,ABSERR,NEVAL,NSUB,IFAIL,VALUES,
     +                  ERRORS,CENTRS,HWIDTS,GREATE,DIR,OLDRES,WORK,G,W,
     +                  RULPTS,CENTER,HWIDTH,X,SCALES,NORMS)
C***BEGIN PROLOGUE DADHRE
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals, I, over a hyper-rectangular
C            region hopefully satisfying for each component of I the
C            following claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            DADHRE repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with  greatest
C            estimated errors until the error request
C            is met or MAXSUB subregions are stored.
C            The regions are divided in two equally sized parts along
C            the direction with greatest absolute fourth divided
C            difference.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            Defines the number of new subregions that are divided
C            in each subdivision step.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINSUB Integer.
C            The computations proceed until there are at least
C            MINSUB subregions in the data structure.
C     MAXSUB Integer.
C            The computations proceed until there are at most
C            MAXSUB subregions in the data structure.
C
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand in the given
C            evaluation point.
C            It must have parameters (NDIM,X,NUMFUN,FUNVLS)
C            Input parameters:
C              NDIM   Integer that defines the dimension of the
C                     integral.
C              X      Real array of dimension NDIM
C                     that defines the evaluation point.
C              NUMFUN Integer that defines the number of
C                     components of I.
C            Output parameter:
C              FUNVLS Real array of dimension NUMFUN
C                     that defines NUMFUN components of the integrand.
C
C     EPSABS Real.
C            Requested absolute error.
C     EPSREL Real.
C            Requested relative error.
C     KEY    Integer.
C            Key to selected local integration rule.
C            KEY = 0 is the default value.
C                  For NDIM = 2 the degree 13 rule is selected.
C                  For NDIM = 3 the degree 11 rule is selected.
C                  For NDIM > 3 the degree  9 rule is selected.
C            KEY = 1 gives the user the 2 dimensional degree 13
C                  integration rule that uses 65 evaluation points.
C            KEY = 2 gives the user the 3 dimensional degree 11
C                  integration rule that uses 127 evaluation points.
C            KEY = 3 gives the user the degree 9 integration rule.
C            KEY = 4 gives the user the degree 7 integration rule.
C                  This is the recommended rule for problems that
C                  require great adaptivity.
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            (In this case the output parameters from DADHRE
C            must not be changed since the last
C            exit from DADHRE.)
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     LENW   Integer.
C            Defines the length of the working array WORK.
C            LENW should be greater or equal to
C            16*MDIV*NUMFUN.
C     WTLENG Integer.
C            The number of weights in the basic integration rule.
C     NSUB   Integer.
C            If RESTAR = 1, then NSUB must specify the number
C            of subregions stored in the previous call to DADHRE.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute accuracies.
C     NEVAL  Integer.
C            Number of function evaluations used by DADHRE.
C     NSUB   Integer.
C            Number of stored subregions.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less
C              subregions processed for all values of K,
C              1 <=  K <=  NUMFUN.
C            IFAIL = 1 if MAXSUB was too small for DADHRE
C              to obtain the required accuracy. In this case DADHRE
C              returns values of RESULT with estimated absolute
C              accuracies ABSERR.
C     VALUES Real array of dimension (NUMFUN,MAXSUB).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,MAXSUB).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,MAXSUB).
C            Used to store the centers of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,MAXSUB).
C            Used to store the half widths of the stored subregions.
C     GREATE Real array of dimension MAXSUB.
C            Used to store the greatest estimated errors in
C            all subregions.
C     DIR    Real array of dimension MAXSUB.
C            DIR is used to store the directions for
C            further subdivision.
C     OLDRES Real array of dimension (NUMFUN,MDIV).
C            Used to store old estimates of the integrals over the
C            subregions.
C     WORK   Real array of dimension LENW.
C            Used  in DRLHRE and DTRHRE.
C     G      Real array of dimension (NDIM,WTLENG,2*MDIV).
C            The fully symmetric sum generators for the rules.
C            G(1,J,1),...,G(NDIM,J,1) are the generators for the
C            points associated with the Jth weights.
C            When MDIV subregions are divided in 2*MDIV
C            subregions, the subregions may be processed on different
C            processors and we must make a copy of the generators
C            for each processor.
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1), ..., W(1,WTLENG) are weights for the basic rule.
C            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights.
C     RULPTS Real array of dimension WTLENG.
C            Work array used in DINHRE.
C     CENTER Real array of dimension NDIM.
C            Work array used in DTRHRE.
C     HWIDTH Real array of dimension NDIM.
C            Work array used in DTRHRE.
C     X      Real array of dimension (NDIM,2*MDIV).
C            Work array used in DRLHRE.
C     SCALES Real array of dimension (3,WTLENG).
C            Work array used by DINHRE and DRLHRE.
C     NORMS  Real array of dimension (3,WTLENG).
C            Work array used by DINHRE and DRLHRE.
C
C***REFERENCES
C
C   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm
C   for numerical integration over an n-dimensional cube, J.Comput.Appl.
C   Math. 2(1976)207-217.
C
C   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006:
C   An adaptive algorithm for numerical integration over an
C   N-dimensional rectangular region,J.Comput.Appl.Math. 6(1980)295-302.
C
C***ROUTINES CALLED DTRHRE,DINHRE,DRLHRE
C***END PROLOGUE DADHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MDIV,MINSUB,MAXSUB,KEY,LENW,RESTAR
      INTEGER NUM,NEVAL,NSUB,IFAIL,WTLENG
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN)
      DOUBLE PRECISION VALUES(NUMFUN,MAXSUB),ERRORS(NUMFUN,MAXSUB)
      DOUBLE PRECISION CENTRS(NDIM,MAXSUB)
      DOUBLE PRECISION HWIDTS(NDIM,MAXSUB)
      DOUBLE PRECISION GREATE(MAXSUB),DIR(MAXSUB)
      DOUBLE PRECISION OLDRES(NUMFUN,MDIV)
      DOUBLE PRECISION WORK(LENW),RULPTS(WTLENG)
      DOUBLE PRECISION G(NDIM,WTLENG,2*MDIV),W(5,WTLENG)
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM,2*MDIV)
      DOUBLE PRECISION SCALES(3,WTLENG),NORMS(3,WTLENG)
C
C   Local variables.
C
C   INTSGN is used to get correct sign on the integral.
C   SBRGNS is the number of stored subregions.
C   NDIV   The number of subregions to be divided in each main step.
C   POINTR Pointer to the position in the data structure where
C          the new subregions are to be stored.
C   DIRECT Direction of subdivision.
C   ERRCOF Heuristic error coeff. defined in DINHRE and used by DRLHRE
C          and DADHRE.
C
      INTEGER I,J,K
      INTEGER INTSGN,SBRGNS
      INTEGER L1
      INTEGER NDIV,POINTR,DIRECT,INDEX
      DOUBLE PRECISION OLDCEN,EST1,EST2,ERRCOF(6)
C
C***FIRST EXECUTABLE STATEMENT DADHRE
C
C   Get the correct sign on the integral.
C
      INTSGN = 1
      DO 10 J = 1,NDIM
          IF (B(J).LT.A(J)) THEN
              INTSGN = - INTSGN
          END IF
10    CONTINUE
C
C   Call DINHRE to compute the weights and abscissas of
C   the function evaluation points.
C
      CALL DINHRE(NDIM,KEY,WTLENG,W,G,ERRCOF,RULPTS,SCALES,NORMS)
C
C   If RESTAR = 1, then this is a restart run.
C
      IF (RESTAR.EQ.1) THEN
          SBRGNS = NSUB
          GO TO 110
      END IF
C
C   Initialize the SBRGNS, CENTRS and HWIDTS.
C
      SBRGNS = 1
      DO 15 J = 1,NDIM
          CENTRS(J,1) = (A(J)+B(J))/2
          HWIDTS(J,1) = ABS(B(J)-A(J))/2
15    CONTINUE
C
C   Initialize RESULT, ABSERR and NEVAL.
C
      DO 20 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
20    CONTINUE
      NEVAL = 0
C
C   Apply DRLHRE over the whole region.
C
      CALL DRLHRE(NDIM,CENTRS(1,1),HWIDTS(1,1),WTLENG,G,W,ERRCOF,NUMFUN,
     +            FUNSUB,SCALES,NORMS,X,WORK,VALUES(1,1),ERRORS(1,1),
     +            DIR(1))
      NEVAL = NEVAL + NUM
C
C   Add the computed values to RESULT and ABSERR.
C
      DO 55 J = 1,NUMFUN
          RESULT(J) = RESULT(J) + VALUES(J,1)
55    CONTINUE
      DO 65 J = 1,NUMFUN
          ABSERR(J) = ABSERR(J) + ERRORS(J,1)
65    CONTINUE
C
C   Store results in heap.
C
      INDEX = 1
      CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,HWIDTS,
     +            GREATE,WORK(1),WORK(NUMFUN+1),CENTER,HWIDTH,DIR)
C
C***End initialisation.
C
C***Begin loop while the error is too great,
C   and SBRGNS+1 is less than MAXSUB.
C
110   IF (SBRGNS+1.LE.MAXSUB) THEN
C
C   If we are allowed to divide further,
C   prepare to apply basic rule over each half of the
C   NDIV subregions with greatest errors.
C   If MAXSUB is great enough, NDIV = MDIV
C
          IF (MDIV.GT.1) THEN
              NDIV = MAXSUB - SBRGNS
              NDIV = MIN(NDIV,MDIV,SBRGNS)
          ELSE
              NDIV = 1
          END IF
C
C   Divide the NDIV subregions in two halves, and compute
C   integral and error over each half.
C
          DO 150 I = 1,NDIV
              POINTR = SBRGNS + NDIV + 1 - I
C
C   Adjust RESULT and ABSERR.
C
              DO 115 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) - VALUES(J,1)
                  ABSERR(J) = ABSERR(J) - ERRORS(J,1)
115           CONTINUE
C
C   Compute first half region.
C
              DO 120 J = 1,NDIM
                  CENTRS(J,POINTR) = CENTRS(J,1)
                  HWIDTS(J,POINTR) = HWIDTS(J,1)
120           CONTINUE
              DIRECT = DIR(1)
              DIR(POINTR) = DIRECT
              HWIDTS(DIRECT,POINTR) = HWIDTS(DIRECT,1)/2
              OLDCEN = CENTRS(DIRECT,1)
              CENTRS(DIRECT,POINTR) = OLDCEN - HWIDTS(DIRECT,POINTR)
C
C   Save the computed values of the integrals.
C
              DO 125 J = 1,NUMFUN
                  OLDRES(J,NDIV-I+1) = VALUES(J,1)
125           CONTINUE
C
C   Adjust the heap.
C
              CALL DTRHRE(1,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
     +                    HWIDTH,DIR)
C
C   Compute second half region.
C
              DO 130 J = 1,NDIM
                  CENTRS(J,POINTR-1) = CENTRS(J,POINTR)
                  HWIDTS(J,POINTR-1) = HWIDTS(J,POINTR)
130           CONTINUE
              CENTRS(DIRECT,POINTR-1) = OLDCEN + HWIDTS(DIRECT,POINTR)
              HWIDTS(DIRECT,POINTR-1) = HWIDTS(DIRECT,POINTR)
              DIR(POINTR-1) = DIRECT
150       CONTINUE
C
C   Make copies of the generators for each processor.
C
          DO 190 I = 2,2*NDIV
              DO 190 J = 1,NDIM
                  DO 190 K = 1,WTLENG
                      G(J,K,I) = G(J,K,1)
190       CONTINUE
C
C   Apply basic rule.
C
Cvd$l cncall
          DO 200 I = 1,2*NDIV
              INDEX = SBRGNS + I
              L1 = 1 + (I-1)*8*NUMFUN
              CALL DRLHRE(NDIM,CENTRS(1,INDEX),HWIDTS(1,INDEX),WTLENG,
     +                    G(1,1,I),W,ERRCOF,NUMFUN,FUNSUB,SCALES,NORMS,
     +                    X(1,I),WORK(L1),VALUES(1,INDEX),
     +                    ERRORS(1,INDEX),DIR(INDEX))
200       CONTINUE
          NEVAL = NEVAL + 2*NDIV*NUM
C
C   Add new contributions to RESULT.
C
          DO 220 I = 1,2*NDIV
              DO 210 J = 1,NUMFUN
                  RESULT(J) = RESULT(J) + VALUES(J,SBRGNS+I)
210           CONTINUE
220       CONTINUE
C
C   Check consistency of results and if necessary adjust
C   the estimated errors.
C
          DO 240 I = 1,NDIV
              GREATE(SBRGNS+2*I-1) = 0
              GREATE(SBRGNS+2*I) = 0
              DO 230 J = 1,NUMFUN
                  EST1 = ABS(OLDRES(J,I)- (VALUES(J,
     +                   SBRGNS+2*I-1)+VALUES(J,SBRGNS+2*I)))
                  EST2 = ERRORS(J,SBRGNS+2*I-1) + ERRORS(J,SBRGNS+2*I)
                  IF (EST2.GT.0) THEN
                      ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)*
     +                  (1+ERRCOF(5)*EST1/EST2)
                      ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)*
     +                                       (1+ERRCOF(5)*EST1/EST2)
                  END IF
                  ERRORS(J,SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1) +
     +                                     ERRCOF(6)*EST1
                  ERRORS(J,SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I) +
     +                                   ERRCOF(6)*EST1
                  IF (ERRORS(J,SBRGNS+2*I-1).GT.
     +                GREATE(SBRGNS+2*I-1)) THEN
                      GREATE(SBRGNS+2*I-1) = ERRORS(J,SBRGNS+2*I-1)
                  END IF
                  IF (ERRORS(J,SBRGNS+2*I).GT.GREATE(SBRGNS+2*I)) THEN
                      GREATE(SBRGNS+2*I) = ERRORS(J,SBRGNS+2*I)
                  END IF
                  ABSERR(J) = ABSERR(J) + ERRORS(J,SBRGNS+2*I-1) +
     +                        ERRORS(J,SBRGNS+2*I)
230           CONTINUE
240       CONTINUE
C
C   Store results in heap.
C
          DO 250 I = 1,2*NDIV
              INDEX = SBRGNS + I
              CALL DTRHRE(2,NDIM,NUMFUN,INDEX,VALUES,ERRORS,CENTRS,
     +                    HWIDTS,GREATE,WORK(1),WORK(NUMFUN+1),CENTER,
     +                    HWIDTH,DIR)
250       CONTINUE
          SBRGNS = SBRGNS + 2*NDIV
C
C   Check for termination.
C
          IF (SBRGNS.LT.MINSUB) THEN
              GO TO 110
          END IF
          DO 255 J = 1,NUMFUN
              IF (ABSERR(J).GT.EPSREL*ABS(RESULT(J)) .AND.
     +            ABSERR(J).GT.EPSABS) THEN
                  GO TO 110
              END IF
255       CONTINUE
          IFAIL = 0
          GO TO 499
C
C   Else we did not succeed with the
C   given value of MAXSUB.
C
      ELSE
          IFAIL = 1
      END IF
C
C   Compute more accurate values of RESULT and ABSERR.
C
499   CONTINUE
      DO 500 J = 1,NUMFUN
          RESULT(J) = 0
          ABSERR(J) = 0
500   CONTINUE
      DO 510 I = 1,SBRGNS
          DO 505 J = 1,NUMFUN
              RESULT(J) = RESULT(J) + VALUES(J,I)
              ABSERR(J) = ABSERR(J) + ERRORS(J,I)
505       CONTINUE
510   CONTINUE
C
C   Compute correct sign on the integral.
C
      DO 600 J = 1,NUMFUN
          RESULT(J) = RESULT(J)*INTSGN
600   CONTINUE
      NSUB = SBRGNS
      RETURN
C
C***END DADHRE
C
      END
