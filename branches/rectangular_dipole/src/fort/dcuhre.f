      SUBROUTINE DCUHRE(NDIM,NUMFUN,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,
     +                  EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,
     +                  WORK)
C***BEGIN PROLOGUE DCUHRE
C***DATE WRITTEN   900116   (YYMMDD)
C***REVISION DATE  900116   (YYMMDD)
C***CATEGORY NO. H2B1A1
C***AUTHOR
C            Jarle Berntsen, The Computing Centre,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544055
C            Email..  jarle@eik.ii.uib.no
C            Terje O. Espelid, Department of Informatics,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, Norway
C            Phone..  47-5-544180
C            Email..  terje@eik.ii.uib.no
C            Alan Genz, Computer Science Department, Washington State
C            University, Pullman, WA 99163-1210, USA
C            Phone.. 509-335-2131
C            Email..  acg@cs2.cs.wsu.edu
C***KEYWORDS automatic multidimensional integrator,
C            n-dimensional hyper-rectangles,
C            general purpose, global adaptive
C***PURPOSE  The routine calculates an approximation to a given
C            vector of definite integrals
C
C      B(1) B(2)     B(NDIM)
C     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1),
C      A(1) A(2)     A(NDIM)  1  2      NUMFUN
C
C       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN.
C              I   I  1  2      NDIM
C
C            hopefully satisfying for each component of I the following
C            claim for accuracy:
C            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K)))
C***DESCRIPTION Computation of integrals over hyper-rectangular
C            regions.
C            DCUHRE is a driver for the integration routine
C            DADHRE, which repeatedly subdivides the region
C            of integration and estimates the integrals and the
C            errors over the subregions with greatest
C            estimated errors until the error request
C            is met or MAXPTS function evaluations have been used.
C
C            For NDIM = 2 the default integration rule is of
C            degree 13 and uses 65 evaluation points.
C            For NDIM = 3 the default integration rule is of
C            degree 11 and uses 127 evaluation points.
C            For NDIM greater then 3 the default integration rule
C            is of degree 9 and uses NUM evaluation points where
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            The degree 9 rule may also be applied for NDIM = 2
C            and NDIM = 3.
C            A rule of degree 7 is available in all dimensions.
C            The number of evaluation
C            points used by the degree 7 rule is
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C
C            When DCUHRE computes estimates to a vector of
C            integrals, all components of the vector are given
C            the same treatment. That is, I(F ) and I(F ) for
C                                            J         K
C            J not equal to K, are estimated with the same
C            subdivision of the region of integration.
C            For integrals with enough similarity, we may save
C            time by applying DCUHRE to all integrands in one call.
C            For integrals that vary continuously as functions of
C            some parameter, the estimates produced by DCUHRE will
C            also vary continuously when the same subdivision is
C            applied to all components. This will generally not be
C            the case when the different components are given
C            separate treatment.
C
C            On the other hand this feature should be used with
C            caution when the different components of the integrals
C            require clearly different subdivisions.
C
C   ON ENTRY
C
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <=  15.
C     NUMFUN Integer.
C            Number of components of the integral.
C     A      Real array of dimension NDIM.
C            Lower limits of integration.
C     B      Real array of dimension NDIM.
C            Upper limits of integration.
C     MINPTS Integer.
C            Minimum number of function evaluations.
C     MAXPTS Integer.
C            Maximum number of function evaluations.
C            The number of function evaluations over each subregion
C            is NUM.
C            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then
C              NUM = 65
C            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then
C              NUM = 127
C            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then
C              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) +
C                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM
C            Elseif (KEY = 4) Then
C              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM
C            MAXPTS >= 3*NUM and MAXPTS >= MINPTS
C            For 3 < NDIM < 13 the minimum values for MAXPTS are:
C             NDIM =    4   5   6    7    8    9    10   11    12
C            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235
C            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299
C     FUNSUB Externally declared subroutine for computing
C            all components of the integrand at the given
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
C     NW     Integer.
C            Defines the length of the working array WORK.
C            Let MAXSUB denote the maximum allowed number of subregions
C            for the given values of MAXPTS, KEY and NDIM.
C            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1
C            NW should be greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1
C            For efficient execution on parallel computers
C            NW should be chosen greater or equal to
C            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1
C            where MDIV is the number of subregions that are divided in
C            each subdivision step.
C            MDIV is default set internally in DCUHRE equal to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C            In this case the only parameters for DCUHRE that may
C            be changed (with respect to the previous call of DCUHRE)
C            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR.
C
C   ON RETURN
C
C     RESULT Real array of dimension NUMFUN.
C            Approximations to all components of the integral.
C     ABSERR Real array of dimension NUMFUN.
C            Estimates of absolute errors.
C     NEVAL  Integer.
C            Number of function evaluations used by DCUHRE.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or
C              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less
C              function evaluations for all values of K,
C              1 <= K <= NUMFUN .
C            IFAIL = 1 if MAXPTS was too small for DCUHRE
C              to obtain the required accuracy. In this case DCUHRE
C              returns values of RESULT with estimated absolute
C              errors ABSERR.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than 15.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN is less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     WORK   Real array of dimension NW.
C            Used as working storage.
C            WORK(NW) = NSUB, the number of subregions in the data
C            structure.
C            Let WRKSUB=(NW-1-17*NUMFUN*MDIV)/(2*NDIM+2*NUMFUN+2)
C            WORK(1),...,WORK(NUMFUN*WRKSUB) contain
C              the estimated components of the integrals over the
C              subregions.
C            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain
C              the estimated errors over the subregions.
C            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB+NDIM*
C              WRKSUB) contain the centers of the subregions.
C            WORK(2*NUMFUN*WRKSUB+NDIM*WRKSUB+1),...,WORK((2*NUMFUN+
C              NDIM)*WRKSUB+NDIM*WRKSUB) contain subregion half widths.
C            WORK(2*NUMFUN*WRKSUB+2*NDIM*WRKSUB+1),...,WORK(2*NUMFUN*
C              WRKSUB+2*NDIM*WRKSUB+WRKSUB) contain the greatest errors
C              in each subregion.
C            WORK((2*NUMFUN+2*NDIM+1)*WRKSUB+1),...,WORK((2*NUMFUN+
C              2*NDIM+1)*WRKSUB+WRKSUB) contain the direction of
C              subdivision in each subregion.
C            WORK(2*(NDIM+NUMFUN+1)*WRKSUB),...,WORK(2*(NDIM+NUMFUN+1)*
C              WRKSUB+ 17*MDIV*NUMFUN) is used as temporary
C              storage in DADHRE.
C
C
C        DCUHRE Example Test Program
C
C
C   DTEST1 is a simple test driver for DCUHRE.
C
C   Output produced on a SUN 3/50.
c
C       DCUHRE TEST RESULTS
C
C    FTEST CALLS = 3549, IFAIL =  0
C   N   ESTIMATED ERROR   INTEGRAL
C   1     0.00000010     0.13850818
C   2     0.00000013     0.06369469
C   3     0.00000874     0.05861748
C   4     0.00000021     0.05407034
C   5     0.00000019     0.05005614
C   6     0.00000009     0.04654608
C
C     PROGRAM DTEST1
C     EXTERNAL FTEST
C     INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW
C     PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1)
C     DOUBLE PRECISION A(NDIM), B(NDIM), WRKSTR(NW)
C     DOUBLE PRECISION ABSEST(NF), FINEST(NF), ABSREQ, RELREQ
C     DO 10 N = 1,NDIM
C        A(N) = 0
C        B(N) = 1
C  10 CONTINUE
C     MINCLS = 0
C     MAXCLS = 10000
C     KEY = 0
C     ABSREQ = 0
C     RELREQ = 1E-3
C     CALL DCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ,
C    * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR)
C     PRINT 9999, NEVAL, IFAIL
C9999 FORMAT (8X, 'DCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4,
C    * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR   INTEGRAL')
C     DO 20 N = 1,NF
C        PRINT 9998, N, ABSEST(N), FINEST(N)
C9998    FORMAT (3X, I2, 2F15.8)
C  20 CONTINUE
C     END
C     SUBROUTINE FTEST(NDIM, Z, NFUN, F)
C     INTEGER N, NDIM, NFUN
C     DOUBLE PRECISION Z(NDIM), F(NFUN), SUM
C     SUM = 0
C     DO 10 N = 1,NDIM
C        SUM = SUM + N*Z(N)**2
C  10 CONTINUE
C     F(1) = EXP(-SUM/2)
C     DO 20 N = 1,NDIM
C        F(N+1) = Z(N)*F(1)
C  20 CONTINUE
C     END
C
C***LONG DESCRIPTION
C
C   The information for each subregion is contained in the
C   data structure WORK.
C   When passed on to DADHRE, WORK is split into eight
C   arrays VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR,
C   OLDRES and WORK.
C   VALUES contains the estimated values of the integrals.
C   ERRORS contains the estimated errors.
C   CENTRS contains the centers of the subregions.
C   HWIDTS contains the half widths of the subregions.
C   GREATE contains the greatest estimated error for each subregion.
C   DIR    contains the directions for further subdivision.
C   OLDRES and WORK are used as work arrays in DADHRE.
C
C   The data structures for the subregions are in DADHRE organized
C   as a heap, and the size of GREATE(I) defines the position of
C   region I in the heap. The heap is maintained by the program
C   DTRHRE.
C
C   The subroutine DADHRE is written for efficient execution on shared
C   memory parallel computer. On a computer with NPROC processors we wil
C   in each subdivision step divide MDIV regions, where MDIV is
C   chosen such that MOD(2*MDIV,NPROC) = 0, in totally 2*MDIV new region
C   Each processor will then compute estimates of the integrals and erro
C   over 2*MDIV/NPROC subregions in each subdivision step.
C   The subroutine for estimating the integral and the error over
C   each subregion, DRLHRE, uses WORK2 as a work array.
C   We must make sure that each processor writes its results to
C   separate parts of the memory, and therefore the sizes of WORK and
C   WORK2 are functions of MDIV.
C   In order to achieve parallel processing of subregions, compiler
C   directives should be placed in front of the DO 200
C   loop in DADHRE on machines like Alliant and CRAY.
C
C***REFERENCES
C   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm
C   for the Approximate Calculation of Multiple Integrals,
C   To be published.
C
C   J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive
C   Multidimensional Integration Routine for a Vector of
C   Integrals, To be published.
C
C***ROUTINES CALLED DCHHRE,DADHRE
C***END PROLOGUE DCUHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN,MINPTS,MAXPTS,KEY,NW,RESTAR
      INTEGER NEVAL,IFAIL
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
      DOUBLE PRECISION RESULT(NUMFUN),ABSERR(NUMFUN),WORK(NW)
C
C   Local variables.
C
C   MDIV   Integer.
C          MDIV is the number of subregions that are divided in
C          each subdivision step in DADHRE.
C          MDIV is chosen default to 1.
C          For efficient execution on parallel computers
C          with NPROC processors MDIV should be set equal to
C          the smallest integer such that MOD(2*MDIV,NPROC) = 0.
C   MAXDIM Integer.
C          The maximum allowed value of NDIM.
C   MAXWT  Integer. The maximum number of weights used by the
C          integration rule.
C   WTLENG Integer.
C          The number of generators used by the selected rule.
C   WORK2  Real work space. The length
C          depends on the parameters MDIV,MAXDIM and MAXWT.
C   MAXSUB Integer.
C          The maximum allowed number of subdivisions
C          for the given values of KEY, NDIM and MAXPTS.
C   MINSUB Integer.
C          The minimum allowed number of subregions for the given
C          values of MINPTS, KEY and NDIM.
C   WRKSUB Integer.
C          The maximum allowed number of subregions as a function
C          of NW, NUMFUN, NDIM and MDIV. This determines the length
C          of the main work arrays.
C   NUM    Integer. The number of integrand evaluations needed
C          over each subregion.
C
      INTEGER MDIV,MAXWT,WTLENG,MAXDIM,LENW2,MAXSUB,MINSUB
      INTEGER NUM,NSUB,LENW,KEYF
      PARAMETER (MDIV=1)
      PARAMETER (MAXDIM=15)
      PARAMETER (MAXWT=14)
      PARAMETER (LENW2=2*MDIV*MAXDIM* (MAXWT+1)+12*MAXWT+2*MAXDIM)
      INTEGER WRKSUB,I1,I2,I3,I4,I5,I6,I7,I8,K1,K2,K3,K4,K5,K6,K7,K8
      DOUBLE PRECISION WORK2(LENW2)
C
C***FIRST EXECUTABLE STATEMENT DCUHRE
C
C   Compute NUM, WTLENG, MAXSUB and MINSUB,
C   and check the input parameters.
C
      CALL DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,EPSABS,
     +            EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,KEYF,
     +            IFAIL,WTLENG)
      WRKSUB = (NW - 1 - 17*MDIV*NUMFUN)/(2*NDIM + 2*NUMFUN + 2)
      IF (IFAIL.NE.0) THEN
          GO TO 999
      END IF
C
C   Split up the work space.
C
      I1 = 1
      I2 = I1 + WRKSUB*NUMFUN
      I3 = I2 + WRKSUB*NUMFUN
      I4 = I3 + WRKSUB*NDIM
      I5 = I4 + WRKSUB*NDIM
      I6 = I5 + WRKSUB
      I7 = I6 + WRKSUB
      I8 = I7 + NUMFUN*MDIV
      K1 = 1
      K2 = K1 + 2*MDIV*WTLENG*NDIM
      K3 = K2 + WTLENG*5
      K4 = K3 + WTLENG
      K5 = K4 + NDIM
      K6 = K5 + NDIM
      K7 = K6 + 2*MDIV*NDIM
      K8 = K7 + 3*WTLENG
C
C   On restart runs the number of subregions from the
C   previous call is assigned to NSUB.
C
      IF (RESTAR.EQ.1) THEN
          NSUB = WORK(NW)
      END IF
C
C   Compute the size of the temporary work space needed in DADHRE.
C
      LENW = 16*MDIV*NUMFUN
      CALL DADHRE(NDIM,NUMFUN,MDIV,A,B,MINSUB,MAXSUB,FUNSUB,EPSABS,
     +            EPSREL,KEYF,RESTAR,NUM,LENW,WTLENG,
     +            RESULT,ABSERR,NEVAL,NSUB,IFAIL,WORK(I1),WORK(I2),
     +            WORK(I3),WORK(I4),WORK(I5),WORK(I6),WORK(I7),WORK(I8),
     +            WORK2(K1),WORK2(K2),WORK2(K3),WORK2(K4),WORK2(K5),
     +            WORK2(K6),WORK2(K7),WORK2(K8))
      WORK(NW) = NSUB
999   RETURN
C
C***END DCUHRE
C
      END
