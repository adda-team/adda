      SUBROUTINE DCHHRE(MAXDIM,NDIM,NUMFUN,MDIV,A,B,MINPTS,MAXPTS,
     +                  EPSABS,EPSREL,KEY,NW,RESTAR,NUM,MAXSUB,MINSUB,
     +                  KEYF,IFAIL,WTLENG)
C***BEGIN PROLOGUE DCHHRE
C***PURPOSE  DCHHRE checks the validity of the
C            input parameters to DCUHRE.
C***DESCRIPTION
C            DCHHRE computes NUM, MAXSUB, MINSUB, KEYF, WTLENG and
C            IFAIL as functions of the input parameters to DCUHRE,
C            and checks the validity of the input parameters to DCUHRE.
C
C   ON ENTRY
C
C     MAXDIM Integer.
C            The maximum allowed number of dimensions.
C     NDIM   Integer.
C            Number of variables. 1 < NDIM <= MAXDIM.
C     NUMFUN Integer.
C            Number of components of the integral.
C     MDIV   Integer.
C            MDIV is the number of subregions that are divided in
C            each subdivision step in DADHRE.
C            MDIV is chosen default to 1.
C            For efficient execution on parallel computers
C            with NPROC processors MDIV should be set equal to
C            the smallest integer such that MOD(2*MDIV,NPROC) = 0.
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
C     RESTAR Integer.
C            If RESTAR = 0, this is the first attempt to compute
C            the integral.
C            If RESTAR = 1, then we restart a previous attempt.
C
C   ON RETURN
C
C     NUM    Integer.
C            The number of function evaluations over each subregion.
C     MAXSUB Integer.
C            The maximum allowed number of subregions for the
C            given values of MAXPTS, KEY and NDIM.
C     MINSUB Integer.
C            The minimum allowed number of subregions for the given
C            values of MINPTS, KEY and NDIM.
C     IFAIL  Integer.
C            IFAIL = 0 for normal exit.
C            IFAIL = 2 if KEY is less than 0 or KEY greater than 4.
C            IFAIL = 3 if NDIM is less than 2 or NDIM greater than
C                      MAXDIM.
C            IFAIL = 4 if KEY = 1 and NDIM not equal to 2.
C            IFAIL = 5 if KEY = 2 and NDIM not equal to 3.
C            IFAIL = 6 if NUMFUN less than 1.
C            IFAIL = 7 if volume of region of integration is zero.
C            IFAIL = 8 if MAXPTS is less than 3*NUM.
C            IFAIL = 9 if MAXPTS is less than MINPTS.
C            IFAIL = 10 if EPSABS < 0 and EPSREL < 0.
C            IFAIL = 11 if NW is too small.
C            IFAIL = 12 if unlegal RESTAR.
C     KEYF   Integer.
C            Key to selected integration rule.
C     WTLENG Integer.
C            The number of generators of the chosen integration rule.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DCHHRE
C
C   Global variables.
C
      INTEGER NDIM,NUMFUN,MDIV,MINPTS,MAXPTS,KEY,NW,MINSUB,MAXSUB
      INTEGER RESTAR,NUM,KEYF,IFAIL,MAXDIM,WTLENG
      DOUBLE PRECISION A(NDIM),B(NDIM),EPSABS,EPSREL
C
C   Local variables.
C
      INTEGER LIMIT,J
C
C***FIRST EXECUTABLE STATEMENT DCHHRE
C
      IFAIL = 0
C
C   Check on legal KEY.
C
      IF (KEY.LT.0 .OR. KEY.GT.4) THEN
          IFAIL = 2
          GO TO 999
      END IF
C
C   Check on legal NDIM.
C
      IF (NDIM.LT.2 .OR. NDIM.GT.MAXDIM) THEN
          IFAIL = 3
          GO TO 999
      END IF
C
C   For KEY = 1, NDIM must be equal to 2.
C
      IF (KEY.EQ.1 .AND. NDIM.NE.2) THEN
          IFAIL = 4
          GO TO 999
      END IF
C
C   For KEY = 2, NDIM must be equal to 3.
C
      IF (KEY.EQ.2 .AND. NDIM.NE.3) THEN
          IFAIL = 5
          GO TO 999
      END IF
C
C   For KEY = 0, we point at the selected integration rule.
C
      IF (KEY.EQ.0) THEN
          IF (NDIM.EQ.2) THEN
              KEYF = 1
          ELSE IF (NDIM.EQ.3) THEN
              KEYF = 2
          ELSE
              KEYF = 3
          ENDIF
      ELSE
          KEYF = KEY
      ENDIF
C
C   Compute NUM and WTLENG as a function of KEYF and NDIM.
C
      IF (KEYF.EQ.1) THEN
          NUM = 65
          WTLENG = 14
      ELSE IF (KEYF.EQ.2) THEN
          NUM = 127
          WTLENG = 13
      ELSE IF (KEYF.EQ.3) THEN
          NUM = 1 + 4*2*NDIM + 2*NDIM* (NDIM-1) + 4*NDIM* (NDIM-1) +
     +          4*NDIM* (NDIM-1)* (NDIM-2)/3 + 2**NDIM
          WTLENG = 9
          IF (NDIM.EQ.2) WTLENG = 8
      ELSE IF (KEYF.EQ.4) THEN
          NUM = 1 + 3*2*NDIM + 2*NDIM* (NDIM-1) + 2**NDIM
          WTLENG = 6
      END IF
C
C   Compute MAXSUB.
C
      MAXSUB = (MAXPTS-NUM)/ (2*NUM) + 1
C
C   Compute MINSUB.
C
      MINSUB = (MINPTS-NUM)/ (2*NUM) + 1
      IF (MOD(MINPTS-NUM,2*NUM).NE.0) THEN
          MINSUB = MINSUB + 1
      END IF
      MINSUB = MAX(2,MINSUB)
C
C   Check on positive NUMFUN.
C
      IF (NUMFUN.LT.1) THEN
          IFAIL = 6
          GO TO 999
      END IF
C
C   Check on legal upper and lower limits of integration.
C
      DO 10 J = 1,NDIM
          IF (A(J)-B(J).EQ.0) THEN
              IFAIL = 7
              GO TO 999
          END IF
10    CONTINUE
C
C   Check on MAXPTS < 3*NUM.
C
      IF (MAXPTS.LT.3*NUM) THEN
          IFAIL = 8
          GO TO 999
      END IF
C
C   Check on MAXPTS >= MINPTS.
C
      IF (MAXPTS.LT.MINPTS) THEN
          IFAIL = 9
          GO TO 999
      END IF
C
C   Check on legal accuracy requests.
C
      IF (EPSABS.LT.0 .AND. EPSREL.LT.0) THEN
          IFAIL = 10
          GO TO 999
      END IF
C
C   Check on big enough double precision workspace.
C
      LIMIT = MAXSUB* (2*NDIM+2*NUMFUN+2) + 17*MDIV*NUMFUN + 1
      IF (NW.LT.LIMIT) THEN
          IFAIL = 11
          GO TO 999
      END IF
C
C    Check on legal RESTAR.
C
      IF (RESTAR.NE.0 .AND. RESTAR.NE.1) THEN
          IFAIL = 12
          GO TO 999
      END IF
999   RETURN
C
C***END DCHHRE
C
      END
