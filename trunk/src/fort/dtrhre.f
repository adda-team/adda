      SUBROUTINE DTRHRE(DVFLAG,NDIM,NUMFUN,SBRGNS,VALUES,ERRORS,CENTRS,
     +                  HWIDTS,GREATE,ERROR,VALUE,CENTER,HWIDTH,DIR)
C***BEGIN PROLOGUE DTRHRE
C***PURPOSE DTRHRE maintains a heap of subregions.
C***DESCRIPTION DTRHRE maintains a heap of subregions.
C            The subregions are ordered according to the size
C            of the greatest error estimates of each subregion(GREATE).
C
C   PARAMETERS
C
C     DVFLAG Integer.
C            If DVFLAG = 1, we remove the subregion with
C            greatest error from the heap.
C            If DVFLAG = 2, we insert a new subregion in the heap.
C     NDIM   Integer.
C            Number of variables.
C     NUMFUN Integer.
C            Number of components of the integral.
C     SBRGNS Integer.
C            Number of subregions in the heap.
C     VALUES Real array of dimension (NUMFUN,SBRGNS).
C            Used to store estimated values of the integrals
C            over the subregions.
C     ERRORS Real array of dimension (NUMFUN,SBRGNS).
C            Used to store the corresponding estimated errors.
C     CENTRS Real array of dimension (NDIM,SBRGNS).
C            Used to store the center limits of the stored subregions.
C     HWIDTS Real array of dimension (NDIM,SBRGNS).
C            Used to store the hwidth limits of the stored subregions.
C     GREATE Real array of dimension SBRGNS.
C            Used to store the greatest estimated errors in
C            all subregions.
C     ERROR  Real array of dimension NUMFUN.
C            Used as intermediate storage for the error of a subregion.
C     VALUE  Real array of dimension NUMFUN.
C            Used as intermediate storage for the estimate
C            of the integral over a subregion.
C     CENTER Real array of dimension NDIM.
C            Used as intermediate storage for the center of
C            the subregion.
C     HWIDTH Real array of dimension NDIM.
C            Used as intermediate storage for the half width of
C            the subregion.
C     DIR    Integer array of dimension SBRGNS.
C            DIR is used to store the directions for
C            further subdivision.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE DTRHRE
C
C   Global variables.
C
      INTEGER DVFLAG,NDIM,NUMFUN,SBRGNS
      DOUBLE PRECISION VALUES(NUMFUN,*),ERRORS(NUMFUN,*)
      DOUBLE PRECISION CENTRS(NDIM,*)
      DOUBLE PRECISION HWIDTS(NDIM,*)
      DOUBLE PRECISION GREATE(*)
      DOUBLE PRECISION ERROR(NUMFUN),VALUE(NUMFUN)
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM)
      DOUBLE PRECISION DIR(*)
C
C   Local variables.
C
C   GREAT  is used as intermediate storage for the greatest error of a
C          subregion.
C   DIRECT is used as intermediate storage for the direction of further
C          subdivision.
C   SUBRGN Position of child/parent subregion in the heap.
C   SUBTMP Position of parent/child subregion in the heap.
C
      INTEGER J,SUBRGN,SUBTMP
      DOUBLE PRECISION GREAT,DIRECT
C
C***FIRST EXECUTABLE STATEMENT DTRHRE
C
C   Save values to be stored in their correct place in the heap.
C
      GREAT = GREATE(SBRGNS)
      DIRECT = DIR(SBRGNS)
      DO 5 J = 1,NUMFUN
          ERROR(J) = ERRORS(J,SBRGNS)
          VALUE(J) = VALUES(J,SBRGNS)
5     CONTINUE
      DO 10 J = 1,NDIM
          CENTER(J) = CENTRS(J,SBRGNS)
          HWIDTH(J) = HWIDTS(J,SBRGNS)
10    CONTINUE
C
C    If DVFLAG = 1, we will remove the region
C    with greatest estimated error from the heap.
C
      IF (DVFLAG.EQ.1) THEN
          SBRGNS = SBRGNS - 1
          SUBRGN = 1
20        SUBTMP = 2*SUBRGN
          IF (SUBTMP.LE.SBRGNS) THEN
              IF (SUBTMP.NE.SBRGNS) THEN
C
C   Find max. of left and right child.
C
                  IF (GREATE(SUBTMP).LT.GREATE(SUBTMP+1)) THEN
                      SUBTMP = SUBTMP + 1
                  END IF
              END IF
C
C   Compare max.child with parent.
C   If parent is max., then done.
C
              IF (GREAT.LT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp up the heap.
C
                  GREATE(SUBRGN) = GREATE(SUBTMP)
                  DO 25 J = 1,NUMFUN
                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
25                CONTINUE
                  DIR(SUBRGN) = DIR(SUBTMP)
                  DO 30 J = 1,NDIM
                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
30                CONTINUE
                  SUBRGN = SUBTMP
                  GO TO 20
              END IF
          END IF
      ELSE IF (DVFLAG.EQ.2) THEN
C
C   If DVFLAG = 2, then insert new region in the heap.
C
          SUBRGN = SBRGNS
40        SUBTMP = SUBRGN/2
          IF (SUBTMP.GE.1) THEN
C
C   Compare max.child with parent.
C   If parent is max, then done.
C
              IF (GREAT.GT.GREATE(SUBTMP)) THEN
C
C   Move the values at position subtmp down the heap.
C
                  GREATE(SUBRGN) = GREATE(SUBTMP)
                  DO 45 J = 1,NUMFUN
                      ERRORS(J,SUBRGN) = ERRORS(J,SUBTMP)
                      VALUES(J,SUBRGN) = VALUES(J,SUBTMP)
45                CONTINUE
                  DIR(SUBRGN) = DIR(SUBTMP)
                  DO 50 J = 1,NDIM
                      CENTRS(J,SUBRGN) = CENTRS(J,SUBTMP)
                      HWIDTS(J,SUBRGN) = HWIDTS(J,SUBTMP)
50                CONTINUE
                  SUBRGN = SUBTMP
                  GO TO 40
              END IF
          END IF
      END IF
C
C    Insert the saved values in their correct places.
C
      IF (SBRGNS.GT.0) THEN
          GREATE(SUBRGN) = GREAT
          DO 55 J = 1,NUMFUN
              ERRORS(J,SUBRGN) = ERROR(J)
              VALUES(J,SUBRGN) = VALUE(J)
55        CONTINUE
          DIR(SUBRGN) = DIRECT
          DO 60 J = 1,NDIM
              CENTRS(J,SUBRGN) = CENTER(J)
              HWIDTS(J,SUBRGN) = HWIDTH(J)
60        CONTINUE
      END IF
C
C***END DTRHRE
C
      RETURN
      END
