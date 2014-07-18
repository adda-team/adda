      SUBROUTINE DFSHRE(NDIM,CENTER,HWIDTH,X,G,NUMFUN,FUNSUB,FULSMS,
     +                  FUNVLS)
C***BEGIN PROLOGUE DFSHRE
C***KEYWORDS fully symmetric sum
C***PURPOSE  To compute fully symmetric basic rule sums
C***AUTHOR   Alan Genz, Computer Science Department, Washington
C            State University, Pullman, WA 99163-1210 USA
C***LAST MODIFICATION 88-04-08
C***DESCRIPTION DFSHRE computes a fully symmetric sum for a vector
C            of integrand values over a hyper-rectangular region.
C            The sum is fully symmetric with respect to the center of
C            the region and is taken over all sign changes and
C            permutations of the generators for the sum.
C
C   ON ENTRY
C
C   NDIM   Integer.
C          Number of variables.
C   CENTER Real array of dimension NDIM.
C          The coordinates for the center of the region.
C   HWIDTH Real Array of dimension NDIM.
C          HWIDTH(I) is half of the width of dimension I of the region.
C   X      Real Array of dimension NDIM.
C          A work array.
C   G      Real Array of dimension NDIM.
C          The generators for the fully symmetric sum. These MUST BE
C          non-negative and non-increasing.
C   NUMFUN Integer.
C          Number of components for the vector integrand.
C   FUNSUB Externally declared subroutine.
C          For computing the components of the integrand at a point X.
C          It must have parameters (NDIM, X, NUMFUN, FUNVLS).
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
C   ON RETURN
C
C   FULSMS Real array of dimension NUMFUN.
C          The values for the fully symmetric sums for each component
C          of the integrand.
C   FUNVLS Real array of dimension NUMFUN.
C          A work array.
C
C***ROUTINES CALLED: FUNSUB
C
C***END PROLOGUE DFSHRE
C
C   Global variables.
C
      EXTERNAL FUNSUB
      INTEGER NDIM,NUMFUN
      DOUBLE PRECISION CENTER(NDIM),HWIDTH(NDIM),X(NDIM),G(NDIM),
     +                 FULSMS(NUMFUN),FUNVLS(NUMFUN)
C
C   Local variables.
C
      INTEGER IXCHNG,LXCHNG,I,J,L
      DOUBLE PRECISION GL,GI
C
C***FIRST EXECUTABLE STATEMENT DFSHRE
C
      DO 10 J = 1,NUMFUN
          FULSMS(J) = 0
10    CONTINUE
C
C     Compute centrally symmetric sum for permutation of G
C
20    DO 30 I = 1,NDIM
          X(I) = CENTER(I) + G(I)*HWIDTH(I)
30    CONTINUE
40    CALL FUNSUB(NDIM,X,NUMFUN,FUNVLS)
      DO 50 J = 1,NUMFUN
          FULSMS(J) = FULSMS(J) + FUNVLS(J)
50    CONTINUE
      DO 60 I = 1,NDIM
          G(I) = - G(I)
          X(I) = CENTER(I) + G(I)*HWIDTH(I)
          IF (G(I).LT.0) GO TO 40
60    CONTINUE
C
C       Find next distinct permutation of G and loop back for next sum.
C       Permutations are generated in reverse lexicographic order.
C
      DO 80 I = 2,NDIM
          IF (G(I-1).GT.G(I)) THEN
              GI = G(I)
              IXCHNG = I - 1
              DO 70 L = 1, (I-1)/2
                  GL = G(L)
                  G(L) = G(I-L)
                  G(I-L) = GL
                  IF (GL.LE.GI) IXCHNG = IXCHNG - 1
                  IF (G(L).GT.GI) LXCHNG = L
70            CONTINUE
              IF (G(IXCHNG).LE.GI) IXCHNG = LXCHNG
              G(I) = G(IXCHNG)
              G(IXCHNG) = GI
              GO TO 20
          END IF
80    CONTINUE
C
C     Restore original order to generators
C
      DO 90 I = 1,NDIM/2
          GI = G(I)
          G(I) = G(NDIM-I+1)
          G(NDIM-I+1) = GI
90    CONTINUE
C
C***END DFSHRE
C
      END
