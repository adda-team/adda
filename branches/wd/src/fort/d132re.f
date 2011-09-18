      SUBROUTINE D132RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D132RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE D132RE computes abscissas and weights of a 2 dimensional
C            integration rule of degree 13.
C            Two null rules of degree 11, one null rule of degree 9
C            and one null rule of degree 7 to be used in error
C            estimation are also computed.
C ***DESCRIPTION D132RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to
C            G and W. The heuristic error coefficients ERRCOF
C            will also be assigned.
C
C
C   ON ENTRY
C
C     WTLENG Integer.
C            The number of weights in each of the rules.
C
C   ON RETURN
C
C     W      Real array of dimension (5,WTLENG).
C            The weights for the basic and null rules.
C            W(1,1),...,W(1,WTLENG) are weights for the basic rule.
C            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights.
C     G      Real array of dimension (NDIM,WTLENG).
C            The fully symmetric sum generators for the rules.
C            G(1,J),...,G(NDIM,J) are the generators for the points
C            associated with the the Jth weights.
C     ERRCOF Real array of dimension 6.
C            Heuristic error coefficients that are used in the
C            error estimation in BASRUL.
C     RULPTS Real array of dimension WTLENG.
C            The number of points produced by each generator.
C***REFERENCES S.Eriksen,
C              Thesis of the degree cand.scient, Dept. of Informatics,
C              Univ. of Bergen,Norway, 1984.
C
C***ROUTINES CALLED-NONE
C***END PROLOGUE D132RE
C
C   Global variables
C
      INTEGER WTLENG
      DOUBLE PRECISION W(5,WTLENG),G(2,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local variables.
C
      INTEGER I,J
      DOUBLE PRECISION DIM2G(16)
      DOUBLE PRECISION DIM2W(14,5)
C
      DATA (DIM2G(I),I=1,16)/0.2517129343453109D+00,
     +     0.7013933644534266D+00,0.9590960631619962D+00,
     +     0.9956010478552127D+00,0.5000000000000000D+00,
     +     0.1594544658297559D+00,0.3808991135940188D+00,
     +     0.6582769255267192D+00,0.8761473165029315D+00,
     +     0.9982431840531980D+00,0.9790222658168462D+00,
     +     0.6492284325645389D+00,0.8727421201131239D+00,
     +     0.3582614645881228D+00,0.5666666666666666D+00,
     +     0.2077777777777778D+00/
C
      DATA (DIM2W(I,1),I=1,14)/0.3379692360134460D-01,
     +     0.9508589607597761D-01,0.1176006468056962D+00,
     +     0.2657774586326950D-01,0.1701441770200640D-01,
     +     0.0000000000000000D+00,0.1626593098637410D-01,
     +     0.1344892658526199D+00,0.1328032165460149D+00,
     +     0.5637474769991870D-01,0.3908279081310500D-02,
     +     0.3012798777432150D-01,0.1030873234689166D+00,
     +     0.6250000000000000D-01/
C
      DATA (DIM2W(I,2),I=1,14)/0.3213775489050763D+00,
     +     - .1767341636743844D+00,0.7347600537466072D-01,
     +     - .3638022004364754D-01,0.2125297922098712D-01,
     +     0.1460984204026913D+00,0.1747613286152099D-01,
     +     0.1444954045641582D+00,0.1307687976001325D-03,
     +     0.5380992313941161D-03,0.1042259576889814D-03,
     +     - .1401152865045733D-02,0.8041788181514763D-02,
     +     - .1420416552759383D+00/
C
      DATA (DIM2W(I,3),I=1,14)/0.3372900883288987D+00,
     +     - .1644903060344491D+00,0.7707849911634622D-01,
     +     - .3804478358506310D-01,0.2223559940380806D-01,
     +     0.1480693879765931D+00,0.4467143702185814D-05,
     +     0.1508944767074130D+00,0.3647200107516215D-04,
     +     0.5777198999013880D-03,0.1041757313688177D-03,
     +     - .1452822267047819D-02,0.8338339968783705D-02,
     +     - .1472796329231960D+00/
C
      DATA (DIM2W(I,4),I=1,14)/ - .8264123822525677D+00,
     +     0.3065838614094360D+00,0.2389292538329435D-02,
     +     - .1343024157997222D+00,0.8833366840533900D-01,
     +     0.0000000000000000D+00,0.9786283074168292D-03,
     +     - .1319227889147519D+00,0.7990012200150630D-02,
     +     0.3391747079760626D-02,0.2294915718283264D-02,
     +     - .1358584986119197D-01,0.4025866859057809D-01,
     +     0.3760268580063992D-02/
C
      DATA (DIM2W(I,5),I=1,14)/0.6539094339575232D+00,
     +     - .2041614154424632D+00, - .1746981515794990D+00,
     +     0.3937939671417803D-01,0.6974520545933992D-02,
     +     0.0000000000000000D+00,0.6667702171778258D-02,
     +     0.5512960621544304D-01,0.5443846381278607D-01,
     +     0.2310903863953934D-01,0.1506937747477189D-01,
     +     - .6057021648901890D-01,0.4225737654686337D-01,
     +     0.2561989142123099D-01/
C
C***FIRST EXECUTABLE STATEMENT D132RE
C
C   Assign values to W.
C
      DO 10 I = 1,14
          DO 10 J = 1,5
              W(J,I) = DIM2W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
      DO 20 I = 1,2
          DO 20 J = 1,14
              G(I,J) = 0
20    CONTINUE
      G(1,2) = DIM2G(1)
      G(1,3) = DIM2G(2)
      G(1,4) = DIM2G(3)
      G(1,5) = DIM2G(4)
      G(1,6) = DIM2G(5)
      G(1,7) = DIM2G(6)
      G(2,7) = G(1,7)
      G(1,8) = DIM2G(7)
      G(2,8) = G(1,8)
      G(1,9) = DIM2G(8)
      G(2,9) = G(1,9)
      G(1,10) = DIM2G(9)
      G(2,10) = G(1,10)
      G(1,11) = DIM2G(10)
      G(2,11) = G(1,11)
      G(1,12) = DIM2G(11)
      G(2,12) = DIM2G(12)
      G(1,13) = DIM2G(13)
      G(2,13) = DIM2G(14)
      G(1,14) = DIM2G(15)
      G(2,14) = DIM2G(16)
C
C   Assign values to RULPTS.
C
      RULPTS(1) = 1
      DO 30 I = 2,11
          RULPTS(I) = 4
30    CONTINUE
      RULPTS(12) = 8
      RULPTS(13) = 8
      RULPTS(14) = 8
C
C   Assign values to ERRCOF.
C
      ERRCOF(1) = 10
      ERRCOF(2) = 10
      ERRCOF(3) = 1.
      ERRCOF(4) = 5.
      ERRCOF(5) = 0.5
      ERRCOF(6) = 0.25
C
C***END D132RE
C
      RETURN
      END
