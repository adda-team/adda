      SUBROUTINE D113RE(WTLENG,W,G,ERRCOF,RULPTS)
C***BEGIN PROLOGUE D113RE
C***AUTHOR   Jarle Berntsen, EDB-senteret,
C            University of Bergen, Thormohlens gt. 55,
C            N-5008 Bergen, NORWAY
C***PURPOSE D113RE computes abscissas and weights of a 3 dimensional
C            integration rule of degree 11.
C            Two null rules of degree 9, one null rule of degree 7
C            and one null rule of degree 5 to be used in error
C            estimation are also computed.
C***DESCRIPTION D113RE will select the correct values of the abscissas
C            and corresponding weights for different
C            integration rules and null rules and assign them to G
C            and W.
C            The heuristic error coefficients ERRCOF
C            will also be computed.
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
C            The number of points used by each generator.
C
C***REFERENCES  J.Berntsen, Cautious adaptive numerical integration
C               over the 3-cube, Reports in Informatics 17, Dept. of
C               Inf.,Univ. of Bergen, Norway, 1985.
C               J.Berntsen and T.O.Espelid, On the construction of
C               higher degree three-dimensional embedded integration
C               rules, SIAM J. Numer. Anal.,Vol. 25,No. 1, pp.222-234,
C               1988.
C***ROUTINES CALLED-NONE
C***END PROLOGUE D113RE
C
C   Global variables.
C
      INTEGER WTLENG
      DOUBLE PRECISION W(5,WTLENG),G(3,WTLENG),ERRCOF(6)
      DOUBLE PRECISION RULPTS(WTLENG)
C
C   Local variables.
C
      INTEGER I,J
      DOUBLE PRECISION DIM3G(14)
      DOUBLE PRECISION DIM3W(13,5)
C
      DATA (DIM3G(I),I=1,14)/0.1900000000000000D+00,
     +     0.5000000000000000D+00,0.7500000000000000D+00,
     +     0.8000000000000000D+00,0.9949999999999999D+00,
     +     0.9987344998351400D+00,0.7793703685672423D+00,
     +     0.9999698993088767D+00,0.7902637224771788D+00,
     +     0.4403396687650737D+00,0.4378478459006862D+00,
     +     0.9549373822794593D+00,0.9661093133630748D+00,
     +     0.4577105877763134D+00/
C
      DATA (DIM3W(I,1),I=1,13)/0.7923078151105734D-02,
     +     0.6797177392788080D-01,0.1086986538805825D-02,
     +     0.1838633662212829D+00,0.3362119777829031D-01,
     +     0.1013751123334062D-01,0.1687648683985235D-02,
     +     0.1346468564512807D+00,0.1750145884600386D-02,
     +     0.7752336383837454D-01,0.2461864902770251D+00,
     +     0.6797944868483039D-01,0.1419962823300713D-01/
C
      DATA (DIM3W(I,2),I=1,13)/0.1715006248224684D+01,
     +     - .3755893815889209D+00,0.1488632145140549D+00,
     +     - .2497046640620823D+00,0.1792501419135204D+00,
     +     0.3446126758973890D-02, - .5140483185555825D-02,
     +     0.6536017839876425D-02, - .6513454939229700D-03,
     +     - .6304672433547204D-02,0.1266959399788263D-01,
     +     - .5454241018647931D-02,0.4826995274768427D-02/
C
      DATA (DIM3W(I,3),I=1,13)/0.1936014978949526D+01,
     +     - .3673449403754268D+00,0.2929778657898176D-01,
     +     - .1151883520260315D+00,0.5086658220872218D-01,
     +     0.4453911087786469D-01, - .2287828257125900D-01,
     +     0.2908926216345833D-01, - .2898884350669207D-02,
     +     - .2805963413307495D-01,0.5638741361145884D-01,
     +     - .2427469611942451D-01,0.2148307034182882D-01/
C
      DATA (DIM3W(I,4),I=1,13)/0.5170828195605760D+00,
     +     0.1445269144914044D-01, - .3601489663995932D+00,
     +     0.3628307003418485D+00,0.7148802650872729D-02,
     +     - .9222852896022966D-01,0.1719339732471725D-01,
     +     - .1021416537460350D+00, - .7504397861080493D-02,
     +     0.1648362537726711D-01,0.5234610158469334D-01,
     +     0.1445432331613066D-01,0.3019236275367777D-02/
C
      DATA (DIM3W(I,5),I=1,13)/0.2054404503818520D+01,
     +     0.1377759988490120D-01, - .5768062917904410D+00,
     +     0.3726835047700328D-01,0.6814878939777219D-02,
     +     0.5723169733851849D-01, - .4493018743811285D-01,
     +     0.2729236573866348D-01,0.3547473950556990D-03,
     +     0.1571366799739551D-01,0.4990099219278567D-01,
     +     0.1377915552666770D-01,0.2878206423099872D-02/
C
C***FIRST EXECUTABLE STATEMENT D113RE
C
C   Assign values to W.
C
      DO 10 I = 1,13
          DO 10 J = 1,5
              W(J,I) = DIM3W(I,J)
10    CONTINUE
C
C   Assign values to G.
C
      DO 20 I = 1,3
          DO 20 J = 1,13
              G(I,J) = 0
20    CONTINUE
      G(1,2) = DIM3G(1)
      G(1,3) = DIM3G(2)
      G(1,4) = DIM3G(3)
      G(1,5) = DIM3G(4)
      G(1,6) = DIM3G(5)
      G(1,7) = DIM3G(6)
      G(2,7) = G(1,7)
      G(1,8) = DIM3G(7)
      G(2,8) = G(1,8)
      G(1,9) = DIM3G(8)
      G(2,9) = G(1,9)
      G(3,9) = G(1,9)
      G(1,10) = DIM3G(9)
      G(2,10) = G(1,10)
      G(3,10) = G(1,10)
      G(1,11) = DIM3G(10)
      G(2,11) = G(1,11)
      G(3,11) = G(1,11)
      G(1,12) = DIM3G(12)
      G(2,12) = DIM3G(11)
      G(3,12) = G(2,12)
      G(1,13) = DIM3G(13)
      G(2,13) = G(1,13)
      G(3,13) = DIM3G(14)
C
C   Assign values to RULPTS.
C
      RULPTS(1) = 1
      RULPTS(2) = 6
      RULPTS(3) = 6
      RULPTS(4) = 6
      RULPTS(5) = 6
      RULPTS(6) = 6
      RULPTS(7) = 12
      RULPTS(8) = 12
      RULPTS(9) = 8
      RULPTS(10) = 8
      RULPTS(11) = 8
      RULPTS(12) = 24
      RULPTS(13) = 24
C
C   Assign values to ERRCOF.
C
      ERRCOF(1) = 4
      ERRCOF(2) = 4.
      ERRCOF(3) = 0.5
      ERRCOF(4) = 3.
      ERRCOF(5) = 0.5
      ERRCOF(6) = 0.25
C
C***END D113RE
C
      RETURN
      END
