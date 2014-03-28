C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      DOUBLE PRECISION FUNCTION NORMDIST(X,MEAN,VAR)
C================================================================================== 
C     PURPOSE: COMPUTES THE NORMAL DISTRIBUTION AT X GIVEN THE MEAN AND THE VARIANCE
C==================================================================================
C            VARIABLE     DESCRIPTION                      DATA TYPE
C            --------     -----------                      --------- 
C
C     INPUT:
C
C            X                                             DOUBLE PRECISION 
C            MFMEAN                                        DOUBLE PRECISION 
C            MFVAR                                         DOUBLE PRECISION
C
C     RETURNS:
C
C            NORMDIST     NORMAL DISTRIBUTION              DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION X,MEAN,VAR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      NORMDIST = (D_ONE/DSQRT(D_TWO*PI*VAR))
     $     * DEXP(-((X-MEAN)**D_TWO)/(D_TWO*VAR))
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      DOUBLE PRECISION FUNCTION CDFNORMDIST(X,MEAN,VAR)
C================================================================================== 
C     PURPOSE: COMPUTES THE CDF OF THE NORMAL DISTRIBUTION AT X GIVEN THE MEAN 
C              AND THE VARIANCE
C==================================================================================
C            VARIABLE     DESCRIPTION                      DATA TYPE
C            --------     -----------                      --------- 
C
C     INPUT:
C
C            X                                             DOUBLE PRECISION 
C            MFMEAN                                        DOUBLE PRECISION 
C            MFVAR                                         DOUBLE PRECISION
C
C     RETURNS:
C
C            CDFNORMDIST  CFD OF THE NORMAL DISTRIBUTION   DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION X,MEAN,VAR
      DOUBLE PRECISION ERF
      EXTERNAL ERF
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      CDFNORMDIST = D_HALF*(D_ONE + ERF((X-MEAN)/DSQRT(D_TWO*VAR)))
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 
C
C     THE FOLLOWING (FUNCTION PSI AND ITS DEPENDENCIES) ARE PART OF THE
C     SLATEC COMMON MATHEMATICAL LIBRARY. THE SOURCE WAS OBTAINED FROM
C     http://netlib.org/
C
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK BETA
      DOUBLE PRECISION FUNCTION BETA (A, B)
C***BEGIN PROLOGUE  BETA
C***PURPOSE  Compute the complete Beta function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7B
C***TYPE      DOUBLE PRECISION (BETA-S, BETA-D, CBETA-C)
C***KEYWORDS  COMPLETE BETA FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C BETA(A,B) calculates the double precision complete beta function
C for double precision arguments A and B.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DGAMLM, DGAMMA, DLBETA, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  BETA
      DOUBLE PRECISION A, B, ALNSML, XMAX, XMIN, DLBETA, DGAMMA, D1MACH
      LOGICAL FIRST
      EXTERNAL DGAMMA
      SAVE XMAX, ALNSML, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BETA
      IF (FIRST) THEN
         CALL DGAMLM (XMIN, XMAX)
         ALNSML = LOG (D1MACH(1))
      ENDIF
      FIRST = .FALSE.
C
      IF (A .LE. 0.D0 .OR. B .LE. 0.D0) CALL XERMSG ('SLATEC', 'BETA',
     +   'BOTH ARGUMENTS MUST BE GT 0', 2, 2)
C
      IF (A+B.LT.XMAX) BETA = DGAMMA(A)*DGAMMA(B)/DGAMMA(A+B)
      IF (A+B.LT.XMAX) RETURN
C
      BETA = DLBETA (A, B)
      IF (BETA.LT.ALNSML) GO TO 20
      BETA = EXP (BETA)
      RETURN
C
 20   BETA = 0.D0
      CALL XERMSG ('SLATEC', 'BETA',
     +   'A AND/OR B SO BIG BETA UNDERFLOWS', 1, 1)
      RETURN
C
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK BETAI
      DOUBLE PRECISION FUNCTION BETAI (PIN, QIN, X)
C***BEGIN PROLOGUE  BETAI
C***PURPOSE  Calculate the incomplete Beta function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7F
C***TYPE      DOUBLE PRECISION (BETAI-S, BETAI-D)
C***KEYWORDS  FNLIB, INCOMPLETE BETA FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C   BETAI calculates the DOUBLE PRECISION incomplete beta function.
C
C   The incomplete beta function ratio is the probability that a
C   random variable from a beta distribution having parameters PIN and
C   QIN will be less than or equal to X.
C
C     -- Input Arguments -- All arguments are DOUBLE PRECISION.
C   X      upper limit of integration.  X must be in (0,1) inclusive.
C   PIN    first beta distribution parameter.  PIN must be .GT. 0.0.
C   QIN    second beta distribution parameter.  QIN must be .GT. 0.0.
C
C***REFERENCES  Nancy E. Bosten and E. L. Battiste, Remark on Algorithm
C                 179, Communications of the ACM 17, 3 (March 1974),
C                 pp. 156.
C***ROUTINES CALLED  D1MACH, DLBETA, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920528  DESCRIPTION and REFERENCES sections revised.  (WRB)
C***END PROLOGUE  BETAI
      DOUBLE PRECISION X, PIN, QIN, ALNEPS, ALNSML, C, EPS, FINSUM, P,
     1  PS, Q, SML, TERM, XB, XI, Y, D1MACH, DLBETA, P1
      LOGICAL FIRST
      SAVE EPS, ALNEPS, SML, ALNSML, FIRST
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  BETAI
      IF (FIRST) THEN
         EPS = D1MACH(3)
         ALNEPS = LOG (EPS)
         SML = D1MACH(1)
         ALNSML = LOG (SML)
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 0.D0 .OR. X .GT. 1.D0) CALL XERMSG ('SLATEC', 'BETAI',
     +   'X IS NOT IN THE RANGE (0,1)', 1, 2)
      IF (PIN .LE. 0.D0 .OR. QIN .LE. 0.D0) CALL XERMSG ('SLATEC',
     +   'BETAI', 'P AND/OR Q IS LE ZERO', 2, 2)
C
      Y = X
      P = PIN
      Q = QIN
      IF (Q.LE.P .AND. X.LT.0.8D0) GO TO 20
      IF (X.LT.0.2D0) GO TO 20
      Y = 1.0D0 - Y
      P = QIN
      Q = PIN
C
 20   IF ((P+Q)*Y/(P+1.D0).LT.EPS) GO TO 80
C
C EVALUATE THE INFINITE SUM FIRST.  TERM WILL EQUAL
C Y**P/BETA(PS,P) * (1.-PS)-SUB-I * Y**I / FAC(I) .
C
      PS = Q - AINT(Q)
      IF (PS.EQ.0.D0) PS = 1.0D0
      XB = P*LOG(Y) - DLBETA(PS,P) - LOG(P)
      BETAI = 0.0D0
      IF (XB.LT.ALNSML) GO TO 40
C
      BETAI = EXP (XB)
      TERM = BETAI*P
      IF (PS.EQ.1.0D0) GO TO 40
      N = MAX (ALNEPS/LOG(Y), 4.0D0)
      DO 30 I=1,N
        XI = I
        TERM = TERM * (XI-PS)*Y/XI
        BETAI = BETAI + TERM/(P+XI)
 30   CONTINUE
C
C NOW EVALUATE THE FINITE SUM, MAYBE.
C
 40   IF (Q.LE.1.0D0) GO TO 70
C
      XB = P*LOG(Y) + Q*LOG(1.0D0-Y) - DLBETA(P,Q) - LOG(Q)
      IB = MAX (XB/ALNSML, 0.0D0)
      TERM = EXP(XB - IB*ALNSML)
      C = 1.0D0/(1.D0-Y)
      P1 = Q*C/(P+Q-1.D0)
C
      FINSUM = 0.0D0
      N = Q
      IF (Q.EQ.DBLE(N)) N = N - 1
      DO 50 I=1,N
        IF (P1.LE.1.0D0 .AND. TERM/EPS.LE.FINSUM) GO TO 60
        XI = I
        TERM = (Q-XI+1.0D0)*C*TERM/(P+Q-XI)
C
        IF (TERM.GT.1.0D0) IB = IB - 1
        IF (TERM.GT.1.0D0) TERM = TERM*SML
C
        IF (IB.EQ.0) FINSUM = FINSUM + TERM
 50   CONTINUE
C
 60   BETAI = BETAI + FINSUM
 70   IF (Y.NE.X .OR. P.NE.PIN) BETAI = 1.0D0 - BETAI
      BETAI = MAX (MIN (BETAI, 1.0D0), 0.0D0)
      RETURN
C
 80   BETAI = 0.0D0
      XB = P*LOG(MAX(Y,SML)) - LOG(P) - DLBETA(P,Q)
      IF (XB.GT.ALNSML .AND. Y.NE.0.0D0) BETAI = EXP(XB)
      IF (Y.NE.X .OR. P.NE.PIN) BETAI = 1.0D0 - BETAI
C
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK DLBETA
      DOUBLE PRECISION FUNCTION DLBETA (A, B)
C***BEGIN PROLOGUE  DLBETA
C***PURPOSE  Compute the natural logarithm of the complete Beta
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7B
C***TYPE      DOUBLE PRECISION (ALBETA-S, DLBETA-D, CLBETA-C)
C***KEYWORDS  FNLIB, LOGARITHM OF THE COMPLETE BETA FUNCTION,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DLBETA(A,B) calculates the double precision natural logarithm of
C the complete beta function for double precision arguments
C A and B.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D9LGMC, DGAMMA, DLNGAM, DLNREL, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DLBETA
      DOUBLE PRECISION A, B, P, Q, CORR, SQ2PIL, D9LGMC, DGAMMA, DLNGAM,
     1  DLNREL
      EXTERNAL DGAMMA
      SAVE SQ2PIL
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
C***FIRST EXECUTABLE STATEMENT  DLBETA
      P = MIN (A, B)
      Q = MAX (A, B)
C
      IF (P .LE. 0.D0) CALL XERMSG ('SLATEC', 'DLBETA',
     +   'BOTH ARGUMENTS MUST BE GT ZERO', 1, 2)
C
      IF (P.GE.10.D0) GO TO 30
      IF (Q.GE.10.D0) GO TO 20
C
C P AND Q ARE SMALL.
C
      DLBETA = LOG (DGAMMA(P) * (DGAMMA(Q)/DGAMMA(P+Q)) )
      RETURN
C
C P IS SMALL, BUT Q IS BIG.
C
 20   CORR = D9LGMC(Q) - D9LGMC(P+Q)
      DLBETA = DLNGAM(P) + CORR + P - P*LOG(P+Q)
     1  + (Q-0.5D0)*DLNREL(-P/(P+Q))
      RETURN
C
C P AND Q ARE BIG.
C
 30   CORR = D9LGMC(P) + D9LGMC(Q) - D9LGMC(P+Q)
      DLBETA = -0.5D0*LOG(Q) + SQ2PIL + CORR + (P-0.5D0)*LOG(P/(P+Q))
     1  + Q*DLNREL(-P/(P+Q))
      RETURN
C
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK DGAMMA
      DOUBLE PRECISION FUNCTION DGAMMA (X)
C***BEGIN PROLOGUE  DGAMMA
C***PURPOSE  Compute the complete Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      DOUBLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DGAMMA(X) calculates the double precision complete Gamma function
C for double precision argument X.
C
C Series for GAM        on the interval  0.          to  1.00000E+00
C                                        with weighted error   5.79E-32
C                                         log weighted error  31.24
C                               significant figures required  30.00
C                                    decimal places required  32.05
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9LGMC, DCSEVL, DGAMLM, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  DGAMMA
      DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX,
     1  XMIN, Y, D9LGMC, DCSEVL, D1MACH
      LOGICAL FIRST
C
      SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
      DATA GAMCS(  1) / +.8571195590 9893314219 2006239994 2 D-2      /
      DATA GAMCS(  2) / +.4415381324 8410067571 9131577165 2 D-2      /
      DATA GAMCS(  3) / +.5685043681 5993633786 3266458878 9 D-1      /
      DATA GAMCS(  4) / -.4219835396 4185605010 1250018662 4 D-2      /
      DATA GAMCS(  5) / +.1326808181 2124602205 8400679635 2 D-2      /
      DATA GAMCS(  6) / -.1893024529 7988804325 2394702388 6 D-3      /
      DATA GAMCS(  7) / +.3606925327 4412452565 7808221722 5 D-4      /
      DATA GAMCS(  8) / -.6056761904 4608642184 8554829036 5 D-5      /
      DATA GAMCS(  9) / +.1055829546 3022833447 3182350909 3 D-5      /
      DATA GAMCS( 10) / -.1811967365 5423840482 9185589116 6 D-6      /
      DATA GAMCS( 11) / +.3117724964 7153222777 9025459316 9 D-7      /
      DATA GAMCS( 12) / -.5354219639 0196871408 7408102434 7 D-8      /
      DATA GAMCS( 13) / +.9193275519 8595889468 8778682594 0 D-9      /
      DATA GAMCS( 14) / -.1577941280 2883397617 6742327395 3 D-9      /
      DATA GAMCS( 15) / +.2707980622 9349545432 6654043308 9 D-10     /
      DATA GAMCS( 16) / -.4646818653 8257301440 8166105893 3 D-11     /
      DATA GAMCS( 17) / +.7973350192 0074196564 6076717535 9 D-12     /
      DATA GAMCS( 18) / -.1368078209 8309160257 9949917230 9 D-12     /
      DATA GAMCS( 19) / +.2347319486 5638006572 3347177168 8 D-13     /
      DATA GAMCS( 20) / -.4027432614 9490669327 6657053469 9 D-14     /
      DATA GAMCS( 21) / +.6910051747 3721009121 3833697525 7 D-15     /
      DATA GAMCS( 22) / -.1185584500 2219929070 5238712619 2 D-15     /
      DATA GAMCS( 23) / +.2034148542 4963739552 0102605193 2 D-16     /
      DATA GAMCS( 24) / -.3490054341 7174058492 7401294910 8 D-17     /
      DATA GAMCS( 25) / +.5987993856 4853055671 3505106602 6 D-18     /
      DATA GAMCS( 26) / -.1027378057 8722280744 9006977843 1 D-18     /
      DATA GAMCS( 27) / +.1762702816 0605298249 4275966074 8 D-19     /
      DATA GAMCS( 28) / -.3024320653 7353062609 5877211204 2 D-20     /
      DATA GAMCS( 29) / +.5188914660 2183978397 1783355050 6 D-21     /
      DATA GAMCS( 30) / -.8902770842 4565766924 4925160106 6 D-22     /
      DATA GAMCS( 31) / +.1527474068 4933426022 7459689130 6 D-22     /
      DATA GAMCS( 32) / -.2620731256 1873629002 5732833279 9 D-23     /
      DATA GAMCS( 33) / +.4496464047 8305386703 3104657066 6 D-24     /
      DATA GAMCS( 34) / -.7714712731 3368779117 0390152533 3 D-25     /
      DATA GAMCS( 35) / +.1323635453 1260440364 8657271466 6 D-25     /
      DATA GAMCS( 36) / -.2270999412 9429288167 0231381333 3 D-26     /
      DATA GAMCS( 37) / +.3896418998 0039914493 2081663999 9 D-27     /
      DATA GAMCS( 38) / -.6685198115 1259533277 9212799999 9 D-28     /
      DATA GAMCS( 39) / +.1146998663 1400243843 4761386666 6 D-28     /
      DATA GAMCS( 40) / -.1967938586 3451346772 9510399999 9 D-29     /
      DATA GAMCS( 41) / +.3376448816 5853380903 3489066666 6 D-30     /
      DATA GAMCS( 42) / -.5793070335 7821357846 2549333333 3 D-31     /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DGAMMA
      IF (FIRST) THEN
         NGAM = INITDS (GAMCS, 42, 0.1*DBLE(D1MACH(3)) )
C
         CALL DGAMLM (XMIN, XMAX)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.10.D0) GO TO 50
C
C COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
C GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
C
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - N
      N = N - 1
      DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
C COMPUTE GAMMA(X) FOR X .LT. 1.0
C
      N = -N
      IF (X .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', 'X IS 0', 4, 2)
      IF (X .LT. 0.0 .AND. X+N-2 .EQ. 0.D0) CALL XERMSG ('SLATEC',
     +   'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
      IF (X .LT. (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X) .LT. DXREL)
     +   CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',
     +   1, 1)
C
      DO 20 I=1,N
        DGAMMA = DGAMMA/(X+I-1 )
 20   CONTINUE
      RETURN
C
C GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
C
 30   DO 40 I=1,N
        DGAMMA = (Y+I) * DGAMMA
 40   CONTINUE
      RETURN
C
C GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
C
 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'X SO BIG GAMMA OVERFLOWS', 3, 2)
C
      DGAMMA = 0.D0
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
      IF (X.LT.XMIN) RETURN
C
      DGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
C
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
     +   'DGAMMA',
     +   'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
C
      SINPIY = SIN (PI*Y)
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA',
     +   'X IS A NEGATIVE INTEGER', 4, 2)
C
      DGAMMA = -PI/(Y*SINPIY*DGAMMA)
C
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK DGAMLM
      SUBROUTINE DGAMLM (XMIN, XMAX)
C***BEGIN PROLOGUE  DGAMLM
C***PURPOSE  Compute the minimum and maximum bounds for the argument in
C            the Gamma function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A, R2
C***TYPE      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Calculate the minimum and maximum legal bounds for X in gamma(X).
C XMIN and XMAX are not the only bounds, but they are the only non-
C trivial ones to calculate.
C
C             Output Arguments --
C XMIN   double precision minimum legal value of X in gamma(X).  Any
C        smaller value of X might result in underflow.
C XMAX   double precision maximum legal value of X in gamma(X).  Any
C        larger value of X might cause overflow.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DGAMLM
      DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD, D1MACH
C***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = LOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML)
     1    / (XMIN*XLN+0.5D0)
        IF (ABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMIN', 1, 2)
C
 20   XMIN = -XMIN + 0.01D0
C
      ALNBIG = LOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG)
     1    / (XMAX*XLN-0.5D0)
        IF (ABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMAX', 2, 2)
C
 40   XMAX = XMAX - 0.01D0
      XMIN = MAX (XMIN, -XMAX+1.D0)
C
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK DLNREL
      DOUBLE PRECISION FUNCTION DLNREL (X)
C***BEGIN PROLOGUE  DLNREL
C***PURPOSE  Evaluate ln(1+X) accurate in the sense of relative error.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4B
C***TYPE      DOUBLE PRECISION (ALNREL-S, DLNREL-D, CLNREL-C)
C***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DLNREL(X) calculates the double precision natural logarithm of
C (1.0+X) for double precision argument X.  This routine should
C be used when X is small and accurate to calculate the logarithm
C accurately (in the relative error sense) in the neighborhood
C of 1.0.
C
C Series for ALNR       on the interval -3.75000E-01 to  3.75000E-01
C                                        with weighted error   6.35E-32
C                                         log weighted error  31.20
C                               significant figures required  30.93
C                                    decimal places required  32.01
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DLNREL
      DOUBLE PRECISION ALNRCS(43), X, XMIN,  DCSEVL, D1MACH
      LOGICAL FIRST
      SAVE ALNRCS, NLNREL, XMIN, FIRST
      DATA ALNRCS(  1) / +.1037869356 2743769800 6862677190 98 D+1     /
      DATA ALNRCS(  2) / -.1336430150 4908918098 7660415531 33 D+0     /
      DATA ALNRCS(  3) / +.1940824913 5520563357 9261993747 50 D-1     /
      DATA ALNRCS(  4) / -.3010755112 7535777690 3765377765 92 D-2     /
      DATA ALNRCS(  5) / +.4869461479 7154850090 4563665091 37 D-3     /
      DATA ALNRCS(  6) / -.8105488189 3175356066 8099430086 22 D-4     /
      DATA ALNRCS(  7) / +.1377884779 9559524782 9382514960 59 D-4     /
      DATA ALNRCS(  8) / -.2380221089 4358970251 3699929149 35 D-5     /
      DATA ALNRCS(  9) / +.4164041621 3865183476 3918599019 89 D-6     /
      DATA ALNRCS( 10) / -.7359582837 8075994984 2668370319 98 D-7     /
      DATA ALNRCS( 11) / +.1311761187 6241674949 1522943450 11 D-7     /
      DATA ALNRCS( 12) / -.2354670931 7742425136 6960923301 75 D-8     /
      DATA ALNRCS( 13) / +.4252277327 6034997775 6380529625 67 D-9     /
      DATA ALNRCS( 14) / -.7719089413 4840796826 1081074933 00 D-10    /
      DATA ALNRCS( 15) / +.1407574648 1359069909 2153564721 91 D-10    /
      DATA ALNRCS( 16) / -.2576907205 8024680627 5370786275 84 D-11    /
      DATA ALNRCS( 17) / +.4734240666 6294421849 1543950059 38 D-12    /
      DATA ALNRCS( 18) / -.8724901267 4742641745 3012632926 75 D-13    /
      DATA ALNRCS( 19) / +.1612461490 2740551465 7398331191 15 D-13    /
      DATA ALNRCS( 20) / -.2987565201 5665773006 7107924168 15 D-14    /
      DATA ALNRCS( 21) / +.5548070120 9082887983 0413216972 79 D-15    /
      DATA ALNRCS( 22) / -.1032461915 8271569595 1413339619 32 D-15    /
      DATA ALNRCS( 23) / +.1925023920 3049851177 8785032448 68 D-16    /
      DATA ALNRCS( 24) / -.3595507346 5265150011 1897078442 66 D-17    /
      DATA ALNRCS( 25) / +.6726454253 7876857892 1945742267 73 D-18    /
      DATA ALNRCS( 26) / -.1260262416 8735219252 0824256375 46 D-18    /
      DATA ALNRCS( 27) / +.2364488440 8606210044 9161589555 19 D-19    /
      DATA ALNRCS( 28) / -.4441937705 0807936898 8783891797 33 D-20    /
      DATA ALNRCS( 29) / +.8354659446 4034259016 2412939946 66 D-21    /
      DATA ALNRCS( 30) / -.1573155941 6479562574 8992535210 66 D-21    /
      DATA ALNRCS( 31) / +.2965312874 0247422686 1543697066 66 D-22    /
      DATA ALNRCS( 32) / -.5594958348 1815947292 1560132266 66 D-23    /
      DATA ALNRCS( 33) / +.1056635426 8835681048 1872841386 66 D-23    /
      DATA ALNRCS( 34) / -.1997248368 0670204548 3149994666 66 D-24    /
      DATA ALNRCS( 35) / +.3778297781 8839361421 0498559999 99 D-25    /
      DATA ALNRCS( 36) / -.7153158688 9081740345 0381653333 33 D-26    /
      DATA ALNRCS( 37) / +.1355248846 3674213646 5020245333 33 D-26    /
      DATA ALNRCS( 38) / -.2569467304 8487567430 0798293333 33 D-27    /
      DATA ALNRCS( 39) / +.4874775606 6216949076 4595199999 99 D-28    /
      DATA ALNRCS( 40) / -.9254211253 0849715321 1323733333 33 D-29    /
      DATA ALNRCS( 41) / +.1757859784 1760239233 2697600000 00 D-29    /
      DATA ALNRCS( 42) / -.3341002667 7731010351 3770666666 66 D-30    /
      DATA ALNRCS( 43) / +.6353393618 0236187354 1802666666 66 D-31    /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DLNREL
      IF (FIRST) THEN
         NLNREL = INITDS (ALNRCS, 43, 0.1*DBLE(D1MACH(3)))
         XMIN = -1.0D0 + SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LE. (-1.D0)) CALL XERMSG ('SLATEC', 'DLNREL', 'X IS LE -1'
     +   , 2, 2)
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DLNREL',
     +   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 1, 1)
C
      IF (ABS(X).LE.0.375D0) DLNREL = X*(1.D0 -
     1  X*DCSEVL (X/.375D0, ALNRCS, NLNREL))
C
      IF (ABS(X).GT.0.375D0) DLNREL = LOG (1.0D0+X)
C
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK DLNGAM
      DOUBLE PRECISION FUNCTION DLNGAM (X)
C***BEGIN PROLOGUE  DLNGAM
C***PURPOSE  Compute the logarithm of the absolute value of the Gamma
C            function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7A
C***TYPE      DOUBLE PRECISION (ALNGAM-S, DLNGAM-D, CLNGAM-C)
C***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DLNGAM(X) calculates the double precision logarithm of the
C absolute value of the Gamma function for double precision
C argument X.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, D9LGMC, DGAMMA, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C***END PROLOGUE  DLNGAM
      DOUBLE PRECISION X, DXREL, PI, SINPIY, SQPI2L, SQ2PIL, XMAX,
     1  Y, DGAMMA, D9LGMC, D1MACH, TEMP
      LOGICAL FIRST
      EXTERNAL DGAMMA
      SAVE SQ2PIL, SQPI2L, PI, XMAX, DXREL, FIRST
      DATA SQ2PIL / 0.9189385332 0467274178 0329736405 62 D0 /
      DATA SQPI2L / +.2257913526 4472743236 3097614947 441 D+0    /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DLNGAM
      IF (FIRST) THEN
         TEMP = 1.D0/LOG(D1MACH(2))
         XMAX = TEMP*D1MACH(2)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS (X)
      IF (Y.GT.10.D0) GO TO 20
C
C LOG (ABS (DGAMMA(X)) ) FOR ABS(X) .LE. 10.0
C
      DLNGAM = LOG (ABS (DGAMMA(X)) )
      RETURN
C
C LOG ( ABS (DGAMMA(X)) ) FOR ABS(X) .GT. 10.0
C
 20   IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DLNGAM',
     +   'ABS(X) SO BIG DLNGAM OVERFLOWS', 2, 2)
C
      IF (X.GT.0.D0) DLNGAM = SQ2PIL + (X-0.5D0)*LOG(X) - X + D9LGMC(Y)
      IF (X.GT.0.D0) RETURN
C
      SINPIY = ABS (SIN(PI*Y))
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DLNGAM',
     +   'X IS A NEGATIVE INTEGER', 3, 2)
C
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC',
     +   'DLNGAM',
     +   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',
     +   1, 1)
C
      DLNGAM = SQPI2L + (X-0.5D0)*LOG(Y) - X - LOG(SINPIY) - D9LGMC(Y)
      RETURN
C
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK D9LGMC
      DOUBLE PRECISION FUNCTION D9LGMC (X)
C***BEGIN PROLOGUE  D9LGMC
C***SUBSIDIARY
C***PURPOSE  Compute the log Gamma correction factor so that
C            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
C            + D9LGMC(X).
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7E
C***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-C)
C***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
C             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C Compute the log gamma correction factor for X .GE. 10. so that
C LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
C
C Series for ALGM       on the interval  0.          to  1.00000E-02
C                                        with weighted error   1.28E-31
C                                         log weighted error  30.89
C                               significant figures required  29.81
C                                    decimal places required  31.48
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900720  Routine changed from user-callable to subsidiary.  (WRB)
C***END PROLOGUE  D9LGMC
      DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, DCSEVL, D1MACH
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS(  1) / +.1666389480 4518632472 0572965082 2 D+0      /
      DATA ALGMCS(  2) / -.1384948176 0675638407 3298605913 5 D-4      /
      DATA ALGMCS(  3) / +.9810825646 9247294261 5717154748 7 D-8      /
      DATA ALGMCS(  4) / -.1809129475 5724941942 6330626671 9 D-10     /
      DATA ALGMCS(  5) / +.6221098041 8926052271 2601554341 6 D-13     /
      DATA ALGMCS(  6) / -.3399615005 4177219443 0333059966 6 D-15     /
      DATA ALGMCS(  7) / +.2683181998 4826987489 5753884666 6 D-17     /
      DATA ALGMCS(  8) / -.2868042435 3346432841 4462239999 9 D-19     /
      DATA ALGMCS(  9) / +.3962837061 0464348036 7930666666 6 D-21     /
      DATA ALGMCS( 10) / -.6831888753 9857668701 1199999999 9 D-23     /
      DATA ALGMCS( 11) / +.1429227355 9424981475 7333333333 3 D-24     /
      DATA ALGMCS( 12) / -.3547598158 1010705471 9999999999 9 D-26     /
      DATA ALGMCS( 13) / +.1025680058 0104709120 0000000000 0 D-27     /
      DATA ALGMCS( 14) / -.3401102254 3167487999 9999999999 9 D-29     /
      DATA ALGMCS( 15) / +.1276642195 6300629333 3333333333 3 D-30     /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (FIRST) THEN
         NALGM = INITDS (ALGMCS, 15, DBLE(D1MACH(3)) )
         XBIG = 1.0D0/SQRT(D1MACH(3))
         XMAX = EXP (MIN(LOG(D1MACH(2)/12.D0), -LOG(12.D0*D1MACH(1))))
      ENDIF
      FIRST = .FALSE.
C
      IF (X .LT. 10.D0) CALL XERMSG ('SLATEC', 'D9LGMC',
     +   'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
C
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS,
     1  NALGM) / X
      RETURN
C
 20   D9LGMC = 0.D0
      CALL XERMSG ('SLATEC', 'D9LGMC', 'X SO BIG D9LGMC UNDERFLOWS', 2,
     +   1)
      RETURN
C
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK ERF
      DOUBLE PRECISION FUNCTION ERF (X)
C***BEGIN PROLOGUE  ERF
C***PURPOSE  Compute the error function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C8A, L5A1E
C***TYPE      DOUBLE PRECISION (ERF-S, ERF-D)
C***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C ERF(X) calculates the double precision error function for double
C precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000E+00
C                                        with weighted error   1.28E-32
C                                         log weighted error  31.89
C                               significant figures required  31.05
C                                    decimal places required  32.55
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, ERFC, INITDS
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900727  Added EXTERNAL statement.  (WRB)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  ERF
      DOUBLE PRECISION X, ERFCS(21), SQEPS, SQRTPI, XBIG, Y, D1MACH,
     1  DCSEVL, ERFC
      LOGICAL FIRST
      EXTERNAL ERFC
      SAVE ERFCS, SQRTPI, NTERF, XBIG, SQEPS, FIRST
      DATA ERFCS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERFCS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERFCS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERFCS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERFCS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERFCS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERFCS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERFCS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERFCS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERFCS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERFCS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERFCS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERFCS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERFCS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERFCS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERFCS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERFCS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERFCS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERFCS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERFCS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERFCS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  ERF
      IF (FIRST) THEN
         NTERF = INITDS (ERFCS, 21, 0.1*DBLE(D1MACH(3)))
         XBIG = SQRT(-LOG(SQRTPI*D1MACH(3)))
         SQEPS = SQRT(2.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y.GT.1.D0) GO TO 20
C
C ERF(X) = 1.0 - ERFC(X)  FOR  -1.0 .LE. X .LE. 1.0
C
      IF (Y.LE.SQEPS) ERF = 2.0D0*X*X/SQRTPI
      IF (Y.GT.SQEPS) ERF = X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,
     1  ERFCS, NTERF))
      RETURN
C
C ERF(X) = 1.0 - ERFC(X) FOR ABS(X) .GT. 1.0
C
 20   IF (Y.LE.XBIG) ERF = SIGN (1.0D0-ERFC(Y), X)
      IF (Y.GT.XBIG) ERF = SIGN (1.0D0, X)
C
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK ERFC
      DOUBLE PRECISION FUNCTION ERFC (X)
C***BEGIN PROLOGUE  ERFC
C***PURPOSE  Compute the complementary error function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C8A, L5A1E
C***TYPE      DOUBLE PRECISION (ERFC-S, ERFC-D)
C***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB,
C             SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C ERFC(X) calculates the double precision complementary error function
C for double precision argument X.
C
C Series for ERF        on the interval  0.          to  1.00000E+00
C                                        with weighted Error   1.28E-32
C                                         log weighted Error  31.89
C                               significant figures required  31.05
C                                    decimal places required  32.55
C
C Series for ERC2       on the interval  2.50000E-01 to  1.00000E+00
C                                        with weighted Error   2.67E-32
C                                         log weighted Error  31.57
C                               significant figures required  30.31
C                                    decimal places required  32.42
C
C Series for ERFC       on the interval  0.          to  2.50000E-01
C                                        with weighted error   1.53E-31
C                                         log weighted error  30.82
C                               significant figures required  29.47
C                                    decimal places required  31.70
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770701  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  ERFC
      DOUBLE PRECISION X, ERFCS(21), ERFCCS(59), ERC2CS(49), SQEPS,
     1  SQRTPI, XMAX, TXMAX, XSML, Y, D1MACH, DCSEVL
      LOGICAL FIRST
      SAVE ERFCS, ERC2CS, ERFCCS, SQRTPI, NTERF,
     1 NTERFC, NTERC2, XSML, XMAX, SQEPS, FIRST
      DATA ERFCS(  1) / -.4904612123 4691808039 9845440333 76 D-1     /
      DATA ERFCS(  2) / -.1422612051 0371364237 8247418996 31 D+0     /
      DATA ERFCS(  3) / +.1003558218 7599795575 7546767129 33 D-1     /
      DATA ERFCS(  4) / -.5768764699 7674847650 8270255091 67 D-3     /
      DATA ERFCS(  5) / +.2741993125 2196061034 4221607914 71 D-4     /
      DATA ERFCS(  6) / -.1104317550 7344507604 1353812959 05 D-5     /
      DATA ERFCS(  7) / +.3848875542 0345036949 9613114981 74 D-7     /
      DATA ERFCS(  8) / -.1180858253 3875466969 6317518015 81 D-8     /
      DATA ERFCS(  9) / +.3233421582 6050909646 4029309533 54 D-10    /
      DATA ERFCS( 10) / -.7991015947 0045487581 6073747085 95 D-12    /
      DATA ERFCS( 11) / +.1799072511 3961455611 9672454866 34 D-13    /
      DATA ERFCS( 12) / -.3718635487 8186926382 3168282094 93 D-15    /
      DATA ERFCS( 13) / +.7103599003 7142529711 6899083946 66 D-17    /
      DATA ERFCS( 14) / -.1261245511 9155225832 4954248533 33 D-18    /
      DATA ERFCS( 15) / +.2091640694 1769294369 1705002666 66 D-20    /
      DATA ERFCS( 16) / -.3253973102 9314072982 3641600000 00 D-22    /
      DATA ERFCS( 17) / +.4766867209 7976748332 3733333333 33 D-24    /
      DATA ERFCS( 18) / -.6598012078 2851343155 1999999999 99 D-26    /
      DATA ERFCS( 19) / +.8655011469 9637626197 3333333333 33 D-28    /
      DATA ERFCS( 20) / -.1078892517 7498064213 3333333333 33 D-29    /
      DATA ERFCS( 21) / +.1281188399 3017002666 6666666666 66 D-31    /
      DATA ERC2CS(  1) / -.6960134660 2309501127 3915082619 7 D-1      /
      DATA ERC2CS(  2) / -.4110133936 2620893489 8221208466 6 D-1      /
      DATA ERC2CS(  3) / +.3914495866 6896268815 6114370524 4 D-2      /
      DATA ERC2CS(  4) / -.4906395650 5489791612 8093545077 4 D-3      /
      DATA ERC2CS(  5) / +.7157479001 3770363807 6089414182 5 D-4      /
      DATA ERC2CS(  6) / -.1153071634 1312328338 0823284791 2 D-4      /
      DATA ERC2CS(  7) / +.1994670590 2019976350 5231486770 9 D-5      /
      DATA ERC2CS(  8) / -.3642666471 5992228739 3611843071 1 D-6      /
      DATA ERC2CS(  9) / +.6944372610 0050125899 3127721463 3 D-7      /
      DATA ERC2CS( 10) / -.1371220902 1043660195 3460514121 0 D-7      /
      DATA ERC2CS( 11) / +.2788389661 0071371319 6386034808 7 D-8      /
      DATA ERC2CS( 12) / -.5814164724 3311615518 6479105031 6 D-9      /
      DATA ERC2CS( 13) / +.1238920491 7527531811 8016881795 0 D-9      /
      DATA ERC2CS( 14) / -.2690639145 3067434323 9042493788 9 D-10     /
      DATA ERC2CS( 15) / +.5942614350 8479109824 4470968384 0 D-11     /
      DATA ERC2CS( 16) / -.1332386735 7581195792 8775442057 0 D-11     /
      DATA ERC2CS( 17) / +.3028046806 1771320171 7369724330 4 D-12     /
      DATA ERC2CS( 18) / -.6966648814 9410325887 9586758895 4 D-13     /
      DATA ERC2CS( 19) / +.1620854541 0539229698 1289322762 8 D-13     /
      DATA ERC2CS( 20) / -.3809934465 2504919998 7691305772 9 D-14     /
      DATA ERC2CS( 21) / +.9040487815 9788311493 6897101297 5 D-15     /
      DATA ERC2CS( 22) / -.2164006195 0896073478 0981204700 3 D-15     /
      DATA ERC2CS( 23) / +.5222102233 9958549846 0798024417 2 D-16     /
      DATA ERC2CS( 24) / -.1269729602 3645553363 7241552778 0 D-16     /
      DATA ERC2CS( 25) / +.3109145504 2761975838 3622741295 1 D-17     /
      DATA ERC2CS( 26) / -.7663762920 3203855240 0956671481 1 D-18     /
      DATA ERC2CS( 27) / +.1900819251 3627452025 3692973329 0 D-18     /
      DATA ERC2CS( 28) / -.4742207279 0690395452 2565599996 5 D-19     /
      DATA ERC2CS( 29) / +.1189649200 0765283828 8068307845 1 D-19     /
      DATA ERC2CS( 30) / -.3000035590 3257802568 4527131306 6 D-20     /
      DATA ERC2CS( 31) / +.7602993453 0432461730 1938527709 8 D-21     /
      DATA ERC2CS( 32) / -.1935909447 6068728815 6981104913 0 D-21     /
      DATA ERC2CS( 33) / +.4951399124 7733378810 0004238677 3 D-22     /
      DATA ERC2CS( 34) / -.1271807481 3363718796 0862198988 8 D-22     /
      DATA ERC2CS( 35) / +.3280049600 4695130433 1584165205 3 D-23     /
      DATA ERC2CS( 36) / -.8492320176 8228965689 2479242239 9 D-24     /
      DATA ERC2CS( 37) / +.2206917892 8075602235 1987998719 9 D-24     /
      DATA ERC2CS( 38) / -.5755617245 6965284983 1281950719 9 D-25     /
      DATA ERC2CS( 39) / +.1506191533 6392342503 5414405119 9 D-25     /
      DATA ERC2CS( 40) / -.3954502959 0187969531 0428569599 9 D-26     /
      DATA ERC2CS( 41) / +.1041529704 1515009799 8464505173 3 D-26     /
      DATA ERC2CS( 42) / -.2751487795 2787650794 5017890133 3 D-27     /
      DATA ERC2CS( 43) / +.7290058205 4975574089 9770368000 0 D-28     /
      DATA ERC2CS( 44) / -.1936939645 9159478040 7750109866 6 D-28     /
      DATA ERC2CS( 45) / +.5160357112 0514872983 7005482666 6 D-29     /
      DATA ERC2CS( 46) / -.1378419322 1930940993 8964480000 0 D-29     /
      DATA ERC2CS( 47) / +.3691326793 1070690422 5109333333 3 D-30     /
      DATA ERC2CS( 48) / -.9909389590 6243654206 5322666666 6 D-31     /
      DATA ERC2CS( 49) / +.2666491705 1953884133 2394666666 6 D-31     /
      DATA ERFCCS(  1) / +.7151793102 0292477450 3697709496 D-1        /
      DATA ERFCCS(  2) / -.2653243433 7606715755 8893386681 D-1        /
      DATA ERFCCS(  3) / +.1711153977 9208558833 2699194606 D-2        /
      DATA ERFCCS(  4) / -.1637516634 5851788416 3746404749 D-3        /
      DATA ERFCCS(  5) / +.1987129350 0552036499 5974806758 D-4        /
      DATA ERFCCS(  6) / -.2843712412 7665550875 0175183152 D-5        /
      DATA ERFCCS(  7) / +.4606161308 9631303696 9379968464 D-6        /
      DATA ERFCCS(  8) / -.8227753025 8792084205 7766536366 D-7        /
      DATA ERFCCS(  9) / +.1592141872 7709011298 9358340826 D-7        /
      DATA ERFCCS( 10) / -.3295071362 2528432148 6631665072 D-8        /
      DATA ERFCCS( 11) / +.7223439760 4005554658 1261153890 D-9        /
      DATA ERFCCS( 12) / -.1664855813 3987295934 4695966886 D-9        /
      DATA ERFCCS( 13) / +.4010392588 2376648207 7671768814 D-10       /
      DATA ERFCCS( 14) / -.1004816214 4257311327 2170176283 D-10       /
      DATA ERFCCS( 15) / +.2608275913 3003338085 9341009439 D-11       /
      DATA ERFCCS( 16) / -.6991110560 4040248655 7697812476 D-12       /
      DATA ERFCCS( 17) / +.1929492333 2617070862 4205749803 D-12       /
      DATA ERFCCS( 18) / -.5470131188 7543310649 0125085271 D-13       /
      DATA ERFCCS( 19) / +.1589663309 7626974483 9084032762 D-13       /
      DATA ERFCCS( 20) / -.4726893980 1975548392 0369584290 D-14       /
      DATA ERFCCS( 21) / +.1435873376 7849847867 2873997840 D-14       /
      DATA ERFCCS( 22) / -.4449510561 8173583941 7250062829 D-15       /
      DATA ERFCCS( 23) / +.1404810884 7682334373 7305537466 D-15       /
      DATA ERFCCS( 24) / -.4513818387 7642108962 5963281623 D-16       /
      DATA ERFCCS( 25) / +.1474521541 0451330778 7018713262 D-16       /
      DATA ERFCCS( 26) / -.4892621406 9457761543 6841552532 D-17       /
      DATA ERFCCS( 27) / +.1647612141 4106467389 5301522827 D-17       /
      DATA ERFCCS( 28) / -.5626817176 3294080929 9928521323 D-18       /
      DATA ERFCCS( 29) / +.1947443382 2320785142 9197867821 D-18       /
      DATA ERFCCS( 30) / -.6826305642 9484207295 6664144723 D-19       /
      DATA ERFCCS( 31) / +.2421988887 2986492401 8301125438 D-19       /
      DATA ERFCCS( 32) / -.8693414133 5030704256 3800861857 D-20       /
      DATA ERFCCS( 33) / +.3155180346 2280855712 2363401262 D-20       /
      DATA ERFCCS( 34) / -.1157372324 0496087426 1239486742 D-20       /
      DATA ERFCCS( 35) / +.4288947161 6056539462 3737097442 D-21       /
      DATA ERFCCS( 36) / -.1605030742 0576168500 5737770964 D-21       /
      DATA ERFCCS( 37) / +.6063298757 4538026449 5069923027 D-22       /
      DATA ERFCCS( 38) / -.2311404251 6979584909 8840801367 D-22       /
      DATA ERFCCS( 39) / +.8888778540 6618855255 4702955697 D-23       /
      DATA ERFCCS( 40) / -.3447260576 6513765223 0718495566 D-23       /
      DATA ERFCCS( 41) / +.1347865460 2069650682 7582774181 D-23       /
      DATA ERFCCS( 42) / -.5311794071 1250217364 5873201807 D-24       /
      DATA ERFCCS( 43) / +.2109341058 6197831682 8954734537 D-24       /
      DATA ERFCCS( 44) / -.8438365587 9237891159 8133256738 D-25       /
      DATA ERFCCS( 45) / +.3399982524 9452089062 7359576337 D-25       /
      DATA ERFCCS( 46) / -.1379452388 0732420900 2238377110 D-25       /
      DATA ERFCCS( 47) / +.5634490311 8332526151 3392634811 D-26       /
      DATA ERFCCS( 48) / -.2316490434 4770654482 3427752700 D-26       /
      DATA ERFCCS( 49) / +.9584462844 6018101526 3158381226 D-27       /
      DATA ERFCCS( 50) / -.3990722880 3301097262 4224850193 D-27       /
      DATA ERFCCS( 51) / +.1672129225 9444773601 7228709669 D-27       /
      DATA ERFCCS( 52) / -.7045991522 7660138563 8803782587 D-28       /
      DATA ERFCCS( 53) / +.2979768402 8642063541 2357989444 D-28       /
      DATA ERFCCS( 54) / -.1262522466 4606192972 2422632994 D-28       /
      DATA ERFCCS( 55) / +.5395438704 5424879398 5299653154 D-29       /
      DATA ERFCCS( 56) / -.2380992882 5314591867 5346190062 D-29       /
      DATA ERFCCS( 57) / +.1099052830 1027615735 9726683750 D-29       /
      DATA ERFCCS( 58) / -.4867713741 6449657273 2518677435 D-30       /
      DATA ERFCCS( 59) / +.1525877264 1103575676 3200828211 D-30       /
      DATA SQRTPI / 1.772453850 9055160272 9816748334 115D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  ERFC
      IF (FIRST) THEN
         ETA = 0.1*DBLE(D1MACH(3))
         NTERF = INITDS (ERFCS, 21, ETA)
         NTERFC = INITDS (ERFCCS, 59, ETA)
         NTERC2 = INITDS (ERC2CS, 49, ETA)
C
         XSML = -SQRT(-LOG(SQRTPI*D1MACH(3)))
         TXMAX = SQRT(-LOG(SQRTPI*D1MACH(1)))
         XMAX = TXMAX - 0.5D0*LOG(TXMAX)/TXMAX - 0.01D0
         SQEPS = SQRT(2.0D0*D1MACH(3))
      ENDIF
      FIRST = .FALSE.
C
      IF (X.GT.XSML) GO TO 20
C
C ERFC(X) = 1.0 - ERF(X)  FOR  X .LT. XSML
C
      ERFC = 2.0D0
      RETURN
C
 20   IF (X.GT.XMAX) GO TO 40
      Y = ABS(X)
      IF (Y.GT.1.0D0) GO TO 30
C
C ERFC(X) = 1.0 - ERF(X)  FOR ABS(X) .LE. 1.0
C
      IF (Y.LT.SQEPS) ERFC = 1.0D0 - 2.0D0*X/SQRTPI
      IF (Y.GE.SQEPS) ERFC = 1.0D0 - X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0,
     1  ERFCS, NTERF))
      RETURN
C
C ERFC(X) = 1.0 - ERF(X)  FOR  1.0 .LT. ABS(X) .LE. XMAX
C
 30   Y = Y*Y
      IF (Y.LE.4.D0) ERFC = EXP(-Y)/ABS(X) * (0.5D0 + DCSEVL (
     1  (8.D0/Y-5.D0)/3.D0, ERC2CS, NTERC2) )
      IF (Y.GT.4.D0) ERFC = EXP(-Y)/ABS(X) * (0.5D0 + DCSEVL (
     1  8.D0/Y-1.D0, ERFCCS, NTERFC) )
      IF (X.LT.0.D0) ERFC = 2.0D0 - ERFC
      RETURN
C
 40   CALL XERMSG ('SLATEC', 'ERFC', 'X SO BIG ERFC UNDERFLOWS', 1, 1)
      ERFC = 0.D0
      RETURN
C
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

*DECK DPSI
      DOUBLE PRECISION FUNCTION DPSI (X)
C***BEGIN PROLOGUE  DPSI
C***PURPOSE  Compute the Psi (or Digamma) function.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C7C
C***TYPE      DOUBLE PRECISION (PSI-S, DPSI-D, CPSI-C)
C***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DPSI calculates the double precision Psi (or Digamma) function for
C double precision argument X.  PSI(X) is the logarithmic derivative
C of the Gamma function of X.
C
C Series for PSI        on the interval  0.          to  1.00000E+00
C                                        with weighted error   5.79E-32
C                                         log weighted error  31.24
C                               significant figures required  30.93
C                                    decimal places required  32.05
C
C
C Series for APSI       on the interval  0.          to  1.00000E-02
C                                        with weighted error   7.75E-33
C                                         log weighted error  32.11
C                               significant figures required  28.88
C                                    decimal places required  32.71
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCOT, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890911  Removed unnecessary intrinsics.  (WRB)
C   890911  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900727  Added EXTERNAL statement.  (WRB)
C   920618  Removed space from variable name.  (RWC, WRB)
C***END PROLOGUE  DPSI
      DOUBLE PRECISION X, PSICS(42), APSICS(16), AUX, DXREL, PI, XBIG,
     1  Y, DCOT, DCSEVL, D1MACH
      LOGICAL FIRST
      EXTERNAL DCOT
      SAVE PSICS, APSICS, PI, NTPSI, NTAPSI, XBIG, DXREL, FIRST
      DATA PSICS(  1) / -.3805708083 5217921520 4376776670 39 D-1     /
      DATA PSICS(  2) / +.4914153930 2938712748 2046996542 77 D+0     /
      DATA PSICS(  3) / -.5681574782 1244730242 8920647340 81 D-1     /
      DATA PSICS(  4) / +.8357821225 9143131362 7756507478 62 D-2     /
      DATA PSICS(  5) / -.1333232857 9943425998 0792741723 93 D-2     /
      DATA PSICS(  6) / +.2203132870 6930824892 8723979795 21 D-3     /
      DATA PSICS(  7) / -.3704023817 8456883592 8890869492 29 D-4     /
      DATA PSICS(  8) / +.6283793654 8549898933 6514187176 90 D-5     /
      DATA PSICS(  9) / -.1071263908 5061849855 2835417470 74 D-5     /
      DATA PSICS( 10) / +.1831283946 5484165805 7315898103 78 D-6     /
      DATA PSICS( 11) / -.3135350936 1808509869 0057797968 85 D-7     /
      DATA PSICS( 12) / +.5372808776 2007766260 4719191436 15 D-8     /
      DATA PSICS( 13) / -.9211681415 9784275717 8806326247 30 D-9     /
      DATA PSICS( 14) / +.1579812652 1481822782 2528840328 23 D-9     /
      DATA PSICS( 15) / -.2709864613 2380443065 4405894097 07 D-10    /
      DATA PSICS( 16) / +.4648722859 9096834872 9473195295 49 D-11    /
      DATA PSICS( 17) / -.7975272563 8303689726 5047977727 37 D-12    /
      DATA PSICS( 18) / +.1368272385 7476992249 2510538928 38 D-12    /
      DATA PSICS( 19) / -.2347515606 0658972717 3206779807 19 D-13    /
      DATA PSICS( 20) / +.4027630715 5603541107 9079250062 81 D-14    /
      DATA PSICS( 21) / -.6910251853 1179037846 5474229747 71 D-15    /
      DATA PSICS( 22) / +.1185604713 8863349552 9291395257 68 D-15    /
      DATA PSICS( 23) / -.2034168961 6261559308 1542104842 23 D-16    /
      DATA PSICS( 24) / +.3490074968 6463043850 3742329323 51 D-17    /
      DATA PSICS( 25) / -.5988014693 4976711003 0110813934 93 D-18    /
      DATA PSICS( 26) / +.1027380162 8080588258 3980057122 13 D-18    /
      DATA PSICS( 27) / -.1762704942 4561071368 3592601053 86 D-19    /
      DATA PSICS( 28) / +.3024322801 8156920457 4540354901 33 D-20    /
      DATA PSICS( 29) / -.5188916830 2092313774 2860888746 66 D-21    /
      DATA PSICS( 30) / +.8902773034 5845713905 0058874879 99 D-22    /
      DATA PSICS( 31) / -.1527474289 9426728392 8949719040 00 D-22    /
      DATA PSICS( 32) / +.2620731479 8962083136 3583180799 99 D-23    /
      DATA PSICS( 33) / -.4496464273 8220696772 5983880533 33 D-24    /
      DATA PSICS( 34) / +.7714712959 6345107028 9193642666 66 D-25    /
      DATA PSICS( 35) / -.1323635476 1887702968 1026389333 33 D-25    /
      DATA PSICS( 36) / +.2270999436 2408300091 2773119999 99 D-26    /
      DATA PSICS( 37) / -.3896419021 5374115954 4913919999 99 D-27    /
      DATA PSICS( 38) / +.6685198138 8855302310 6798933333 33 D-28    /
      DATA PSICS( 39) / -.1146998665 4920864872 5299199999 99 D-28    /
      DATA PSICS( 40) / +.1967938588 6541405920 5154133333 33 D-29    /
      DATA PSICS( 41) / -.3376448818 9750979801 9072000000 00 D-30    /
      DATA PSICS( 42) / +.5793070319 3214159246 6773333333 33 D-31    /
      DATA APSICS(  1) / -.8327107910 6929076017 4456932269 D-3        /
      DATA APSICS(  2) / -.4162518421 9273935282 1627121990 D-3        /
      DATA APSICS(  3) / +.1034315609 7874129117 4463193961 D-6        /
      DATA APSICS(  4) / -.1214681841 3590415298 7299556365 D-9        /
      DATA APSICS(  5) / +.3113694319 9835615552 1240278178 D-12       /
      DATA APSICS(  6) / -.1364613371 9317704177 6516100945 D-14       /
      DATA APSICS(  7) / +.9020517513 1541656513 0837974000 D-17       /
      DATA APSICS(  8) / -.8315429974 2159146482 9933635466 D-19       /
      DATA APSICS(  9) / +.1012242570 7390725418 8479482666 D-20       /
      DATA APSICS( 10) / -.1562702494 3562250762 0478933333 D-22       /
      DATA APSICS( 11) / +.2965427168 0890389613 3226666666 D-24       /
      DATA APSICS( 12) / -.6746868867 6570216374 1866666666 D-26       /
      DATA APSICS( 13) / +.1803453116 9718990421 3333333333 D-27       /
      DATA APSICS( 14) / -.5569016182 4598360746 6666666666 D-29       /
      DATA APSICS( 15) / +.1958679226 0773625173 3333333333 D-30       /
      DATA APSICS( 16) / -.7751958925 2333568000 0000000000 D-32       /
      DATA PI / 3.1415926535 8979323846 2643383279 50 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DPSI
      IF (FIRST) THEN
         NTPSI = INITDS (PSICS, 42, 0.1*DBLE(D1MACH(3)) )
         NTAPSI = INITDS (APSICS, 16, 0.1*DBLE(D1MACH(3)) )
C
         XBIG = 1.0D0/SQRT(D1MACH(3))
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
C
      IF (Y.GT.10.0D0) GO TO 50
C
C DPSI(X) FOR ABS(X) .LE. 2
C
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - N
      N = N - 1
      DPSI = DCSEVL (2.D0*Y-1.D0, PSICS, NTPSI)
      IF (N.EQ.0) RETURN
C
      IF (N.GT.0) GO TO 30
C
      N = -N
      IF (X .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DPSI', 'X IS 0', 2, 2)
      IF (X .LT. 0.D0 .AND. X+N-2 .EQ. 0.D0) CALL XERMSG ('SLATEC',
     +   'DPSI', 'X IS A NEGATIVE INTEGER', 3, 2)
      IF (X .LT. (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X) .LT. DXREL)
     +   CALL XERMSG ('SLATEC', 'DPSI',
     +   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',
     +   1, 1)
C
      DO 20 I=1,N
        DPSI = DPSI - 1.D0/(X+I-1)
 20   CONTINUE
      RETURN
C
C DPSI(X) FOR X .GE. 2.0 AND X .LE. 10.0
C
 30   DO 40 I=1,N
        DPSI = DPSI + 1.0D0/(Y+I)
 40   CONTINUE
      RETURN
C
C DPSI(X) FOR ABS(X) .GT. 10.0
C
 50   AUX = 0.D0
      IF (Y.LT.XBIG) AUX = DCSEVL (2.D0*(10.D0/Y)**2-1.D0, APSICS,
     1  NTAPSI)
C
      IF (X.LT.0.D0) DPSI = LOG(ABS(X)) - 0.5D0/X + AUX
     1  - PI*DCOT(PI*X)
      IF (X.GT.0.D0) DPSI = LOG(X) - 0.5D0/X + AUX
      RETURN
C
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK DCOT
      DOUBLE PRECISION FUNCTION DCOT (X)
C***BEGIN PROLOGUE  DCOT
C***PURPOSE  Compute the cotangent.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C4A
C***TYPE      DOUBLE PRECISION (COT-S, DCOT-D, CCOT-C)
C***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C DCOT(X) calculates the double precision trigonometric cotangent
C for double precision argument X.  X is in units of radians.
C
C Series for COT        on the interval  0.          to  6.25000E-02
C                                        with weighted error   5.52E-34
C                                         log weighted error  33.26
C                               significant figures required  32.34
C                                    decimal places required  33.85
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920618  Removed space from variable names.  (RWC, WRB)
C***END PROLOGUE  DCOT
      DOUBLE PRECISION X, COTCS(15), AINTY, AINTY2, PI2REC, SQEPS,
     1  XMAX, XMIN, XSML, Y, YREM, PRODBG, DCSEVL, D1MACH
      LOGICAL FIRST
      SAVE COTCS, PI2REC, NTERMS, XMAX, XSML, XMIN, SQEPS, FIRST
      DATA COTCS(  1) / +.2402591609 8295630250 9553617744 970 D+0    /
      DATA COTCS(  2) / -.1653303160 1500227845 4746025255 758 D-1    /
      DATA COTCS(  3) / -.4299839193 1724018935 6476228239 895 D-4    /
      DATA COTCS(  4) / -.1592832233 2754104602 3490851122 445 D-6    /
      DATA COTCS(  5) / -.6191093135 1293487258 8620579343 187 D-9    /
      DATA COTCS(  6) / -.2430197415 0726460433 1702590579 575 D-11   /
      DATA COTCS(  7) / -.9560936758 8000809842 7062083100 000 D-14   /
      DATA COTCS(  8) / -.3763537981 9458058041 6291539706 666 D-16   /
      DATA COTCS(  9) / -.1481665746 4674657885 2176794666 666 D-18   /
      DATA COTCS( 10) / -.5833356589 0366657947 7984000000 000 D-21   /
      DATA COTCS( 11) / -.2296626469 6464577392 8533333333 333 D-23   /
      DATA COTCS( 12) / -.9041970573 0748332671 9999999999 999 D-26   /
      DATA COTCS( 13) / -.3559885519 2060006400 0000000000 000 D-28   /
      DATA COTCS( 14) / -.1401551398 2429866666 6666666666 666 D-30   /
      DATA COTCS( 15) / -.5518004368 7253333333 3333333333 333 D-33   /
      DATA PI2REC / .01161977236 7581343075 5350534900 57 D0 /
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DCOT
      IF (FIRST) THEN
         NTERMS = INITDS (COTCS, 15, 0.1*DBLE(D1MACH(3)) )
         XMAX = 1.0D0/D1MACH(4)
         XSML = SQRT(3.0D0*D1MACH(3))
         XMIN = EXP (MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + 0.01D0)
         SQEPS = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
C
      Y = ABS(X)
      IF (Y .LT. XMIN) CALL XERMSG ('SLATEC', 'DCOT',
     +   'ABS(X) IS ZERO OR SO SMALL DCOT OVERFLOWS', 2, 2)
      IF (Y .GT. XMAX) CALL XERMSG ('SLATEC', 'DCOT',
     +   'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 3, 2)
C
C CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC)
C = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z
C = AINT(.625*Y) + AINT(Z) + REM(Z)
C
      AINTY = AINT (Y)
      YREM = Y - AINTY
      PRODBG = 0.625D0*AINTY
      AINTY = AINT (PRODBG)
      Y = (PRODBG-AINTY) + 0.625D0*YREM + PI2REC*Y
      AINTY2 = AINT (Y)
      AINTY = AINTY + AINTY2
      Y = Y - AINTY2
C
      IFN = MOD (AINTY, 2.0D0)
      IF (IFN.EQ.1) Y = 1.0D0 - Y
C
      IF (ABS(X) .GT. 0.5D0 .AND. Y .LT. ABS(X)*SQEPS) CALL XERMSG
     +   ('SLATEC', 'DCOT',
     +   'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI ' //
     +   '(N.NE.0)', 1, 1)
C
      IF (Y.GT.0.25D0) GO TO 20
      DCOT = 1.0D0/X
      IF (Y.GT.XSML) DCOT = (0.5D0 + DCSEVL (32.0D0*Y*Y-1.D0, COTCS,
     1  NTERMS)) / Y
      GO TO 40
C
 20   IF (Y.GT.0.5D0) GO TO 30
      DCOT = (0.5D0 + DCSEVL (8.D0*Y*Y-1.D0, COTCS, NTERMS))/(0.5D0*Y)
      DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
      GO TO 40
C
 30   DCOT = (0.5D0 + DCSEVL (2.D0*Y*Y-1.D0, COTCS, NTERMS))/(.25D0*Y)
      DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
      DCOT = (DCOT*DCOT-1.D0)*0.5D0/DCOT
C
 40   IF (X.NE.0.D0) DCOT = SIGN (DCOT, X)
      IF (IFN.EQ.1) DCOT = -DCOT
C
      RETURN
      END


C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK DCSEVL
      DOUBLE PRECISION FUNCTION DCSEVL (X, CS, N)
C***BEGIN PROLOGUE  DCSEVL
C***PURPOSE  Evaluate a Chebyshev series.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
C***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Evaluate the N-term Chebyshev series CS at X.  Adapted from
C  a method presented in the paper by Broucke referenced below.
C
C       Input Arguments --
C  X    value at which the series is to be evaluated.
C  CS   array of N terms of a Chebyshev series.  In evaluating
C       CS, only half the first coefficient is summed.
C  N    number of terms in array CS.
C
C***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
C                 Chebyshev series, Algorithm 446, Communications of
C                 the A.C.M. 16, (1973) pp. 254-256.
C               L. Fox and I. B. Parker, Chebyshev Polynomials in
C                 Numerical Analysis, Oxford University Press, 1968,
C                 page 56.
C***ROUTINES CALLED  D1MACH, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770401  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900329  Prologued revised extensively and code rewritten to allow
C           X to be slightly outside interval (-1,+1).  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DCSEVL
      DOUBLE PRECISION B0, B1, B2, CS(*), ONEPL, TWOX, X, D1MACH
      LOGICAL FIRST
      SAVE FIRST, ONEPL
      DATA FIRST /.TRUE./
C***FIRST EXECUTABLE STATEMENT  DCSEVL
      IF (FIRST) ONEPL = 1.0D0 + D1MACH(4)
      FIRST = .FALSE.
      IF (N .LT. 1) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'NUMBER OF TERMS .LE. 0', 2, 2)
      IF (N .GT. 1000) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'NUMBER OF TERMS .GT. 1000', 3, 2)
      IF (ABS(X) .GT. ONEPL) CALL XERMSG ('SLATEC', 'DCSEVL',
     +   'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
C
      B1 = 0.0D0
      B0 = 0.0D0
      TWOX = 2.0D0*X
      DO 10 I = 1,N
         B2 = B1
         B1 = B0
         NI = N + 1 - I
         B0 = TWOX*B1 - B2 + CS(NI)
   10 CONTINUE
C
      DCSEVL = 0.5D0*(B0-B2)
C
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK INITDS
      FUNCTION INITDS (OS, NOS, ETA)
C***BEGIN PROLOGUE  INITDS
C***PURPOSE  Determine the number of terms needed in an orthogonal
C            polynomial series so that it meets a specified accuracy.
C***LIBRARY   SLATEC (FNLIB)
C***CATEGORY  C3A2
C***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
C***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
C             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
C***AUTHOR  Fullerton, W., (LANL)
C***DESCRIPTION
C
C  Initialize the orthogonal series, represented by the array OS, so
C  that INITDS is the number of terms needed to insure the error is no
C  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
C  machine precision.
C
C             Input Arguments --
C   OS     double precision array of NOS coefficients in an orthogonal
C          series.
C   NOS    number of coefficients in OS.
C   ETA    single precision scalar containing requested accuracy of
C          series.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   770601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891115  Modified error message.  (WRB)
C   891115  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  INITDS
      DOUBLE PRECISION OS(*)
C***FIRST EXECUTABLE STATEMENT  INITDS
      IF (NOS .LT. 1) CALL XERMSG ('SLATEC', 'INITDS',
     +   'Number of coefficients is less than 1', 2, 1)
C
      ERR = 0.
      DO 10 II = 1,NOS
        I = NOS + 1 - II
        ERR = ERR + DABS(DBLE(OS(I))) !ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
C
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITDS',
     +   'Chebyshev series too short for specified accuracy', 1, 1)
      INITDS = I
C
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      DOUBLE PRECISION FUNCTION ERFINV(X)
      IMPLICIT NONE
      DOUBLE PRECISION X,XX,P,Q,MEAN,SD,BOUND
      INTEGER WHICH,STATUS

      IF((X.LT.-1.0D0).OR.(X.GT.1.0D0)) THEN
         WRITE(*,*)'ERROR IN FUNCTION ERFINV(X)'
         WRITE(*,*)'X MUST BE IN ]-1,+1['
         IF(X.LT.-1.0D0) WRITE(*,*)'IN THIS CALL X < -1'
         IF(X.GT.1.0D0)  WRITE(*,*)'IN THIS CALL X > +1'
         WRITE(*,*)'COMPUTATIONS ABORTED'
         STOP
      ELSEIF(X.EQ.1.0D0) THEN
         ERFINV = HUGE(X)
      ELSEIF(X.EQ.-1.0D0) THEN
         ERFINV = -HUGE(X)
      ELSE
         WHICH = 2
         P = (X+1.0D0)/2.0D0
         Q = 1.0D0-P
         XX = X
         SD = 1.0D0
         MEAN = 0.0D0
         CALL CDFNOR(WHICH,P,Q,XX,MEAN,SD,STATUS,BOUND)     
         ERFINV = XX/DSQRT(2.0D0)
      ENDIF
      RETURN
      END


C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE CDFNOR(WHICH,P,Q,X,MEAN,SD,STATUS,BOUND)
C**********************************************************************
C
C      SUBROUTINE CDFNOR( WHICH, P, Q, X, MEAN, SD, STATUS, BOUND )
C               CUMULATIVE DISTRIBUTION FUNCTION
C               NORMAL DISTRIBUTION
C
C
C                              FUNCTION
C
C
C     CALCULATES ANY ONE PARAMETER OF THE NORMAL
C     DISTRIBUTION GIVEN VALUES FOR THE OTHERS.
C
C
C                              ARGUMENTS
C
C
C     WHICH  --> INTEGER INDICATING  WHICH OF THE  NEXT  PARAMETER
C     VALUES IS TO BE CALCULATED USING VALUES  OF THE OTHERS.
C     LEGAL RANGE: 1..4
C               IWHICH = 1 : CALCULATE P AND Q FROM X,MEAN AND SD
C               IWHICH = 2 : CALCULATE X FROM P,Q,MEAN AND SD
C               IWHICH = 3 : CALCULATE MEAN FROM P,Q,X AND SD
C               IWHICH = 4 : CALCULATE SD FROM P,Q,X AND MEAN
C                    INTEGER WHICH
C
C     P <--> THE INTEGRAL FROM -INFINITY TO X OF THE NORMAL DENSITY.
C            INPUT RANGE: (0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            INPUT RANGE: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     X < --> UPPER LIMIT OF INTEGRATION OF THE NORMAL-DENSITY.
C             INPUT RANGE: ( -INFINITY, +INFINITY)
C                    DOUBLE PRECISION X
C
C     MEAN <--> THE MEAN OF THE NORMAL DENSITY.
C               INPUT RANGE: (-INFINITY, +INFINITY)
C                    DOUBLE PRECISION MEAN
C
C     SD <--> STANDARD DEVIATION OF THE NORMAL DENSITY.
C             INPUT RANGE: (0, +INFINITY).
C                    DOUBLE PRECISION SD
C
C     STATUS <-- 0 IF CALCULATION COMPLETED CORRECTLY
C               -I IF INPUT PARAMETER NUMBER I IS OUT OF RANGE
C                1 IF ANSWER APPEARS TO BE LOWER THAN LOWEST
C                  SEARCH BOUND
C                2 IF ANSWER APPEARS TO BE HIGHER THAN GREATEST
C                  SEARCH BOUND
C                3 IF P + Q .NE. 1
C                    INTEGER STATUS
C
C     BOUND <-- UNDEFINED IF STATUS IS 0
C
C               BOUND EXCEEDED BY PARAMETER NUMBER I IF STATUS
C               IS NEGATIVE.
C
C               LOWER SEARCH BOUND IF STATUS IS 1.
C
C               UPPER SEARCH BOUND IF STATUS IS 2.
C
C
C                              METHOD
C
C
C
C
C     A SLIGHTLY MODIFIED VERSION OF ANORM FROM
C
C     CODY, W.D. (1993). "ALGORITHM 715: SPECFUN - A PORTABEL FORTRAN
C     PACKAGE OF SPECIAL FUNCTION ROUTINES AND TEST DRIVERS"
C     ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE. 19, 22-32.
C
C     IS USED TO CALULATE THE  CUMULATIVE STANDARD NORMAL DISTRIBUTION.
C
C     THE RATIONAL FUNCTIONS FROM PAGES  90-95  OF KENNEDY AND GENTLE,
C     STATISTICAL  COMPUTING,  MARCEL  DEKKER, NY,  1980 ARE  USED  AS
C     STARTING VALUES TO NEWTON'S ITERATIONS WHICH COMPUTE THE INVERSE
C     STANDARD NORMAL.  THEREFORE NO  SEARCHES  ARE NECESSARY FOR  ANY
C     PARAMETER.
C
C     FOR X < -15, THE ASYMPTOTIC EXPANSION FOR THE NORMAL IS USED  AS
C     THE STARTING VALUE IN FINDING THE INVERSE STANDARD NORMAL.
C     THIS IS FORMULA 26.2.12 OF ABRAMOWITZ AND STEGUN.
C
C
C                              NOTE
C
C
C      THE NORMAL DENSITY IS PROPORTIONAL TO
C      EXP( - 0.5 * (( X - MEAN)/SD)**2)
C
C
C**********************************************************************
C     .. PARAMETERS ..
C     ..
C     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION BOUND,MEAN,P,Q,SD,X
      INTEGER STATUS,WHICH
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION Z,PQ
C     ..
C     .. EXTERNAL FUNCTIONS ..

      DOUBLE PRECISION DINVNR,SPMPAR
      EXTERNAL DINVNR,SPMPAR
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL CUMNOR
C     ..
C     .. EXECUTABLE STATEMENTS ..
C
C     CHECK ARGUMENTS
C
      STATUS = 0
      IF (.NOT. ((WHICH.LT.1).OR. (WHICH.GT.4))) GO TO 30
      IF (.NOT. (WHICH.LT.1)) GO TO 10
      BOUND = 1.0D0
      GO TO 20

   10 BOUND = 4.0D0
   20 STATUS = -1
      RETURN

   30 IF (WHICH.EQ.1) GO TO 70
C
C     P
C
      IF (.NOT. ((P.LE.0.0D0).OR. (P.GT.1.0D0))) GO TO 60
      IF (.NOT. (P.LE.0.0D0)) GO TO 40
      BOUND = 0.0D0
      GO TO 50

   40 BOUND = 1.0D0
   50 STATUS = -2
      RETURN

   60 CONTINUE
   70 IF (WHICH.EQ.1) GO TO 110
C
C     Q
C
      IF (.NOT. ((Q.LE.0.0D0).OR. (Q.GT.1.0D0))) GO TO 100
      IF (.NOT. (Q.LE.0.0D0)) GO TO 80
      BOUND = 0.0D0
      GO TO 90

   80 BOUND = 1.0D0
   90 STATUS = -3
      RETURN

  100 CONTINUE
  110 IF (WHICH.EQ.1) GO TO 150
C
C     P + Q
C
      PQ = P + Q
      IF (.NOT. (ABS(((PQ)-0.5D0)-0.5D0).GT.
     +    (3.0D0*SPMPAR(1)))) GO TO 140
      IF (.NOT. (PQ.LT.0.0D0)) GO TO 120
      BOUND = 0.0D0
      GO TO 130

  120 BOUND = 1.0D0
  130 STATUS = 3
      RETURN

  140 CONTINUE
  150 IF (WHICH.EQ.4) GO TO 170
C
C     SD
C
      IF (.NOT. (SD.LE.0.0D0)) GO TO 160
      BOUND = 0.0D0
      STATUS = -6
      RETURN

  160 CONTINUE
C
C     CALCULATE ANSWERS
C
  170 IF ((1).EQ. (WHICH)) THEN
C
C     COMPUTING P
C
          Z = (X-MEAN)/SD
          CALL CUMNOR(Z,P,Q)

      ELSE IF ((2).EQ. (WHICH)) THEN
C
C     COMPUTING X
C
          Z = DINVNR(P,Q)
          X = SD*Z + MEAN

      ELSE IF ((3).EQ. (WHICH)) THEN
C
C     COMPUTING THE MEAN
C
          Z = DINVNR(P,Q)
          MEAN = X - SD*Z

      ELSE IF ((4).EQ. (WHICH)) THEN
C
C     COMPUTING SD
C
          Z = DINVNR(P,Q)
          SD = (X-MEAN)/Z
      END IF

      RETURN

      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION SPMPAR(I)
C-----------------------------------------------------------------------
C
C     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
C     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
C
C        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
C
C        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
C
C        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
C
C-----------------------------------------------------------------------
C     WRITTEN BY
C        ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN VIRGINIA
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
C     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
C     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
C-----------------------------------------------------------------------
C     .. SCALAR ARGUMENTS ..
      INTEGER I
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION B,BINV,BM1,ONE,W,Z
      INTEGER EMAX,EMIN,IBETA,M
C     ..
C     .. EXTERNAL FUNCTIONS ..
      INTEGER IPMPAR
      EXTERNAL IPMPAR
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DBLE
C     ..
C     .. EXECUTABLE STATEMENTS ..
C
      IF (I.GT.1) GO TO 10
      B = IPMPAR(4)
      M = IPMPAR(8)
      SPMPAR = B** (1-M)
      RETURN
C
   10 IF (I.GT.2) GO TO 20
      B = IPMPAR(4)
      EMIN = IPMPAR(9)
      ONE = DBLE(1)
      BINV = ONE/B
      W = B** (EMIN+2)
      SPMPAR = ((W*BINV)*BINV)*BINV
      RETURN
C
   20 IBETA = IPMPAR(4)
      M = IPMPAR(8)
      EMAX = IPMPAR(10)
C
      B = IBETA
      BM1 = IBETA - 1
      ONE = DBLE(1)
      Z = B** (M-1)
      W = ((Z-ONE)*B+BM1)/ (B*Z)
C
      Z = B** (EMAX-2)
      SPMPAR = ((W*Z)*B)*B
      RETURN

      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE CUMNOR(ARG,RESULT,CCUM)
C**********************************************************************
C
C     SUBROUINE CUMNOR(X,RESULT,CCUM)
C
C
C                              FUNCTION
C
C
C     COMPUTES THE CUMULATIVE  OF    THE  NORMAL   DISTRIBUTION,   I.E.,
C     THE INTEGRAL FROM -INFINITY TO X OF
C          (1/SQRT(2*PI)) EXP(-U*U/2) DU
C
C     X --> UPPER LIMIT OF INTEGRATION.
C                                        X IS DOUBLE PRECISION
C
C     RESULT <-- CUMULATIVE NORMAL DISTRIBUTION.
C                                        RESULT IS DOUBLE PRECISION
C
C     CCUM <-- COMPLIMENT OF CUMULATIVE NORMAL DISTRIBUTION.
C                                        CCUM IS DOUBLE PRECISION
C
C
C     RENAMING OF FUNCTION ANORM FROM:
C
C     CODY, W.D. (1993). "ALGORITHM 715: SPECFUN - A PORTABEL FORTRAN
C     PACKAGE OF SPECIAL FUNCTION ROUTINES AND TEST DRIVERS"
C     ACM TRANSACTIONS ON MATHEMATICAL SOFTWARE. 19, 22-32.
C
C     WITH SLIGHT MODIFICATIONS TO RETURN CCUM AND TO DEAL WITH
C     MACHINE CONSTANTS.
C
C**********************************************************************
C
C
C ORIGINAL COMMENTS:
C------------------------------------------------------------------
C
C THIS FUNCTION EVALUATES THE NORMAL DISTRIBUTION FUNCTION:
C
C                              / X
C                     1       |       -T*T/2
C          P(X) = ----------- |      E       DT
C                 SQRT(2 PI)  |
C                             /-OO
C
C   THE MAIN COMPUTATION EVALUATES NEAR-MINIMAX APPROXIMATIONS
C   DERIVED FROM THOSE IN "RATIONAL CHEBYSHEV APPROXIMATIONS FOR
C   THE ERROR FUNCTION" BY W. J. CODY, MATH. COMP., 1969, 631-637.
C   THIS TRANSPORTABLE PROGRAM USES RATIONAL FUNCTIONS THAT
C   THEORETICALLY APPROXIMATE THE NORMAL DISTRIBUTION FUNCTION TO
C   AT LEAST 18 SIGNIFICANT DECIMAL DIGITS.  THE ACCURACY ACHIEVED
C   DEPENDS ON THE ARITHMETIC SYSTEM, THE COMPILER, THE INTRINSIC
C   FUNCTIONS, AND PROPER SELECTION OF THE MACHINE-DEPENDENT
C   CONSTANTS.
C
C*******************************************************************
C*******************************************************************
C
C EXPLANATION OF MACHINE-DEPENDENT CONSTANTS.
C
C   MIN   = SMALLEST MACHINE REPRESENTABLE NUMBER.
C
C   EPS   = ARGUMENT BELOW WHICH ANORM(X) MAY BE REPRESENTED BY
C           0.5  AND ABOVE WHICH  X*X  WILL NOT UNDERFLOW.
C           A CONSERVATIVE VALUE IS THE LARGEST MACHINE NUMBER X
C           SUCH THAT   1.0 + X = 1.0   TO MACHINE PRECISION.
C*******************************************************************
C*******************************************************************
C
C ERROR RETURNS
C
C  THE PROGRAM RETURNS  ANORM = 0     FOR  ARG .LE. XLOW.
C
C
C INTRINSIC FUNCTIONS REQUIRED ARE:
C
C     ABS, AINT, EXP
C
C
C  AUTHOR: W. J. CODY
C          MATHEMATICS AND COMPUTER SCIENCE DIVISION
C          ARGONNE NATIONAL LABORATORY
C          ARGONNE, IL 60439
C
C  LATEST MODIFICATION: MARCH 15, 1992
C
C------------------------------------------------------------------
      INTEGER I
      DOUBLE PRECISION A,ARG,B,C,D,DEL,EPS,HALF,P,ONE,Q,RESULT,SIXTEN,
     +                 TEMP,SQRPI,THRSH,ROOT32,X,XDEN,XNUM,Y,XSQ,ZERO,
     +                 MIN,CCUM
      DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
C------------------------------------------------------------------
C  EXTERNAL FUNCTION
C------------------------------------------------------------------
      DOUBLE PRECISION SPMPAR
      EXTERNAL SPMPAR
C------------------------------------------------------------------
C  MATHEMATICAL CONSTANTS
C
C  SQRPI = 1 / SQRT(2*PI), ROOT32 = SQRT(32), AND
C  THRSH IS THE ARGUMENT FOR WHICH ANORM = 0.75.
C------------------------------------------------------------------
      DATA ONE,HALF,ZERO,SIXTEN/1.0D0,0.5D0,0.0D0,1.60D0/,
     +     SQRPI/3.9894228040143267794D-1/,THRSH/0.66291D0/,
     +     ROOT32/5.656854248D0/
C------------------------------------------------------------------
C  COEFFICIENTS FOR APPROXIMATION IN FIRST INTERVAL
C------------------------------------------------------------------
      DATA A/2.2352520354606839287D00,1.6102823106855587881D02,
     +     1.0676894854603709582D03,1.8154981253343561249D04,
     +     6.5682337918207449113D-2/
      DATA B/4.7202581904688241870D01,9.7609855173777669322D02,
     +     1.0260932208618978205D04,4.5507789335026729956D04/
C------------------------------------------------------------------
C  COEFFICIENTS FOR APPROXIMATION IN SECOND INTERVAL
C------------------------------------------------------------------
      DATA C/3.9894151208813466764D-1,8.8831497943883759412D00,
     +     9.3506656132177855979D01,5.9727027639480026226D02,
     +     2.4945375852903726711D03,6.8481904505362823326D03,
     +     1.1602651437647350124D04,9.8427148383839780218D03,
     +     1.0765576773720192317D-8/
      DATA D/2.2266688044328115691D01,2.3538790178262499861D02,
     +     1.5193775994075548050D03,6.4855582982667607550D03,
     +     1.8615571640885098091D04,3.4900952721145977266D04,
     +     3.8912003286093271411D04,1.9685429676859990727D04/
C------------------------------------------------------------------
C  COEFFICIENTS FOR APPROXIMATION IN THIRD INTERVAL
C------------------------------------------------------------------
      DATA P/2.1589853405795699D-1,1.274011611602473639D-1,
     +     2.2235277870649807D-2,1.421619193227893466D-3,
     +     2.9112874951168792D-5,2.307344176494017303D-2/
      DATA Q/1.28426009614491121D00,4.68238212480865118D-1,
     +     6.59881378689285515D-2,3.78239633202758244D-3,
     +     7.29751555083966205D-5/
C------------------------------------------------------------------
C  MACHINE DEPENDENT CONSTANTS
C------------------------------------------------------------------
      EPS = SPMPAR(1)*0.5D0
      MIN = SPMPAR(2)
C------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y.LE.THRSH) THEN
C------------------------------------------------------------------
C  EVALUATE  ANORM  FOR  |X| <= 0.66291
C------------------------------------------------------------------
          XSQ = ZERO
          IF (Y.GT.EPS) XSQ = X*X
          XNUM = A(5)*XSQ
          XDEN = XSQ
          DO 10 I = 1,3
              XNUM = (XNUM+A(I))*XSQ
              XDEN = (XDEN+B(I))*XSQ
   10     CONTINUE
          RESULT = X* (XNUM+A(4))/ (XDEN+B(4))
          TEMP = RESULT
          RESULT = HALF + TEMP
          CCUM = HALF - TEMP
C------------------------------------------------------------------
C  EVALUATE  ANORM  FOR 0.66291 <= |X| <= SQRT(32)
C------------------------------------------------------------------
      ELSE IF (Y.LE.ROOT32) THEN
          XNUM = C(9)*Y
          XDEN = Y
          DO 20 I = 1,7
              XNUM = (XNUM+C(I))*Y
              XDEN = (XDEN+D(I))*Y
   20     CONTINUE
          RESULT = (XNUM+C(8))/ (XDEN+D(8))
          XSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-XSQ)* (Y+XSQ)
          RESULT = EXP(-XSQ*XSQ*HALF)*EXP(-DEL*HALF)*RESULT
          CCUM = ONE - RESULT
          IF (X.GT.ZERO) THEN
              TEMP = RESULT
              RESULT = CCUM
              CCUM = TEMP
          END IF
C------------------------------------------------------------------
C  EVALUATE  ANORM  FOR |X| > SQRT(32)
C------------------------------------------------------------------
      ELSE
          RESULT = ZERO
          XSQ = ONE/ (X*X)
          XNUM = P(6)*XSQ
          XDEN = XSQ
          DO 30 I = 1,4
              XNUM = (XNUM+P(I))*XSQ
              XDEN = (XDEN+Q(I))*XSQ
   30     CONTINUE
          RESULT = XSQ* (XNUM+P(5))/ (XDEN+Q(5))
          RESULT = (SQRPI-RESULT)/Y
          XSQ = AINT(X*SIXTEN)/SIXTEN
          DEL = (X-XSQ)* (X+XSQ)
          RESULT = EXP(-XSQ*XSQ*HALF)*EXP(-DEL*HALF)*RESULT
          CCUM = ONE - RESULT
          IF (X.GT.ZERO) THEN
              TEMP = RESULT
              RESULT = CCUM
              CCUM = TEMP
          END IF

      END IF

      IF (RESULT.LT.MIN) RESULT = 0.0D0
      IF (CCUM.LT.MIN) CCUM = 0.0D0
C------------------------------------------------------------------
C  FIX UP FOR NEGATIVE ARGUMENT, ERF, ETC.
C------------------------------------------------------------------
C----------LAST CARD OF ANORM ----------
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION DINVNR(P,Q)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DINVNR(P,Q)
C     DOUBLE PRECISION NORMAL DISTRIBUTION INVERSE
C
C
C                              FUNCTION
C
C
C     RETURNS X  SUCH THAT CUMNOR(X)  =   P,  I.E., THE  INTEGRAL FROM -
C     INFINITY TO X OF (1/SQRT(2*PI)) EXP(-U*U/2) DU IS P
C
C
C                              ARGUMENTS
C
C
C     P --> THE PROBABILITY WHOSE NORMAL DEVIATE IS SOUGHT.
C                    P IS DOUBLE PRECISION
C
C     Q --> 1-P
C                    P IS DOUBLE PRECISION
C
C
C                              METHOD
C
C
C     THE  RATIONAL   FUNCTION   ON  PAGE 95    OF KENNEDY  AND  GENTLE,
C     STATISTICAL COMPUTING, MARCEL DEKKER, NY , 1980 IS USED AS A START
C     VALUE FOR THE NEWTON METHOD OF FINDING ROOTS.
C
C
C                              NOTE
C
C
C     IF P OR Q .LT. MACHINE EPS RETURNS +/- DINVNR(EPS)
C
C**********************************************************************
C     .. PARAMETERS ..
      INTEGER MAXIT
      PARAMETER (MAXIT=100)
      DOUBLE PRECISION EPS
      PARAMETER (EPS=1.0D-13)
      DOUBLE PRECISION R2PI
      PARAMETER (R2PI=0.3989422804014326D0)
      DOUBLE PRECISION NHALF
      PARAMETER (NHALF=-0.5D0)
C     ..
C     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION P,Q
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION STRTX,XCUR,CUM,CCUM,PP,DX
      INTEGER I
      LOGICAL QPORQ
C     ..
C     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION STVALN
      EXTERNAL STVALN
C     ..
C     .. EXTERNAL SUBROUTINES ..
      EXTERNAL CUMNOR
C     ..
C     .. STATEMENT FUNCTIONS ..
      DOUBLE PRECISION DENNOR,X

      DENNOR(X) = R2PI*EXP(NHALF*X*X)
C     ..
C     .. EXECUTABLE STATEMENTS ..
C
C     FIND MINIMUM OF P AND Q
C
      QPORQ = P .LE. Q
      IF (.NOT. (QPORQ)) GO TO 10
      PP = P
      GO TO 20

   10 PP = Q
C
C     INITIALIZATION STEP
C
   20 STRTX = STVALN(PP)
      XCUR = STRTX
C
C     NEWTON INTERATIONS
C
      DO 30,I = 1,MAXIT
          CALL CUMNOR(XCUR,CUM,CCUM)
          DX = (CUM-PP)/DENNOR(XCUR)
          XCUR = XCUR - DX
          IF (ABS(DX/XCUR).LT.EPS) GO TO 40
   30 CONTINUE
      DINVNR = STRTX
C
C     IF WE GET HERE, NEWTON HAS FAILED
C
      IF (.NOT.QPORQ) DINVNR = -DINVNR
      RETURN
C
C     IF WE GET HERE, NEWTON HAS SUCCEDED
C
   40 DINVNR = XCUR
      IF (.NOT.QPORQ) DINVNR = -DINVNR
      RETURN

      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      INTEGER FUNCTION IPMPAR(I)
C-----------------------------------------------------------------------
C
C     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
C     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
C     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
C
C  INTEGERS.
C
C     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
C
C               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
C
C     IPMPAR(1) = A, THE BASE.
C
C     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
C
C     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
C     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
C     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
C
C               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
C
C               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
C               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
C
C     IPMPAR(4) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
C
C-----------------------------------------------------------------------
C
C     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED, ACTIVATE
C     THE DATA STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM
C     COLUMN 1. (ALL THE OTHER DATA STATEMENTS SHOULD HAVE C IN
C     COLUMN 1.)
C
C-----------------------------------------------------------------------
C
C     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
C     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
C     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
C     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
C
C-----------------------------------------------------------------------
C     .. SCALAR ARGUMENTS ..
      INTEGER I
C     ..
C     .. LOCAL ARRAYS ..
      INTEGER IMACH(10)
C     ..
C     .. DATA STATEMENTS ..
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C     DATA IMACH( 1) /   2 /
C     DATA IMACH( 2) /  31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /  16 /
C     DATA IMACH( 5) /   6 /
C     DATA IMACH( 6) / -64 /
C     DATA IMACH( 7) /  63 /
C     DATA IMACH( 8) /  14 /
C     DATA IMACH( 9) / -64 /
C     DATA IMACH(10) /  63 /
C
C     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
C     PC 7300, AND AT&T 6300.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   33 /
C     DATA IMACH( 3) / 8589934591 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -256 /
C     DATA IMACH( 7) /  255 /
C     DATA IMACH( 8) /   60 /
C     DATA IMACH( 9) / -256 /
C     DATA IMACH(10) /  255 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   39 /
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /    8 /
C     DATA IMACH( 5) /   13 /
C     DATA IMACH( 6) /  -50 /
C     DATA IMACH( 7) /   76 /
C     DATA IMACH( 8) /   26 /
C     DATA IMACH( 9) /  -50 /
C     DATA IMACH(10) /   76 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C     DATA IMACH( 1) /      2 /
C     DATA IMACH( 2) /     39 /
C     DATA IMACH( 3) / 549755813887 /
C     DATA IMACH( 4) /      8 /
C     DATA IMACH( 5) /     13 /
C     DATA IMACH( 6) /    -50 /
C     DATA IMACH( 7) /     76 /
C     DATA IMACH( 8) /     26 /
C     DATA IMACH( 9) / -32754 /
C     DATA IMACH(10) /  32780 /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
C     ARITHMETIC (NOS OPERATING SYSTEM).
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   48 /
C     DATA IMACH( 3) / 281474976710655 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   48 /
C     DATA IMACH( 6) / -974 /
C     DATA IMACH( 7) / 1070 /
C     DATA IMACH( 8) /   95 /
C     DATA IMACH( 9) / -926 /
C     DATA IMACH(10) / 1070 /
C
C     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
C     ARITHMETIC (NOS/VE OPERATING SYSTEM).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    48 /
C     DATA IMACH( 6) / -4096 /
C     DATA IMACH( 7) /  4095 /
C     DATA IMACH( 8) /    96 /
C     DATA IMACH( 9) / -4096 /
C     DATA IMACH(10) /  4095 /
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    63 /
C     DATA IMACH( 3) / 9223372036854775807 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    47 /
C     DATA IMACH( 6) / -8189 /
C     DATA IMACH( 7) /  8190 /
C     DATA IMACH( 8) /    94 /
C     DATA IMACH( 9) / -8099 /
C     DATA IMACH(10) /  8190 /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   15 /
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /   16 /
C     DATA IMACH( 5) /    6 /
C     DATA IMACH( 6) /  -64 /
C     DATA IMACH( 7) /   63 /
C     DATA IMACH( 8) /   14 /
C     DATA IMACH( 9) /  -64 /
C     DATA IMACH(10) /   63 /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   23 /
C     DATA IMACH( 3) / 8388607 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   23 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   38 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
C     AND DPS 8/70 SERIES.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   63 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   15 /
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   23 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   39 /
C     DATA IMACH( 9) / -128 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   15 /
C     DATA IMACH( 3) / 32767 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   23 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   55 /
C     DATA IMACH( 9) / -128 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE HP 9000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -126 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
C     5/7/9 AND THE SEL SYSTEMS 85/86.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /   16 /
C     DATA IMACH( 5) /    6 /
C     DATA IMACH( 6) /  -64 /
C     DATA IMACH( 7) /   63 /
C     DATA IMACH( 8) /   14 /
C     DATA IMACH( 9) /  -64 /
C     DATA IMACH(10) /   63 /
C
C     MACHINE CONSTANTS FOR THE IBM PC.
C
C      DATA IMACH(1)/2/
C      DATA IMACH(2)/31/
C      DATA IMACH(3)/2147483647/
C      DATA IMACH(4)/2/
C      DATA IMACH(5)/24/
C      DATA IMACH(6)/-125/
C      DATA IMACH(7)/128/
C      DATA IMACH(8)/53/
C      DATA IMACH(9)/-1021/
C      DATA IMACH(10)/1024/
C
C     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
C     MACFORTRAN II.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   56 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   54 /
C     DATA IMACH( 9) / -101 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   62 /
C     DATA IMACH( 9) / -128 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   56 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
C     SERIES (MIPS R3000 PROCESSOR).
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
      DATA IMACH( 1) /     2 /
      DATA IMACH( 2) /    31 /
      DATA IMACH( 3) / 2147483647 /
      DATA IMACH( 4) /     2 /
      DATA IMACH( 5) /    24 /
      DATA IMACH( 6) /  -125 /
      DATA IMACH( 7) /   128 /
      DATA IMACH( 8) /    53 /
      DATA IMACH( 9) / -1021 /
      DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   35 /
C     DATA IMACH( 3) / 34359738367 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   27 /
C     DATA IMACH( 6) / -128 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   60 /
C     DATA IMACH( 9) /-1024 /
C     DATA IMACH(10) / 1023 /
C
C     MACHINE CONSTANTS FOR THE VAX 11/780.
C
C     DATA IMACH( 1) /    2 /
C     DATA IMACH( 2) /   31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /    2 /
C     DATA IMACH( 5) /   24 /
C     DATA IMACH( 6) / -127 /
C     DATA IMACH( 7) /  127 /
C     DATA IMACH( 8) /   56 /
C     DATA IMACH( 9) / -127 /
C     DATA IMACH(10) /  127 /
C
      IPMPAR = IMACH(I)
      RETURN

      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION STVALN(P)
C
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION STVALN(P)
C                    STARTING VALUE FOR NETON-RAPHON
C                CALCULATION OF NORMAL DISTRIBUTION INVERSE
C
C
C                              FUNCTION
C
C
C     RETURNS X  SUCH THAT CUMNOR(X)  =   P,  I.E., THE  INTEGRAL FROM -
C     INFINITY TO X OF (1/SQRT(2*PI)) EXP(-U*U/2) DU IS P
C
C
C                              ARGUMENTS
C
C
C     P --> THE PROBABILITY WHOSE NORMAL DEVIATE IS SOUGHT.
C                    P IS DOUBLE PRECISION
C
C
C                              METHOD
C
C
C     THE  RATIONAL   FUNCTION   ON  PAGE 95    OF KENNEDY  AND  GENTLE,
C     STATISTICAL COMPUTING, MARCEL DEKKER, NY , 1980.
C
C**********************************************************************
C
C     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION P
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION SIGN,Y,Z
C     ..
C     .. LOCAL ARRAYS ..
      DOUBLE PRECISION XDEN(5),XNUM(5)
C     ..
C     .. EXTERNAL FUNCTIONS ..
      DOUBLE PRECISION DEVLPL
      EXTERNAL DEVLPL
C     ..
C     .. INTRINSIC FUNCTIONS ..
      INTRINSIC DBLE,LOG,SQRT
C     ..
C     .. DATA STATEMENTS ..
      DATA XNUM/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,
     +     -0.204231210245D-1,-0.453642210148D-4/
      DATA XDEN/0.993484626060D-1,0.588581570495D0,0.531103462366D0,
     +     0.103537752850D0,0.38560700634D-2/
C     ..
C     .. EXECUTABLE STATEMENTS ..
      IF (.NOT. (P.LE.0.5D0)) GO TO 10
      SIGN = -1.0D0
      Z = P
      GO TO 20

   10 SIGN = 1.0D0
      Z = 1.0D0 - P
   20 Y = SQRT(-2.0D0*LOG(Z))
      STVALN = Y + DEVLPL(XNUM,5,Y)/DEVLPL(XDEN,5,Y)
      STVALN = SIGN*STVALN
      RETURN

      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)
C              DOUBLE PRECISION EVALUATE A POLYNOMIAL AT X
C
C
C                              FUNCTION
C
C
C     RETURNS
C          A(1) + A(2)*X + ... + A(N)*X**(N-1)
C
C
C                              ARGUMENTS
C
C
C     A --> ARRAY OF COEFFICIENTS OF THE POLYNOMIAL.
C                                        A IS DOUBLE PRECISION(N)
C
C     N --> LENGTH OF A, ALSO DEGREE OF POLYNOMIAL - 1.
C                                        N IS INTEGER
C
C     X --> POINT AT WHICH THE POLYNOMIAL IS TO BE EVALUATED.
C                                        X IS DOUBLE PRECISION
C
C**********************************************************************
C
C     .. SCALAR ARGUMENTS ..
      DOUBLE PRECISION X
      INTEGER N
C     ..
C     .. ARRAY ARGUMENTS ..
      DOUBLE PRECISION A(N)
C     ..
C     .. LOCAL SCALARS ..
      DOUBLE PRECISION TERM
      INTEGER I
C     ..
C     .. EXECUTABLE STATEMENTS ..
      TERM = A(N)
      DO 10,I = N - 1,1,-1
          TERM = A(I) + TERM*X
   10 CONTINUE
      DEVLPL = TERM
      RETURN

      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================


C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

!DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
      IMPLICIT NONE
      INTEGER :: I
      DOUBLE PRECISION :: B, X
!***BEGIN PROLOGUE  D1MACH
!***PURPOSE  Return floating point machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      SINGLE PRECISION (D1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        A = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of A above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH(3) = B**(-T), the smallest relative spacing.
!   D1MACH(4) = B**(1-T), the largest relative spacing.
!   D1MACH(5) = LOG10(B)
!
!   Assume single precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
!   EMIN .LE. E .LE. EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(11) = T, the number of base-B digits.
!   I1MACH(12) = EMIN, the smallest exponent E.
!   I1MACH(13) = EMAX, the largest exponent E.
!
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   960329  Modified for Fortran 90 (BE after suggestions by EHG)      
!***END PROLOGUE  D1MACH
!      
      X = 1.0D0
      B = RADIX(X)
      SELECT CASE (I)
        CASE (1)
          D1MACH = B**(MINEXPONENT(X)-1) ! the smallest positive magnitude.
        CASE (2)
          D1MACH = HUGE(X)               ! the largest magnitude.
        CASE (3)
          D1MACH = B**(-DIGITS(X))       ! the smallest relative spacing.
        CASE (4)
          D1MACH = B**(1-DIGITS(X))      ! the largest relative spacing.
        CASE (5)
          D1MACH = LOG10(B)
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN D1MACH - I OUT OF BOUNDS')
          STOP
      END SELECT
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

!DECK I1MACH
      INTEGER FUNCTION I1MACH (I)
      IMPLICIT NONE
      INTEGER :: I
      REAL :: X
      DOUBLE PRECISION :: XX
!***BEGIN PROLOGUE  I1MACH
!***PURPOSE  Return integer machine dependent constants.
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      INTEGER (I1MACH-I)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   I1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument and can be referenced as follows:
!
!        K = I1MACH(I)
!
!   where I=1,...,16.  The (output) value of K above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   I/O unit numbers:
!     I1MACH( 1) = the standard input unit.
!     I1MACH( 2) = the standard output unit.
!     I1MACH( 3) = the standard punch unit.
!     I1MACH( 4) = the standard error message unit.
!
!   Words:
!     I1MACH( 5) = the number of bits per integer storage unit.
!     I1MACH( 6) = the number of characters per integer storage unit.
!
!   Integers:
!     assume integers are represented in the S-digit, base-A form
!
!                sign ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
!
!                where 0 .LE. X(I) .LT. A for I=0,...,S-1.
!     I1MACH( 7) = A, the base.
!     I1MACH( 8) = S, the number of base-A digits.
!     I1MACH( 9) = A**S - 1, the largest magnitude.
!
!   Floating-Point Numbers:
!     Assume floating-point numbers are represented in the T-digit,
!     base-B form
!                sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!                where 0 .LE. X(I) .LT. B for I=1,...,T,
!                0 .LT. X(1), and EMIN .LE. E .LE. EMAX.
!     I1MACH(10) = B, the base.
!
!   Single-Precision:
!     I1MACH(11) = T, the number of base-B digits.
!     I1MACH(12) = EMIN, the smallest exponent E.
!     I1MACH(13) = EMAX, the largest exponent E.
!
!   Double-Precision:
!     I1MACH(14) = T, the number of base-B digits.
!     I1MACH(15) = EMIN, the smallest exponent E.
!     I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of I1MACH(1) - I1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   960411  Modified for Fortran 90 (BE after suggestions by EHG).   
!   980727  Modified value of I1MACH(6) (BE after suggestion by EHG).   
!***END PROLOGUE  I1MACH
!
      X  = 1.0      
      XX = 1.0D0

      SELECT CASE (I)
        CASE (1)
          I1MACH = 5 ! Input unit
        CASE (2)
          I1MACH = 6 ! Output unit
        CASE (3)
          I1MACH = 0 ! Punch unit is no longer used
        CASE (4)
c         I1MACH = 0 ! Error message unit
          I1MACH = 6 ! Error message unit	! cfd
        CASE (5)
          I1MACH = BIT_SIZE(I)
        CASE (6)
          I1MACH = 4            ! Characters per integer is hopefully no
                                ! longer used. 
                                ! If it is used it has to be set manually.
                                ! The value 4 is correct on IEEE-machines.
        CASE (7)
          I1MACH = RADIX(1)
        CASE (8)
          I1MACH = BIT_SIZE(I) - 1
        CASE (9)
          I1MACH = HUGE(1)
        CASE (10)
          I1MACH = RADIX(X)
        CASE (11)
          I1MACH = DIGITS(X)
        CASE (12)
          I1MACH = MINEXPONENT(X)
        CASE (13)
          I1MACH = MAXEXPONENT(X)
        CASE (14)
          I1MACH = DIGITS(XX)
        CASE (15)
          I1MACH = MINEXPONENT(XX)
        CASE (16)
          I1MACH = MAXEXPONENT(XX) 
        CASE DEFAULT
          WRITE (*, FMT = 9000)
 9000     FORMAT ('1ERROR    1 IN I1MACH - I OUT OF BOUNDS')
          STOP
        END SELECT
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C***BEGIN PROLOGUE  XERMSG
C***PURPOSE  Process error messages for SLATEC and other libraries.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERMSG-A)
C***KEYWORDS  ERROR MESSAGE, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C   XERMSG processes a diagnostic message in a manner determined by the
C   value of LEVEL and the current value of the library error control
C   flag, KONTRL.  See subroutine XSETF for details.
C
C    LIBRAR   A character constant (or character variable) with the name
C             of the library.  This will be 'SLATEC' for the SLATEC
C             Common Math Library.  The error handling package is
C             general enough to be used by many libraries
C             simultaneously, so it is desirable for the routine that
C             detects and reports an error to identify the library name
C             as well as the routine name.
C
C    SUBROU   A character constant (or character variable) with the name
C             of the routine that detected the error.  Usually it is the
C             name of the routine that is calling XERMSG.  There are
C             some instances where a user callable library routine calls
C             lower level subsidiary routines where the error is
C             detected.  In such cases it may be more informative to
C             supply the name of the routine the user called rather than
C             the name of the subsidiary routine that detected the
C             error.
C
C    MESSG    A character constant (or character variable) with the text
C             of the error or warning message.  In the example below,
C             the message is a character constant that contains a
C             generic message.
C
C                   CALL XERMSG ('SLATEC', 'MMPY',
C                  *'THE ORDER OF THE MATRIX EXCEEDS THE ROW DIMENSION',
C                  *3, 1)
C
C             It is possible (and is sometimes desirable) to generate a
C             specific message--e.g., one that contains actual numeric
C             values.  Specific numeric values can be converted into
C             character strings using formatted WRITE statements into
C             character variables.  This is called standard Fortran
C             internal file I/O and is exemplified in the first three
C             lines of the following example.  You can also catenate
C             substrings of characters to construct the error message.
C             Here is an example showing the use of both writing to
C             an internal file and catenating character strings.
C
C                   CHARACTER*5 CHARN, CHARL
C                   WRITE (CHARN,10) N
C                   WRITE (CHARL,10) LDA
C                10 FORMAT(I5)
C                   CALL XERMSG ('SLATEC', 'MMPY', 'THE ORDER'//CHARN//
C                  *   ' OF THE MATRIX EXCEEDS ITS ROW DIMENSION OF'//
C                  *   CHARL, 3, 1)
C
C             There are two subtleties worth mentioning.  One is that
C             the // for character catenation is used to construct the
C             error message so that no single character constant is
C             continued to the next line.  This avoids confusion as to
C             whether there are trailing blanks at the end of the line.
C             The second is that by catenating the parts of the message
C             as an actual argument rather than encoding the entire
C             message into one large character variable, we avoid
C             having to know how long the message will be in order to
C             declare an adequate length for that large character
C             variable.  XERMSG calls XERPRN to print the message using
C             multiple lines if necessary.  If the message is very long,
C             XERPRN will break it into pieces of 72 characters (as
C             requested by XERMSG) for printing on multiple lines.
C             Also, XERMSG asks XERPRN to prefix each line with ' *  '
C             so that the total line length could be 76 characters.
C             Note also that XERPRN scans the error message backwards
C             to ignore trailing blanks.  Another feature is that
C             the substring '$$' is treated as a new line sentinel
C             by XERPRN.  If you want to construct a multiline
C             message without having to count out multiples of 72
C             characters, just use '$$' as a separator.  '$$'
C             obviously must occur within 72 characters of the
C             start of each line to have its intended effect since
C             XERPRN is asked to wrap around at 72 characters in
C             addition to looking for '$$'.
C
C    NERR     An integer value that is chosen by the library routine's
C             author.  It must be in the range -99 to 999 (three
C             printable digits).  Each distinct error should have its
C             own error number.  These error numbers should be described
C             in the machine readable documentation for the routine.
C             The error numbers need be unique only within each routine,
C             so it is reasonable for each routine to start enumerating
C             errors from 1 and proceeding to the next integer.
C
C    LEVEL    An integer value in the range 0 to 2 that indicates the
C             level (severity) of the error.  Their meanings are
C
C            -1  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.  An attempt is made to only print this
C                message once.
C
C             0  A warning message.  This is used if it is not clear
C                that there really is an error, but the user's attention
C                may be needed.
C
C             1  A recoverable error.  This is used even if the error is
C                so serious that the routine cannot return any useful
C                answer.  If the user has told the error package to
C                return after recoverable errors, then XERMSG will
C                return to the Library routine which can then return to
C                the user's routine.  The user may also permit the error
C                package to terminate the program upon encountering a
C                recoverable error.
C
C             2  A fatal error.  XERMSG will not return to its caller
C                after it receives a fatal error.  This level should
C                hardly ever be used; it is much better to allow the
C                user a chance to recover.  An example of one of the few
C                cases in which it is permissible to declare a level 2
C                error is a reverse communication Library routine that
C                is likely to be called repeatedly until it integrates
C                across some interval.  If there is a serious error in
C                the input such that another step cannot be taken and
C                the Library routine is called again without the input
C                error having been corrected by the caller, the Library
C                routine will probably be called forever with improper
C                input.  In this case, it is reasonable to declare the
C                error to be fatal.
C
C    Each of the arguments to XERMSG is input; none will be modified by
C    XERMSG.  A routine may make multiple calls to XERMSG with warning
C    level messages; however, after a call to XERMSG with a recoverable
C    error, the routine should return to the user.  Do not try to call
C    XERMSG with a second recoverable error after the first recoverable
C    error because the error package saves the error number.  The user
C    can retrieve this error number by calling another entry point in
C    the error handling package and then clear the error number when
C    recovering from the error.  Calling XERMSG in succession causes the
C    old error number to be overwritten by the latest error number.
C    This is considered harmless for error numbers associated with
C    warning messages but must not be done for error numbers of serious
C    errors.  After a call to XERMSG with a recoverable error, the user
C    must be given a chance to call NUMXER or XERCLR to retrieve or
C    clear the error number.
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  FDUMP, J4SAVE, XERCNT, XERHLT, XERPRN, XERSVE
C***REVISION HISTORY  (YYMMDD)
C   880101  DATE WRITTEN
C   880621  REVISED AS DIRECTED AT SLATEC CML MEETING OF FEBRUARY 1988.
C           THERE ARE TWO BASIC CHANGES.
C           1.  A NEW ROUTINE, XERPRN, IS USED INSTEAD OF XERPRT TO
C               PRINT MESSAGES.  THIS ROUTINE WILL BREAK LONG MESSAGES
C               INTO PIECES FOR PRINTING ON MULTIPLE LINES.  '$$' IS
C               ACCEPTED AS A NEW LINE SENTINEL.  A PREFIX CAN BE
C               ADDED TO EACH LINE TO BE PRINTED.  XERMSG USES EITHER
C               ' ***' OR ' *  ' AND LONG MESSAGES ARE BROKEN EVERY
C               72 CHARACTERS (AT MOST) SO THAT THE MAXIMUM LINE
C               LENGTH OUTPUT CAN NOW BE AS GREAT AS 76.
C           2.  THE TEXT OF ALL MESSAGES IS NOW IN UPPER CASE SINCE THE
C               FORTRAN STANDARD DOCUMENT DOES NOT ADMIT THE EXISTENCE
C               OF LOWER CASE.
C   880708  REVISED AFTER THE SLATEC CML MEETING OF JUNE 29 AND 30.
C           THE PRINCIPAL CHANGES ARE
C           1.  CLARIFY COMMENTS IN THE PROLOGUES
C           2.  RENAME XRPRNT TO XERPRN
C           3.  REWORK HANDLING OF '$$' IN XERPRN TO HANDLE BLANK LINES
C               SIMILAR TO THE WAY FORMAT STATEMENTS HANDLE THE /
C               CHARACTER FOR NEW RECORDS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           CLEAN UP THE CODING.
C   890721  REVISED TO USE NEW FEATURE IN XERPRN TO COUNT CHARACTERS IN
C           PREFIX.
C   891013  REVISED TO CORRECT COMMENTS.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Changed test on NERR to be -9999999 < NERR < 99999999, but
C           NERR .ne. 0, and on LEVEL to be -2 < LEVEL < 3.  Added
C           LEVEL=-1 logic, changed calls to XERSAV to XERSVE, and
C           XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERMSG
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8 XLIBR, XSUBR
      CHARACTER*72  TEMP
      CHARACTER*20  LFIRST
C***FIRST EXECUTABLE STATEMENT  XERMSG
      LKNTRL = J4SAVE (2, 0, .FALSE.)
      MAXMES = J4SAVE (4, 0, .FALSE.)
C
C       LKNTRL IS A LOCAL COPY OF THE CONTROL FLAG KONTRL.
C       MAXMES IS THE MAXIMUM NUMBER OF TIMES ANY PARTICULAR MESSAGE
C          SHOULD BE PRINTED.
C
C       WE PRINT A FATAL ERROR MESSAGE AND TERMINATE FOR AN ERROR IN
C          CALLING XERMSG.  THE ERROR NUMBER SHOULD BE POSITIVE,
C          AND THE LEVEL SHOULD BE BETWEEN 0 AND 2.
C
      IF (NERR.LT.-9999999 .OR. NERR.GT.99999999 .OR. NERR.EQ.0 .OR.
     *   LEVEL.LT.-1 .OR. LEVEL.GT.2) THEN
         CALL XERPRN (' ***', -1, 'FATAL ERROR IN...$$ ' //
     *      'XERMSG -- INVALID ERROR NUMBER OR LEVEL$$ '//
     *      'JOB ABORT DUE TO FATAL ERROR.', 72)
         CALL XERSVE (' ', ' ', ' ', 0, 0, 0, KDUMMY)
         CALL XERHLT (' ***XERMSG -- INVALID INPUT')
         RETURN
      ENDIF
C
C       RECORD THE MESSAGE.
C
      I = J4SAVE (1, NERR, .TRUE.)
      CALL XERSVE (LIBRAR, SUBROU, MESSG, 1, NERR, LEVEL, KOUNT)
C
C       HANDLE PRINT-ONCE WARNING MESSAGES.
C
      IF (LEVEL.EQ.-1 .AND. KOUNT.GT.1) RETURN
C
C       ALLOW TEMPORARY USER OVERRIDE OF THE CONTROL FLAG.
C
      XLIBR  = LIBRAR
      XSUBR  = SUBROU
      LFIRST = MESSG
      LERR   = NERR
      LLEVEL = LEVEL
      CALL XERCNT (XLIBR, XSUBR, LFIRST, LERR, LLEVEL, LKNTRL)
C
      LKNTRL = MAX(-2, MIN(2,LKNTRL))
      MKNTRL = ABS(LKNTRL)
C
C       SKIP PRINTING IF THE CONTROL FLAG VALUE AS RESET IN XERCNT IS
C       ZERO AND THE ERROR IS NOT FATAL.
C
      IF (LEVEL.LT.2 .AND. LKNTRL.EQ.0) GO TO 30
      IF (LEVEL.EQ.0 .AND. KOUNT.GT.MAXMES) GO TO 30
      IF (LEVEL.EQ.1 .AND. KOUNT.GT.MAXMES .AND. MKNTRL.EQ.1) GO TO 30
      IF (LEVEL.EQ.2 .AND. KOUNT.GT.MAX(1,MAXMES)) GO TO 30
C
C       ANNOUNCE THE NAMES OF THE LIBRARY AND SUBROUTINE BY BUILDING A
C       MESSAGE IN CHARACTER VARIABLE TEMP (NOT EXCEEDING 66 CHARACTERS)
C       AND SENDING IT OUT VIA XERPRN.  PRINT ONLY IF CONTROL FLAG
C       IS NOT ZERO.
C
      IF (LKNTRL .NE. 0) THEN
         TEMP(1:21) = 'MESSAGE FROM ROUTINE '
         I = MIN(LEN(SUBROU), 16)
         TEMP(22:21+I) = SUBROU(1:I)
         TEMP(22+I:33+I) = ' IN LIBRARY '
         LTEMP = 33 + I
         I = MIN(LEN(LIBRAR), 16)
         TEMP(LTEMP+1:LTEMP+I) = LIBRAR (1:I)
         TEMP(LTEMP+I+1:LTEMP+I+1) = '.'
         LTEMP = LTEMP + I + 1
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       IF LKNTRL IS POSITIVE, PRINT AN INTRODUCTORY LINE BEFORE
C       PRINTING THE MESSAGE.  THE INTRODUCTORY LINE TELLS THE CHOICE
C       FROM EACH OF THE FOLLOWING THREE OPTIONS.
C       1.  LEVEL OF THE MESSAGE
C              'INFORMATIVE MESSAGE'
C              'POTENTIALLY RECOVERABLE ERROR'
C              'FATAL ERROR'
C       2.  WHETHER CONTROL FLAG WILL ALLOW PROGRAM TO CONTINUE
C              'PROG CONTINUES'
C              'PROG ABORTED'
C       3.  WHETHER OR NOT A TRACEBACK WAS REQUESTED.  (THE TRACEBACK
C           MAY NOT BE IMPLEMENTED AT SOME SITES, SO THIS ONLY TELLS
C           WHAT WAS REQUESTED, NOT WHAT WAS DELIVERED.)
C              'TRACEBACK REQUESTED'
C              'TRACEBACK NOT REQUESTED'
C       NOTICE THAT THE LINE INCLUDING FOUR PREFIX CHARACTERS WILL NOT
C       EXCEED 74 CHARACTERS.
C       WE SKIP THE NEXT BLOCK IF THE INTRODUCTORY LINE IS NOT NEEDED.
C
      IF (LKNTRL .GT. 0) THEN
C
C       THE FIRST PART OF THE MESSAGE TELLS ABOUT THE LEVEL.
C
         IF (LEVEL .LE. 0) THEN
            TEMP(1:20) = 'INFORMATIVE MESSAGE,'
            LTEMP = 20
         ELSEIF (LEVEL .EQ. 1) THEN
            TEMP(1:30) = 'POTENTIALLY RECOVERABLE ERROR,'
            LTEMP = 30
         ELSE
            TEMP(1:12) = 'FATAL ERROR,'
            LTEMP = 12
         ENDIF
C
C       THEN WHETHER THE PROGRAM WILL CONTINUE.
C
         IF ((MKNTRL.EQ.2 .AND. LEVEL.GE.1) .OR.
     *       (MKNTRL.EQ.1 .AND. LEVEL.EQ.2)) THEN
            TEMP(LTEMP+1:LTEMP+14) = ' PROG ABORTED,'
            LTEMP = LTEMP + 14
         ELSE
            TEMP(LTEMP+1:LTEMP+16) = ' PROG CONTINUES,'
            LTEMP = LTEMP + 16
         ENDIF
C
C       FINALLY TELL WHETHER THERE SHOULD BE A TRACEBACK.
C
         IF (LKNTRL .GT. 0) THEN
            TEMP(LTEMP+1:LTEMP+20) = ' TRACEBACK REQUESTED'
            LTEMP = LTEMP + 20
         ELSE
            TEMP(LTEMP+1:LTEMP+24) = ' TRACEBACK NOT REQUESTED'
            LTEMP = LTEMP + 24
         ENDIF
         CALL XERPRN (' ***', -1, TEMP(1:LTEMP), 72)
      ENDIF
C
C       NOW SEND OUT THE MESSAGE.
C
      CALL XERPRN (' *  ', -1, MESSG, 72)
C
C       IF LKNTRL IS POSITIVE, WRITE THE ERROR NUMBER AND REQUEST A
C          TRACEBACK.
C
      IF (LKNTRL .GT. 0) THEN
         WRITE (TEMP, '(''ERROR NUMBER = '', I8)') NERR
         DO 10 I=16,22
            IF (TEMP(I:I) .NE. ' ') GO TO 20
   10    CONTINUE
C
   20    CALL XERPRN (' *  ', -1, TEMP(1:15) // TEMP(I:23), 72)
         CALL FDUMP
      ENDIF
C
C       IF LKNTRL IS NOT ZERO, PRINT A BLANK LINE AND AN END OF MESSAGE.
C
      IF (LKNTRL .NE. 0) THEN
         CALL XERPRN (' *  ', -1, ' ', 72)
         CALL XERPRN (' ***', -1, 'END OF MESSAGE', 72)
         CALL XERPRN ('    ',  0, ' ', 72)
      ENDIF
C
C       IF THE ERROR IS NOT FATAL OR THE ERROR IS RECOVERABLE AND THE
C       CONTROL FLAG IS SET FOR RECOVERY, THEN RETURN.
C
   30 IF (LEVEL.LE.0 .OR. (LEVEL.EQ.1 .AND. MKNTRL.LE.1)) RETURN
C
C       THE PROGRAM WILL BE STOPPED DUE TO AN UNRECOVERED ERROR OR A
C       FATAL ERROR.  PRINT THE REASON FOR THE ABORT AND THE ERROR
C       SUMMARY IF THE CONTROL FLAG AND THE MAXIMUM ERROR COUNT PERMIT.
C
      IF (LKNTRL.GT.0 .AND. KOUNT.LT.MAX(1,MAXMES)) THEN
         IF (LEVEL .EQ. 1) THEN
            CALL XERPRN
     *         (' ***', -1, 'JOB ABORT DUE TO UNRECOVERED ERROR.', 72)
         ELSE
            CALL XERPRN(' ***', -1, 'JOB ABORT DUE TO FATAL ERROR.', 72)
         ENDIF
         CALL XERSVE (' ', ' ', ' ', -1, 0, 0, KDUMMY)
         CALL XERHLT (' ')
      ELSE
         CALL XERHLT (MESSG)
      ENDIF
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK XERPRN
      SUBROUTINE XERPRN (PREFIX, NPREF, MESSG, NWRAP)
C***BEGIN PROLOGUE  XERPRN
C***SUBSIDIARY
C***PURPOSE  Print error messages processed by XERMSG.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERPRN-A)
C***KEYWORDS  ERROR MESSAGES, PRINTING, XERROR
C***AUTHOR  Fong, Kirby, (NMFECC at LLNL)
C***DESCRIPTION
C
C This routine sends one or more lines to each of the (up to five)
C logical units to which error messages are to be sent.  This routine
C is called several times by XERMSG, sometimes with a single line to
C print and sometimes with a (potentially very long) message that may
C wrap around into multiple lines.
C
C PREFIX  Input argument of type CHARACTER.  This argument contains
C         characters to be put at the beginning of each line before
C         the body of the message.  No more than 16 characters of
C         PREFIX will be used.
C
C NPREF   Input argument of type INTEGER.  This argument is the number
C         of characters to use from PREFIX.  If it is negative, the
C         intrinsic function LEN is used to determine its length.  If
C         it is zero, PREFIX is not used.  If it exceeds 16 or if
C         LEN(PREFIX) exceeds 16, only the first 16 characters will be
C         used.  If NPREF is positive and the length of PREFIX is less
C         than NPREF, a copy of PREFIX extended with blanks to length
C         NPREF will be used.
C
C MESSG   Input argument of type CHARACTER.  This is the text of a
C         message to be printed.  If it is a long message, it will be
C         broken into pieces for printing on multiple lines.  Each line
C         will start with the appropriate prefix and be followed by a
C         piece of the message.  NWRAP is the number of characters per
C         piece; that is, after each NWRAP characters, we break and
C         start a new line.  In addition the characters '$$' embedded
C         in MESSG are a sentinel for a new line.  The counting of
C         characters up to NWRAP starts over for each new line.  The
C         value of NWRAP typically used by XERMSG is 72 since many
C         older error messages in the SLATEC Library are laid out to
C         rely on wrap-around every 72 characters.
C
C NWRAP   Input argument of type INTEGER.  This gives the maximum size
C         piece into which to break MESSG for printing on multiple
C         lines.  An embedded '$$' ends a line, and the count restarts
C         at the following character.  If a line break does not occur
C         on a blank (it would split a word) that word is moved to the
C         next line.  Values of NWRAP less than 16 will be treated as
C         16.  Values of NWRAP greater than 132 will be treated as 132.
C         The actual line length will be NPREF + NWRAP after NPREF has
C         been adjusted to fall between 0 and 16 and NWRAP has been
C         adjusted to fall between 16 and 132.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   880621  DATE WRITTEN
C   880708  REVISED AFTER THE SLATEC CML SUBCOMMITTEE MEETING OF
C           JUNE 29 AND 30 TO CHANGE THE NAME TO XERPRN AND TO REWORK
C           THE HANDLING OF THE NEW LINE SENTINEL TO BEHAVE LIKE THE
C           SLASH CHARACTER IN FORMAT STATEMENTS.
C   890706  REVISED WITH THE HELP OF FRED FRITSCH AND REG CLEMENS TO
C           STREAMLINE THE CODING AND FIX A BUG THAT CAUSED EXTRA BLANK
C           LINES TO BE PRINTED.
C   890721  REVISED TO ADD A NEW FEATURE.  A NEGATIVE VALUE OF NPREF
C           CAUSES LEN(PREFIX) TO BE USED AS THE LENGTH.
C   891013  REVISED TO CORRECT ERROR IN CALCULATING PREFIX LENGTH.
C   891214  Prologue converted to Version 4.0 format.  (WRB)
C   900510  Added code to break messages between words.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERPRN
      CHARACTER*(*) PREFIX, MESSG
      INTEGER NPREF, NWRAP
      CHARACTER*148 CBUFF
      INTEGER IU(5), NUNIT
      CHARACTER*2 NEWLIN
      PARAMETER (NEWLIN = '$$')
C***FIRST EXECUTABLE STATEMENT  XERPRN
      CALL XGETUA(IU,NUNIT)
C
C       A ZERO VALUE FOR A LOGICAL UNIT NUMBER MEANS TO USE THE STANDARD
C       ERROR MESSAGE UNIT INSTEAD.  I1MACH(4) RETRIEVES THE STANDARD
C       ERROR MESSAGE UNIT.
C
      N = I1MACH(4)
      DO 10 I=1,NUNIT
         IF (IU(I) .EQ. 0) IU(I) = N
   10 CONTINUE
C
C       LPREF IS THE LENGTH OF THE PREFIX.  THE PREFIX IS PLACED AT THE
C       BEGINNING OF CBUFF, THE CHARACTER BUFFER, AND KEPT THERE DURING
C       THE REST OF THIS ROUTINE.
C
      IF ( NPREF .LT. 0 ) THEN
         LPREF = LEN(PREFIX)
      ELSE
         LPREF = NPREF
      ENDIF
      LPREF = MIN(16, LPREF)
      IF (LPREF .NE. 0) CBUFF(1:LPREF) = PREFIX
C
C       LWRAP IS THE MAXIMUM NUMBER OF CHARACTERS WE WANT TO TAKE AT ONE
C       TIME FROM MESSG TO PRINT ON ONE LINE.
C
      LWRAP = MAX(16, MIN(132, NWRAP))
C
C       SET LENMSG TO THE LENGTH OF MESSG, IGNORE ANY TRAILING BLANKS.
C
      LENMSG = LEN(MESSG)
      N = LENMSG
      DO 20 I=1,N
         IF (MESSG(LENMSG:LENMSG) .NE. ' ') GO TO 30
         LENMSG = LENMSG - 1
   20 CONTINUE
   30 CONTINUE
C
C       IF THE MESSAGE IS ALL BLANKS, THEN PRINT ONE BLANK LINE.
C
      IF (LENMSG .EQ. 0) THEN
         CBUFF(LPREF+1:LPREF+1) = ' '
         DO 40 I=1,NUNIT
            WRITE(IU(I), '(A)') CBUFF(1:LPREF+1)
   40    CONTINUE
         RETURN
      ENDIF
C
C       SET NEXTC TO THE POSITION IN MESSG WHERE THE NEXT SUBSTRING
C       STARTS.  FROM THIS POSITION WE SCAN FOR THE NEW LINE SENTINEL.
C       WHEN NEXTC EXCEEDS LENMSG, THERE IS NO MORE TO PRINT.
C       WE LOOP BACK TO LABEL 50 UNTIL ALL PIECES HAVE BEEN PRINTED.
C
C       WE LOOK FOR THE NEXT OCCURRENCE OF THE NEW LINE SENTINEL.  THE
C       INDEX INTRINSIC FUNCTION RETURNS ZERO IF THERE IS NO OCCURRENCE
C       OR IF THE LENGTH OF THE FIRST ARGUMENT IS LESS THAN THE LENGTH
C       OF THE SECOND ARGUMENT.
C
C       THERE ARE SEVERAL CASES WHICH SHOULD BE CHECKED FOR IN THE
C       FOLLOWING ORDER.  WE ARE ATTEMPTING TO SET LPIECE TO THE NUMBER
C       OF CHARACTERS THAT SHOULD BE TAKEN FROM MESSG STARTING AT
C       POSITION NEXTC.
C
C       LPIECE .EQ. 0   THE NEW LINE SENTINEL DOES NOT OCCUR IN THE
C                       REMAINDER OF THE CHARACTER STRING.  LPIECE
C                       SHOULD BE SET TO LWRAP OR LENMSG+1-NEXTC,
C                       WHICHEVER IS LESS.
C
C       LPIECE .EQ. 1   THE NEW LINE SENTINEL STARTS AT MESSG(NEXTC:
C                       NEXTC).  LPIECE IS EFFECTIVELY ZERO, AND WE
C                       PRINT NOTHING TO AVOID PRODUCING UNNECESSARY
C                       BLANK LINES.  THIS TAKES CARE OF THE SITUATION
C                       WHERE THE LIBRARY ROUTINE HAS A MESSAGE OF
C                       EXACTLY 72 CHARACTERS FOLLOWED BY A NEW LINE
C                       SENTINEL FOLLOWED BY MORE CHARACTERS.  NEXTC
C                       SHOULD BE INCREMENTED BY 2.
C
C       LPIECE .GT. LWRAP+1  REDUCE LPIECE TO LWRAP.
C
C       ELSE            THIS LAST CASE MEANS 2 .LE. LPIECE .LE. LWRAP+1
C                       RESET LPIECE = LPIECE-1.  NOTE THAT THIS
C                       PROPERLY HANDLES THE END CASE WHERE LPIECE .EQ.
C                       LWRAP+1.  THAT IS, THE SENTINEL FALLS EXACTLY
C                       AT THE END OF A LINE.
C
      NEXTC = 1
   50 LPIECE = INDEX(MESSG(NEXTC:LENMSG), NEWLIN)
      IF (LPIECE .EQ. 0) THEN
C
C       THERE WAS NO NEW LINE SENTINEL FOUND.
C
         IDELTA = 0
         LPIECE = MIN(LWRAP, LENMSG+1-NEXTC)
         IF (LPIECE .LT. LENMSG+1-NEXTC) THEN
            DO 52 I=LPIECE+1,2,-1
               IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
                  LPIECE = I-1
                  IDELTA = 1
                  GOTO 54
               ENDIF
   52       CONTINUE
         ENDIF
   54    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSEIF (LPIECE .EQ. 1) THEN
C
C       WE HAVE A NEW LINE SENTINEL AT MESSG(NEXTC:NEXTC+1).
C       DON'T PRINT A BLANK LINE.
C
         NEXTC = NEXTC + 2
         GO TO 50
      ELSEIF (LPIECE .GT. LWRAP+1) THEN
C
C       LPIECE SHOULD BE SET DOWN TO LWRAP.
C
         IDELTA = 0
         LPIECE = LWRAP
         DO 56 I=LPIECE+1,2,-1
            IF (MESSG(NEXTC+I-1:NEXTC+I-1) .EQ. ' ') THEN
               LPIECE = I-1
               IDELTA = 1
               GOTO 58
            ENDIF
   56    CONTINUE
   58    CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC = NEXTC + LPIECE + IDELTA
      ELSE
C
C       IF WE ARRIVE HERE, IT MEANS 2 .LE. LPIECE .LE. LWRAP+1.
C       WE SHOULD DECREMENT LPIECE BY ONE.
C
         LPIECE = LPIECE - 1
         CBUFF(LPREF+1:LPREF+LPIECE) = MESSG(NEXTC:NEXTC+LPIECE-1)
         NEXTC  = NEXTC + LPIECE + 2
      ENDIF
C
C       PRINT
C
      DO 60 I=1,NUNIT
         WRITE(IU(I), '(A)') CBUFF(1:LPREF+LPIECE)
   60 CONTINUE
C
      IF (NEXTC .LE. LENMSG) GO TO 50
      RETURN
      END


C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK XERSVE
      SUBROUTINE XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL,
     +   ICOUNT)
C***BEGIN PROLOGUE  XERSVE
C***SUBSIDIARY
C***PURPOSE  Record that an error has occurred.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3
C***TYPE      ALL (XERSVE-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C *Usage:
C
C        INTEGER  KFLAG, NERR, LEVEL, ICOUNT
C        CHARACTER * (len) LIBRAR, SUBROU, MESSG
C
C        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT)
C
C *Arguments:
C
C        LIBRAR :IN    is the library that the message is from.
C        SUBROU :IN    is the subroutine that the message is from.
C        MESSG  :IN    is the message to be saved.
C        KFLAG  :IN    indicates the action to be performed.
C                      when KFLAG > 0, the message in MESSG is saved.
C                      when KFLAG=0 the tables will be dumped and
C                      cleared.
C                      when KFLAG < 0, the tables will be dumped and
C                      not cleared.
C        NERR   :IN    is the error number.
C        LEVEL  :IN    is the error severity.
C        ICOUNT :OUT   the number of times this message has been seen,
C                      or zero if the table has overflowed and does not
C                      contain this message specifically.  When KFLAG=0,
C                      ICOUNT will not be altered.
C
C *Description:
C
C   Record that this error occurred and possibly dump and clear the
C   tables.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  I1MACH, XGETUA
C***REVISION HISTORY  (YYMMDD)
C   800319  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900413  Routine modified to remove reference to KFLAG.  (WRB)
C   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling
C           sequence, use IF-THEN-ELSE, make number of saved entries
C           easily changeable, changed routine name from XERSAV to
C           XERSVE.  (RWC)
C   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERSVE
      PARAMETER (LENTAB=10)
      INTEGER LUN(5)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      CHARACTER*8  LIBTAB(LENTAB), SUBTAB(LENTAB), LIB, SUB
      CHARACTER*20 MESTAB(LENTAB), MES
      DIMENSION NERTAB(LENTAB), LEVTAB(LENTAB), KOUNT(LENTAB)
      SAVE LIBTAB, SUBTAB, MESTAB, NERTAB, LEVTAB, KOUNT, KOUNTX, NMSG
      DATA KOUNTX/0/, NMSG/0/
C***FIRST EXECUTABLE STATEMENT  XERSVE
C
      IF (KFLAG.LE.0) THEN
C
C        Dump the table.
C
         IF (NMSG.EQ.0) RETURN
C
C        Print to each unit.
C
         CALL XGETUA (LUN, NUNIT)
         DO 20 KUNIT = 1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
C
C           Print the table header.
C
            WRITE (IUNIT,9000)
C
C           Print body of table.
C
            DO 10 I = 1,NMSG
               WRITE (IUNIT,9010) LIBTAB(I), SUBTAB(I), MESTAB(I),
     *            NERTAB(I),LEVTAB(I),KOUNT(I)
   10       CONTINUE
C
C           Print number of other errors.
C
            IF (KOUNTX.NE.0) WRITE (IUNIT,9020) KOUNTX
            WRITE (IUNIT,9030)
   20    CONTINUE
C
C        Clear the error tables.
C
         IF (KFLAG.EQ.0) THEN
            NMSG = 0
            KOUNTX = 0
         ENDIF
      ELSE
C
C        PROCESS A MESSAGE...
C        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG,
C        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL.
C
         LIB = LIBRAR
         SUB = SUBROU
         MES = MESSG
         DO 30 I = 1,NMSG
            IF (LIB.EQ.LIBTAB(I) .AND. SUB.EQ.SUBTAB(I) .AND.
     *         MES.EQ.MESTAB(I) .AND. NERR.EQ.NERTAB(I) .AND.
     *         LEVEL.EQ.LEVTAB(I)) THEN
                  KOUNT(I) = KOUNT(I) + 1
                  ICOUNT = KOUNT(I)
                  RETURN
            ENDIF
   30    CONTINUE
C
         IF (NMSG.LT.LENTAB) THEN
C
C           Empty slot found for new message.
C
            NMSG = NMSG + 1
            LIBTAB(I) = LIB
            SUBTAB(I) = SUB
            MESTAB(I) = MES
            NERTAB(I) = NERR
            LEVTAB(I) = LEVEL
            KOUNT (I) = 1
            ICOUNT    = 1
         ELSE
C
C           Table is full.
C
            KOUNTX = KOUNTX+1
            ICOUNT = 0
         ENDIF
      ENDIF
      RETURN
C
C     Formats.
C
 9000 FORMAT ('0          ERROR MESSAGE SUMMARY' /
     +   ' LIBRARY    SUBROUTINE MESSAGE START             NERR',
     +   '     LEVEL     COUNT')
 9010 FORMAT (1X,A,3X,A,3X,A,3I10)
 9020 FORMAT ('0OTHER ERRORS NOT INDIVIDUALLY TABULATED = ', I10)
 9030 FORMAT (1X)
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK XERHLT
      SUBROUTINE XERHLT (MESSG)
C***BEGIN PROLOGUE  XERHLT
C***SUBSIDIARY
C***PURPOSE  Abort program execution and print error message.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERHLT-A)
C***KEYWORDS  ABORT PROGRAM EXECUTION, ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        ***Note*** machine dependent routine
C        XERHLT aborts the execution of the program.
C        The error message causing the abort is given in the calling
C        sequence, in case one needs it for printing on a dayfile,
C        for example.
C
C     Description of Parameters
C        MESSG is as in XERMSG.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to delete length of character
C           and changed routine name from XERABT to XERHLT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERHLT
      CHARACTER*(*) MESSG
C***FIRST EXECUTABLE STATEMENT  XERHLT
      STOP
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

*DECK XERCNT
      SUBROUTINE XERCNT (LIBRAR, SUBROU, MESSG, NERR, LEVEL, KONTRL)
C***BEGIN PROLOGUE  XERCNT
C***SUBSIDIARY
C***PURPOSE  Allow user control over handling of errors.
C***LIBRARY   SLATEC (XERROR)
C***CATEGORY  R3C
C***TYPE      ALL (XERCNT-A)
C***KEYWORDS  ERROR, XERROR
C***AUTHOR  Jones, R. E., (SNLA)
C***DESCRIPTION
C
C     Abstract
C        Allows user control over handling of individual errors.
C        Just after each message is recorded, but before it is
C        processed any further (i.e., before it is printed or
C        a decision to abort is made), a call is made to XERCNT.
C        If the user has provided his own version of XERCNT, he
C        can then override the value of KONTROL used in processing
C        this message by redefining its value.
C        KONTRL may be set to any value from -2 to 2.
C        The meanings for KONTRL are the same as in XSETF, except
C        that the value of KONTRL changes only for this message.
C        If KONTRL is set to a value outside the range from -2 to 2,
C        it will be moved back into that range.
C
C     Description of Parameters
C
C      --Input--
C        LIBRAR - the library that the routine is in.
C        SUBROU - the subroutine that XERMSG is being called from
C        MESSG  - the first 20 characters of the error message.
C        NERR   - same as in the call to XERMSG.
C        LEVEL  - same as in the call to XERMSG.
C        KONTRL - the current value of the control flag as set
C                 by a call to XSETF.
C
C      --Output--
C        KONTRL - the new value of KONTRL.  If KONTRL is not
C                 defined, it will remain at its original value.
C                 This changed value of control affects only
C                 the current occurrence of the current message.
C
C***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC
C                 Error-handling Package, SAND82-0800, Sandia
C                 Laboratories, 1982.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   790801  DATE WRITTEN
C   861211  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900206  Routine changed from user-callable to subsidiary.  (WRB)
C   900510  Changed calling sequence to include LIBRARY and SUBROUTINE
C           names, changed routine name from XERCTL to XERCNT.  (RWC)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  XERCNT
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
C***FIRST EXECUTABLE STATEMENT  XERCNT
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
C
C     ABSTRACT
C        J4SAVE SAVES AND RECALLS SEVERAL GLOBAL VARIABLES NEEDED
C        BY THE LIBRARY ERROR HANDLING ROUTINES.
C
C     DESCRIPTION OF PARAMETERS
C      --INPUT--
C        IWHICH - INDEX OF ITEM DESIRED.
C                 = 1 REFERS TO CURRENT ERROR NUMBER.
C                 = 2 REFERS TO CURRENT ERROR CONTROL FLAG.
C                 = 3 REFERS TO CURRENT UNIT NUMBER TO WHICH ERROR
C                     MESSAGES ARE TO BE SENT.  (0 MEANS USE STANDARD.)
C                 = 4 REFERS TO THE MAXIMUM NUMBER OF TIMES ANY
C                     MESSAGE IS TO BE PRINTED (AS SET BY XERMAX).
C                 = 5 REFERS TO THE TOTAL NUMBER OF UNITS TO WHICH
C                     EACH ERROR MESSAGE IS TO BE WRITTEN.
C                 = 6 REFERS TO THE 2ND UNIT FOR ERROR MESSAGES
C                 = 7 REFERS TO THE 3RD UNIT FOR ERROR MESSAGES
C                 = 8 REFERS TO THE 4TH UNIT FOR ERROR MESSAGES
C                 = 9 REFERS TO THE 5TH UNIT FOR ERROR MESSAGES
C        IVALUE - THE VALUE TO BE SET FOR THE IWHICH-TH PARAMETER,
C                 IF ISET IS .TRUE. .
C        ISET   - IF ISET=.TRUE., THE IWHICH-TH PARAMETER WILL BE
C                 GIVEN THE VALUE, IVALUE.  IF ISET=.FALSE., THE
C                 IWHICH-TH PARAMETER WILL BE UNCHANGED, AND IVALUE
C                 IS A DUMMY PARAMETER.
C      --OUTPUT--
C        THE (OLD) VALUE OF THE IWHICH-TH PARAMETER WILL BE RETURNED
C        IN THE FUNCTION VALUE, J4SAVE.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     ADAPTED FROM BELL LABORATORIES PORT LIBRARY ERROR HANDLER
C     LATEST REVISION ---  23 MAY 1979
C
      LOGICAL ISET
      INTEGER IPARAM(9)
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,1,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END


C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      SUBROUTINE FDUMP
C***BEGIN PROLOGUE  FDUMP
C***DATE WRITTEN   790801   (YYMMDD)
C***REVISION DATE  820801   (YYMMDD)
C***CATEGORY NO.  Z
C***KEYWORDS  ERROR,XERROR PACKAGE
C***AUTHOR  JONES, R. E., (SNLA)
C***PURPOSE  SYMBOLIC DUMP (SHOULD BE LOCALLY WRITTEN).
C***DESCRIPTION
C        ***NOTE*** MACHINE DEPENDENT ROUTINE
C        FDUMP IS INTENDED TO BE REPLACED BY A LOCALLY WRITTEN
C        VERSION WHICH PRODUCES A SYMBOLIC DUMP.  FAILING THIS,
C        IT SHOULD BE REPLACED BY A VERSION WHICH PRINTS THE
C        SUBPROGRAM NESTING LIST.  NOTE THAT THIS DUMP MUST BE
C        PRINTED ON EACH OF UP TO FIVE FILES, AS INDICATED BY THE
C        XGETUA ROUTINE.  SEE XSETUA AND XGETUA FOR DETAILS.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C     LATEST REVISION ---  23 MAY 1979
C***ROUTINES CALLED  (NONE)
C***END PROLOGUE  FDUMP
C***FIRST EXECUTABLE STATEMENT  FDUMP
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      SUBROUTINE XGETUA(IUNIT,N)
C
C     ABSTRACT
C        XGETUA MAY BE CALLED TO DETERMINE THE UNIT NUMBER OR NUMBERS
C        TO WHICH ERROR MESSAGES ARE BEING SENT.
C        THESE UNIT NUMBERS MAY HAVE BEEN SET BY A CALL TO XSETUN,
C        OR A CALL TO XSETUA, OR MAY BE A DEFAULT VALUE.
C
C     DESCRIPTION OF PARAMETERS
C      --OUTPUT--
C        IUNIT - AN ARRAY OF ONE TO FIVE UNIT NUMBERS, DEPENDING
C                ON THE VALUE OF N.  A VALUE OF ZERO REFERS TO THE
C                DEFAULT UNIT, AS DEFINED BY THE I1MACH MACHINE
C                CONSTANT ROUTINE.  ONLY IUNIT(1),...,IUNIT(N) ARE
C                DEFINED BY XGETUA.  THE VALUES OF IUNIT(N+1),...,
C                IUNIT(5) ARE NOT DEFINED (FOR N.LT.5) OR ALTERED
C                IN ANY WAY BY XGETUA.
C        N     - THE NUMBER OF UNITS TO WHICH COPIES OF THE
C                ERROR MESSAGES ARE BEING SENT.  N WILL BE IN THE
C                RANGE FROM 1 TO 5.
C
C     WRITTEN BY RON JONES, WITH SLATEC COMMON MATH LIBRARY SUBCOMMITTEE
C
      DIMENSION IUNIT(5)
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNIT(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END







