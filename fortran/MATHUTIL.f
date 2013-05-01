C Copyright (C) 2013 Ahmad El Sayed
C 
C  This file is part of PMFpack.
C 
C  PMFpack is free software: you can redistribute it and/or modify
C  it under the terms of the GNU Lesser General Public License as published by
C  the Free Software Foundation, either version 3 of the License, or
C  (at your option) any later version.
C 
C  PMFpack is distributed in the hope that it will be useful,
C  but WITHOUT ANY WARRANTY; without even the implied warranty of
C  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C  GNU Lesser General Public License for more details.
C 
C  You should have received a copy of the GNU Lesser General Public License
C  along with PMFpack. If not, see <http://www.gnu.org/licenses/>.
C 
C 
      DOUBLE PRECISION FUNCTION ERFINV(Y)
C     COMPUTES THE INVERSE OF THE ERROR FUNCTION
      IMPLICIT NONE 
      DOUBLE PRECISION X,Y,TOL,ZERO,LOW,HIGH,TRIALX,TRIALY
      PARAMETER(TOL =1.0D-15,ZERO = +0.0D0)
      DOUBLE PRECISION ERF
      EXTERNAL ERF

      LOW = -10.0D0
      HIGH = -LOW

      IF(Y.GE.1.0D0)THEN
C         ERFINV = -DLOG(ZERO)   !TO GENERATE +INF
         ERFINV = 1.0/0.0
         GOTO 666
      ENDIF
      
      IF(Y.LE.-1.0D0)THEN
C         ERFINV = DLOG(ZERO)    !TO GENERATE -INF
         ERFINV = -1.0/0.0
         GOTO 666
      ENDIF

      DO WHILE((HIGH-LOW).GT.TOL)
         TRIALX = (LOW+HIGH)/2.0D0
         TRIALY = ERF(TRIALX)
         IF(TRIALY.LT.Y)THEN
            LOW = TRIALX
         ELSE
            HIGH = TRIALX
         ENDIF
      ENDDO
      
      ERFINV = (LOW+HIGH)/2.0D0
     
 666  RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

c$$$      DOUBLE PRECISION FUNCTION ERF(X)
c$$$C     COMPUTES THE ERROR FUNCTION
c$$$C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
c$$$C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
c$$$C             CAMBRIDGE UNIVERSITY PRESS, 1992
c$$$
c$$$C     USES GAMMP
c$$$C     RETURNS THE ERROR FUNCTION ERF(X).
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION X
c$$$      DOUBLE PRECISION GAMMP
c$$$      EXTERNAL GAMMP
c$$$
c$$$      IF(X.LT.0.0D0)THEN
c$$$         ERF=-GAMMP(0.5D0,X**2.0D0)
c$$$      ELSE
c$$$         ERF=GAMMP(0.5D0,X**2.0D0)
c$$$      ENDIF
c$$$      RETURN
c$$$      END

      DOUBLE PRECISION FUNCTION ERF(X)
C     COMPUTES THE ERROR FUNCTION
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      IMPLICIT NONE
      DOUBLE PRECISION ERFC,X
      DOUBLE PRECISION T,Z
      Z = DABS(X)
      T = 1.0D0/(1.0D0+0.5D0*Z)
      ERFC = T*DEXP(-Z*Z-1.26551223D0
     $     +T*(1.00002368D0+T*(0.37409196D0+
     $     T*(0.09678418D0+T*(-0.18628806D0
     $     +T*(0.27886807D0+T*(-1.13520398D0+
     $     T*(1.48851587D0+T*(-0.82215223D0
     $     +T*0.17087277D0)))))))))
      IF (X.LT.0.0D0) ERFC = 2.0D0-ERFC
      ERF = 1.0D0-ERFC
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION FUNCTION BETA(Z,W)
C     USES GAMMLN
C     RETURNS THE VALUE OF THE BETA FUNCTION B(Z, W).
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      IMPLICIT NONE
      DOUBLE PRECISION W,Z
      DOUBLE PRECISION  GAMMLN
      EXTERNAL GAMMLN
      BETA=DEXP(GAMMLN(Z)+GAMMLN(W)-GAMMLN(Z+W))
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION FUNCTION GAMMP(A,X)
C     USES GCF,GSER
C     RETURNS THE INCOMPLETE GAMMA FUNCTION P (A, X).
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      IMPLICIT NONE
      DOUBLE PRECISION A,X
      DOUBLE PRECISION GAMMCF,GAMSER,GLN

      IF(X.LT.0.0D0.OR.A.LE.0.0D0)PAUSE 'BAD ARGUMENTS IN GAMMP'
      IF(X.LT.A+1.0D0)THEN      !USE THE SERIES REPRESENTATION.
         CALL GSER(GAMSER,A,X,GLN)
         GAMMP=GAMSER
      ELSE                      !USE THE CONTINUED FRACTION REPRESENTATION
         CALL GCF(GAMMCF,A,X,GLN)
         GAMMP=1.0D0-GAMMCF     !AND TAKE ITS COMPLEMENT.
      ENDIF
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  
      
      SUBROUTINE GCF(GAMMCF,A,X,GLN)
C     USES GAMMLN
C     RETURNS THE INCOMPLETE GAMMA FUNCTION Q(A, X) EVALUATED BY ITS CONTINUED FRACTION REPRE-
C     SENTATION AS GAMMCF. ALSO RETURNS LN GAMMA(A) AS GLN.
C     PARAMETERS: ITMAX IS THE MAXIMUM ALLOWED NUMBER OF ITERATIONS; EPS IS THE RELATIVE ACCU-
C     RACY; FPMIN IS A NUMBER NEAR THE SMALLEST REPRESENTABLE FLOATING-POINT NUMBER.
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      IMPLICIT NONE
      INTEGER ITMAX
      DOUBLE PRECISION A,GAMMCF,GLN,X,EPS,FPMIN
      PARAMETER (ITMAX=10000,EPS=1.0D-15,FPMIN=1.0D-30)
      INTEGER I
      DOUBLE PRECISION AN,B,C,D,DEL,H,GAMMLN
      EXTERNAL GAMMLN

      GLN=GAMMLN(A)
      B=X+1.0D0-A                  ! SET UP FOR EVALUATING CONTINUED FRACTION BY MODIﬁED
      C=1.0D0/FPMIN                ! LENT'S METHOD WITH B0 = 0.
      D=1.0D0/B
      H=D
      DO  I=1,ITMAX             !ITERATE TO CONVERGENCE.
         AN=-I*(I-A)
         B=B+2.0D0
         D=AN*D+B
         IF(DABS(D).LT.FPMIN)D=FPMIN
         C=B+AN/C
         IF(DABS(C).LT.FPMIN)C=FPMIN
         D=1.0D0/D
         DEL=D*C
         
         H=H*DEL
         IF(DABS(DEL-1.0D0).LT.EPS)GOTO 1
      ENDDO 
      PAUSE 'A TOO LARGE, ITMAX TOO SMALL IN GCF'
 1    GAMMCF=DEXP(-X+A*DLOG(X)-GLN)*H      !PUT FACTORS IN FRONT.
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      SUBROUTINE GSER(GAMSER,A,X,GLN)
C     USES GAMMLN
C     RETURNS THE INCOMPLETE GAMMA FUNCTION P (A, X) EVALUATED BY ITS SERIES 
C     REPRESENTATION AS GAMSER. ALSO RETURNS LN GAMMA(A) AS GLN.
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      IMPLICIT NONE
      INTEGER ITMAX
      DOUBLE PRECISION A,GAMSER,GLN,X,EPS
      PARAMETER (ITMAX=10000,EPS=1.0D-15)
      INTEGER N
      DOUBLE PRECISION AP,DEL,SUM,GAMMLN
      EXTERNAL GAMMLN

      GLN=GAMMLN(A)
      IF(X.LE.0.0D0)THEN
         IF(X.LT.0.0D0)PAUSE 'X < 0 IN GSER'
         GAMSER=0.0D0
         RETURN
      ENDIF
      AP=A
      SUM=1.0D0/A
      DEL=SUM
      DO  N=1,ITMAX
         AP=AP+1.0D0
         DEL=DEL*X/AP
         SUM=SUM+DEL
         IF(DABS(DEL).LT.DABS(SUM)*EPS)GOTO 1
      ENDDO
      PAUSE 'A TOO LARGE, ITMAX TOO SMALL IN GSER'
 1    GAMSER=SUM*DEXP(-X+A*DLOG(X)-GLN)
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================       

      DOUBLE PRECISION FUNCTION GAMMLN(XX)
C     RETURNS THE VALUE LN[GAMMA(XX)] FOR XX > 0.
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      IMPLICIT NONE
      DOUBLE PRECISION XX
      INTEGER J
      DOUBLE PRECISION SER,STP,TMP,X,Y,COF(6)
      SAVE COF,STP
      DATA COF,STP/76.18009172947146D0,
     $     -86.50532032941677D0,
     $     24.01409824083091D0,
     $     -1.231739572450155D0,
     $     0.1208650973866179D-2,
     $     -0.5395239384953D-5,
     $     2.5066282746310005D0/

      X=XX
      Y=X
      TMP=X+5.5D0
      TMP=(X+0.5D0)*DLOG(TMP)-TMP
      SER=1.000000000190015D0
      DO J=1,6
         Y=Y+1.D0
         SER=SER+COF(J)/Y
      ENDDO
      GAMMLN=TMP+DLOG(STP*SER/X)
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION FUNCTION BETAI(A,B,X)
C     RETURNS THE INCOMPLETE BETA FUNCTION I_X(A, B) 
C     USES BETACF,GAMMLN
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      DOUBLE PRECISION  A,B,X
      DOUBLE PRECISION  BT,BETACF,GAMMLN

      IF(X.LT.0.0D0.OR.X.GT.1.0D0)PAUSE 'BAD ARGUMENT X IN BETAI'
      IF(X.EQ.0.0D0.OR.X.EQ.1.0D0)THEN
         BT = 0.0D0
      ELSE  ! FACTORS IN FRONT OF THE CONTINUED FRACTION.
         BT = DEXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     $        +A*DLOG(X)+B*DLOG(1.0D0-X))
      ENDIF
      IF(X.LT.(A+1.0D0)/(A+B+2.0D0)) THEN ! USE CONTINUED FRACTION DIRECTLY.
         BETAI = BT*BETACF(A,B,X)/A
         RETURN
      ELSE
         BETAI = 1.0D0-BT*BETACF(B,A,1.0D0-X)/B ! USE CONTINUED FRACTION AFTER MAKING 
                                                ! THE SYMMETRY TRANSFORMATION
         RETURN                             
      ENDIF
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================    

      DOUBLE PRECISION FUNCTION BETACF(A,B,X)
C     USED BY BETAI: EVALUATES CONTINUED FRACTION FOR INCOMPLETE BETA 
C     FUNCTION BY MODIFIED  LENTZ'S METHOD
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      INTEGER MAXIT
      DOUBLE PRECISION A,B,X,EPS,FPMIN
      PARAMETER (MAXIT=10000,EPS=3.0D-7,FPMIN=1.0D-30)
      INTEGER M,M2
      DOUBLE PRECISION AA,C,D,DEL,H,QAB,QAM,QAP 

      QAB=A+B                   ! THESE Q'S WILL BE USED IN FACTORS THAT 
                                ! OCCUR IN THE COEFFICIENTS
      QAP=A+1.0D0                                     
      QAM=A-1.0D0
      C=1.0D0                   ! FIRST STEP OF LENTZ'S METHOD
      D=1.0D0-QAB*X/QAP
      IF(DABS(D).LT.FPMIN)D=FPMIN
      D=1.0D0/D
      H=D
      DO M=1,MAXIT
         M2=2*M
         AA=M*(B-M)*X/((QAM+M2)*(A+M2))
         D=1.0D0+AA*D                         ! ONE STEP (THE EVEN ONE) OF THE RECURRENCE
         IF(DABS(D).LT.FPMIN)D=FPMIN
         C=1.0D0+AA/C
         IF(DABS(C).LT.FPMIN)C=FPMIN
         D=1.0D0/D
         H=H*D*C
         AA=-(A+M)*(QAB+M)*X/((A+M2)*(QAP+M2))
         D=1.0D0+AA*D                         ! NEXT STEP OF THE RECURRENCE (THE ODD ONE).
         IF(DABS(D).LT.FPMIN)D=FPMIN
         C=1.0D0+AA/C
         IF(DABS(C).LT.FPMIN)C=FPMIN
         D=1.0D0/D
         DEL=D*C
         H=H*DEL
         IF(DABS(DEL-1.0D0).LT.EPS) GOTO 1     !ARE WE DONE?
      ENDDO
      !PAUSE 'A OR B TOO BIG, OR MAXIT TOO SMALL IN BETACF'
 1    BETACF=H
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION TRAP(F,R,N)
C     COMPUTES THE INTEGRAL OF A DISCRETE FUNCTION USING TRAPEZOIDAL INTEGRATION
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION F(N),R(N)	
      TRAP = 0.0D0
      DO I = 1,N-1
         TRAP = TRAP + (R(I+1)-R(I)) * (F(I)+F(I+1))
      ENDDO
      TRAP = TRAP * 0.50D0

      RETURN
      END
    
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      LOGICAL FUNCTION ISNAN(A)
C     TEST FOR NAN OCCURRENCE
      DOUBLE PRECISION A
      IF(DABS(A).LE.HUGE(A))THEN 
         ISNAN = .FALSE.
      ELSE
         ISNAN = .TRUE.
      END IF
      RETURN
      END
     
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 
 
      LOGICAL FUNCTION ISINF(A)
C     TEST FOR INF OCCURRENCE 
      DOUBLE PRECISION A
      IF ((A+1.0D0).EQ.A) THEN
         ISINF = .TRUE.
      ELSE
         ISINF = .FALSE.
      END IF
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
         NTPSI = INITDS (PSICS, 42, 0.1*REAL(D1MACH(3)) )
         NTAPSI = INITDS (APSICS, 16, 0.1*REAL(D1MACH(3)) )
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
         NTERMS = INITDS (COTCS, 15, 0.1*REAL(D1MACH(3)) )
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
        ERR = ERR + ABS(REAL(OS(I)))
        IF (ERR.GT.ETA) GO TO 20
   10 CONTINUE
C
   20 IF (I .EQ. NOS) CALL XERMSG ('SLATEC', 'INITDS',
     +   'Chebyshev series too short for specified accuracy', 1, 1)
      INITDS = I
C
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
