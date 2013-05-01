C===============================================================================
C===============================================================================
C===============================================================================

      DOUBLE PRECISION FUNCTION ZRIDDR(FUNC,X1,X2,XACC)
      
C     USING RIDDERS' METHOD, RETURN THE ROOT OF A FUNCTION FUNC KNOWN 
C     TO LIE BETWEEN X1 AND X2. THE ROOT, RETURNED AS ZRIDDR, WILL BE 
C     REFINED TO AN APPROXIMATE ACCURACY XACC.
      
      INTEGER MAXIT
      DOUBLE PRECISION X1,X2,XACC,FUNC,UNUSED
      PARAMETER (MAXIT=5000,UNUSED=-1.11D30)
      EXTERNAL FUNC
      INTEGER J
      DOUBLE PRECISION FH,FL,FM,FNEW,S,XH,XL,XM,XNEW
      
      FL=FUNC(X1)
      FH=FUNC(X2)
      IF((FL.GT.0.0D0.AND.FH.LT.0.0D0)
     $     .OR.(FL.LT.0.D0.AND.FH.GT.0.0D0))THEN
         XL=X1
         XH=X2
         ZRIDDR=UNUSED                              
         DO J=1,MAXIT                              
            XM=0.5D0*(XL+XH)
            FM=FUNC(XM)                           
            S=DSQRT(FM**2-FL*FH)                       
            IF(S.EQ.0.0D0) RETURN
            XNEW=XM+(XM-XL)*(DSIGN(1.0D0,FL-FH)*FM/S)     
            IF (DABS(XNEW-ZRIDDR).LE.XACC) RETURN
            ZRIDDR=XNEW
            FNEW=FUNC(ZRIDDR)                     
            IF (FNEW.EQ.0.0D0) RETURN                 
            IF(DSIGN(FM,FNEW).NE.FM) THEN         
               XL=XM                             
               FL=FM
               XH=ZRIDDR
               FH=FNEW
            ELSE IF(DSIGN(FL,FNEW).NE.FL) THEN
               XH=ZRIDDR
               FH=FNEW
            ELSE IF(DSIGN(FH,FNEW).NE.FH) THEN
               XL=ZRIDDR
               FL=FNEW
            ELSE
               PAUSE 'NEVER GET HERE IN ZRIDDR'
            ENDIF
            IF(DABS(XH-XL).LE.XACC) RETURN
         ENDDO
         PAUSE 'ZRIDDR EXCEED MAXIMUM ITERATIONS'
      ELSE IF (FL.EQ.0.0D0) THEN
         ZRIDDR=X1
      ELSE IF (FH.EQ.0.0D0) THEN
         ZRIDDR=X2
      ELSE
         PAUSE 'ROOT MUST BE BRACKETED IN ZRIDDR'
      ENDIF
      RETURN
      END

C===============================================================================
C===============================================================================
C===============================================================================

      DOUBLE PRECISION FUNCTION ZBRENT(FUNC,X1,X2,TOL)

C     USING BRENT'S METHOD, FIND THE ROOT OF A FUNCTION FUNC KNOWN TO LIE 
C     BETWEEN X1 AND X2. THE ROOT, RETURNED AS ZBRENT, WILL BE REFINED UNTIL 
C     ITS ACCURACY IS TOL. PARAMETERS: MAXIMUM ALLOWED NUMBER OF ITERATIONS, 
C     AND MACHINE FLOATING-POINT PRECISION.
      
      INTEGER ITMAX
      DOUBLE PRECISION TOL,X1,X2,FUNC,EPS
      EXTERNAL FUNC
      INTEGER ITER
      DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
      PARAMETER (ITMAX=5000,EPS = 1.0D-15)
      
      A=X1
      B=X2
      FA=FUNC(A)
      FB=FUNC(B)
      IF((FA.GT.0.0D0.AND.FB.GT.0.0D0).OR.(FA.LT.0.0D0.AND.FB.LT.0.0D0))
     $     PAUSE 'ROOT MUST BE BRACKETED FOR ZBRENT'
      C=B
      FC=FB
      DO  ITER=1,ITMAX
         IF((FB.GT.0.0D0.AND.FC.GT.0.0D0)
     $        .OR.(FB.LT.0.0D0.AND.FC.LT.0.0D0))THEN
            C=A   
            FC=FA
            D=B-A
            E=D
         ENDIF
         IF(DABS(FC).LT.DABS(FB)) THEN
            A=B
            B=C
            C=A
            FA=FB
            FB=FC
            FC=FA
         ENDIF
         TOL1=2.0D0*EPS*ABS(B)+0.5D0*TOL     
         XM=0.5D0*(C-B)
         IF(DABS(XM).LE.TOL1 .OR. FB.EQ.0.0D0)THEN
            ZBRENT=B
            RETURN
         ENDIF
         IF(DABS(E).GE.TOL1 .AND. DABS(FA).GT.DABS(FB)) THEN
            S=FB/FA                      
            IF(A.EQ.C) THEN
               P=2.0D0*XM*S
               Q=1.0D0-S
            ELSE
               Q=FA/FC
               R=FB/FC
               P=S*(2.0D0*XM*Q*(Q-R)-(B-A)*(R-1.0D0))
               Q=(Q-1.0D0)*(R-1.0D0)*(S-1.0D0)
               
            ENDIF
            IF(P.GT.0.0D0) Q=-Q 
            P=DABS(P)
            IF(2.0D0*P .LT. DMIN1(3.0D0*XM*Q-DABS(TOL1*Q),DABS(E*Q))) 
     $           THEN
               E=D                   
               D=P/Q
            ELSE
               D=XM                 
               E=D
            ENDIF
         ELSE                           
            D=XM
            E=D
         ENDIF
         A=B                            
         FA=FB
         IF(DABS(D) .GT. TOL1) THEN       
            B=B+D
         ELSE
            B=B+DSIGN(TOL1,XM)
         ENDIF
         FB=FUNC(B)
      ENDDO
      PAUSE 'ZBRENT EXCEEDING MAXIMUM ITERATIONS'
      ZBRENT=B
      RETURN
      END

C===============================================================================
C===============================================================================
C===============================================================================

      SUBROUTINE ZBRAC(FUNC,X1,X2,SUCCES)

C     GIVEN A FUNCTION FUNC AND AN INITIAL GUESSED RANGE X1 TO X2, THE ROUTINE 
C     EXPANDS THE RANGE GEOMETRICALLY UNTIL A ROOT IS BRACKETED BY THE RETURNED
C     VALUES X1 AND X2 (IN WHICH CASE SUCCES RETURNS AS .TRUE.) OR UNTIL THE 
C     RANGE BECOMES UNACCEPTABLY LARGE(IN WHICH CASE SUCCES RETURNS AS .FALSE.).

      INTEGER NTRY
      DOUBLE PRECISION X1,X2,FUNC,FACTOR
      EXTERNAL FUNC
      INTEGER J
      REAL F1,F2
      LOGICAL SUCCES
      PARAMETER (FACTOR=1.6D0,NTRY=500)

      IF(X1.EQ.X2) PAUSE 'YOU HAVE TO GUESS AN INITIAL RANGE IN ZBRAC'
      F1=FUNC(X1)
      F2=FUNC(X2)
      SUCCES=.TRUE.
      DO J=1,NTRY
         IF(F1*F2.LT.0.0D0) RETURN
         IF(ABS(F1).LT.ABS(F2)) THEN
            X1=X1+FACTOR*(X1-X2)
            F1=FUNC(X1)
         ELSE
            X2=X2+FACTOR*(X2-X1)
            F2=FUNC(X2)
         ENDIF
      ENDDO
      SUCCES=.FALSE.
      RETURN
      END

C===============================================================================
C===============================================================================
C===============================================================================

      SUBROUTINE ZBRAK(FX,X1,X2,N,XB1,XB2,NB)

C     GIVEN A FUNCTION FX DEFINED ON THE INTERVAL FROM X1-X2 SUBDIVIDE THE 
C     INTERVAL INTO N EQUALLY SPACED SEGMENTS, AND SEARCH FOR ZERO CROSSINGS
C     OF THE FUNCTION. NB IS INPUT AS THE MAXIMUM NUMBER OF ROOTS SOUGHT, 
C     AND IS RESET TO THE NUMBER OF BRACKETING PAIRS XB1(1:NB),
C     XB2(1:NB) THAT ARE FOUND.

      INTEGER N,NB
      DOUBLE PRECISION X1,X2,XB1(NB),XB2(NB),FX
      EXTERNAL FX
      INTEGER I,NBB
      DOUBLE PRECISION DX,FC,FP,X

      NBB=0
      X=X1
      DX=(X2-X1)/DBLE(N)
      FP=FX(X)
      DO I=1,N                     
         X=X+DX
         FC=FX(X)
         IF(FC*FP.LE.0.0D0) THEN       
            NBB=NBB+1
            XB1(NBB)=X-DX
            XB2(NBB)=X
            IF(NBB.EQ.NB)GOTO 1
         ENDIF
         FP=FC
      ENDDO
 1    CONTINUE
      NB=NBB
      RETURN
      END

C===============================================================================
C===============================================================================
C===============================================================================

      SUBROUTINE ZEROIN(F,B,C,RE,AE,IFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER IFLAG
      EXTERNAL F
C
C     SANDIA MATHEMATICAL PROGRAM LIBRARY
C     APPLIED MATHEMATICS DIVISION 2646
C     SANDIA LABORATORIES
C     ALBUQUERQUE, NEW MEXICO  87185
C     CONTROL DATA 6600/7600  VERSION 8.1  AUGUST 1980
C                   *************************
C                   *       ISSUED BY       *
C                   *  SANDIA LABORATORIES, *
C                   *   A PRIME CONTRACTOR  *
C                   ********     TO THE     *
C                          *  UNITED STATES *
C                          *   DEPARTMENT   *
C                          *       OF       *
C                          *     ENERGY     *
C      *********************  ---NOTICE---  *********************
C      *THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED*
C      *  BY THE UNITED STATES GOVERNMENT.  NEITHER THE UNITED  *
C      *   STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,   *
C      *               NOR ANY OF THEIR EMPLOYEES,              *
C      * NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR *
C      * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR  *
C      * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE  *
C      *          **********    ACCURACY,   **********          *
C      *          *        *  COMPLETENESS  *        *          *
C      *          *        *  OR USEFULNESS *        *          *
C      *          *        *     OF ANY     *        *          *
C      *          *        *  INFORMATION,  *        *          *
C      *          *        *   APPARATUS,   *        *          *
C      *       ****        *     PRODUCT    *        ****       *
C      *       *           *   OR PROCESS   *           *       *
C      *       *           *   DISCLOSED,   *           *       *
C      *       *           *  OR REPRESENTS *           *       *
C      *       *          **    THAT ITS    **          *       *
C      *       *          **  USE WOULD NOT **          *       *
C      *********          **    INFRINGE    **          *********
C                         **    PRIVATELY   **
C                         **      OWNED     **
C                         **     RIGHTS.    **
C                         **                **
C                         **                **
C                         **                **
C                         ********************
C
C     BASED ON A METHOD BY T J DEKKER
C     WRITTEN BY L F SHAMPINE AND H A WATTS
C     MODIFIED FOR THE MATH LIBRARY BY C B BAILEY
C
C     ABSTRACT
C        ZEROIN SEARCHES FOR A ZERO OF A FUNCTION F(X) BETWEEN
C        THE GIVEN VALUES B AND C UNTIL THE WIDTH OF THE INTERVAL
C        (B,C) HAS COLLAPSED TO WITHIN A TOLERANCE SPECIFIED BY
C        THE STOPPING CRITERION, ABS(B-C) .LE. 2.*(RW*ABS(B)+AE).
C        THE METHOD USED IS AN EFFICIENT COMBINATION OF BISECTION AND
C        THE SECANT RULE.  IN ORDER TO INSURE THAT ZEROIN WILL CONVERGE
C        TO A ZERO, THE USER SHOULD PICK VALUES FOR B AND C AT WHICH
C        THE FUNCTION DIFFERS IN SIGN.
C
C     DESCRIPTION OF ARGUMENTS
C     F,B,C,RE AND AE ARE INPUT PARAMETERS
C     B,C AND IFLAG ARE OUTPUT PARAMETERS
C        F     - NAME OF THE REAL VALUED EXTERNAL FUNCTION.  THIS NAME
C                MUST BE IN AN EXTERNAL STATEMENT IN THE CALLING
C                PROGRAM.  F MUST BE A FUNCTION OF ONE REAL ARGUMENT.
C        B     - ONE END OF THE INTERVAL (B,C).  THE VALUE RETURNED FOR
C                B USUALLY IS THE BETTER APPROXIMATION TO A ZERO OF F.
C        C     - THE OTHER END OF THE INTERVAL (B,C)
C        RE    - RELATIVE ERROR USED FOR RW IN THE STOPPING CRITERION.
C                IF THE REQUESTED RE IS LESS THAN MACHINE PRECISION,
C                THEN RW IS SET TO APPROXIMATELY MACHINE PRECISION.
C        AE    - ABSOLUTE ERROR USED IN THE STOPPING CRITERION.  IF THE
C                GIVEN INTERVAL (B,C) CONTAINS THE ORIGIN, THEN A
C                NONZERO VALUE SHOULD BE CHOSEN FOR AE.
C        IFLAG - A STATUS CODE.  USER MUST CHECK IFLAG AFTER EACH CALL.
C                CONTROL RETURNS TO THE USER FROM ZEROIN IN ALL CASES.
C                XERROR DOES NOT PROCESS DIAGNOSTICS IN THESE CASES.
C                 1 B IS WITHIN THE REQUESTED TOLERANCE OF A ZERO.
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE, THE FUNCTION CHANGES SIGN IN (B,C), AND
C                   F(X) DECREASED IN MAGNITUDE AS (B,C) COLLAPSED.
C                 2 F(B) = 0.  HOWEVER, THE INTERVAL (B,C) MAY NOT HAVE
C                   COLLAPSED TO THE REQUESTED TOLERANCE.
C                 3 B MAY BE NEAR A SINGULAR POINT OF F(X).
C                   THE INTERVAL (B,C) COLLAPSED TO THE REQUESTED
C                   TOLERANCE AND THE FUNCTION CHANGES SIGN IN (B,C) BUT
C                   F(X) INCREASED IN MAGNITUDE AS (B,C) COLLAPSED,I.E.
C                     ABS(F(B OUT)) .GT. MAX(ABS(F(B IN)),ABS(F(C IN)))
C                 4 NO CHANGE IN SIGN OF F(X) WAS FOUND ALTHOUGH THE
C                   INTERVAL (B,C) COLLAPSED TO THE REQUESTED TOLERANCE.
C                   THE USER MUST EXAMINE THIS CASE AND DECIDE WHETHER
C                   B IS NEAR A LOCAL MINIMUM OF F(X), OR B IS NEAR A
C                   ZERO OF EVEN MULTIPLICITY, OR NEITHER OF THESE.
C                 5 TOO MANY (.GT. 500) FUNCTION EVALUATIONS USED.
C
C     REFERENCES
C       1.  L F SHAMPINE AND H A WATTS, ZEROIN, A ROOT-SOLVING CODE,
C           SC-TM-70-631, SEPT 1970.
C       2.  T J DEKKER, FINDING A ZERO BY MEANS OF SUCCESSIVE LINEAR
C           INTERPOLATION, *CONSTRUCTIVE ASPECTS OF THE FUNDAMENTAL
C           THEOREM OF ALGEBRA*, EDITED BY B DEJON AND P HENRICI, 1969.

      ER = 2.0D0 * D1MACH(4)
C
C     INITIALIZE
      RW=DMAX1(RE,ER)
      AW=DMAX1(AE,0.0D0)
      IC=0
      ACBS=DABS(B-C)
      A=C
      T=A
      FA=F(T)
      T=B
      FB=F(T)
      FC=FA
      KOUNT=2
      FX=DMAX1(DABS(FB),DABS(FC))
C
    1 IF (DABS(FC) .GE. DABS(FB)) GO TO 2
C     PERFORM INTERCHANGE
      A=B
      FA=FB
      B=C
      FB=FC
      C=A
      FC=FA
C
    2 IF (FB .EQ. 0.0D0) GO TO 11
      CMB=0.5D0*(C-B)
      ACMB=DABS(CMB)
      TOL=RW*DABS(B)+AW
C
C     TEST STOPPING CRITERION
      IF (ACMB .LE. TOL) GO TO 10
C
C     CALCULATE NEW ITERATE IMPLICITLY AS B+P/Q
C     WHERE WE ARRANGE P .GE. 0.
C     THE IMPLICIT FORM IS USED TO PREVENT OVERFLOW.
      P=(B-A)*FB
      Q=FA-FB
      IF (P .GE. 0.0D0) GO TO 3
      P=-P
      Q=-Q
C
C     UPDATE A AND CHECK FOR SATISFACTORY REDUCTION
C     IN THE SIZE OF OUR BOUNDING INTERVAL.
    3 A=B
      FA=FB
      IC=IC+1
      IF (IC .LT. 4) GO TO 4
      IF (8.0D0*ACMB .GE. ACBS) GO TO 6
      IC=0
      ACBS=ACMB
C
C     TEST FOR TOO SMALL A CHANGE
    4 IF (P .GT. DABS(Q)*TOL) GO TO 5
C
C     INCREMENT BY TOLERANCE
      B=B+DSIGN(TOL,CMB)
      GO TO 7
C
C     ROOT OUGHT TO BE BETWEEN B AND (C+B)/2.0D0
    5 IF (P .GE. CMB*Q) GO TO 6
C
C     INTERPOLATE
      B=B+P/Q
      GO TO 7
C
    6 B=0.5D0*(C+B)
C     BISECT
C
C     HAVE COMPLETED COMPUTATION FOR NEW ITERATE B
    7 T=B
      FB=F(T)
      IF (FB .EQ. 0.0D0) GO TO 11
C
C     DECIDE WHETHER NEXT STEP IS INTERPOLATION OR EXTRAPOLATION
      IF (DSIGN(1.0D0,FB) .NE. DSIGN(1.0D0,FC)) GO TO 8
      C=A
      FC=FA
    8 KOUNT=KOUNT+1
      IF (KOUNT .GT. 5000) GO TO 15 !500) GO TO 15
      GO TO 1
C
C
C     FINISHED. PROCESS RESULTS FOR PROPER SETTING OF IFLAG
C
   10 IF (DSIGN(1.0D0,FB) .EQ. DSIGN(1.0D0,FC)) GO TO 13
      IF (DABS(FB) .GT. FX) GO TO 12
      IFLAG = 1
      RETURN
   11 IFLAG = 2
      RETURN
   12 IFLAG = 3
      RETURN
   13 IFLAG = 4
      RETURN
   15 IFLAG = 5
      RETURN
      END

C*KIN change to D1MACH from D1MACH.....2nd APRIL 92
      DOUBLE PRECISION FUNCTION D1MACH (IDUM)
      INTEGER IDUM
C-----------------------------------------------------------------------
C THIS ROUTINE COMPUTES THE UNIT ROUNDOFF OF THE MACHINE IN DOUBLE
C PRECISION.  THIS IS DEFINED AS THE SMALLEST POSITIVE MACHINE NUMBER
C U SUCH THAT  1.0D0 + U .NE. 1.0D0 (IN DOUBLE PRECISION).
C-----------------------------------------------------------------------
      DOUBLE PRECISION U, COMP
      U = 1.0D0
 10   U = U*0.5D0
      COMP = 1.0D0 + U
      IF (COMP .NE. 1.0D0) GO TO 10
      D1MACH = U*2.0D0
      RETURN
C----------------------- END OF FUNCTION D1MACH ------------------------
      END


