CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CLUSTERNEARZERO(BETA,XMIN,XMAX,NX,X)
C     TRANSFORMATION 1 IN COMPUTATIONAL FLUID MECHANICS AND HEAT TRANSFER
C     (TANNEHILL AND ANDERSON, 1997)
C     CLUSTERS POINTS NEAR X = MIN AS BETA APPROACHES 1.0
C     PARAMETERS: 1 < BETA < INF
      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION X(NX),XU(NX)
      DOUBLE PRECISION XMIN,XMAX
      DOUBLE PRECISION BETA,H,FACT
      INTEGER I
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO = 0.0D0,ONE = 1.0D0)
      CALL UNIFORMGRID(ZERO,ONE,NX,XU)
      H = XMAX-XMIN
      DO I = 1,NX
         FACT = ((BETA+1.0D0)/(BETA-1.0D0))**(1.0D0-XU(I))
         X(I) = H*((BETA+1.0D0)-(BETA-1.0D0)*FACT)/(FACT+1.0D0)
      ENDDO
      X(1) = XMIN
      X(NX) = XMAX
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CLUSTERATBOUNDARIES(ALPHA,BETA,XMIN,XMAX,NX,X)
C     TRANSFORMATION 2 IN COMPUTATIONAL FLUID MECHANICS AND HEAT TRANSFER
C     (TANNEHILL AND ANDERSON, 1997)     
C     PARAMETERS: 1 < BETA < INF
C                 ALPHA = 0   -> REFINEMENT NEAR X = XMAX
C                       = 0.5 -> EQUAL REFINEMENT NEAR X = XMIN AND X = XMAX
C                       > 0.5 -> MORE POINTS NEAR X = XMIN
C                       < 0.5 -> MORE POINTS NEAR X = XMAX
      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION X(NX),XU(NX)
      DOUBLE PRECISION XMIN,XMAX
      DOUBLE PRECISION ALPHA,BETA,H,FACT
      INTEGER I
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO = 0.0D0,ONE = 1.0D0)
      CALL UNIFORMGRID(ZERO,ONE,NX,XU)
      H = XMAX-XMIN
      DO I = 1,NX
         FACT = ((BETA+1.0D0)/(BETA-1.0D0))
     $        **((XU(I)-ALPHA)/(1.0D0-ALPHA))
         X(I) = H*((BETA+2.0D0*ALPHA)*FACT-BETA+2.0D0*ALPHA)
     $        /((2.0D0*ALPHA+1.0D0)*(1.0D0+FACT))
      ENDDO
      X(1) = XMIN
      X(NX) = XMAX
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE CLUSTERATINTERIORPT(XC,TAU,XMIN,XMAX,NX,X)
C     TRANSFORMATION 1 IN COMPUTATIONAL FLUID MECHANICS AND HEAT TRANSFER
C     (TANNEHILL AND ANDERSON, 1997)
C     LARGER TAU LEADS TO MORE CLUSTERING ABOUT XC
      IMPLICIT NONE
      INTEGER NX
      DOUBLE PRECISION X(NX),XU(NX)
      DOUBLE PRECISION XC,XMIN,XMAX
      DOUBLE PRECISION TAU,B,H
      INTEGER I
      DOUBLE PRECISION SINHF
      EXTERNAL SINHF
      DOUBLE PRECISION ZERO,ONE
      PARAMETER(ZERO = 0.0D0,ONE = 1.0D0)
      CALL UNIFORMGRID(ZERO,ONE,NX,XU)
      H = XMAX-XMIN
      B = (1.0D0/(2.0D0*TAU))
     $     *DLOG((1.0D0 + (DEXP(TAU)-1.0D0)*(XC/H))
     $     /(1.0D0 + (DEXP(-TAU)-1.0D0)*(XC/H)))
      DO I = 1,NX
         X(I) = XC*(1.0D0 + SINHF(TAU*(XU(I)-B))/SINHF(TAU*B))
      ENDDO
      X(1) = XMIN
      X(NX) = XMAX
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE UNIFORMGRID(XMIN,XMAX,NX,X)
      IMPLICIT NONE
      INTEGER IX,NX
      DOUBLE PRECISION X(NX),XMIN,XMAX,DELTA
      DELTA = (XMAX-XMIN)/DBLE(NX-1)
      X(1) = XMIN
      DO IX = 2,NX
         X(IX) = X(IX-1)+DELTA
      ENDDO
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE GRIDEXP(GRIDMIN,GRIDMAX,MS,NGRID,GRID,DGRID)
C     THIS SUBROUTINE CREATES AN EXPONENTIAL GRID DISTRIBUTION
C     BY MAPPING THE LINEAR SPACE [0 1] GIVEN A MAPPING STRENGTH MS
C     IF MS.EQ.0, THIS ROUTINE RETURNS A UNIFORM GRID    
      IMPLICIT NONE
      INTEGER NGRID,I
      DOUBLE PRECISION GRID(NGRID),DGRID(NGRID-1)
      DOUBLE PRECISION MS	!MAP STRENGTH
      DOUBLE PRECISION GRIDMIN,GRIDMAX,GRIDLENGTH
      DOUBLE PRECISION LINSPACE(NGRID-1),DLINSPACE,SPACING
      
      IF(MS.LT.0.0D0) THEN
         WRITE(*,*)'MAP STRENGTH (MS) MUST BE POSITIVE.'
         WRITE(*,*)'SOLBVER HAS STOPPED.'
         STOP
      ENDIF
      DLINSPACE  = 1.0D0/DBLE(NGRID-1)
      GRIDLENGTH = GRIDMAX-GRIDMIN
      SPACING  = GRIDLENGTH/DBLE(NGRID-1)
      DO I = 1,NGRID-1
         IF(MS.EQ.0) THEN
            GRID(I) = GRIDMIN + DBLE(I-1)*SPACING
         ELSE
            LINSPACE(I) = GRIDMIN+DBLE(I-1)*DLINSPACE
            GRID(I) = GRIDMIN + GRIDLENGTH
     $           *(DEXP(MS*LINSPACE(I))-1.0D0)/(DEXP(MS)-1.0D0)
         ENDIF
      ENDDO
      GRID(NGRID) = GRIDMAX
      DO I = 1,NGRID-1
         DGRID(I) = GRID(I+1)-GRID(I)
      ENDDO
      RETURN 
      END     

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
     
      DOUBLE PRECISION FUNCTION SINHF(X)
      IMPLICIT NONE
      DOUBLE PRECISION X
      SINHF = (DEXP(X)-DEXP(-X))/2.0D0
      RETURN
      END 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
