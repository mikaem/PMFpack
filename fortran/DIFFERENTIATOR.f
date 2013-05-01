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
      DOUBLE PRECISION FUNCTION DFRIDR(FUNC,X,H,ERR)
C     USES FUNC
C     RETURNS THE DERIVATIVE OF A FUNCTION FUNC AT A POINT X BY RIDDERS' 
C     METHOD OF POLYNOMIAL EXTRAPOLATION. THE VALUE H IS INPUT AS AN ESTIMATED 
C     INITIAL STEPSIZE; IT NEED NOT BE SMALL, BUT RATHER SHOULD BE AN INCREMENT 
C     IN X OVER WHICH FUNC CHANGES SUBSTANTIALLY. AN ESTIMATE OF THE ERROR IN 
C     THE DERIVATIVE IS RETURNED AS ERR.
C     PARAMETERS: STEPSIZE IS DECREASED BY CON AT EACH ITERATION. MAX SIZE OF 
C     TABLEAU IS SET BY NTAB. RETURN WHEN ERROR IS SAFE WORSE THAN THE BEST SO FAR.
C     SOURCE: W.H. PRESS, B.P. FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING,
C             NUMERICAL RECIPES IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING,
C             CAMBRIDGE UNIVERSITY PRESS, 1992
      IMPLICIT NONE
      INTEGER NTAB
      DOUBLE PRECISION ERR,H,X,FUNC,CON,CON2,BIG,SAFE
      PARAMETER(CON=1.2D0,CON2=CON*CON,BIG=1.0D+30,NTAB=40,SAFE=2.0D0) !CON WAS 1.15D0 NTAB WAS 100
      EXTERNAL FUNC
      INTEGER I,J
      DOUBLE PRECISION ERRT,FAC,HH,A(NTAB,NTAB)
      
      IF(H.EQ.0.0D0) PAUSE 'H MUST BE NONZERO IN DFRIDR'
      HH=H
      A(1,1)=(FUNC(X+HH)-FUNC(X-HH))/(2.0D0*HH)
      ERR=BIG
      DO I=2,NTAB            
         HH=HH/CON              
         A(1,I)=(FUNC(X+HH)-FUNC(X-HH))/(2.0D0*HH) 
         FAC=CON2
         DO J=2,I            
            A(J,I)=(A(J-1,I)*FAC-A(J-1,I-1))/(FAC-1.0D0)      
            FAC=CON2*FAC
            ERRT=DMAX1(DABS(A(J,I)-A(J-1,I)),DABS(A(J,I)-A(J-1,I-1)))
            
            IF (ERRT.LE.ERR) THEN               
               ERR=ERRT
               DFRIDR=A(J,I)
            ENDIF
         ENDDO
         IF(DABS(A(I,I)-A(I-1,I-1)).GE.SAFE*ERR) RETURN
      ENDDO
      RETURN 
      END



C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 
C
C     THE FOLLOWING FUNCTIONS, D2FRIDR AND D2FRIDRMIXED, ARE MODIFIED VERSIONS  
C     OF THE FUNCTION DFRIDR PROVIDED ABOVE. THE MODIFICATIONS EXTEND THE APPLICA-
C     BILITY OF RIDDER'S ALGORITH (AS SUPPLIED IN [W.H. PRESS, B.P. FLANNERY, S.A. 
C     TEUKOLSKY, W.T. VETTERLING, NUMERICAL RECIPES IN FORTRAN: THE ART OF 
C     SCIENTIFIC COMPUTING, CAMBRIDGE UNIVERSITY PRESS, 1992]) TO SECOND ORDER
C     DERIVATIVES (D2FRIDR) AND MIXED DERIVATIVES (D2FRIDRMIXED).
C
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION D2FRIDR(FUNC,X,H,ERR)
  
C     USES FUNC
C     RETURNS THE SECOND DERIVATIVE OF A FUNCTION FUNC AT A POINT X BY RIDDERS' 
C     METHOD OF POLYNOMIAL EXTRAPOLATION. THE VALUE H IS INPUT AS AN ESTIMATED 
C     INITIAL STEPSIZE; IT NEED NOT BE SMALL, BUT RATHER SHOULD BE AN INCREMENT 
C     IN X OVER WHICH FUNC CHANGES SUBSTANTIALLY. AN ESTIMATE OF THE ERROR IN 
C     THE DERIVATIVE IS RETURNED AS ERR.
C     PARAMETERS: STEPSIZE IS DECREASED BY CON AT EACH ITERATION. MAX SIZE OF 
C     TABLEAU IS SET BY NTAB. RETURN WHEN ERROR IS SAFE WORSE THAN THE BEST SO FAR.

      IMPLICIT NONE
      INTEGER NTAB
      DOUBLE PRECISION ERR,H,X,FUNC,CON,CON2,BIG,SAFE
      PARAMETER(CON=1.2D0,CON2=CON*CON,BIG=1.0D+30,NTAB=40,SAFE=2.0D0) 
      EXTERNAL FUNC
      INTEGER I,J
      DOUBLE PRECISION ERRT,FAC,HH,A(NTAB,NTAB)
      
      IF(H.EQ.0.0D0) PAUSE 'H MUST BE NONZERO IN DFRIDR'
      HH=H
      A(1,1)=(FUNC(X+HH)-2.0D0*FUNC(X)+FUNC(X-HH))/(HH**2.0D0)
      ERR=BIG
      DO I=2,NTAB            
         HH=HH/CON              
         A(1,I)=(FUNC(X+HH)-2.0D0*FUNC(X)+FUNC(X-HH))/(HH**2.0D0) 
         FAC=CON2
         DO J=2,I            
            A(J,I)=(A(J-1,I)*FAC-A(J-1,I-1))/(FAC-1.0D0)      
            FAC=CON2*FAC
            ERRT=DMAX1(DABS(A(J,I)-A(J-1,I)),DABS(A(J,I)-A(J-1,I-1)))
            
            IF (ERRT.LE.ERR) THEN               
               ERR=ERRT
               D2FRIDR=A(J,I)
            ENDIF
         ENDDO
         IF(DABS(A(I,I)-A(I-1,I-1)).GE.SAFE*ERR) RETURN
      ENDDO
      RETURN 
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION D2FRIDRMIXED(FUNC,X,Y,HX,HY,ERR)
 
C     USES FUNC
C     RETURNS THE SECOND DERIVATIVE OF A FUNCTION FUNC AT A POINT X BY RIDDERS' 
C     METHOD OF POLYNOMIAL EXTRAPOLATION. THE VALUE H IS INPUT AS AN ESTIMATED 
C     INITIAL STEPSIZE; IT NEED NOT BE SMALL, BUT RATHER SHOULD BE AN INCREMENT 
C     IN X OVER WHICH FUNC CHANGES SUBSTANTIALLY. AN ESTIMATE OF THE ERROR IN 
C     THE DERIVATIVE IS RETURNED AS ERR.
C     PARAMETERS: STEPSIZE IS DECREASED BY CON AT EACH ITERATION. MAX SIZE OF 
C     TABLEAU IS SET BY NTAB. RETURN WHEN ERROR IS SAFE WORSE THAN THE BEST SO FAR.

      IMPLICIT NONE
      INTEGER NTAB
      DOUBLE PRECISION ERR,HX,HY,X,Y,FUNC,CON,CON2,BIG,SAFE
      PARAMETER(CON=1.2D0,CON2=CON*CON,BIG=1.0D+30,NTAB=40,SAFE=2.0D0)
      EXTERNAL FUNC
      INTEGER I,J
      DOUBLE PRECISION ERRT,FAC,HHX,HHY,A(NTAB,NTAB)
      
      IF(HX.EQ.0.0D0) PAUSE 'HX MUST BE NONZERO IN DFRIDR' 
      IF(HY.EQ.0.0D0) PAUSE 'HY MUST BE NONZERO IN DFRIDR'
      HHX=HX
      HHY=HY
      A(1,1)=(FUNC(X+HHX,Y+HHY)-FUNC(X+HHX,Y-HHY)
     $     -FUNC(X-HHX,Y+HHY)+FUNC(X-HHX,Y-HHY))/(4*HHX*HHY)
      ERR=BIG
      DO I=2,NTAB            
         HHX=HHX/CON 
         HHY=HHY/CON 
         A(1,I)=(FUNC(X+HHX,Y+HHY)-FUNC(X+HHX,Y-HHY)
     $     -FUNC(X-HHX,Y+HHY)+FUNC(X-HHX,Y-HHY))/(4.0D0*HHX*HHY)
         FAC=CON2
         DO J=2,I            
            A(J,I)=(A(J-1,I)*FAC-A(J-1,I-1))/(FAC-1.0D0)      
            FAC=CON2*FAC
            ERRT=DMAX1(DABS(A(J,I)-A(J-1,I)),DABS(A(J,I)-A(J-1,I-1)))
            
            IF (ERRT.LE.ERR) THEN               
               ERR=ERRT
               D2FRIDRMIXED=A(J,I)
            ENDIF
         ENDDO
         IF(DABS(A(I,I)-A(I-1,I-1)).GE.SAFE*ERR) RETURN
      ENDDO
      RETURN 
      END
