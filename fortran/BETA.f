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
C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
C================================================================================== 
C     PURPOSE: COMPUTES THE BETA PROBABILITY DENSITY FUNCTION
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            ETA       MIXTURE FRACTION GRID          DOUBLE PRECISION (ARRAY)
C            NETA      SIZE OF ETA AND PDF            INTEGER
C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION
C
C     OUTPUT:
C
C            PDF       PROBABILITY DENSITY FUNCTION   DOUBLE PRECISION (ARRAY)
C================================================================================== 
      IMPLICIT NONE
      INCLUDE 'formats.h'
      INTEGER IETA,NETA
      DOUBLE PRECISION MFMEAN,MFVAR,ETA(NETA),PDF(NETA)
      DOUBLE PRECISION G,V,W,B
      DOUBLE PRECISION BETA
      EXTERNAL BETA
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE

      CALL CHECKPARMS(MFMEAN,MFVAR)
      
      G = (MFMEAN*(D_ONE-MFMEAN)/MFVAR) - D_ONE
      V = MFMEAN*G
      W = (D_ONE-MFMEAN)*G
      B = BETA(V,W)
      DO IETA = 1,NETA
         PDF(IETA) = (ETA(IETA)**(V-D_ONE))
     $        *((D_ONE-ETA(IETA))**(W-D_ONE))/B
      ENDDO
      ! A TREATMENT IS STILL REQUIRED FOR SMALL PROBABILITIES AND AT THE BOUNDARIES!!
      
      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'BETA: PROBABILITY DENSITY FUNCTION:'
         WRITE(*,50)'============================================='
         WRITE(*,100) 'V =',V
         WRITE(*,100) 'W =',W
         WRITE(*,100) 'B =',B
         WRITE(*,50)'============================================='
         WRITE(*,200)'INDEX','ETA','PDF'
         WRITE(*,50)'============================================='
         DO IETA = 1,NETA
            WRITE(*,300) IETA, ETA(IETA), PDF(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,50)' '
      ENDIF

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE BETA_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,VEL,CV)
C========================================================================================== 
C     PURPOSE: COMPUTES THE CONDITIONAL VELOCITY USING THE BETA-PDF
C==========================================================================================
C            VARIABLE   DESCRIPTION                     DATA TYPE
C            --------   -----------                     --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION GRI            DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA AND CSDRI           INTEGER
C            MFMEAN     MIXTURE FRACTION MEAN           DOUBLE PRECISION 
C            MFVAR      MIXTURE FRACTION VARIANCE       DOUBLE PRECISION
C            MFMEANGRAD GRAD. OF MIX. FRAC. MEAN        DOUBLE PRECISION (ARRAY, SIZE=3) 
C            MFVARGRAD  GRAD. OF MIX. FRAC. VARIANCE    DOUBLE PRECISION (ARRAY, SIZE=3)
C            DT         TURBULENT DIFFUSIVITY           DOUBLE PRECISION 
C            VEL        MEAN VELOCITY VECTOR            DOUBLE PRECISION (ARRAY, SIZE=3)
C
C     OUTPUT:
C
C            CV         COND. VELOCITY                  DOUBLE PRECISION (ARRAY, SIZE=NETA)
C========================================================================================== 
      IMPLICIT NONE
      INCLUDE 'formats.h'
      INTEGER I,IETA,NETA
      DOUBLE PRECISION MFMEAN,MFVAR,DT
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3),VEL(3)
      DOUBLE PRECISION ETA(NETA),CV(3,NETA)
      DOUBLE PRECISION G,V,W,B,PDF(NETA)
      DOUBLE PRECISION V_VAR,V_MEAN,W_VAR,W_MEAN,PDF_V,PDF_W,C1,C2
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION BETA,DPSI
      EXTERNAL BETA,DPSI
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE

      CALL CHECKPARMS(MFMEAN,MFVAR)

      G = (MFMEAN*(D_ONE-MFMEAN)/MFVAR) - D_ONE
      V = MFMEAN*G
      W = (D_ONE-MFMEAN)*G
      B = BETA(V,W)
      V_MEAN = (D_TWO*MFMEAN - D_THREE*(MFMEAN**D_TWO))/MFVAR - D_ONE
      W_MEAN = ((D_ONE-MFMEAN) * (D_ONE - D_THREE*MFMEAN))/MFVAR + D_ONE
      V_VAR = -((MFMEAN**D_TWO)*(D_ONE-MFMEAN))/(MFVAR**D_TWO)
      W_VAR = -(MFMEAN*((D_ONE-MFMEAN)**D_TWO))/(MFVAR**D_TWO)
      C1 = -DPSI(V) + DPSI(V+W)
      C2 = -DPSI(W) + DPSI(V+W)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

      DO IETA = 1,NETA
         PDF_V = PDF(IETA)*(DLOG(ETA(IETA)) + C1)
         PDF_W = PDF(IETA)*(DLOG(D_ONE-ETA(IETA)) + C2)
         DO I = 1,3
            CV(I,IETA) =  VEL(I) - (DT/ PDF(IETA))
     $           *((PDF_V*V_MEAN + PDF_W*W_MEAN)*MFMEANGRAD(I)
     $           + (PDF_V*V_VAR + PDF_W*W_VAR)*MFVARGRAD(I))
         ENDDO
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'BETA: CONDITIONAL VELOCITY:'
         WRITE(*,50)'============================================='
         WRITE(*,400)'INDEX','ETA','CV_X','CV_Y','CV_Z'
         DO IETA = 1,NETA
            WRITE(*,500) IETA,ETA(IETA),CV(1,IETA),CV(2,IETA),CV(3,IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,50)' '
      ENDIF

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE BETA_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [HOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE BETA-PDF
C==========================================================================================
C            VARIABLE   DESCRIPTION                     DATA TYPE
C            --------   -----------                     --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION GRI            DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA AND CSDRI           INTEGER
C            MFMEAN     MIXTURE FRACTION MEAN           DOUBLE PRECISION 
C            MFVAR      MIXTURE FRACTION VARIANCE       DOUBLE PRECISION
C            CHI        MEAN SCALAR DISSIPATION RATE    DOUBLE PRECISION 
C
C     OUTPUT:
C
C            CSDRH      COND. SCALAR DISSIPATION RATE   DOUBLE PRECISION (ARRAY, SIZE=NETA)
C==========================================================================================
      IMPLICIT NONE
      INCLUDE 'formats.h'
      INTEGER I,IETA,NETA
      DOUBLE PRECISION ETA(NETA),PDF(NETA),CSDRH(NETA)
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,DT
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3)
      DOUBLE PRECISION DIIDV(NETA),ERR_DIIDV(NETA),MAXERR_DIIDV
      DOUBLE PRECISION IM,H,H_NEW,ERR,ERR_OLD    
      LOGICAL VERBOSE_TEMP

      DOUBLE PRECISION DFRIDR,II_FIXEDMFMEAN
      EXTERNAL DFRIDR,II_FIXEDMFMEAN

      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION ETAPASS
      COMMON/EPASS/ETAPASS

      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION DIFF_MAXERR,DIFF_FACT
      COMMON/DIFFERENTIATORVARS/DIFF_MAXERR,DIFF_FACT
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      
      CALL CHECKPARMS(MFMEAN,MFVAR)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      IM = MFVAR + MFMEAN**D_TWO
      H  = DMIN1(MFVAR, MFMEAN - IM)
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     + DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     - DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
C      H  = DMIN1(1.0D-3, H * 0.2D0)

C      H = DMIN1(1.0D-4, H*0.05D0)

      MFMEANPASS = MFMEAN
      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         ERR = HUGE(1.0D0)
         H_NEW = H
         DO WHILE(ERR.GT.DIFF_MAXERR)
            ERR_OLD = ERR
            DIIDV(IETA) = DFRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
            IF(ERR.GT.ERR_OLD) THEN
               H_NEW = H_NEW*DIFF_FACT
               DIIDV(IETA) = DFRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
               EXIT
            ELSE
               H_NEW = H_NEW/DIFF_FACT
            ENDIF
         ENDDO
         ERR_DIIDV(IETA) = ERR
      ENDDO
     
C     BOUNDARY VALUES ARE KNOWN
      CSDRH(1)    = D_ZERO
      CSDRH(NETA) = D_ZERO
C     COMPUTE CSDR AT INTERNAL GRID POINTS     
      DO IETA = 2,NETA-1
         CSDRH(IETA) = (D_TWO/PDF(IETA))*CHI*DIIDV(IETA)                                         
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'BETA: COND. SCAL. DISS. RATE (HOMOGENEOUS):'
         WRITE(*,50)'============================================='
         WRITE(*,250)'INDEX','ETA','DIIDV','ERROR'
         DO IETA = 1,NETA
            WRITE(*,350) IETA, ETA(IETA), DIIDV(IETA), ERR_DIIDV(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,200)'INDEX','ETA','CSDRH'
         DO IETA = 1,NETA
            WRITE(*,300) IETA, ETA(IETA), CSDRH(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,50)' '
      ENDIF

      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE BETA_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,
     $     MFVARGRAD,DT,CHI,CSDRI)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [INHOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE BETA-PDF
C==========================================================================================
C            VARIABLE   DESCRIPTION                     DATA TYPE
C            --------   -----------                     --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION GRI            DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA AND CSDRI           INTEGER
C            MFMEAN     MIXTURE FRACTION MEAN           DOUBLE PRECISION 
C            MFVAR      MIXTURE FRACTION VARIANCE       DOUBLE PRECISION
C            MFMEANGRAD GRAD. OF MIX. FRAC. MEAN        DOUBLE PRECISION (ARRAY, SIZE=3) 
C            MFVARGRAD  GRAD. OF MIX. FRAC. VARIANCE    DOUBLE PRECISION (ARRAY, SIZE=3)
C            DT         TURBULENT DIFFUSIVITY           DOUBLE PRECISION 
C            CHI        MEAN SCALAR DISSIPATION RATE    DOUBLE PRECISION 
C
C     OUTPUT:
C
C            CSDRI      COND. SCALAR DISSIPATION RATE   DOUBLE PRECISION (ARRAY, SIZE=NETA)
C==========================================================================================
      IMPLICIT NONE
      INCLUDE 'formats.h'
      INTEGER I,IETA,NETA
      DOUBLE PRECISION ETA(NETA),PDF(NETA),CSDRI(NETA)
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,DT
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3)
      DOUBLE PRECISION PROD1,PROD2,PROD3
      DOUBLE PRECISION IM,H,H_NEW,ERR,ERR_OLD,H_MEAN,H_VAR
      DOUBLE PRECISION H_MEAN_NEW, H_VAR_NEW
      DOUBLE PRECISION DIIDV(NETA),ERR_DIIDV(NETA)
      DOUBLE PRECISION DIIDM(NETA),ERR_DIIDM(NETA)
      DOUBLE PRECISION D2IIDV2(NETA),ERR_D2IIDV2(NETA)
      DOUBLE PRECISION D2IIDM2(NETA),ERR_D2IIDM2(NETA)
      DOUBLE PRECISION D2IIDMDV(NETA),ERR_D2IIDMDV(NETA)
      
      LOGICAL VERBOSE_TEMP

      DOUBLE PRECISION DFRIDR,D2FRIDR ,D2FRIDRMIXED
      EXTERNAL DFRIDR,D2FRIDR ,D2FRIDRMIXED

      DOUBLE PRECISION II_FIXEDMFVAR,II_FIXEDMFMEAN,II_MIXEDDERIV
      EXTERNAL II_FIXEDMFVAR,II_FIXEDMFMEAN,II_MIXEDDERIV

      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION ETAPASS
      COMMON/EPASS/ETAPASS

      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION DIFF_MAXERR,DIFF_FACT
      COMMON/DIFFERENTIATORVARS/DIFF_MAXERR,DIFF_FACT
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      
      CALL CHECKPARMS(MFMEAN,MFVAR)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      IM = MFVAR + MFMEAN**D_TWO
      H  = DMIN1(MFVAR, MFMEAN - IM)
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     + DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     - DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
      !H  = DMIN1(1.0D-3, H*0.2D0)    

      H_MEAN = 2.0D0*H !H
      H_VAR  = H/4.0D0 !DMIN1(1.0D-4, H*0.05D0)
      !H_VAR  = DMIN1(1.0D-4, H*2.0D0)

      MFMEANPASS = MFMEAN
      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         ERR = HUGE(1.0D0)
         H_NEW = H_VAR
         DO WHILE(ERR.GT.DIFF_MAXERR)
            ERR_OLD = ERR
            DIIDV(IETA) = DFRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
            IF(ERR.GT.ERR_OLD) THEN
               H_NEW = H_NEW*DIFF_FACT
                DIIDV(IETA) = DFRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
               EXIT
            ELSE
               H_NEW = H_NEW/DIFF_FACT
            ENDIF
         ENDDO
         ERR_DIIDV(IETA) = ERR
      ENDDO

      MFMEANPASS = MFMEAN
      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         ERR = HUGE(1.0D0)
         H_NEW = H_VAR
         DO WHILE(ERR.GT.DIFF_MAXERR)
            ERR_OLD = ERR
            D2IIDV2(IETA) = D2FRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
            IF(ERR.GT.ERR_OLD) THEN
               H_NEW = H_NEW*DIFF_FACT
                D2IIDV2(IETA) = D2FRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
               EXIT
            ELSE
               H_NEW = H_NEW/DIFF_FACT
            ENDIF
         ENDDO
         ERR_D2IIDV2(IETA) = ERR
      ENDDO

      MFVARPASS = MFVAR
      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         ERR = HUGE(1.0D0)
         H_NEW = H_MEAN
         DO WHILE(ERR.GT.DIFF_MAXERR)
            ERR_OLD = ERR
            D2IIDM2(IETA) = D2FRIDR(II_FIXEDMFVAR,MFMEAN,H_NEW,ERR)
            IF(ERR.GT.ERR_OLD) THEN
               H_NEW = H_NEW*DIFF_FACT
                D2IIDM2(IETA) = D2FRIDR(II_FIXEDMFVAR,MFMEAN,H_NEW,ERR)
               EXIT
            ELSE
               H_NEW = H_NEW/DIFF_FACT
            ENDIF
         ENDDO
         ERR_D2IIDM2(IETA) = ERR
      ENDDO

      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         ERR = HUGE(1.0D0)
         H_MEAN_NEW = H_MEAN
         H_VAR_NEW  = H_VAR
         DO WHILE(ERR.GT.DIFF_MAXERR)
            ERR_OLD = ERR
            D2IIDMDV(IETA) = D2FRIDRMIXED(II_MIXEDDERIV,MFMEAN,MFVAR,
     $           H_MEAN_NEW,H_VAR_NEW,ERR)
            IF(ERR.GT.ERR_OLD) THEN
               H_MEAN_NEW = H_MEAN_NEW*DIFF_FACT
               H_VAR_NEW  = H_VAR_NEW*DIFF_FACT
               D2IIDMDV(IETA) = D2FRIDRMIXED(II_MIXEDDERIV,MFMEAN,MFVAR,
     $              H_MEAN_NEW,H_VAR_NEW,ERR)
               EXIT
            ELSE
               H_MEAN_NEW = H_MEAN_NEW/DIFF_FACT
               H_VAR_NEW  = H_VAR_NEW/DIFF_FACT
            ENDIF
         ENDDO
         ERR_D2IIDMDV(IETA) = ERR
      ENDDO

C==================================================

c$$$      IM = MFVAR + MFMEAN**D_TWO
c$$$      H  = DMIN1(MFVAR, MFMEAN - IM)
c$$$      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
c$$$     $     + DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
c$$$      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
c$$$     $     - DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
c$$$      !H  = DMIN1(1.0D-3, H * 0.2D0)   

c$$$      MFMEANPASS = MFMEAN
c$$$      DO IETA = 2,NETA-1
c$$$         ETAPASS = ETA(IETA)
c$$$         ERR = HUGE(1.0D0)
c$$$         H_NEW = H
c$$$         DO WHILE(ERR.GT.DIFF_MAXERR)
c$$$            ERR_OLD = ERR
c$$$            DIIDV(IETA) = DFRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
c$$$            IF(ERR.GT.ERR_OLD) THEN
c$$$               H_NEW = H_NEW*DIFF_FACT
c$$$                DIIDV(IETA) = DFRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
c$$$               EXIT
c$$$            ELSE
c$$$               H_NEW = H_NEW/DIFF_FACT
c$$$            ENDIF
c$$$         ENDDO
c$$$         ERR_DIIDV(IETA) = ERR
c$$$      ENDDO
c$$$
c$$$      MFMEANPASS = MFMEAN
c$$$      DO IETA = 2,NETA-1
c$$$         ETAPASS = ETA(IETA)
c$$$         ERR = HUGE(1.0D0)
c$$$         H_NEW = H
c$$$         DO WHILE(ERR.GT.DIFF_MAXERR)
c$$$            ERR_OLD = ERR
c$$$            D2IIDV2(IETA) = D2FRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
c$$$            IF(ERR.GT.ERR_OLD) THEN
c$$$               H_NEW = H_NEW*DIFF_FACT
c$$$                D2IIDV2(IETA) = D2FRIDR(II_FIXEDMFMEAN,MFVAR,H_NEW,ERR)
c$$$               EXIT
c$$$            ELSE
c$$$               H_NEW = H_NEW/DIFF_FACT
c$$$            ENDIF
c$$$         ENDDO
c$$$         ERR_D2IIDV2(IETA) = ERR
c$$$      ENDDO
c$$$
c$$$      MFVARPASS = MFVAR
c$$$      DO IETA = 2,NETA-1
c$$$         ETAPASS = ETA(IETA)
c$$$         ERR = HUGE(1.0D0)
c$$$         H_NEW = H
c$$$         DO WHILE(ERR.GT.DIFF_MAXERR)
c$$$            ERR_OLD = ERR
c$$$            D2IIDM2(IETA) = D2FRIDR(II_FIXEDMFVAR,MFMEAN,H_NEW,ERR)
c$$$            IF(ERR.GT.ERR_OLD) THEN
c$$$               H_NEW = H_NEW*DIFF_FACT
c$$$                D2IIDM2(IETA) = D2FRIDR(II_FIXEDMFVAR,MFMEAN,H_NEW,ERR)
c$$$               EXIT
c$$$            ELSE
c$$$               H_NEW = H_NEW/DIFF_FACT
c$$$            ENDIF
c$$$         ENDDO
c$$$         ERR_D2IIDM2(IETA) = ERR
c$$$      ENDDO
c$$$
c$$$      DO IETA = 2,NETA-1
c$$$         ETAPASS = ETA(IETA)
c$$$         ERR = HUGE(1.0D0)
c$$$         H_NEW = H
c$$$         DO WHILE(ERR.GT.DIFF_MAXERR)
c$$$            ERR_OLD = ERR
c$$$            D2IIDMDV(IETA) = D2FRIDRMIXED(II_MIXEDDERIV,MFMEAN,MFVAR,
c$$$     $           H_NEW,H_NEW,ERR)
c$$$            IF(ERR.GT.ERR_OLD) THEN
c$$$               H_NEW = H_NEW*DIFF_FACT
c$$$               D2IIDMDV(IETA) = D2FRIDRMIXED(II_MIXEDDERIV,MFMEAN,MFVAR,
c$$$     $              H_NEW,H_NEW,ERR)
c$$$               EXIT
c$$$            ELSE
c$$$               H_NEW = H_NEW/DIFF_FACT
c$$$            ENDIF
c$$$         ENDDO
c$$$         ERR_D2IIDMDV(IETA) = ERR
c$$$      ENDDO
C==================================================

C     COMPUTE THE SCALAR PRODUCTS
      PROD1 = D_ZERO            ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. MEAN
      PROD2 = D_ZERO            ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. VARIANCE
      PROD3 = D_ZERO            ! SCALAR PRODUCT OF MIX FRAC. MEAN AND MIX FRAC VARIANCE GRADS
      DO I = 1,3
         PROD1 = PROD1 + MFMEANGRAD(I)**D_TWO 
         PROD2 = PROD2 + MFVARGRAD(I)**D_TWO
         PROD3 = PROD3 + MFMEANGRAD(I)*MFVARGRAD(I)
      ENDDO
      
C     BOUNDARY VALUES ARE KNOWN
      CSDRI(1)    = D_ZERO
      CSDRI(NETA) = D_ZERO
C     COMPUTE CSDR AT INTERNAL GRID POINTS     
      DO IETA = 2,NETA-1

         CSDRI(IETA) = (D_TWO/PDF(IETA))
     $        *(-DIIDV(IETA)* (-CHI + D_TWO*DT*PROD1)
     $        + DT*(PROD2* D2IIDV2(IETA) + PROD1*D2IIDM2(IETA)
     $        + D_TWO*PROD3*D2IIDMDV(IETA)))
                                         
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'BETA: COND. SCAL. DISS. RATE (HOMOGENEOUS):'
         WRITE(*,50)'============================================='
         WRITE(*,250)'INDEX','ETA','DIIDV','ERROR'
         DO IETA = 1,NETA
            WRITE(*,350) IETA, ETA(IETA), DIIDV(IETA), ERR_DIIDV(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,250)'INDEX','ETA','D2IIDV2','ERROR'
         DO IETA = 1,NETA
            WRITE(*,350) IETA, ETA(IETA),D2IIDV2(IETA),ERR_D2IIDV2(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,250)'INDEX','ETA','D2IIDM2','ERROR'
         DO IETA = 1,NETA
            WRITE(*,350) IETA, ETA(IETA),D2IIDM2(IETA),ERR_D2IIDM2(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,250)'INDEX','ETA','D2IIDMDV','ERROR'
         DO IETA = 1,NETA
            WRITE(*,350)IETA,ETA(IETA),D2IIDMDV(IETA),ERR_D2IIDMDV(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,200)'INDEX','ETA','CSDRH'
         DO IETA = 1,NETA
            WRITE(*,300) IETA, ETA(IETA), CSDRI(IETA)
         ENDDO
         WRITE(*,50)'============================================='
         WRITE(*,50)' '
      ENDIF

      RETURN 
      END


C===========================================================================
C===========================================================================
C===========================================================================

      DOUBLE PRECISION FUNCTION II_FIXEDMFMEAN(MFVAR)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST- AND 
C              SECOND-ORDER PARTIAL DERIVATIVES OF II)ETA) WITH RESPECT TO THE 
C              [MIXTURE FRACTION VARIANCE]
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION MFVAR
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION BETAI,BETA
      EXTERNAL BETAI,BETA
      DOUBLE PRECISION B1,B2,G,P
      
      G = (MFMEANPASS*(1.0D0-MFMEANPASS)/MFVAR) - 1.0D0
      B1 = MFMEANPASS*G
      B2 = (1.0D0-MFMEANPASS)*G
      P = (ETAPASS**B1)*((1.0D0-ETAPASS)**B2)/BETA(B1+1.0D0,B2+1.0D0)
      II_FIXEDMFMEAN=(ETAPASS-MFMEANPASS)*BETAI(B1,B2,ETAPASS)+ MFVAR*P

      RETURN
      END

C===========================================================================
C===========================================================================
C===========================================================================

      DOUBLE PRECISION FUNCTION II_FIXEDMFVAR(MFMEAN)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE SECOND-ORDER 
C              PARTIAL DERIVATIVE OF II(ETA) WITH RESPECT TO THE 
C              [MIXTURE FRACTION MEAN]
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION MFMEAN
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION BETAI,BETA
      EXTERNAL BETAI,BETA
      DOUBLE PRECISION B1,B2,G,P
      
      G = (MFMEAN*(1.0D0-MFMEAN)/MFVARPASS) - 1.0D0
      B1 = MFMEAN*G
      B2 = (1.0D0-MFMEAN)*G
      P = (ETAPASS**B1)*((1.0D0-ETAPASS)**B2)/BETA(B1+1.0D0,B2+1.0D0)
      II_FIXEDMFVAR=(ETAPASS-MFMEAN)*BETAI(B1,B2,ETAPASS)+ MFVARPASS*P

      RETURN
      END

C===========================================================================
C===========================================================================
C===========================================================================

      DOUBLE PRECISION FUNCTION II_MIXEDDERIV(MFMEAN,MFVAR)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE MIXED 
C              DERIVATIVE OF II(ETA) WITH RESPECT TO THE [MIXTURE FRACTION 
C              MEAN AND THE MIXTURE FRACTION VARIANCE]
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION 
C==================================================================================
      IMPLICIT NONE
      DOUBLE PRECISION MFMEAN
      DOUBLE PRECISION MFVAR
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION BETAI,BETA
      EXTERNAL BETAI,BETA
      DOUBLE PRECISION B1,B2,G,P

      G = (MFMEAN*(1.0D0-MFMEAN)/MFVAR) - 1.0D0
      B1 = MFMEAN*G
      B2 = (1.0D0-MFMEAN)*G
      P = (ETAPASS**B1)*((1.0D0-ETAPASS)**B2)/BETA(B1+1.0D0,B2+1.0D0)
      II_MIXEDDERIV = (ETAPASS-MFMEAN)*BETAI(B1,B2,ETAPASS) + MFVAR*P

      RETURN
      END
