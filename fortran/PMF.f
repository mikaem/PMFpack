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
C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE INITIALISE
C================================================================================== 
C     PURPOSE: SPECIFIES SOME OPTIONS AND SETS SOLVER TOLERANCES FOR 
C              THE INTEGRATOR AND DIFFERENTIATOR(MUST BE CALLED ONCE        
C               PRIOR TO CALLING ANY OF THE BELOW SUBROUTINES AND BEFORE 
C              BEFORE LOOPING OVER ALL THE GRID POINTS OF THE PHYSICAL DOMAIN
C==================================================================================
      IMPLICIT NONE

      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE

      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR

      INTEGER ROOTF_METH
      COMMON/ROOTFINDERVARS1/ROOTF_METH
      DOUBLE PRECISION ROOTF_RE,ROOTF_AE
      COMMON/ROOTFINDERVARS2/ROOTF_RE,ROOTF_AE

      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM

      DOUBLE PRECISION DIFF_MAXERR,DIFF_FACT
      COMMON/DIFFERENTIATORVARS/DIFF_MAXERR,DIFF_FACT

C     VERBOSITY: [.FALSE.]OFF / [.TRUE.]ON
C
      VERBOSE = .TRUE. !.FALSE.

C     CONSTANTS
C
      PI     =  4.0D0*DATAN(1.0D0)
      D_ZERO  = 0.0D0
      D_ONE  =  1.0D0
      D_HALF =  0.5D0
      D_TWO  =  2.0D0
      D_THREE = 3.0D0
      I_TWO   = 2
      I_FOUR  = 4

C     BISECTION CONTROL (REQUIRED FOR SUBROUTINE ZEROIN)
C
      ROOTF_METH = 3          ! [1] RIDDER'S METHOD [2] BRENT'S METHOD/ [3] DEKKER'S METHOD
      ROOTF_RE = 1.0D-9
      ROOTF_AE = 1.0D-9

C     INTEGRATOR CONTROL (REQUIRED FOR QUADPACK)
C
      INTEG_LIM = 5000         ! MAXIMUM NUMBER OF SUBINTERVALS IN THE 
                               ! PARTITION OF THE GIVEN INTEGRATION INTERVAL
      INTEG_EPSABS = 1.0D-12   ! ABSOLUTE ACCURACY REQUESTED
      INTEG_EPSREL = 1.0D-10   ! RELATIVE ACCURACY REQUESTED

C     DIFFERENTIATOR CONTROL (REQUIRED FOR RIDDER'S ALGORITHM)
C
      DIFF_MAXERR = 1.0D-10    ! TARGET ERROR
                               ! THE WAY THE DIFFERENTIATOR IS SETUP DOES NOT NECESSARILY
                               ! RETURN AN ERROR SMALLER OR EQUAL THAN DIFF_MAXERR. THE ERROR
                               ! DEPENDS THE INITIAL STEP SIZE OF THE DIFFERENTIATOR, H (SEE THE
                               ! DESCRIPTION OF THE FACTOR DIFF_FACT BELOW). H IS MULTIPLIED BY 
                               ! DIFF_FACT IN SUCCESSIVIE ITERATIONS TO DECREASE THE ERROR UP TO
                               ! THE POINT WHERE THE ERROR STARTS INCREASING (H BECOMES TOO LARGE).
                               ! AT THIS POINT THE PREVIOUS VALUE OF H IS EMPLOYED TO COMPUTE THE 
                               ! DERIVATIVE AND AN ERROR LARGER THAN DIFF_FACT IS RETURNED.
      DIFF_FACT   = 10.0D0     ! FACTOR BY WHICH THE INITIAL STEP SIZE OF THE DIFFERENTIATOR, H,
                               ! IS MULTIPLIED IN SUCCESSIVE  ITERATIONS. THE INITIAL STEP SIZE 
                               ! IS CHOSEN FOLLOWING REF [2] AS MIN(A,B,C,D) WHERE
                               ! A = MFVAR
                               ! B = MFMEAN*(1-MFMEAN) -MFVAR
                               ! C = ABS((1-2*MFMEAN-SQRT(1+4*MFVAR))/2)
                               ! D = ABS((1-2*MFMEAN-SQRT(1-4*MFVAR))/2)

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CHECKPARMS(MFMEAN,MFVAR)
C================================================================================== 
C     PURPOSE: CHECKS THE VALIDITY OF THE VALUES OF THE MIXTURE FRACTION MEAN,
C              THE MIXTURE FRACTION AND VARIANCE, AND THE INTENSITY OF SEGREGATION
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
      INCLUDE 'formats.h'
      DOUBLE PRECISION MFMEAN,MFVAR,ISEG
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE

      ISEG = MFVAR/(MFMEAN*(D_ONE-MFMEAN))

      IF(VERBOSE) THEN
         WRITE(*,50)'============================================='
         WRITE(*,50)'PARAMETERS CHECK:'
         WRITE(*,50)'============================================='
         WRITE(*,100)'MIXTURE FRACTION MEAN     =', MFMEAN
         WRITE(*,100)'MIXTURE FRACTION VARIANCE =', MFVAR
         WRITE(*,100)'INTENSITY OF SEGREGATION  =', ISEG
         WRITE(*,50)' '
      ENDIF
      IF((MFMEAN.LT.D_ZERO).OR.(MFMEAN.GT.D_ONE)) THEN
         WRITE(*,50)'============================================='
         WRITE(*,50)'WRONG INPUT:'
         WRITE(*,50)'MIXTURE FRACTION MEAN IS NOT BETWEEN 0 AND 1.'
         WRITE(*,50)'COMPUTATIONS ABORTED.'
         WRITE(*,50)'============================================='
         STOP
      ELSE
         IF(VERBOSE)
     $        WRITE(*,50)'MIXTUR FRACTION MEAN      -> VALUE IS VALID.'
      ENDIF

      IF(MFVAR.LT.D_ZERO) THEN
         WRITE(*,50)'============================================='
         WRITE(*,50)'WRONG INPUT:'
         WRITE(*,50)'MIXTURE FRACTION VARIANCE IS LESS THAN 0.'
         WRITE(*,50)'COMPUTATIONS ABORTED.'
         WRITE(*,50)'============================================='
         STOP
      ELSE
         IF(VERBOSE)
     $        WRITE(*,50)'MIXTUR FRACTION VARAINCE  -> VALUE IS VALID.'   
      ENDIF

      IF((ISEG.LT.D_ZERO).OR.(ISEG.GT.D_ONE)) THEN
         WRITE(*,50)'WRONG INPUT:'
         WRITE(*,50)'INTENSITY OF SEGREGATION IS NOT BETWEEN 0 AND 1.'
         WRITE(*,50)'COMPUTATIONS ABORTED.'
         WRITE(*,50)'============================================='
         STOP
      ELSE
         IF(VERBOSE)
     $        WRITE(*,50)'INTENSITY OF SEGREGATION  -> VALUE IS VALID.'
      ENDIF
      IF(VERBOSE)
     $     WRITE(*,50)'============================================='

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE PMF_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
C================================================================================== 
C     PURPOSE: COMPUTES THE PMF PROBABILITY DENSITY FUNCTION
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
      DOUBLE PRECISION MFMEAN,MFVAR,ALPHA,TAU,SIGMA2
      DOUBLE PRECISION ETA(NETA),PDF(NETA),E(NETA),PHI(NETA)
      DOUBLE PRECISION MFMEANPASS,MFVARPASS,ALPHAPASS
      DOUBLE PRECISION ERFINV
      LOGICAL ISNAN,ISINF
      EXTERNAL ERFINV,ISNAN,ISINF
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      CALL CHECKPARMS(MFMEAN,MFVAR)

      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)

      CALL FIND_TAU(MFMEAN,MFVAR,TAU)
      
      SIGMA2 = D_ONE - D_TWO*TAU

      DO IETA = 1,NETA
         E(IETA) = ERFINV(D_TWO*ETA(IETA)-D_ONE)
         PHI(IETA) = ALPHA + D_TWO*DSQRT(TAU)*E(IETA) 
         PDF(IETA) = DSQRT(D_TWO*TAU/SIGMA2)
     $        * DEXP(E(IETA)**D_TWO - (PHI(IETA)**D_TWO)/(D_TWO*SIGMA2))
         IF(ISNAN(PDF(IETA))) PDF(IETA) = D_ZERO
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF: PROBABILITY DENSITY FUNCTION:'
         WRITE(*,50)'============================================='
         WRITE(*,100) 'ALPHA   =',ALPHA
         WRITE(*,150) 'TAU     =',TAU, 
     $        '-> MOD. OF ABS. ERR. OF INTEGRATION =',INTEG_ABSERR
         WRITE(*,100) 'SIGMA^2 =',SIGMA2
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

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      SUBROUTINE PMF_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,VEL,CV)
C========================================================================================== 
C     PURPOSE: COMPUTES THE CONDITIONAL VELOCITY USING THE PMF-PDF
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
      INTEGER IETA,NETA,I
      DOUBLE PRECISION MFMEAN,MFVAR,ALPHA,TAU,SIGMA2,DT
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3),VEL(3)
      DOUBLE PRECISION ETA(NETA),CV(3,NETA),E(NETA),PHI(NETA)
      DOUBLE PRECISION ALPHA_M,TAU_M,TAU_V
      DOUBLE PRECISION ERR1,ERR2,H
      DOUBLE PRECISION ERR_OLD, H_NEW,H_MEAN,H_VAR
      DOUBLE PRECISION ERFINV,DFRIDR,TAU_FUNC,TAU_FIXEDMFVAR,
     $     TAU_FIXEDMFMEAN
      EXTERNAL ERFINV,TAU_FUNC,DFRIDR,TAU_FIXEDMFVAR,
     $     TAU_FIXEDMFMEAN  
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION MFMEANPASS_DERIVS,MFVARPASS_DERIVS,ALPHAPASS
      COMMON/COMVARSMEAN/MFMEANPASS_DERIVS
      COMMON/COMVARSVAR/MFVARPASS_DERIVS
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION DIFF_MAXERR,DIFF_FACT
      COMMON/DIFFERENTIATORVARS/DIFF_MAXERR,DIFF_FACT

      DOUBLE PRECISION IM
   
      CALL CHECKPARMS(MFMEAN,MFVAR)

      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)

      ALPHA_M = -DSQRT(D_TWO*PI)*DEXP((ALPHA**D_TWO)/D_TWO)

      CALL FIND_TAU(MFMEAN,MFVAR,TAU)

      SIGMA2 = D_ONE - D_TWO*TAU


C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      IM = MFVAR + MFMEAN**D_TWO
      H  = DMIN1(MFVAR, MFMEAN - IM)
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     + DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     - DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
C      H  = DMIN1(1.0D-3, H * 0.2D0)
      
      H_MEAN = 2.0D0*H          !H
      H_VAR  = DMIN1(1.0D-4, H*0.05D0)
      !H_VAR  = DMIN1(1.0D-4, H*2.0D0)

      ERR1 = HUGE(1.0D0)
      H_NEW = H_MEAN
      MFVARPASS_DERIVS = MFVAR
      DO WHILE(ERR1.GT.DIFF_MAXERR)
         ERR_OLD = ERR1
         TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
         IF(ERR1.GT.ERR_OLD) THEN
            H_NEW = H_NEW*DIFF_FACT
            TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
            EXIT
         ELSE
            H_NEW = H_NEW/DIFF_FACT
         ENDIF
      ENDDO

      ERR2 = HUGE(1.0D0)
      H_NEW = H_VAR
      MFMEANPASS_DERIVS = MFMEAN
      DO WHILE(ERR2.GT.DIFF_MAXERR)
         ERR_OLD = ERR2
         TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
         IF(ERR2.GT.ERR_OLD) THEN
            H_NEW = H_NEW*DIFF_FACT
            TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
            EXIT
         ELSE
            H_NEW = H_NEW/DIFF_FACT
         ENDIF
      ENDDO
    
c$$$C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
c$$$      IM = MFVAR + MFMEAN**D_TWO
c$$$      H  = DMIN1(MFVAR, MFMEAN - IM)
c$$$      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
c$$$     $     + DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
c$$$      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
c$$$     $     - DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
c$$$C      H  = DMIN1(1.0D-3, H * 0.2D0)
c$$$
c$$$      ERR1 = HUGE(1.0D0)
c$$$      H_NEW = H
c$$$      MFVARPASS_DERIVS = MFVAR
c$$$      DO WHILE(ERR1.GT.DIFF_MAXERR)
c$$$         ERR_OLD = ERR1
c$$$         TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
c$$$         IF(ERR1.GT.ERR_OLD) THEN
c$$$            H_NEW = H_NEW*DIFF_FACT
c$$$            TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
c$$$            EXIT
c$$$         ELSE
c$$$            H_NEW = H_NEW/DIFF_FACT
c$$$         ENDIF
c$$$      ENDDO
c$$$
c$$$      ERR2 = HUGE(1.0D0)
c$$$      H_NEW = H
c$$$      MFMEANPASS_DERIVS = MFMEAN
c$$$      DO WHILE(ERR2.GT.DIFF_MAXERR)
c$$$         ERR_OLD = ERR2
c$$$         TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
c$$$         IF(ERR2.GT.ERR_OLD) THEN
c$$$            H_NEW = H_NEW*DIFF_FACT
c$$$            TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
c$$$            EXIT
c$$$         ELSE
c$$$            H_NEW = H_NEW/DIFF_FACT
c$$$         ENDIF
c$$$      ENDDO

      DO IETA = 1,NETA
         E(IETA) = ERFINV(D_TWO*ETA(IETA)-D_ONE)
         PHI(IETA) = ALPHA + D_TWO*DSQRT(TAU)*E(IETA) 
         DO I = 1,3
            CV(I,IETA) = VEL(I) + 
     $           (DT/SIGMA2)
     $           *(MFMEANGRAD(I)*ALPHA_M*PHI(IETA)
     $           -(D_ONE/(D_TWO*TAU))*(TAU_M*MFMEANGRAD(I)
     $           + TAU_V*MFVARGRAD(I))
     $           * (D_ONE + ALPHA*PHI(IETA) 
     $           - (PHI(IETA)**D_TWO)/SIGMA2))
         ENDDO
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF: CONDITIONAL VELOCITY:'
         WRITE(*,50)'============================================='
         WRITE(*,100) 'ALPHA   =',ALPHA
         WRITE(*,100) 'ALPHA_M =',ALPHA_M
         WRITE(*,150) 'TAU     =',TAU, 
     $        '-> MOD. OF ABS. ERR. OF INTEGRATION =',INTEG_ABSERR
         WRITE(*,100) 'SIGMA^2 =',SIGMA2
         WRITE(*,100) 'SIGMA^2 =',SIGMA2
         WRITE(*,150) 'TAU_M   =',TAU_M,  '-> APPROX. ERR. =',ERR1
         WRITE(*,150) 'TAU_V   =',TAU_V,  '-> APPROX. ERR. =',ERR2
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

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      SUBROUTINE PMF_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [HOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE PMF-PDF
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
      INTEGER IETA,NETA,I
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,ALPHA,TAU,SIGMA2
      DOUBLE PRECISION ETA(NETA),CSDRH(NETA),E(NETA),PHI(NETA)  
      DOUBLE PRECISION ERFINV
      EXTERNAL ERFINV
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      CALL CHECKPARMS(MFMEAN,MFVAR)

      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)

      CALL FIND_TAU(MFMEAN,MFVAR,TAU)

      SIGMA2 = D_ONE - D_TWO*TAU

C     BOUNDARY VALUES ARE KNOWN
      CSDRH(1) = D_ZERO
      CSDRH(NETA) = D_ZERO   
C     COMPUTE CSDR AT INTERNAL GRID POINTS
      DO IETA = 2,NETA-1
         E(IETA) = ERFINV(D_TWO*ETA(IETA)-D_ONE)
         CSDRH(IETA) = CHI * DSQRT((D_ONE-TAU)/TAU)
     $        * DEXP(-2*(E(IETA)**D_TWO)
     $        + (ALPHA**D_TWO)/(SIGMA2+D_ONE))
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF: COND. SCAL. DISS. RATE (HOMOGENEOUS):'
         WRITE(*,50)'============================================='
         WRITE(*,100) 'ALPHA   =',ALPHA
         WRITE(*,150) 'TAU     =',TAU, 
     $        '-> MOD. OF ABS. ERR. OF INTEGRATION =',INTEG_ABSERR
         WRITE(*,100) 'SIGMA^2 =',SIGMA2
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

      SUBROUTINE PMF_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,CHI,CSDRI)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [INHOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE PMF-PDF
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
      INTEGER IETA,NETA,I
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,DT,ALPHA,TAU,SIGMA2,
     $     ALPHA_M
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3)
      DOUBLE PRECISION ETA(NETA),CSDRH(NETA),CSDRI(NETA),
     $     E(NETA),PHI(NETA),A(NETA),B(NETA),C(NETA)
      DOUBLE PRECISION TAU_M,TAU_V,TAU_MM,TAU_VV
      DOUBLE PRECISION PROD1,PROD2,PROD3
      DOUBLE PRECISION T1,T2,T3,T4,T5,T6,FACT
      DOUBLE PRECISION H,IM,ERR1,ERR2,ERR3,ERR4,ERR_TEMP,INT_FACT
      DOUBLE PRECISION ERR_OLD, H_NEW,H_MEAN,H_VAR
      DOUBLE PRECISION ERFINV,DFRIDR,D2FRIDR,
     $     TAU_FIXEDMFVAR,TAU_FIXEDMFMEAN
      EXTERNAL ERFINV,DFRIDR,D2FRIDR,
     $     TAU_FIXEDMFVAR, TAU_FIXEDMFMEAN  
      DOUBLE PRECISION MFMEANPASS_DERIVS,MFVARPASS_DERIVS,ALPHAPASS
      COMMON/COMVARSMEAN/MFMEANPASS_DERIVS
      COMMON/COMVARSVAR/MFVARPASS_DERIVS
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION DIFF_MAXERR,DIFF_FACT
      COMMON/DIFFERENTIATORVARS/DIFF_MAXERR,DIFF_FACT


      CALL CHECKPARMS(MFMEAN,MFVAR)

      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)

      ALPHA_M = -DSQRT(D_TWO*PI)*DEXP((ALPHA**D_TWO)/D_TWO)

      CALL FIND_TAU(MFMEAN,MFVAR,TAU)

      SIGMA2 = D_ONE - D_TWO*TAU

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      IM = MFVAR + MFMEAN**D_TWO
      H  = DMIN1(MFVAR, MFMEAN - IM)
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     + DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
     $     - DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
C      H  = DMIN1(1.0D-3, H * 0.2D0)
      
      H_MEAN = 2.0D0*H          !H
      H_VAR  = DMIN1(1.0D-4, H*0.05D0)
      !H_VAR  = DMIN1(1.0D-4, H*2.0D0)

      ERR1 = HUGE(1.0D0)
      H_NEW = H_MEAN
      MFVARPASS_DERIVS = MFVAR
      DO WHILE(ERR1.GT.DIFF_MAXERR)
         ERR_OLD = ERR1
         TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
         IF(ERR1.GT.ERR_OLD) THEN
            H_NEW = H_NEW*DIFF_FACT
            TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
            EXIT
         ELSE
            H_NEW = H_NEW/DIFF_FACT
         ENDIF
      ENDDO

      ERR2 = HUGE(1.0D0)
      H_NEW = H_VAR
      MFMEANPASS_DERIVS = MFMEAN
      DO WHILE(ERR2.GT.DIFF_MAXERR)
         ERR_OLD = ERR2
         TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
         IF(ERR2.GT.ERR_OLD) THEN
            H_NEW = H_NEW*DIFF_FACT
            TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
            EXIT
         ELSE
            H_NEW = H_NEW/DIFF_FACT
         ENDIF
      ENDDO

      ERR3 = HUGE(1.0D0)  
      H_NEW = H_MEAN
      MFVARPASS_DERIVS = MFVAR
      DO WHILE(ERR3.GT.DIFF_MAXERR)
         ERR_OLD = ERR3
         TAU_MM  = D2FRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR3)
         IF(ERR3.GT.ERR_OLD) THEN
            H_NEW = H_NEW*DIFF_FACT
            TAU_MM  = D2FRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR3)
            EXIT
         ELSE
            H_NEW = H_NEW/DIFF_FACT
         ENDIF
      ENDDO

      ERR4 = HUGE(1.0D0)
      H_NEW = H_VAR
      MFMEANPASS_DERIVS = MFMEAN
      DO WHILE(ERR4.GT.DIFF_MAXERR)
         ERR_OLD = ERR4
         TAU_VV  = D2FRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR4)
         IF(ERR4.GT.ERR_OLD) THEN
            H_NEW = H_NEW*DIFF_FACT
            TAU_VV  = D2FRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR4)
            EXIT
         ELSE
            H _NEW= H_NEW/DIFF_FACT
         ENDIF
      ENDDO
      
c$$$C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
c$$$      IM = MFVAR + MFMEAN**D_TWO
c$$$      H  = DMIN1(MFVAR, MFMEAN - IM)
c$$$      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
c$$$     $     + DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
c$$$      H  = DMIN1(H, DABS(((D_ONE - D_TWO*MFMEAN) 
c$$$     $     - DSQRT(D_ONE - 4.0D0*MFVAR))/D_TWO))
c$$$C      H  = DMIN1(1.0D-3, H * 0.2D0)
c$$$
c$$$      ERR1 = HUGE(1.0D0)
c$$$      H_NEW = H
c$$$      MFVARPASS_DERIVS = MFVAR
c$$$      DO WHILE(ERR1.GT.DIFF_MAXERR)
c$$$         ERR_OLD = ERR1
c$$$         TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
c$$$         IF(ERR1.GT.ERR_OLD) THEN
c$$$            H_NEW = H_NEW*DIFF_FACT
c$$$            TAU_M  = DFRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR1)
c$$$            EXIT
c$$$         ELSE
c$$$            H_NEW = H_NEW/DIFF_FACT
c$$$         ENDIF
c$$$      ENDDO
c$$$
c$$$      ERR2 = HUGE(1.0D0)
c$$$      H_NEW = H
c$$$      MFMEANPASS_DERIVS = MFMEAN
c$$$      DO WHILE(ERR2.GT.DIFF_MAXERR)
c$$$         ERR_OLD = ERR2
c$$$         TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
c$$$         IF(ERR2.GT.ERR_OLD) THEN
c$$$            H_NEW = H_NEW*DIFF_FACT
c$$$            TAU_V  = DFRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR2)
c$$$            EXIT
c$$$         ELSE
c$$$            H_NEW = H_NEW/DIFF_FACT
c$$$         ENDIF
c$$$      ENDDO
c$$$
c$$$      ERR3 = HUGE(1.0D0)  
c$$$      H_NEW = H
c$$$      MFVARPASS_DERIVS = MFVAR
c$$$      DO WHILE(ERR3.GT.DIFF_MAXERR)
c$$$         ERR_OLD = ERR3
c$$$         TAU_MM  = D2FRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR3)
c$$$         IF(ERR3.GT.ERR_OLD) THEN
c$$$            H_NEW = H_NEW*DIFF_FACT
c$$$            TAU_MM  = D2FRIDR(TAU_FIXEDMFVAR,MFMEAN,H_NEW,ERR3)
c$$$            EXIT
c$$$         ELSE
c$$$            H_NEW = H_NEW/DIFF_FACT
c$$$         ENDIF
c$$$      ENDDO
c$$$
c$$$      ERR4 = HUGE(1.0D0)
c$$$      H_NEW = H
c$$$      MFMEANPASS_DERIVS = MFMEAN
c$$$      DO WHILE(ERR4.GT.DIFF_MAXERR)
c$$$         ERR_OLD = ERR4
c$$$         TAU_VV  = D2FRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR4)
c$$$         IF(ERR4.GT.ERR_OLD) THEN
c$$$            H_NEW = H_NEW*DIFF_FACT
c$$$            TAU_VV  = D2FRIDR(TAU_FIXEDMFMEAN,MFVAR,H_NEW,ERR4)
c$$$            EXIT
c$$$         ELSE
c$$$            H _NEW= H_NEW/DIFF_FACT
c$$$         ENDIF
c$$$      ENDDO

C     BOUNDARY VALUES ARE KNOWN
      CSDRI(1)    = D_ZERO
      CSDRI(NETA) = D_ZERO
      CSDRH(1)    = D_ZERO
      CSDRH(NETA) = D_ZERO
C     COMPUTE CSDR AT INTERNAL GRID POINTS
      DO IETA = 2,NETA-1
         E(IETA) = ERFINV(D_TWO*ETA(IETA)-D_ONE)

         PHI(IETA) = ALPHA + D_TWO*DSQRT(TAU)*E(IETA) 

         A(IETA) = ALPHA/(D_ONE-TAU) - PHI(IETA)/SIGMA2

         B(IETA) = (ALPHA**D_TWO)/(D_TWO*((D_ONE-TAU)**D_TWO))
     $        - (D_ONE/SIGMA2)*(PHI(IETA)*E(IETA)/DSQRT(TAU)
     $        + (PHI(IETA)**D_TWO)/SIGMA2
     $        - D_ONE/(D_ONE+SIGMA2))

         C(IETA) = DSQRT(D_TWO/PI) * (ALPHA_M
     $        + TAU_M*PHI(IETA)/(TAU*SIGMA2))
     $        * DEXP(-D_TWO*(E(IETA)**D_TWO) + (ALPHA**D_TWO)/D_TWO)

         CSDRH(IETA) = CHI * DSQRT((D_ONE-TAU)/TAU)
     $        * DEXP(-2*(E(IETA)**D_TWO)
     $        + (ALPHA**D_TWO)/(SIGMA2+D_ONE))

         PROD1 = D_ZERO ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. MEAN
         PROD2 = D_ZERO ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. VARIANCE
         PROD3 = D_ZERO ! SCALAR PRODUCT OF MIX FRAC. MEAN AND MIX FRAC VARIANCE GRADS
         DO I = 1,3
            PROD1 = PROD1 + MFMEANGRAD(I)**D_TWO 
            PROD2 = PROD2 + MFVARGRAD(I)**D_TWO
            PROD3 = PROD3 + MFMEANGRAD(I)*MFVARGRAD(I)
         ENDDO

         T1 = PROD1*(D_TWO + ((TAU_M**D_TWO)*TAU_VV/(TAU_V**D_THREE))
     $        - TAU_MM/TAU_V)
C         T1 = PROD1*(D_ZERO + ((TAU_M**D_TWO)*TAU_VV/(TAU_V**D_THREE) )
C     $        - TAU_MM/TAU_V)   ! FOR INTEGER MOMENTS
         T2 = D_TWO*PROD3*ALPHA_M*A(IETA)
         T3 = PROD2*TAU_V
         T4 = D_TWO*PROD3*TAU_M
         T5 = PROD1*(TAU_M**D_TWO)/TAU_V
         T6 = PROD1*C(IETA)
         FACT = D_TWO

         CSDRI(IETA) = (CHI - FACT*DT*(T1-T2-(T3+T4+T5)*B(IETA)))
     $        *(CSDRH(IETA)/CHI) - FACT*DT*T6

         !IF(CSDRI(IETA).LT.0.0D0) CSDRI(IETA) = D_ZERO

      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF: COND. SCAL. DISS. RATE (INHOMOGENEOUS):'
         WRITE(*,50)'============================================='
         WRITE(*,100) 'ALPHA   =',ALPHA
         WRITE(*,100) 'ALPHA_M =',ALPHA_M
         WRITE(*,150) 'TAU     =',TAU, 
     $        '-> MOD. OF ABS. ERR. OF INTEGRATION =',INTEG_ABSERR
         WRITE(*,100) 'SIGMA^2 =',SIGMA2
         WRITE(*,150) 'TAU_M   =',TAU_M,  '-> APPROX. ERR. =',ERR1
         WRITE(*,150) 'TAU_V   =',TAU_V,  '-> APPROX. ERR. =',ERR2
         WRITE(*,150) 'TAU_MM  =',TAU_MM, '-> APPROX. ERR. =',ERR3
         WRITE(*,150) 'TAU_VV  =',TAU_VV, '-> APPROX. ERR. =',ERR4
         WRITE(*,50)'============================================='
         WRITE(*,600)'INDEX','ETA','CSDRH','CSDRI'
         DO IETA = 1,NETA
            WRITE(*,700) IETA,ETA(IETA),CSDRH(IETA),CSDRI(IETA)
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
C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================

      DOUBLE PRECISION FUNCTION TAU_FIXEDMFVAR(MFMEAN)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST- AND 
C              SECOND-ORDER PARTIAL DERIVATIVES OF TAU WITH RESPECT TO THE 
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
      DOUBLE PRECISION MFMEAN,TAU
      DOUBLE PRECISION MFVARPASS_DERIVS
      COMMON/COMVARSVAR/MFVARPASS_DERIVS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      CALL FIND_TAU(MFMEAN,MFVARPASS_DERIVS,TAU)
      TAU_FIXEDMFVAR = TAU

      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION TAU_FIXEDMFMEAN(MFVAR)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST- AND 
C              SECOND-ORDER PARTIAL DERIVATIVES OF TAU WITH RESPECT TO THE 
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
      DOUBLE PRECISION MFVAR,TAU
      DOUBLE PRECISION MFMEANPASS_DERIVS
      COMMON/COMVARSMEAN/MFMEANPASS_DERIVS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      CALL FIND_TAU(MFMEANPASS_DERIVS,MFVAR,TAU)
      TAU_FIXEDMFMEAN = TAU

      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE FIND_TAU(MFMEAN,MFVAR,TAU)
C================================================================================== 
C     PURPOSE: COMPUTES THE VALUE OF THE PARAMETER TAU
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION
C
C     OUTPUT:
C
C            TAU       THE PARAMETER TAU              DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      INCLUDE 'formats.h'
      DOUBLE PRECISION TAU_FEXT,ERFINV,ZRIDDR,ZBRENT
      EXTERNAL TAU_FEXT,ERFINV,ZRIDDR,ZBRENT
      DOUBLE PRECISION TAU,MFMEAN,MFVAR
      DOUBLE PRECISION TAU_MIN, TAU_MAX
      INTEGER IFLAG
      DOUBLE PRECISION MFMEANPASS,MFVARPASS
      COMMON/COMVARSMEANVAR/MFMEANPASS,MFVARPASS
      DOUBLE PRECISION ALPHAPASS
      COMMON/COMVARSALPHA/ALPHAPASS

      INTEGER ROOTF_METH
      COMMON/ROOTFINDERVARS1/ROOTF_METH
      DOUBLE PRECISION ROOTF_RE,ROOTF_AE
      COMMON/ROOTFINDERVARS2/ROOTF_RE,ROOTF_AE
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      INTEGER N,NBMAX
      PARAMETER(N = 10,NBMAX = 1)
      DOUBLE PRECISION TAU1(NBMAX),TAU2(NBMAX)
      INTEGER I,NB
      DOUBLE PRECISION ROOTS,TAUR
      LOGICAL SUCCESS
      DOUBLE PRECISION TAU_MIN_NEW,TAU_MAX_NEW

      MFMEANPASS = MFMEAN
      MFVARPASS  = MFVAR
      ALPHAPASS = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)
      
C     TAU VARIES BETWEEN 0 AND 0.5
      TAU_MIN = D_ZERO
      TAU_MAX = D_HALF

      IF(ROOTF_METH.EQ.1) THEN
         CALL ZBRAK(TAU_FEXT,TAU_MIN,TAU_MAX,N,TAU1,TAU2,NB)
         IF(NB.GT.NBMAX) THEN
            WRITE(*,50)'============================================='
            WRITE(*,50)'MORE THAN ONE ROOT WAS DETECTED WHILE'
            WRITE(*,50)'COMPUTING TAU'
            WRITE(*,50)'COMPUTATIONS ABORTED.'
            WRITE(*,50)'============================================='
            STOP
         ENDIF
         TAU = ZRIDDR(TAU_FEXT,TAU1(NB),TAU2(NB),ROOTF_AE)
      ELSEIF(ROOTF_METH.EQ.2) THEN
         CALL ZBRAK(TAU_FEXT,TAU_MIN,TAU_MAX,N,TAU1,TAU2,NB)
         IF(NB.GT.NBMAX) THEN
            WRITE(*,50)'============================================='
            WRITE(*,50)'MORE THAN ONE ROOT WAS DETECTED WHILE'
            WRITE(*,50)'COMPUTING TAU.'
            WRITE(*,50)'COMPUTATIONS ABORTED.'
            WRITE(*,50)'============================================='
            STOP
         ENDIF
         TAU = ZBRENT(TAU_FEXT,TAU1(NB),TAU2(NB),ROOTF_AE)
      ELSEIF(ROOTF_METH.EQ.3) THEN   
         CALL ZEROIN(TAU_FEXT,TAU_MIN,TAU_MAX,ROOTF_RE,ROOTF_AE,IFLAG)
         TAU = TAU_MIN 
      ELSE
         WRITE(*,50)'============================================='
         WRITE(*,50)'WRONG INPUT:'
         WRITE(*,50)'INVALID ROOTFINDING METHOD (CHECK ROOTF_METH)'
         WRITE(*,50)'COMPUTATIONS ABORTED.'
         WRITE(*,50)'============================================='
         STOP
      ENDIF
            
      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION TAU_FEXT(TAU)
C================================================================================== 
C     PURPOSE: COMPUTES THE RIGHT-HAND SIDE OF EQ.(19) IN [1]
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            TAU       THE PARAMETER TAU              DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM
      DOUBLE PRECISION TAU,INTEGRAL
      DOUBLE PRECISION ALPHAPASS,MFMEANPASS,MFVARPASS,TAUPASS
      COMMON/COMVARSMEANVAR/MFMEANPASS,MFVARPASS
      COMMON/COMVARSALPHA/ALPHAPASS
      COMMON/COMVARSTAU/TAUPASS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
    
      DOUBLE PRECISION INTEGRAND
      EXTERNAL INTEGRAND
      DOUBLE PRECISION BOUND
      INTEGER INF,NEVAL,IER,LIMIT,LENW,LAST
      DOUBLE PRECISION EPSABS,EPSREL,RESULT
      INTEGER IWORK(INTEG_LIM)
      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      
      INF     = I_TWO
      LENW    = I_FOUR*INTEG_LIM
      EPSABS  = INTEG_EPSABS
      EPSREL  = INTEG_EPSREL
      TAUPASS = TAU
      INTEG_ABSERR = D_ZERO
      CALL DQAGI(INTEGRAND,BOUND,INF,EPSABS,EPSREL,RESULT,INTEG_ABSERR,
     $     NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
      INTEGRAL = RESULT
      TAU_FEXT = MFMEANPASS**D_TWO + MFVARPASS - INTEGRAL

      RETURN 
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION FUNCTION INTEGRAND(PHI)
C================================================================================== 
C     PURPOSE: COMPUTES THE INTEGRAND OF THE SECOND TERM ON THE RIGHT-HAND SIDE 
C              OF EQ.(19) IN [1]
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            PHI       SAMPLE SPACE VARIABLE OF       DOUBLE PRECISION
C                      THE REFERENCE FIELD PSI            
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION PHI,X,R,SIGMA2
      DOUBLE PRECISION ERF
      EXTERNAL ERF
      DOUBLE PRECISION ALPHAPASS,TAUPASS
      COMMON/COMVARSALPHA/ALPHAPASS
      COMMON/COMVARSTAU/TAUPASS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      SIGMA2 = D_ONE-D_TWO*TAUPASS
      R = (D_ONE/DSQRT(D_TWO*PI*SIGMA2))
     $     * DEXP(-(PHI**D_TWO)/(D_TWO*SIGMA2))
      X = D_HALF*(D_ONE + ERF((PHI-ALPHAPASS)/(D_TWO*DSQRT(TAUPASS))))
      INTEGRAND = (X**D_TWO)*R

      RETURN 
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 
