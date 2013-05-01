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
      PROGRAM PMF_DRIVER
C================================================================================== 
C     PURPOSE: DRIVER PROGRAM FOR PMFfortPack
C              ALL VARIABLES ARE DESCRIBED BELOW
C              ALL SUBROUTINES ARE DESCRIBED IN PMF.f
C==================================================================================
      IMPLICIT NONE
      INCLUDE 'formats.h'
      
C     INPUT VARIABLES (SET BELOW)
      INTEGER 
     $     IETA,          ! MIXTURE FRACTION SPACE CONTER
     $     NETA           ! NUMBER OF GRID POINTS IN MIXTURE FRACTION SPACE
      PARAMETER(NETA = 300)
C
      DOUBLE PRECISION
     $     MFMEAN,        ! MIXTURE FRACTION MEAN
     $     MFVAR,         ! MIXTURE FRACTION VARIANCE
     $     MFMEANGRAD(3), ! GRADIENT OF THE MIXTURE FRACTION MEAN (COMPONENTS: 1=X, 2=Y, 3=Z)
     $     MFVARGRAD(3),  ! GRADIENT OF THE MIXTURE FRACTION VARIANCE (COMPONENTS: 1=X, 2=Y, 3=Z)
     $     DT,            ! TURBULENT DIFFUSIVITY
     $     VEL(3),        ! VELOCITY VECTOR (COMPONENTS: 1=X, 2=Y, 3=Z)
     $     CHI            ! MEAN SCALAR DISSIPATION RATE
C
      DOUBLE PRECISION ETA(NETA)! MIXTURE FRACTION

C     PMF ARRAYS
      DOUBLE PRECISION
     $     PDF_PMF(NETA),     ! PROBABILITY DENSITY FUNCTION
     $     CV_PMF(3,NETA),    ! CONDITIONAL VELOCITY
     $     CSDRH_PMF(NETA),   ! CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS VERSION]]
     $     CSDRI_PMF(NETA)    ! CONDITIONAL SCALAR DISSIPATION RATE [[INHOMOGENEOUS VERSION]]

C     BETA ARRAYS
      DOUBLE PRECISION 
     $     PDF_BETA(NETA),    ! PROBABILITY DENSITY FUNCTION
     $     CV_BETA(3,NETA),   ! CONDITIONAL VELOCITY
     $     CSDRH_BETA(NETA),  ! CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS VERSION]]
     $     CSDRI_BETA(NETA)   ! CONDITIONAL SCALAR DISSIPATION RATE [[INHOMOGENEOUS VERSION]]

      DOUBLE PRECISION  INTEGRAND(NETA),INTEGRAL, RELERR ! A UTILITY ARRAY AND VARIABLES
      DOUBLE PRECISION TSTART,TFINISH,CPUT_PMF,CPUT_BETA

      DOUBLE PRECISION TRAP
      EXTERNAL TRAP

C     SET NEEDED VARIABLES
C
C     VALUES TAKEN FROM THE DEMO IN REF [2]. 
C     I THINK THE MOMENTS ARE INTEGER. CSDR CODE WRITTEN FOR CENTRAL MOMENTS.
C     MODIFICATION IS EASY.
C
c$$$      MFMEAN = 0.70837758D0
c$$$      MFVAR = 0.194360748D0*MFMEAN*(1.0D0-MFMEAN)
c$$$      MFMEANGRAD(1) = -0.3496497D0
c$$$      MFMEANGRAD(2) = 5.50D0
c$$$      MFMEANGRAD(3) = 0.0D0
c$$$      MFVARGRAD(1) = -0.4294126D0
c$$$      MFVARGRAD(2) = 10.5D0
c$$$      MFVARGRAD(3) = 0.0D0
c$$$      DT = 1.0D0
c$$$      VEL(1) = 30.25D0 !SET TO ZERO IN THE DEMO IN REF [2]
c$$$      VEL(2) = 2.60D0  !SET TO ZERO IN THE DEMO IN REF [2]
c$$$      VEL(3) = 0.0D0
c$$$      CHI = 1.003770D0     

C     VALUES TAKEN FROM THE CMC CALCULATIONS OF A LIFTED H2 FLAME
C
      MFMEAN = 0.3520381153D0
      MFVAR = 0.01624656655D0
      MFMEANGRAD(1) = 3.943452664847875D0
      MFMEANGRAD(2) = -1.198706210914406D+2
      MFMEANGRAD(3) = 0.0D0
      MFVARGRAD(1) = -0.051799869317208D0
      MFVARGRAD(2) = -2.482959669128132D0
      MFVARGRAD(3) = 0.0D0
      DT = 0.002229837701000 
      VEL(1) = 38.651126860D0
      VEL(2) = 0.5267138481D0
      VEL(3) = 0.0D0
      CHI = 1.908199463D+2    
     
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                           *****SETUP*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

C     GENERATE UNIFORM MIXTURE FRACTION GRID.
C     OTHER NON-UNIFORM GRID GENERATION TECHNIQUES ARE AVAILABLE IN GRID.f
      CALL UNIFORMGRID(0.0D0,1.0D0,NETA,ETA)

C     INITIALISE
      CALL INITIALISE

C     CHECK THE VALIDITY OF THE INPUT PARAMETERS (MIX. FRAC. MEAN NAD VARIANCE)
      CALL CHECKPARMS(MFMEAN,MFVAR)

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                       *****PMF APPROACH*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      CALL CPU_TIME(TSTART)

C     COMPUTE THE PMF-PDF
      CALL PMF_PDF(ETA,NETA,MFMEAN,MFVAR,PDF_PMF)

C     COMPUTE THE PMF-CV
      CALL PMF_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,DT,
     $     VEL,CV_PMF)

C     COMPUTE THE [[HOMOGENEOUS]] PMF-CSDR
      CALL PMF_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH_PMF)

C     COMPUTE THE [[INHOMOGENEOUS]] PMF-CSDR
      CALL PMF_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,CHI,CSDRI_PMF)

      CALL CPU_TIME(TFINISH)
      
      CPUT_PMF = TFINISH - TSTART

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                       *****BETA APPROACH*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      CALL CPU_TIME(TSTART)

C     COMPUTE THE BETA-PDF
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF_BETA)

C     COMPUTE THE BETA-CV
      CALL BETA_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,DT,
     $     VEL,CV_BETA)

C     COMPUTE THE [[HOMOGENEOUS]] BETA-CSDR
      CALL BETA_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH_BETA)

C     COMPUTE THE [[INHOMOGENEOUS]] BETA-CSDR
      CALL BETA_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,CHI,CSDRI_BETA)

      CALL CPU_TIME(TFINISH) 

      CPUT_BETA = TFINISH -TSTART


C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                       *****VALIDATION*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
      WRITE(*,*)' '
      WRITE(*,*)'============================================='
      WRITE(*,*)' '
      WRITE(*,*)'COMPARE THE MEAN VALUE OF CHI TO THE INTEGRAL'  
      WRITE(*,*)'1'
      WRITE(*,*)'| <X|ETA>P(ETA)dETA'
      WRITE(*,*)'0'
      WRITE(*,*)'NOTE: THE ERRORS CALCULATED BELOW CAN BE LARGE'
      WRITE(*,*)'      WHEN A COARSE ETA GRID IS EMPLOYED'
      WRITE(*,*)'      (SIMPLE TRAPEZOIDAL INTEGRATION IS USED)'
      WRITE(*,*)' '
      WRITE(*,'(A,F6.4)')'CHI MEAN          =',CHI
C
      DO IETA = 1,NETA
         INTEGRAND(IETA) = CSDRH_PMF(IETA)*PDF_PMF(IETA)
      ENDDO
      INTEGRAL = TRAP(INTEGRAND,ETA,NETA)
      RELERR   = 100.0D0*(CHI-INTEGRAL)/CHI
      WRITE(*,'(A,F6.4,A,F6.4)')'PMF_CHIPDF_INT_H  =', INTEGRAL,
     $     ' | REL ERR(%) =',RELERR
C
      DO IETA = 1,NETA
         INTEGRAND(IETA) = CSDRI_PMF(IETA)*PDF_PMF(IETA)
      ENDDO
      INTEGRAL = TRAP(INTEGRAND,ETA,NETA)
      RELERR   = 100.0D0*(CHI-INTEGRAL)/CHI
      WRITE(*,'(A,F6.4,A,F6.4)')'PMF_CHIPDF_INT_I  =', INTEGRAL,
     $     ' | REL ERR(%) =',RELERR
C
      DO IETA = 1,NETA
         INTEGRAND(IETA) = CSDRH_BETA(IETA)*PDF_BETA(IETA)
      ENDDO
      INTEGRAL = TRAP(INTEGRAND,ETA,NETA)
      RELERR   = 100.0D0*(CHI-INTEGRAL)/CHI
      WRITE(*,'(A,F6.4,A,F6.4)')'BETA_CHIPDF_INT_H =', INTEGRAL,
     $     ' | REL ERR(%) =',RELERR
C
      DO IETA = 1,NETA
         INTEGRAND(IETA) = CSDRI_BETA(IETA)*PDF_BETA(IETA)
      ENDDO
      INTEGRAL = TRAP(INTEGRAND,ETA,NETA)
      RELERR   = 100.0D0*(CHI-INTEGRAL)/CHI
      WRITE(*,'(A,F6.4,A,F6.4)')'BETA_CHIPDF_INT_I =', INTEGRAL,
     $     ' | REL ERR(%) =',RELERR
C
      WRITE(*,*)' '
      WRITE(*,*)'============================================='
      WRITE(*,*)' '
      WRITE(*,*)'CPU TIME STATISTICS:'
      WRITE(*,'(1X,A,F6.4,A)')'PMF  CPU TIME              = ', 
     $     CPUT_PMF,' s'
      WRITE(*,'(1X,A,F6.4,A)')'BETA CPU TIME              = ', 
     $     CPUT_BETA,' s'
      WRITE(*,'(1X,A,F6.4,A)')'BETA-TO-PMF CPU TIME RATIO = ', 
     $     CPUT_BETA/CPUT_PMF
      WRITE(*,*)' '
      WRITE(*,*)'============================================='
C     WRITE RESULTS TO FILES
      WRITE(*,50) ' ================================'
      WRITE(*,50) '|       REULTS WRITTEN TO:       |'
      WRITE(*,50) ' ================================'
      WRITE(*,50) '|   QUANTITY   |      FILE       |'
      WRITE(*,50) ' ================================'
      WRITE(*,50) '|    ETA       | ../eta.dat      |'
      WRITE(*,50) ' ================================'
      WRITE(*,50) '|             PMF                |'
      WRITE(*,50) ' ================================'
      WRITE(*,50) '|    PDF       | ../pmfPDF.dat   |'
      WRITE(*,50) '|    CV        | ../pmfCV.dat    |'
      WRITE(*,50) '|  CSDR (H)    | ../pmfCSDRH.dat |'
      WRITE(*,50) '|  CSDR (I)    | ../pmfCSDRI.dat |'
      WRITE(*,50) ' ================================'
      WRITE(*,50) '|            BETA                |'
      WRITE(*,50) ' ================================'
      WRITE(*,50) '|    PDF       | ../betaPDF.dat  |'
      WRITE(*,50) '|    CV        | ../betaCV.dat   |'
      WRITE(*,50) '|  CSDR (H)    | ../betaCSDRH.dat|'
      WRITE(*,50) '|  CSDR (I)    | ../betaCSDRI.dat|'
      WRITE(*,50) ' ================================'
      WRITE(*,50) 'THE RESULTS CAN BE PLOTTED BY RUNNING ../PMF.m'
      OPEN(1,FILE = 'eta.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(2,FILE = 'pmfPDF.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(3,FILE = 'pmfCV.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(4,FILE = 'pmfCSDRH.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(5,FILE = 'pmfCSDRI.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(6,FILE = 'betaPDF.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(7,FILE = 'betaCV.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(8,FILE = 'betaCSDRH.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      OPEN(9,FILE = 'betaCSDRI.dat',FORM = 'FORMATTED',STATUS='UNKNOWN')
      DO IETA = 1,NETA
         WRITE(1,110) ETA(IETA)
         WRITE(2,110) PDF_PMF(IETA)
         WRITE(3,120) CV_PMF(1,IETA),CV_PMF(2,IETA),CV_PMF(3,IETA)
         WRITE(4,110) CSDRH_PMF(IETA)
         WRITE(5,110) CSDRI_PMF(IETA)
         WRITE(6,110) PDF_BETA(IETA)
         WRITE(7,120) CV_BETA(1,IETA),CV_BETA(2,IETA),CV_BETA(3,IETA)
         WRITE(8,110) CSDRH_BETA(IETA)
         WRITE(9,110) CSDRI_BETA(IETA)
      ENDDO
      CLOSE(1)
      CLOSE(2)
      CLOSE(3)
      CLOSE(4)
      CLOSE(5)
      CLOSE(6)
      CLOSE(7)
      CLOSE(8)
      CLOSE(9)

      STOP
      END

      
