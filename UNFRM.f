C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE DELTA_PDF(ETA,NETA,MFMEAN,PDF)
C================================================================================== 
C     PURPOSE: COMPUTES THE DIRCA DELTA FUNCTION GIVEN THE MIXTURE FRACTION MEAN.
C              THE PDF IS ZERO AT EVERY GRID POINT IN MIXTURE FRACTION SPACE 
C              EXCEPT AT THE GRID POINT CLOSEST TO THE MEAN WHERE A TRIANGULAR 
C              FUNCTION IS USED TO REPRESENT THE DELTA FUNCTION. THE MAGNITUDE OF
C              THE PDF AT THIS GRID POINT IS SET SUCH THAT THE INTEGRAL OF THE 
C              RESULTING PDF OVER MIXTURE FRACTION SPACE IS UNITY.
C
C     NOTE:    THE MAGNITUDE OF THE DELTA FUNCTION DEPENDS ON THE RESOLUTION OF 
C              THE ETA GRID, HOWEVER, THIS WILL NOT AFFECT THE COMPUTATIONS OF 
C              THE UNCONDITIONAL AVERAGES (INTEGRATION OF THE CONDITIONAL AVERAGES 
C              MULTIPLIED BY DIRAC_PDF OVER ETA SPACE) SINCE THE SAME ETA GRID IS 
C              EMPLOYED IN THE CALCULATION OF BOTH.
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            ETA       MIXTURE FRACTION GRID          DOUBLE PRECISION (ARRAY)
C            NETA      SIZE OF ETA AND PDF            INTEGER
C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C
C     OUTPUT:
C
C            PDF       PROBABILITY DENSITY FUNCTION   DOUBLE PRECISION (ARRAY)
C================================================================================== 
      IMPLICIT NONE      
      INTEGER IETA,NETA,IMFMEAN
      DOUBLE PRECISION MFMEAN,ETA(NETA),PDF(NETA),PMFMEAN
      INTEGER IMFMEAN_MN,IMFMEAN_PL,ETA_PL,ETA_MN,DIFF_PL,DIFF_MN
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

C     FIND THE LOCATION OF MFMEAN IN ETA
      CALL LOCATE(ETA,NETA,MFMEAN,IMFMEAN)
      
C     CALCULATE THE MAGNITUDE OF THE DELTA FUNCTION. SPECIAL ATTENTION 
C     IS NEEDED IF THE DELTA FUNCTION IS LOCATED AT THE BOUNDARIES OF ETA
      IF(IMFMEAN.GT.NETA-1) THEN
         IMFMEAN = NETA
         PMFMEAN = 2.0D0/(ETA(IMFMEAN)-ETA(IMFMEAN-1))
      ELSEIF(IMFMEAN.LT.2) THEN
         IMFMEAN = 1
         PMFMEAN = D_TWO/(ETA(IMFMEAN+1)-ETA(IMFMEAN))
      ELSE
         PMFMEAN = D_TWO/(ETA(IMFMEAN+1)-ETA(IMFMEAN-1))
      ENDIF     

C     SET THE PDF TO ZERO OVER ETA SPACE
      DO IETA = 1,NETA
         PDF(IETA) = D_ZERO
      ENDDO 
C     OVERRIDE THE PDF AT THE GRID POINT WHERE THE DELTA FUNCTION IS LOCATED
      PDF(IMFMEAN) = PMFMEAN

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE UNFRM_CSDR(ETA,NETA,CHI,CSDR)
C================================================================================== 
C     PURPOSE: SETS A UNIFORM DISTRIBUTION FOR THE CONDITIONAL SCALAR DISSIPATION 
C              RATE IN MIXTURE FRACTION SPACE. THE CONDITIONAL SCALAR DISSIPATION 
C              RATE IS SET EQUAL TO ITS MEAN (FAVRE-AVERAGED) VALUE AT EVERY GRID 
C              POINT IN ETA. 
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            ETA       MIXTURE FRACTION GRID          DOUBLE PRECISION (ARRAY)
C            NETA      SIZE OF ETA AND PDF            INTEGER
C            CHI       MEAN SCALAR DISSIPATYION RATE  DOUBLE PRECISION 
C
C     OUTPUT:
C
C            CSDR      UNIFORMLY DISTRIBUTED          DOUBLE PRECISION (ARRAY)
C                      SCALAR DISSIPATION RATE
C================================================================================== 
      IMPLICIT NONE      
      INTEGER IETA,NETA
      DOUBLE PRECISION ETA(NETA),CHI,CSDR(NETA)
      
      DO IETA = 1,NETA
         CSDR(IETA) = CHI
      ENDDO
     
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE UNFRM_CV(ETA,NETA,VEL,CV)
C================================================================================== 
C     PURPOSE: SETS A UNIFORM DISTRIBUTION FOR THE CONDITIONAL VELOCITY IN MIXTURE 
C              FRACTION SPACE. THE CONDITIONAL VELOCITY COMPONENETS ARE SET EQUAL 
C              TO THEIR RESPECTIVE UNCONDITIONAL (FAVRE-AVERAGED) VALUES AT EVERY 
C              GRID POINT IN ETA. 
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            ETA       MIXTURE FRACTION GRID          DOUBLE PRECISION (ARRAY)
C            NETA      SIZE OF ETA AND PDF            INTEGER
C            VAR       UNCONDITIONAL VELOCITY         DOUBLE PRECISION 
C                      OR SCALAR DISSIPATYION RATE 
C
C     OUTPUT:
C
C            CONDVAR   UNIFORMLY DISTRIBUTED          DOUBLE PRECISION (ARRAY)
C                      VELOCITY OR SCALAR 
C                      DISSIPATION RATE
C================================================================================== 
      IMPLICIT NONE      
      INTEGER IETA,NETA,I
      DOUBLE PRECISION ETA(NETA),VEL(3),CV(3,NETA)
      
      DO I = 1,3
         DO IETA = 1,NETA
            CV(I,IETA) = VEL(I)
         ENDDO
      ENDDO
     
      RETURN
      END

