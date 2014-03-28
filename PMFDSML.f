C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      SUBROUTINE PMFDSML_PDF(ETA,NETA,ETAM,C,MFMEAN,MFVAR,PDF)
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
      DOUBLE PRECISION MFMEAN,MFVAR,TAU,SIGMA2,PSI(NETA),ETAM,C,C1,C2
      DOUBLE PRECISION ETA(NETA),PDF(NETA),E(NETA)
      DOUBLE PRECISION G,V,W
      DOUBLE PRECISION ETA_TEMP1(NETA-1),PDF_TEMP1(NETA-1)
      DOUBLE PRECISION PDF_AREA,PDF_LA,PDF_RA,PDF_SCALEFACT
      
      DOUBLE PRECISION ERFINV,TRAP,CDFNORMDIST
      EXTERNAL ERFINV,TRAP,CDFNORMDIST
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS

      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM

      
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE,
     $     D_FOUR,D_EIGHT
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS2/D_FOUR,D_EIGHT 
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION ALPHA1,ALPHA2
      DOUBLE PRECISION EXP_NUM,EXP_DEN_1,EXP_DEN_2
      DOUBLE PRECISION TERM1,TERM2,X
      DOUBLE PRECISION INVERSEX
      EXTERNAL INVERSEX

C      CALL CHECKPARMS(MFMEAN,MFVAR)

      C1 = C
      C2 = (MFMEAN+C1-D_ONE)/(ETAM-D_ONE)
      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
      CALL FIND_TAU_DSML(MFMEAN,MFVAR,ETAM,ALPHA1,ALPHA2,TAU) 
      SIGMA2 = D_ONE - D_TWO*TAU  

      DO IETA = 2,NETA-1

         PSI(IETA) = INVERSEX(ETA(IETA),ETAM,ALPHA1,ALPHA2,TAU)

         EXP_NUM = DEXP(-(PSI(IETA)**D_TWO)/(D_TWO*(D_ONE-D_TWO*TAU)))
         EXP_DEN_1 = DEXP(-((PSI(IETA)-ALPHA1)**D_TWO)/(D_FOUR*TAU))
         EXP_DEN_2 = DEXP(-((PSI(IETA)-ALPHA2)**D_TWO)/(D_FOUR*TAU))

         PDF(IETA) = DSQRT(D_TWO*TAU/SIGMA2)
     $        * EXP_NUM/(ETAM*EXP_DEN_1 + (D_ONE-ETAM)*EXP_DEN_2)

         IF(PROGRESS) CALL PROGRESSPC(IETA,NETA-1)

      ENDDO

C     DETERMINE THE SHAPE OF THE PDF USING THE FIRST AND LAST TWO INTERIOR POINTS OF THE
C     PDF ARRAY IN ORDER TO SET BOUNDARY VALUES. THE PRECEDURE MAKES USE OF THE FACT
C     THAT R(PHI) IS THE CUMULATIVE DISTRIBUTION FUNCTION OF P(ETA).
C     R(PHI) IS THE CUMULATIVE NORMAL DISTRIBUTION WITH ZERO MEAN AND SIGMA2 VARIANCE
C     AND PHI IS RELATED TO ETA BY: PHI = ALPHA + 2*DSQRT(TAU)*ERFINV(2*ETA-1). 
C     R(PHI) IS USED TO COMPUTE THE AREAS UNDER THE PDF EXTENDING FROM ETA(1) TO ETA(2)
C     AND FROM ETA(NETA-1) TO ETA(NETA). THESE AREARS ARE THEN USED TO SET THE BOUNDARY 
C     VALUES OF THE PDF, PDF(1) AND PDF(NETA).
C
C     USE R(PHI) TO FIND THE AREA BETWEEN THE LEFT BOUNADRY AND THE FIRST INTERIOR GRID POINT
         PDF_LA = CDFNORMDIST(PSI(2),D_ZERO,SIGMA2)
         PDF(1) = D_TWO*PDF_LA/(ETA(2)-ETA(1)) - PDF(2)
         PDF(1) = DMAX1(PDF(1),D_ZERO)
C     USE R(PHI) TO FIND THE AREA BETWEEN THE LEFT BOUNADRY AND THE LAST INTERIOR GRID POINT
         PDF_RA = D_ONE - CDFNORMDIST(PSI(NETA-1),D_ZERO,SIGMA2)
         PDF(NETA) = D_TWO*PDF_RA/(ETA(NETA)-ETA(NETA-1)) - PDF(NETA-1)
         PDF(NETA) = DMAX1(PDF(NETA),D_ZERO)
C
      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF-T: PROBABILITY DENSITY FUNCTION:'
         WRITE(*,50)'============================================='
         WRITE(*,100) 'ALPHA1  =',ALPHA1
         WRITE(*,100) 'ALPHA2  =',ALPHA2
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

      SUBROUTINE PMFDSML_CV(ETA,NETA,ETAM,C,MFMEAN,MFVAR,
     $     MFMEANGRAD,MFVARGRAD,DT,VEL,CV)
C========================================================================================== 
C     PURPOSE: COMPUTES THE CONDITIONAL VELOCITY USING THE PMF-PDF
C==========================================================================================
C            VARIABLE   DESCRIPTION                     DATA TYPE
C            --------   -----------                     --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION GRID           DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA                     INTEGER
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
      DOUBLE PRECISION ETAM,C,MFMEAN,MFVAR,DT
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3),VEL(3)
      DOUBLE PRECISION ETA(NETA),CV(3,NETA)
      DOUBLE PRECISION PDF(NETA)
      DOUBLE PRECISION DPDM,DPDV,ERR_DPDM,ERR_DPDV 

      LOGICAL VERBOSE,PROGRESS,VERBOSE_TEMP,PROGRESS_TEMP
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINFM,LIMSUPM
      DOUBLE PRECISION LIMINFV,LIMSUPV

      DOUBLE PRECISION H_NEW,ERR,ERR_OLD,H_MEAN,H_VAR 
      DOUBLE PRECISION PMFDSML_PDF_M,PMFDSML_PDF_V
      EXTERNAL PMFDSML_PDF_M,PMFDSML_PDF_V

      DOUBLE PRECISION PMFCV_ETA,PMFCV_ETAM,PMFCV_C,
     $     PMFCV_MEAN,PMFCV_VAR
      COMMON/PMFCVETA/PMFCV_ETA
      COMMON/PMFCVETAM/PMFCV_ETAM
      COMMON/PMFCVC/PMFCV_C
      COMMON/PMFCVMEAN/PMFCV_MEAN
      COMMON/PMFCVVAR/PMFCV_VAR

      INTEGER DIFF_NCHEBCOEFS_TEMP

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP  = VERBOSE
      PROGRESS_TEMP = PROGRESS
      VERBOSE  = .FALSE.
      PROGRESS = .FALSE.
      CALL PMFDSML_PDF(ETA,NETA,ETAM,C,MFMEAN,MFVAR,PDF)
      VERBOSE  = VERBOSE_TEMP
      PROGRESS = PROGRESS_TEMP

      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
      LIMINFV = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUPV = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      LIMINFM = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUPM = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)

      PMFCV_ETAM = ETAM
      PMFCV_C = C
      PMFCV_MEAN = MFMEAN
      PMFCV_VAR = MFVAR

      DO IETA = 2,NETA-1
C
         PMFCV_ETA = ETA(IETA)
C
C     ADJUST NUMBER OF COEFS IF NECESSARY
         IF(PDF(IETA).LE.DIFF_PDFMIN) THEN
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
         ELSE
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
         ENDIF        
C         WRITE(*,*)'PMF-T CV',IETA,PDF(IETA),DIFF_NCHEBCOEFS
C     
C     DPDV       
         CALL CHEBDERIV(LIMINFV,LIMSUPV,DIFF_NCHEBCOEFS_TEMP,
     $        PMFDSML_PDF_V,MFVAR,1,DPDV)
C     
C     DPDM
         CALL CHEBDERIV(LIMINFM,LIMSUPM,DIFF_NCHEBCOEFS_TEMP,
     $        PMFDSML_PDF_M,MFMEAN,1,DPDM)       
C
         DO I = 1,3
            CV(I,IETA) =  VEL(I) - (DT/PDF(IETA))
     $           *(DPDM*MFMEANGRAD(I) + DPDV*MFVARGRAD(I))
         ENDDO
C
         IF(PROGRESS) CALL PROGRESSPC(IETA,NETA-1)
C
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF-T: CONDITIONAL VELOCITY:'
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

      SUBROUTINE PMFDSML_CSDR_H(ETA,NETA,ETAM,C,MFMEAN,MFVAR,CHI,CSDRH)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [HOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE PMF-PDF
C==========================================================================================
C            VARIABLE   DESCRIPTION                     DATA TYPE
C            --------   -----------                     --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION GRID           DOUBLE PRECISION (ARRAY,SIZE=NETA)
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
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,TAU,SIGMA2
      DOUBLE PRECISION ETA(NETA),CSDRH(NETA),E,PSI,ETAM,C,C1,C2
      DOUBLE PRECISION ERFINV
      EXTERNAL ERFINV
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE,
     $     D_FOUR,D_EIGHT
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS2/D_FOUR,D_EIGHT 
      DOUBLE PRECISION ALPHA1,ALPHA2
      DOUBLE PRECISION EXP_NUM_1,EXP_NUM_2,EXP_DEN_1,
     $     EXP_DEN_2,EXP_DEN_3
      DOUBLE PRECISION TERM1,TERM2,X
      DOUBLE PRECISION INVERSEX
      EXTERNAL INVERSEX
      
C      CALL CHECKPARMS(MFMEAN,MFVAR)

      C1 = C
      C2 = (MFMEAN+C1-D_ONE)/(ETAM-D_ONE)
      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
      CALL FIND_TAU_DSML(MFMEAN,MFVAR,ETAM,ALPHA1,ALPHA2,TAU)
      SIGMA2 = D_ONE - D_TWO*TAU   

C     BOUNDARY VALUES ARE KNOWN
      CSDRH(1) = D_ZERO
      CSDRH(NETA) = D_ZERO   
C     COMPUTE CSDR AT INTERNAL GRID POINTS
      DO IETA = 2,NETA-1
C
         PSI = INVERSEX(ETA(IETA),ETAM,ALPHA1,ALPHA2,TAU)
C
         EXP_NUM_1 = DEXP(-((PSI-ALPHA1)**D_TWO)/(D_FOUR*TAU))
         EXP_NUM_2 = DEXP(-((PSI-ALPHA2)**D_TWO)/(D_FOUR*TAU))
         EXP_DEN_1 = DEXP((ALPHA1**D_TWO)/(D_TWO*(TAU-D_ONE)))
         EXP_DEN_2 = DEXP((ALPHA1**D_TWO + ALPHA2**D_TWO
     $        + D_TWO*ALPHA1*ALPHA2*(D_TWO*TAU-D_ONE))
     $        /(D_EIGHT*TAU*(TAU-D_ONE)))
         EXP_DEN_3 = DEXP((ALPHA2**D_TWO)/(D_TWO*(TAU-D_ONE)))
C
         CSDRH(IETA) = CHI * DSQRT((D_ONE-TAU)/TAU)
     $        * ((ETAM*EXP_NUM_1 + (D_ONE-ETAM)*EXP_NUM_2)**D_TWO)
     $        / ((ETAM**D_TWO)*EXP_DEN_1 
     $        + D_TWO*ETAM*(D_ONE-ETAM)*EXP_DEN_2
     $        + ((D_ONE-ETAM)**D_TWO)*EXP_DEN_3)
C
         IF(PROGRESS) CALL PROGRESSPC(IETA,NETA-1)
C
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF-T: COND. SCAL. DISS. RATE (HOMOGENEOUS):'
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

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      SUBROUTINE PMFDSML_CSDR_I(ETA,NETA,ETAM,C,MFMEAN,MFVAR,
     $     MFMEANGRAD,MFVARGRAD,DT,CHI,CSDRI)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [HOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE PMF-PDF
C==========================================================================================
C            VARIABLE   DESCRIPTION                     DATA TYPE
C            --------   -----------                     --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION GRID           DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA AND CSDRI           INTEGER
C            MFMEAN     MIXTURE FRACTION MEAN           DOUBLE PRECISION 
C            MFVAR      MIXTURE FRACTION VARIANCE       DOUBLE PRECISION
C            CHI        MEAN SCALAR DISSIPATION RATE    DOUBLE PRECISION 
C
C     OUTPUT:
C
C            CSDRI      COND. SCALAR DISSIPATION RATE   DOUBLE PRECISION (ARRAY, SIZE=NETA)
C========================================================================================== 
      IMPLICIT NONE
      INCLUDE 'formats.h'
      INTEGER IETA,NETA,I
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,TAU,SIGMA2
      DOUBLE PRECISION DT,MFMEANGRAD(3),MFVARGRAD(3)
      DOUBLE PRECISION ETA(NETA),CSDRH(NETA),CSDRI(NETA),PDF(NETA)
      DOUBLE PRECISION E,PSI,ETAM,C,C1,C2,PDFDED,TERM
      DOUBLE PRECISION II_MM,II_VV,II_MV,P
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      LOGICAL VERBOSE,PROGRESS,VERBOSE_TEMP,PROGRESS_TEMP
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE,
     $     D_FOUR,D_EIGHT
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS2/D_FOUR,D_EIGHT 
      DOUBLE PRECISION ALPHA1,ALPHA2
      DOUBLE PRECISION EXP_NUM_1,EXP_NUM_2,EXP_DEN_1,
     $     EXP_DEN_2,EXP_DEN_3
      DOUBLE PRECISION TERM1,TERM2,X
      DOUBLE PRECISION INVERSEX,ERFINV
      EXTERNAL INVERSEX,ERFINV

      DOUBLE PRECISION PRODM,PRODV,PRODMV
      DOUBLE PRECISION H,H_NEW,H_MEAN,H_VAR,H_MEAN_NEW,H_VAR_NEW
      DOUBLE PRECISION ERR,ERR_OLD
      DOUBLE PRECISION DIIDV(NETA)
      DOUBLE PRECISION D2IIDV2(NETA)
      DOUBLE PRECISION D2IIDM2(NETA)
      DOUBLE PRECISION D2IIDMDV(NETA)

      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINF,LIMSUP
      DOUBLE PRECISION LIMINFM,LIMSUPM
      DOUBLE PRECISION LIMINFV,LIMSUPV

      DOUBLE PRECISION IIPMFPDF_M,IIPMFPDF_V,IIPMFPDF_MV,IIPMFPDF_VM
      EXTERNAL IIPMFPDF_M,IIPMFPDF_V,IIPMFPDF_MV,IIPMFPDF_VM

      DOUBLE PRECISION PMFMFMEANPASS
      COMMON/MFMEANBLOK/PMFMFMEANPASS
      DOUBLE PRECISION PMFMFVARPASS
      COMMON/MFVARBLOK/PMFMFVARPASS
      DOUBLE PRECISION PMFCPASS
      COMMON/C1BLOK/PMFCPASS
      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
      DOUBLE PRECISION PMFPDFPASS
      COMMON/PMFPDFBLOK/PMFPDFPASS

      INTEGER DIFF_NCHEBCOEFS_TEMP
      DOUBLE PRECISION DIFF_MAXERR_TEMP

       
C      CALL CHECKPARMS(MFMEAN,MFVAR)

      VERBOSE_TEMP  = VERBOSE
      PROGRESS_TEMP = PROGRESS
      VERBOSE  = .FALSE.
      PROGRESS = .FALSE.
      CALL PMFDSML_PDF(ETA,NETA,ETAM,C,MFMEAN,MFVAR,PDF)
      VERBOSE  = VERBOSE_TEMP
      PROGRESS = PROGRESS_TEMP

C     PARTIAL DERIVATIVES OF II(ETA)

      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
      LIMINFV = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUPV = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      LIMINFM = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUPM = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)

      PMFETAMPASS = ETAM
      PMFCPASS = C

      DO IETA = 2,NETA-1
C
         PMFETAPASS = ETA(IETA)
C
C     ADJUST NUMBER OF COEFS IF NECESSARY
         IF(PDF(IETA).LE.DIFF_PDFMIN) THEN
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
         ELSE
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
         ENDIF
C         WRITE(*,*)'PMF-T CSDR-I',IETA,PDF(IETA),DIFF_NCHEBCOEFS_TEMP
C     
C     D2IIDV2
         PMFMFMEANPASS = MFMEAN        
         CALL CHEBDERIV(LIMINFV,LIMSUPV,DIFF_NCHEBCOEFS_TEMP,IIPMFPDF_V,
     $        MFVAR,2,D2IIDV2(IETA))
C     
C     D2IIDM2
         PMFMFVARPASS = MFVAR
         CALL CHEBDERIV(LIMINFM,LIMSUPM,DIFF_NCHEBCOEFS_TEMP,IIPMFPDF_M,
     $        MFMEAN,2,D2IIDM2(IETA))  
C     
C     D2IIDMDV
         PMFMFVARPASS = MFVAR
         PMFPDFPASS = PDF(IETA)
         CALL CHEBDERIV(LIMINFM,LIMSUPM,DIFF_NCHEBCOEFS_TEMP,
     $        IIPMFPDF_MV,MFMEAN,1,D2IIDMDV(IETA))
C     
         IF(PROGRESS) CALL PROGRESSPC(IETA,NETA-1)
C
      ENDDO    

C     COMPUTE THE SCALAR PRODUCTS
      PRODM  = D_ZERO            ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. MEAN
      PRODV  = D_ZERO            ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. VARIANCE
      PRODMV = D_ZERO            ! SCALAR PRODUCT OF MIX FRAC. MEAN AND MIX FRAC VARIANCE GRADS
      DO I = 1,3
         PRODM  = PRODM + MFMEANGRAD(I)**D_TWO 
         PRODV  = PRODV + MFVARGRAD(I)**D_TWO
         PRODMV = PRODMV + MFMEANGRAD(I)*MFVARGRAD(I)
      ENDDO

      C1 = C
      C2 = (MFMEAN+C1-D_ONE)/(ETAM-D_ONE)
      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
      CALL FIND_TAU_DSML(MFMEAN,MFVAR,ETAM,ALPHA1,ALPHA2,TAU)
      SIGMA2 = D_ONE - D_TWO*TAU   

C     BOUNDARY VALUES ARE KNOWN
      CSDRI(1) = D_ZERO
      CSDRI(NETA) = D_ZERO   
C     COMPUTE CSDR AT INTERNAL GRID POINTS
      DO IETA = 2,NETA-1
C       
         PSI = INVERSEX(ETA(IETA),ETAM,ALPHA1,ALPHA2,TAU)
C
         EXP_NUM_1 = DEXP(-((PSI-ALPHA1)**D_TWO)/(D_FOUR*TAU))
         EXP_NUM_2 = DEXP(-((PSI-ALPHA2)**D_TWO)/(D_FOUR*TAU))
         EXP_DEN_1 = DEXP((ALPHA1**D_TWO)/(D_TWO*(TAU-D_ONE)))
         EXP_DEN_2 = DEXP((ALPHA1**D_TWO + ALPHA2**D_TWO
     $        + D_TWO*ALPHA1*ALPHA2*(D_TWO*TAU-D_ONE))
     $        /(D_EIGHT*TAU*(TAU-D_ONE)))
         EXP_DEN_3 = DEXP((ALPHA2**D_TWO)/(D_TWO*(TAU-D_ONE)))
C
         CSDRI(IETA) = (CHI - D_TWO*DT*PRODM) * DSQRT((D_ONE-TAU)/TAU)
     $        * ((ETAM*EXP_NUM_1 + (D_ONE-ETAM)*EXP_NUM_2)**D_TWO)
     $        / ((ETAM**D_TWO)*EXP_DEN_1 
     $        + D_TWO*ETAM*(D_ONE-ETAM)*EXP_DEN_2
     $        + ((D_ONE-ETAM)**D_TWO)*EXP_DEN_3)
     $        + D_TWO*(DT/PDF(IETA))*(PRODM*D2IIDM2(IETA) 
     $        + PRODV*D2IIDV2(IETA) +  D_TWO*PRODMV*D2IIDMDV(IETA))
C        
      ENDDO

c$$$      OPEN(1,FILE='res/PMFIIDERIVS.dat',FORM='FORMATTED',
c$$$     $     STATUS='UNKNOWN')
c$$$      DO IETA = 1,NETA
c$$$         WRITE(1,120) D2IIDM2(IETA),D2IIDV2(IETA),D2IIDMDV(IETA)
c$$$      ENDDO
c$$$      CLOSE(1)

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF-T: COND. SCAL. DISS. RATE (HOMOGENEOUS):'
         WRITE(*,50)'============================================='
         WRITE(*,200)'INDEX','ETA','CSDRI'
         DO IETA = 1,NETA
            WRITE(*,300) IETA, ETA(IETA), CSDRI(IETA)
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

      DOUBLE PRECISION FUNCTION INVERSEX(ETA,ETAM,ALPHA1,ALPHA2,TAU)
      IMPLICIT NONE 
      DOUBLE PRECISION ETA,ETAM,ALPHA1,ALPHA2,TAU
      DOUBLE PRECISION XFUNC
      EXTERNAL XFUNC
      DOUBLE PRECISION LOW,HIGH
      LOGICAL SUCCES
      INTEGER IFLAG
      DOUBLE PRECISION ROOTF_RE,ROOTF_AE
      COMMON/ROOTFINDERVARS/ROOTF_RE,ROOTF_AE
      DOUBLE PRECISION ETAPASS,ETAMPASS,ALPHA1PASS,ALPHA2PASS,TAUPASS
      COMMON/COMVARSETA/ETAPASS
      COMMON/COMVARSETAM/ETAMPASS
      COMMON/COMVARSALPHA12/ALPHA1PASS,ALPHA2PASS
      COMMON/COMVARSTAU/TAUPASS

      ETAPASS = ETA
      ETAMPASS = ETAM
      ALPHA1PASS = ALPHA1
      ALPHA2PASS = ALPHA2
      TAUPASS = TAU

C      LOW =  -1.0D+40
C      HIGH = -LOW
      
      LOW =  -1.0D0 
      HIGH = -LOW 
      CALL BRACKET(XFUNC,LOW,HIGH,5.0D0)   
      CALL ZEROIN(XFUNC,LOW,HIGH,ROOTF_RE,ROOTF_AE,IFLAG)
      INVERSEX = LOW
      RETURN
      END      

      DOUBLE PRECISION FUNCTION XFUNC(PSI)
      IMPLICIT NONE
      DOUBLE PRECISION PSI,TERM1,TERM2
      DOUBLE PRECISION ERF
      EXTERNAL ERF
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION ETAPASS,ETAMPASS,ALPHA1PASS,ALPHA2PASS,TAUPASS
      COMMON/COMVARSETA/ETAPASS
      COMMON/COMVARSETAM/ETAMPASS
      COMMON/COMVARSALPHA12/ALPHA1PASS,ALPHA2PASS
      COMMON/COMVARSTAU/TAUPASS
      TERM1 = ETAMPASS*ERF((PSI-ALPHA1PASS)/(D_TWO*DSQRT(TAUPASS)))
      TERM2 = (D_ONE-ETAMPASS)*ERF((PSI-ALPHA2PASS)
     $     /(D_TWO*DSQRT(TAUPASS)))
      XFUNC = D_HALF*(D_ONE+TERM1+TERM2) - ETAPASS
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE FIND_TAU_DSML(MFMEAN,MFVAR,ETAM,ALPHA1,ALPHA2,TAU)
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
      DOUBLE PRECISION TAUFUNC_DSML,ERFINV,ZRIDDR,ZBRENT
      EXTERNAL TAUFUNC_DSML,ERFINV,ZRIDDR,ZBRENT
      DOUBLE PRECISION TAU,MFMEAN,MFVAR,ETAM,ALPHA1,ALPHA2
      DOUBLE PRECISION TAU_MIN, TAU_MAX
      INTEGER IFLAG
            
      DOUBLE PRECISION ROOTF_RE,ROOTF_AE
      COMMON/ROOTFINDERVARS/ROOTF_RE,ROOTF_AE
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      DOUBLE PRECISION MFMEANPASS,MFVARPASS,ETAMPASS,
     $     ALPHA1PASS,ALPHA2PASS
      COMMON/COMVARSMEANVAR/MFMEANPASS,MFVARPASS
      COMMON/COMVARSETAM/ETAMPASS
      COMMON/COMVARSALPHA12/ALPHA1PASS,ALPHA2PASS

      MFMEANPASS = MFMEAN
      MFVARPASS  = MFVAR
      ETAMPASS   = ETAM
      ALPHA1PASS = ALPHA1
      ALPHA2PASS = ALPHA2
            
C     TAU VARIES BETWEEN 0 AND 0.5
      TAU_MIN = D_ZERO
      TAU_MAX = D_HALF
      CALL ZEROIN(TAUFUNC_DSML,TAU_MIN,TAU_MAX,ROOTF_RE,ROOTF_AE,IFLAG)
      TAU = TAU_MIN 
            
      RETURN 
      END

c$$$C==========================================================================================
c$$$C==========================================================================================
c$$$C==========================================================================================
c$$$C========================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION TAUFUNC_DSML(TAU)
c$$$C================================================================================== 
c$$$C     PURPOSE: COMPUTES THE RIGHT-HAND SIDE OF EQ.(19) IN [1]
c$$$C==================================================================================
c$$$C            VARIABLE  DESCRIPTION                    DATA TYPE
c$$$C            --------  -----------                    --------- 
c$$$C
c$$$C     INPUT:
c$$$C
c$$$C            TAU       THE PARAMETER TAU              DOUBLE PRECISION
c$$$C================================================================================== 
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
c$$$      INTEGER INTEG_LIM
c$$$      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
c$$$      COMMON/INTEGRATORVARS2/INTEG_LIM
c$$$      DOUBLE PRECISION TAU,INTEGRAL
c$$$      
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      INTEGER I_TWO,I_FOUR
c$$$      COMMON/INTCONSTANTS/I_TWO,I_FOUR
c$$$    
c$$$      DOUBLE PRECISION TAUFUNC_INT_DSML
c$$$      EXTERNAL TAUFUNC_INT_DSML
c$$$      DOUBLE PRECISION BOUND
c$$$      INTEGER INF,NEVAL,IER,LIMIT,LENW,LAST
c$$$      DOUBLE PRECISION EPSABS,EPSREL,RESULT
c$$$      INTEGER IWORK(INTEG_LIM)
c$$$      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
c$$$      DOUBLE PRECISION INTEG_ABSERR
c$$$      COMMON/INTEGRATORABSERR/INTEG_ABSERR
c$$$
c$$$      DOUBLE PRECISION MFMEANPASS,MFVARPASS,TAUPASS
c$$$      COMMON/COMVARSMEANVAR/MFMEANPASS,MFVARPASS
c$$$      COMMON/COMVARSTAU/TAUPASS
c$$$
c$$$      DOUBLE PRECISION ABSERR,RESABS,RESASC
c$$$           
c$$$      INF     = I_TWO
c$$$      LENW    = I_FOUR*INTEG_LIM
c$$$      EPSABS  = INTEG_EPSABS
c$$$      EPSREL  = INTEG_EPSREL
c$$$      TAUPASS = TAU
c$$$      INTEG_ABSERR = D_ZERO
c$$$      CALL DQAGI(TAUFUNC_INT_DSML,BOUND,INF,EPSABS,EPSREL,RESULT,
c$$$     $     INTEG_ABSERR,NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
c$$$      INTEGRAL = RESULT
c$$$      TAUFUNC_DSML = MFMEANPASS**D_TWO + MFVARPASS - INTEGRAL
c$$$
c$$$CCCCCCCCCCCCCCCCCCCCCC
c$$$
c$$$c$$$      DOUBLE PRECISION D,P
c$$$c$$$      INF = 1
c$$$c$$$      NEVAL = 50
c$$$c$$$      EPSREL = 1.0D-5 
c$$$c$$$      D = PI/4.0D0 !PI !PI/2.0D0
c$$$c$$$      P = 2.0D0 !0.0D0 !1.0D0
c$$$c$$$      TAUPASS = TAU
c$$$c$$$      CALL INTHP(BOUND,BOUND,D,TAUFUNC_INT_DSML,NEVAL,P,
c$$$c$$$     $     EPSREL,INF,RESULT)
c$$$c$$$      IF(INF.GE.10) THEN 
c$$$c$$$         WRITE(*,*) 'ERROR, INF =',INF
c$$$c$$$         STOP
c$$$c$$$      ENDIF
c$$$c$$$      INTEGRAL = RESULT
c$$$c$$$      TAUFUNC_DSML = MFMEANPASS**D_TWO + MFVARPASS - INTEGRAL
c$$$
c$$$CCCCCCCCCCCCCCCCCCCCCC
c$$$
c$$$      RETURN 
c$$$      END
c$$$


C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION TAUFUNC_DSML(TAU)
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
      INTEGER INTEG_LIM,INTEG_METH,INTEG_ORD
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM
      COMMON/INTEGRATORVARS3/INTEG_METH,INTEG_ORD

      DOUBLE PRECISION TAU,SIGMA2,INTEGRAL
      
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
    
      DOUBLE PRECISION TAUFUNC_INT_DSML_GH,TAUFUNC_INT_DSML
      EXTERNAL TAUFUNC_INT_DSML_GH,TAUFUNC_INT_DSML
      DOUBLE PRECISION BOUND
      INTEGER INF,NEVAL,IER,LIMIT,LENW,LAST
      DOUBLE PRECISION EPSABS,EPSREL,RESULT
      INTEGER IWORK(INTEG_LIM)
      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR

      DOUBLE PRECISION MFMEANPASS,MFVARPASS,TAUPASS
      COMMON/COMVARSMEANVAR/MFMEANPASS,MFVARPASS
      COMMON/COMVARSTAU/TAUPASS


c$$$      INTEGER INTEG_METH
c$$$      INTEGER ORD
c$$$      PARAMETER(ORD = 16 INTEG_METH = 1) !20) !30)
c$$$      DOUBLE PRECISION X(ORD),W(ORD)

      INTEGER I,ORD
      DOUBLE PRECISION A,B,ALPHA,X(INTEG_ORD),W(INTEG_ORD)


      IF(INTEG_METH.EQ.1) THEN
C
         INF     = I_TWO
         LENW    = I_FOUR*INTEG_LIM
         EPSABS  = INTEG_EPSABS
         EPSREL  = INTEG_EPSREL
         TAUPASS = TAU
         INTEG_ABSERR = D_ZERO
         CALL DQAGI(TAUFUNC_INT_DSML,BOUND,INF,EPSABS,EPSREL,INTEGRAL,
     $        INTEG_ABSERR,NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
C
      ELSEIF(INTEG_METH.EQ.2) THEN
C
         SIGMA2 = D_ONE - D_TWO*TAU  
         A      = D_ZERO
         B      = D_ONE/(D_TWO*SIGMA2)
         ALPHA  = D_ZERO
         ORD    = INTEG_ORD
         CALL SGHNR(A,B,ALPHA,ORD,X,W)     
         TAUPASS  = TAU
         INTEGRAL = 0.0D0
         DO I = 1,ORD
            INTEGRAL = INTEGRAL + W(I)*TAUFUNC_INT_DSML_GH(X(I))
         ENDDO
C
      ELSE
C
         WRITE(*,*)'INVALID INTEGRATION METHOD'
         WRITE(*,*)'CHECK THE VALUE OF INTEG_METH IN INIT.f'
         WRITE(*,*)'COMPUTATIONS ABORTED'
         STOP
C
      ENDIF
C     
      TAUFUNC_DSML = MFMEANPASS**D_TWO + MFVARPASS - INTEGRAL

      RETURN 
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION FUNCTION TAUFUNC_INT_DSML_GH(PSI)
C================================================================================== 
C     PURPOSE: COMPUTES THE INTEGRAND OF THE SECOND TERM ON THE RIGHT-HAND SIDE 
C              OF EQ.(19) IN [1]
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            PSI       SAMPLE SPACE VARIABLE OF       DOUBLE PRECISION
C                      THE REFERENCE FIELD PSI            
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION PSI,X,R,SIGMA2
      DOUBLE PRECISION ERF
      EXTERNAL ERF
 
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION TERM1,TERM2

      DOUBLE PRECISION ETAMPASS,ALPHA1PASS,ALPHA2PASS,TAUPASS
      COMMON/COMVARSETAM/ETAMPASS
      COMMON/COMVARSALPHA12/ALPHA1PASS,ALPHA2PASS
      COMMON/COMVARSTAU/TAUPASS
      
      TERM1 = ETAMPASS*ERF((PSI-ALPHA1PASS)/(D_TWO*DSQRT(TAUPASS)))
      TERM2 = (D_ONE-ETAMPASS)*ERF((PSI-ALPHA2PASS)
     $     /(D_TWO*DSQRT(TAUPASS)))
      X = D_HALF*(D_ONE+TERM1+TERM2)

      TAUFUNC_INT_DSML_GH = X**D_TWO

      RETURN 
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION FUNCTION TAUFUNC_INT_DSML(PSI)
C================================================================================== 
C     PURPOSE: COMPUTES THE INTEGRAND OF THE SECOND TERM ON THE RIGHT-HAND SIDE 
C              OF EQ.(19) IN [1]
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            PSI       SAMPLE SPACE VARIABLE OF       DOUBLE PRECISION
C                      THE REFERENCE FIELD PSI            
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION PSI,X,R,SIGMA2
      DOUBLE PRECISION ERF
      EXTERNAL ERF
 
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION TERM1,TERM2

      DOUBLE PRECISION ETAMPASS,ALPHA1PASS,ALPHA2PASS,TAUPASS
      COMMON/COMVARSETAM/ETAMPASS
      COMMON/COMVARSALPHA12/ALPHA1PASS,ALPHA2PASS
      COMMON/COMVARSTAU/TAUPASS
      
      TERM1 = ETAMPASS*ERF((PSI-ALPHA1PASS)/(D_TWO*DSQRT(TAUPASS)))
      TERM2 = (D_ONE-ETAMPASS)*ERF((PSI-ALPHA2PASS)
     $     /(D_TWO*DSQRT(TAUPASS)))
      X = D_HALF*(D_ONE+TERM1+TERM2)

      SIGMA2 = D_ONE-D_TWO*TAUPASS
      R = (D_ONE/DSQRT(D_TWO*PI*SIGMA2))
     $     * DEXP(-(PSI**D_TWO)/(D_TWO*SIGMA2))

      TAUFUNC_INT_DSML = (X**D_TWO)*R

      RETURN 
      END


C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION IIPMFPDF_V(MFVAR)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST- AND 
C              SECOND-ORDER PARTIAL DERIVATIVES OF II(ETA) WITH RESPECT TO THE 
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
      DOUBLE PRECISION PMFMFMEANPASS
      COMMON/MFMEANBLOK/PMFMFMEANPASS
      DOUBLE PRECISION PMFCPASS
      COMMON/C1BLOK/PMFCPASS
      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
     $     INT_SIGMA2
      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
     $     INT_SIGMA2
      DOUBLE PRECISION C1,C2,ALPHA1,ALPHA2,TAU,SIGMA2
      DOUBLE PRECISION CDF,PHI,X,TERM1,TERM2
      DOUBLE PRECISION ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
      EXTERNAL ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR

      DOUBLE PRECISION BOUND
      INTEGER INF,NEVAL,IER,LIMIT,LENW,LAST
      DOUBLE PRECISION EPSABS,EPSREL,RESULT
      INTEGER IWORK(INTEG_LIM)
      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM

      C1 = PMFCPASS
      C2 = (PMFMFMEANPASS+C1-D_ONE)/(PMFETAMPASS-D_ONE)
      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
      CALL FIND_TAU_DSML(PMFMFMEANPASS,MFVAR,PMFETAMPASS,
     $     ALPHA1,ALPHA2,TAU)
      SIGMA2 = D_ONE - D_TWO*TAU 
      PHI = INVERSEX(PMFETAPASS,PMFETAMPASS,ALPHA1,ALPHA2,TAU)
      TERM1 = PMFETAMPASS*ERF((PHI-ALPHA1)/(D_TWO*DSQRT(TAU)))
      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-ALPHA2)/(D_TWO*DSQRT(TAU)))
      X = D_HALF*(D_ONE+TERM1+TERM2)
      CDF = CDFNORMDIST(PHI,D_ZERO,SIGMA2)

      INF     = -1
      LENW    = I_FOUR*INTEG_LIM
      EPSABS  = INTEG_EPSABS
      EPSREL  = INTEG_EPSREL
      INTEG_ABSERR = D_ZERO
      INT_TAU    = TAU
      INT_ALPHA1 = ALPHA1
      INT_ALPHA2 = ALPHA2
      INT_MFM    = PMFMFMEANPASS
      INT_MFV    = MFVAR
      INT_SIGMA2 = SIGMA2
      CALL DQAGI(IIEXTFUNC,PHI,INF,EPSABS,EPSREL,RESULT,
     $     INTEG_ABSERR,NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
      
      IIPMFPDF_V = X*CDF-RESULT

      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION IIPMFPDF_M(MFMEAN)
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
      DOUBLE PRECISION PMFMFVARPASS
      COMMON/MFVARBLOK/PMFMFVARPASS
      DOUBLE PRECISION PMFCPASS
      COMMON/C1BLOK/PMFCPASS
      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
     $     INT_SIGMA2
      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
     $     INT_SIGMA2
      DOUBLE PRECISION C1,C2,ALPHA1,ALPHA2,TAU,SIGMA2
      DOUBLE PRECISION CDF,PHI,X,TERM1,TERM2
      DOUBLE PRECISION ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
      EXTERNAL ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR

      DOUBLE PRECISION BOUND
      INTEGER INF,NEVAL,IER,LIMIT,LENW,LAST
      DOUBLE PRECISION EPSABS,EPSREL,RESULT
      INTEGER IWORK(INTEG_LIM)
      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM

      C1 = PMFCPASS
      C2 = (MFMEAN+C1-D_ONE)/(PMFETAMPASS-D_ONE)
      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
      CALL FIND_TAU_DSML(MFMEAN,PMFMFVARPASS,PMFETAMPASS,
     $     ALPHA1,ALPHA2,TAU)
      SIGMA2 = D_ONE - D_TWO*TAU 
      PHI = INVERSEX(PMFETAPASS,PMFETAMPASS,ALPHA1,ALPHA2,TAU)
      TERM1 = PMFETAMPASS*ERF((PHI-ALPHA1)/(D_TWO*DSQRT(TAU)))
      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-ALPHA2)/(D_TWO*DSQRT(TAU)))
      X = D_HALF*(D_ONE+TERM1+TERM2)
      CDF = CDFNORMDIST(PHI,D_ZERO,SIGMA2)

      INF     = -1
      LENW    = I_FOUR*INTEG_LIM
      EPSABS  = INTEG_EPSABS
      EPSREL  = INTEG_EPSREL
      INTEG_ABSERR = D_ZERO
      INT_TAU    = TAU
      INT_ALPHA1 = ALPHA1
      INT_ALPHA2 = ALPHA2
      INT_MFM    = MFMEAN
      INT_MFV    = PMFMFVARPASS
      INT_SIGMA2 = SIGMA2
      CALL DQAGI(IIEXTFUNC,PHI,INF,EPSABS,EPSREL,RESULT,
     $     INTEG_ABSERR,NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
      
      IIPMFPDF_M = X*CDF-RESULT

      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION IIPMFPDF_MV(MFMEAN)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE MIXED 
C              DERIVATIVE OF II(ETA) WITH RESPECT TO THE [MIXTURE FRACTION 
C              MEAN AND THE MIXTURE FRACTION VARIANCE].
C              THIS FUNCTION IS INVOKED WHEN DIFFERENTIATION IS PERFORMED
C              USING RIDDERS' METHOD
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
      DOUBLE PRECISION PMFMFMEANPASS
      COMMON/MFMEANBLOK/PMFMFMEANPASS
      DOUBLE PRECISION PMFMFVARPASS
      COMMON/MFVARBLOK/PMFMFVARPASS

      DOUBLE PRECISION H_MEAN,H_VAR,H_NEW
      DOUBLE PRECISION ERR,ERR_OLD
      DOUBLE PRECISION IIPMFPDF_V
      EXTERNAL IIPMFPDF_V
      DOUBLE PRECISION LIMINF,LIMSUP
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX

      INTEGER DIFF_NCHEBCOEFS_TEMP
      DOUBLE PRECISION DIFF_MAXERR_TEMP
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
      DOUBLE PRECISION PMFPDFPASS
      COMMON/PMFPDFBLOK/PMFPDFPASS
      DOUBLE PRECISION PDF

      MFVAR = PMFMFVARPASS
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)

C     DIIDV
      PMFMFMEANPASS = MFMEAN
      IF(PMFPDFPASS.LE.DIFF_PDFMIN) THEN
         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
      ELSE
         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
      ENDIF
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIPMFPDF_V,
     $     MFVAR,1,IIPMFPDF_MV)
    
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION IIPMFPDF_VM(MFVAR)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE MIXED 
C              DERIVATIVE OF II(ETA) WITH RESPECT TO THE [MIXTURE FRACTION 
C              MEAN AND THE MIXTURE FRACTION VARIANCE].
C              THIS FUNCTION IS INVOKED WHEN DIFFERENTIATION IS PERFORMED
C              USING RIDDERS' METHOD
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
      DOUBLE PRECISION PMFMFMEANPASS
      COMMON/MFMEANBLOK/PMFMFMEANPASS
      DOUBLE PRECISION PMFMFVARPASS
      COMMON/MFVARBLOK/PMFMFVARPASS

      DOUBLE PRECISION H_MEAN,H_VAR,H_NEW
      DOUBLE PRECISION ERR,ERR_OLD
      DOUBLE PRECISION IIPMFPDF_M
      EXTERNAL IIPMFPDF_M
      DOUBLE PRECISION LIMINF,LIMSUP
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX

      DOUBLE PRECISION DIFF_NCHEBCOEFS_TEMP
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
      DOUBLE PRECISION PMFPDFPASS
      COMMON/PMFPDFBLOK/PMFPDFPASS
      DOUBLE PRECISION PDF

      MFMEAN = PMFMFMEANPASS
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)

C     DIIDM
      PMFMFVARPASS = MFVAR
      IF(PMFPDFPASS.LE.DIFF_PDFMIN) THEN
         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
      ELSE
         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
      ENDIF
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIPMFPDF_M,
     $     MFMEAN,1,IIPMFPDF_VM)

      RETURN
      END
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION IIEXTFUNC(PHI)
      IMPLICIT NONE
      DOUBLE PRECISION PHI
      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
     $     INT_SIGMA2
      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
     $     INT_SIGMA2
      DOUBLE PRECISION R,X,TERM1,TERM2
      DOUBLE PRECISION ERF,NORMDIST
      EXTERNAL ERF,NORMDIST
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      TERM1 = PMFETAMPASS*ERF((PHI-INT_ALPHA1)
     $     /(D_TWO*DSQRT(INT_TAU)))
      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-INT_ALPHA2)
     $     /(D_TWO*DSQRT(INT_TAU)))
      X = D_HALF*(D_ONE+TERM1+TERM2)
      R = NORMDIST(PHI,D_ZERO,INT_SIGMA2)

      IIEXTFUNC = X*R

      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION  FUNCTION PMFDSML_PDF_M(MFMEAN)
C================================================================================== 
C     PURPOSE: COMPUTES THE PMF PROBABILITY DENSITY FUNCTION AT A GIVEN ETA VALUE:
C              VARIABLE MEAN, FIXED VARIANCE
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            ETA       MIXTURE FRACTION VALUE         DOUBLE PRECISION
C            NETA      SIZE OF ETA AND PDF            INTEGER
C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION
C
C     OUTPUT:
C
C            PDF       PROBABILITY DENSITY FUNCTION   DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      INCLUDE 'formats.h'
      DOUBLE PRECISION MFMEAN,MFVAR,TAU,SIGMA2,PSI,ETAM,C,C1,C2
      DOUBLE PRECISION ETA,PDF,E
      DOUBLE PRECISION PDF_AREA,PDF_LA,PDF_RA,PDF_SCALEFACT
      DOUBLE PRECISION ALPHA1,ALPHA2
      DOUBLE PRECISION EXP_NUM,EXP_DEN_1,EXP_DEN_2
      DOUBLE PRECISION TERM1,TERM2,X
      DOUBLE PRECISION ERFINV,INVERSEX
      EXTERNAL ERFINV,INVERSEX
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE,
     $     D_FOUR,D_EIGHT
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS2/D_FOUR,D_EIGHT 
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR

      DOUBLE PRECISION PMFCV_ETA,PMFCV_ETAM,PMFCV_C,PMFCV_VAR
      COMMON/PMFCVETA/PMFCV_ETA
      COMMON/PMFCVETAM/PMFCV_ETAM
      COMMON/PMFCVC/PMFCV_C
      COMMON/PMFCVVAR/PMFCV_VAR

      ETA = PMFCV_ETA
      ETAM = PMFCV_ETAM
      C = PMFCV_C
      MFVAR = PMFCV_VAR

      C1 = C
      C2 = (MFMEAN+C1-D_ONE)/(ETAM-D_ONE)
      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
      CALL FIND_TAU_DSML(MFMEAN,MFVAR,ETAM,ALPHA1,ALPHA2,TAU) 
      SIGMA2 = D_ONE - D_TWO*TAU  
      PSI = INVERSEX(ETA,ETAM,ALPHA1,ALPHA2,TAU)
      EXP_NUM = DEXP(-(PSI**D_TWO)/(D_TWO*(D_ONE-D_TWO*TAU)))
      EXP_DEN_1 = DEXP(-((PSI-ALPHA1)**D_TWO)/(D_FOUR*TAU))
      EXP_DEN_2 = DEXP(-((PSI-ALPHA2)**D_TWO)/(D_FOUR*TAU))
      PMFDSML_PDF_M = DSQRT(D_TWO*TAU/SIGMA2)
     $     * EXP_NUM/(ETAM*EXP_DEN_1 + (D_ONE-ETAM)*EXP_DEN_2)

      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION  FUNCTION PMFDSML_PDF_V(MFVAR)
C================================================================================== 
C     PURPOSE: COMPUTES THE PMF PROBABILITY DENSITY FUNCTION AT A GIVEN ETA VALUE:
C              FIXED MEAN, VARIABLE VARIANCE
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            ETA       MIXTURE FRACTION VALUE         DOUBLE PRECISION
C            NETA      SIZE OF ETA AND PDF            INTEGER
C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION
C
C     OUTPUT:
C
C            PDF       PROBABILITY DENSITY FUNCTION   DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      INCLUDE 'formats.h'
      DOUBLE PRECISION MFMEAN,MFVAR,TAU,SIGMA2,PSI,ETAM,C,C1,C2
      DOUBLE PRECISION ETA,PDF,E
      DOUBLE PRECISION PDF_AREA,PDF_LA,PDF_RA,PDF_SCALEFACT
      DOUBLE PRECISION ALPHA1,ALPHA2
      DOUBLE PRECISION EXP_NUM,EXP_DEN_1,EXP_DEN_2
      DOUBLE PRECISION TERM1,TERM2,X
      DOUBLE PRECISION ERFINV,INVERSEX
      EXTERNAL ERFINV,INVERSEX
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE,
     $     D_FOUR,D_EIGHT
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS2/D_FOUR,D_EIGHT 
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR

      DOUBLE PRECISION PMFCV_ETA,PMFCV_ETAM,PMFCV_C,PMFCV_MEAN
      COMMON/PMFCVETA/PMFCV_ETA
      COMMON/PMFCVETAM/PMFCV_ETAM
      COMMON/PMFCVC/PMFCV_C
      COMMON/PMFCVMEAN/PMFCV_MEAN

      ETA = PMFCV_ETA
      ETAM = PMFCV_ETAM
      C = PMFCV_C
      MFMEAN = PMFCV_MEAN

      C1 = C
      C2 = (MFMEAN+C1-D_ONE)/(ETAM-D_ONE)
      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
      CALL FIND_TAU_DSML(MFMEAN,MFVAR,ETAM,ALPHA1,ALPHA2,TAU) 
      SIGMA2 = D_ONE - D_TWO*TAU  
      PSI = INVERSEX(ETA,ETAM,ALPHA1,ALPHA2,TAU)
      EXP_NUM = DEXP(-(PSI**D_TWO)/(D_TWO*(D_ONE-D_TWO*TAU)))
      EXP_DEN_1 = DEXP(-((PSI-ALPHA1)**D_TWO)/(D_FOUR*TAU))
      EXP_DEN_2 = DEXP(-((PSI-ALPHA2)**D_TWO)/(D_FOUR*TAU))
      PMFDSML_PDF_V = DSQRT(D_TWO*TAU/SIGMA2)
     $     * EXP_NUM/(ETAM*EXP_DEN_1 + (D_ONE-ETAM)*EXP_DEN_2)

      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================



c$$$C==================================================================================
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION IIPMFPDF_V(MFVAR)
c$$$C================================================================================== 
c$$$C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST- AND 
c$$$C              SECOND-ORDER PARTIAL DERIVATIVES OF II(ETA) WITH RESPECT TO THE 
c$$$C              [MIXTURE FRACTION VARIANCE]
c$$$C==================================================================================
c$$$C            VARIABLE  DESCRIPTION                    DATA TYPE
c$$$C            --------  -----------                    --------- 
c$$$C
c$$$C     INPUT:
c$$$C
c$$$C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION
c$$$C================================================================================== 
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION MFVAR
c$$$      DOUBLE PRECISION PMFMFMEANPASS
c$$$      COMMON/MFMEANBLOK/PMFMFMEANPASS
c$$$      DOUBLE PRECISION PMFCPASS
c$$$      COMMON/C1BLOK/PMFCPASS
c$$$      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
c$$$      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
c$$$      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2,INT_PHI
c$$$      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2,INT_PHI
c$$$      DOUBLE PRECISION C1,C2,ALPHA1,ALPHA2,TAU,SIGMA2
c$$$      DOUBLE PRECISION CDF,PHI,X,TERM1,TERM2
c$$$      DOUBLE PRECISION ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
c$$$      EXTERNAL ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      INTEGER I_TWO,I_FOUR
c$$$      COMMON/INTCONSTANTS/I_TWO,I_FOUR
c$$$
c$$$      DOUBLE PRECISION BOUND
c$$$      INTEGER INF,NEVAL,IER,LIMIT,LENW,LAST
c$$$      DOUBLE PRECISION EPSABS,EPSREL,RESULT
c$$$      INTEGER IWORK(INTEG_LIM)
c$$$      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
c$$$      DOUBLE PRECISION INTEG_ABSERR
c$$$      COMMON/INTEGRATORABSERR/INTEG_ABSERR
c$$$      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
c$$$      INTEGER INTEG_LIM
c$$$      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
c$$$      COMMON/INTEGRATORVARS2/INTEG_LIM
c$$$
c$$$      INTEGER ORD
c$$$      PARAMETER(ORD = 40) !16) !20) !30)
c$$$      DOUBLE PRECISION XGH(ORD),WGH(ORD)
c$$$      DOUBLE PRECISION XGL(ORD),WGL(ORD)
c$$$
c$$$      DOUBLE PRECISION A
c$$$      DOUBLE PRECISION ALPHA
c$$$      DOUBLE PRECISION B
c$$$      DOUBLE PRECISION BETA
c$$$      INTEGER I,SCALE
c$$$      DOUBLE PRECISION INTEGRAL1,INTEGRAL2
c$$$      DOUBLE PRECISION IIEXTFUNC_GH,IIEXTFUNC_GL
c$$$      EXTERNAL IIEXTFUNC_GH,IIEXTFUNC_GL
c$$$
c$$$      C1 = PMFCPASS
c$$$      C2 = (PMFMFMEANPASS+C1-D_ONE)/(PMFETAMPASS-D_ONE)
c$$$      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
c$$$      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
c$$$      CALL FIND_TAU_DSML(PMFMFMEANPASS,MFVAR,PMFETAMPASS,
c$$$     $     ALPHA1,ALPHA2,TAU)
c$$$      SIGMA2 = D_ONE - D_TWO*TAU 
c$$$      PHI = INVERSEX(PMFETAPASS,PMFETAMPASS,ALPHA1,ALPHA2,TAU)
c$$$      TERM1 = PMFETAMPASS*ERF((PHI-ALPHA1)/(D_TWO*DSQRT(TAU)))
c$$$      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-ALPHA2)/(D_TWO*DSQRT(TAU)))
c$$$      X = D_HALF*(D_ONE+TERM1+TERM2)
c$$$      CDF = CDFNORMDIST(PHI,D_ZERO,SIGMA2)
c$$$
c$$$      INF     = -1
c$$$      LENW    = I_FOUR*INTEG_LIM
c$$$      EPSABS  = INTEG_EPSABS
c$$$      EPSREL  = INTEG_EPSREL
c$$$      INTEG_ABSERR = D_ZERO
c$$$      INT_TAU    = TAU
c$$$      INT_ALPHA1 = ALPHA1
c$$$      INT_ALPHA2 = ALPHA2
c$$$      INT_MFM    = PMFMFMEANPASS
c$$$      INT_MFV    = MFVAR
c$$$      INT_SIGMA2 = SIGMA2
c$$$      INT_PHI    = PHI
c$$$C      CALL DQAGI(IIEXTFUNC,PHI,INF,EPSABS,EPSREL,RESULT,
c$$$C     $     INTEG_ABSERR,NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
c$$$
c$$$      A = D_ZERO
c$$$      B = D_ONE/(D_TWO*SIGMA2)
c$$$      ALPHA = D_ZERO
c$$$      CALL SGHNR(A,B,ALPHA,ORD,XGH,WGH)     
c$$$      INTEGRAL1 = 0.0D0
c$$$      DO I = 1,ORD
c$$$         INTEGRAL1 = INTEGRAL1 + WGH(I)*IIEXTFUNC_GH(XGH(I))
c$$$      ENDDO
c$$$
c$$$      A = PHI !D_ZERO
c$$$      B = D_ONE
c$$$      ALPHA = D_ZERO
c$$$      CALL UGLNR(A,B,ALPHA,ORD,XGL,WGL)     
c$$$      INTEGRAL2 = 0.0D0
c$$$      DO I = 1,ORD
c$$$         INTEGRAL2 = INTEGRAL2 + WGL(I)*IIEXTFUNC_GL(XGL(I))
c$$$      ENDDO
c$$$
c$$$      RESULT = INTEGRAL1-INTEGRAL2
c$$$
c$$$      IIPMFPDF_V = X*CDF-RESULT
c$$$
c$$$      RETURN
c$$$      END
c$$$
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION IIPMFPDF_M(MFMEAN)
c$$$C================================================================================== 
c$$$C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE SECOND-ORDER 
c$$$C              PARTIAL DERIVATIVE OF II(ETA) WITH RESPECT TO THE 
c$$$C              [MIXTURE FRACTION MEAN]
c$$$C==================================================================================
c$$$C            VARIABLE  DESCRIPTION                    DATA TYPE
c$$$C            --------  -----------                    --------- 
c$$$C
c$$$C     INPUT:
c$$$C
c$$$C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
c$$$C================================================================================== 
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION MFMEAN
c$$$      DOUBLE PRECISION PMFMFVARPASS
c$$$      COMMON/MFVARBLOK/PMFMFVARPASS
c$$$      DOUBLE PRECISION PMFCPASS
c$$$      COMMON/C1BLOK/PMFCPASS
c$$$      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
c$$$      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
c$$$      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2
c$$$      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2
c$$$      DOUBLE PRECISION C1,C2,ALPHA1,ALPHA2,TAU,SIGMA2
c$$$      DOUBLE PRECISION CDF,PHI,X,TERM1,TERM2
c$$$      DOUBLE PRECISION ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
c$$$      EXTERNAL ERFINV,ERF,INVERSEX,CDFNORMDIST,IIEXTFUNC
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      INTEGER I_TWO,I_FOUR
c$$$      COMMON/INTCONSTANTS/I_TWO,I_FOUR
c$$$
c$$$      DOUBLE PRECISION BOUND
c$$$      INTEGER INF,NEVAL,IER,LIMIT,LENW,LAST
c$$$      DOUBLE PRECISION EPSABS,EPSREL,RESULT
c$$$      INTEGER IWORK(INTEG_LIM)
c$$$      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
c$$$      DOUBLE PRECISION INTEG_ABSERR
c$$$      COMMON/INTEGRATORABSERR/INTEG_ABSERR
c$$$      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
c$$$      INTEGER INTEG_LIM
c$$$      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
c$$$      COMMON/INTEGRATORVARS2/INTEG_LIM
c$$$
c$$$      C1 = PMFCPASS
c$$$      C2 = (MFMEAN+C1-D_ONE)/(PMFETAMPASS-D_ONE)
c$$$      ALPHA1 = DSQRT(D_TWO)*ERFINV(D_TWO*C1 - D_ONE)
c$$$      ALPHA2 = DSQRT(D_TWO)*ERFINV(D_TWO*(C1+C2) - D_ONE)
c$$$      CALL FIND_TAU_DSML(MFMEAN,PMFMFVARPASS,PMFETAMPASS,
c$$$     $     ALPHA1,ALPHA2,TAU)
c$$$      SIGMA2 = D_ONE - D_TWO*TAU 
c$$$      PHI = INVERSEX(PMFETAPASS,PMFETAMPASS,ALPHA1,ALPHA2,TAU)
c$$$      TERM1 = PMFETAMPASS*ERF((PHI-ALPHA1)/(D_TWO*DSQRT(TAU)))
c$$$      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-ALPHA2)/(D_TWO*DSQRT(TAU)))
c$$$      X = D_HALF*(D_ONE+TERM1+TERM2)
c$$$      CDF = CDFNORMDIST(PHI,D_ZERO,SIGMA2)
c$$$
c$$$      INF     = -1
c$$$      LENW    = I_FOUR*INTEG_LIM
c$$$      EPSABS  = INTEG_EPSABS
c$$$      EPSREL  = INTEG_EPSREL
c$$$      INTEG_ABSERR = D_ZERO
c$$$      INT_TAU    = TAU
c$$$      INT_ALPHA1 = ALPHA1
c$$$      INT_ALPHA2 = ALPHA2
c$$$      INT_MFM    = MFMEAN
c$$$      INT_MFV    = PMFMFVARPASS
c$$$      INT_SIGMA2 = SIGMA2
c$$$      CALL DQAGI(IIEXTFUNC,PHI,INF,EPSABS,EPSREL,RESULT,
c$$$     $     INTEG_ABSERR,NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
c$$$      
c$$$      IIPMFPDF_M = X*CDF-RESULT
c$$$
c$$$      RETURN
c$$$      END
c$$$
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION IIPMFPDF_MV(MFMEAN)
c$$$C================================================================================== 
c$$$C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE MIXED 
c$$$C              DERIVATIVE OF II(ETA) WITH RESPECT TO THE [MIXTURE FRACTION 
c$$$C              MEAN AND THE MIXTURE FRACTION VARIANCE].
c$$$C              THIS FUNCTION IS INVOKED WHEN DIFFERENTIATION IS PERFORMED
c$$$C              USING RIDDERS' METHOD
c$$$C==================================================================================
c$$$C            VARIABLE  DESCRIPTION                    DATA TYPE
c$$$C            --------  -----------                    --------- 
c$$$C
c$$$C     INPUT:
c$$$C
c$$$C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
c$$$C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION 
c$$$C==================================================================================
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION MFMEAN
c$$$      DOUBLE PRECISION MFVAR
c$$$      DOUBLE PRECISION PMFMFMEANPASS
c$$$      COMMON/MFMEANBLOK/PMFMFMEANPASS
c$$$      DOUBLE PRECISION PMFMFVARPASS
c$$$      COMMON/MFVARBLOK/PMFMFVARPASS
c$$$
c$$$      DOUBLE PRECISION H_MEAN,H_VAR,H_NEW
c$$$      DOUBLE PRECISION ERR,ERR_OLD
c$$$      DOUBLE PRECISION IIPMFPDF_V
c$$$      EXTERNAL IIPMFPDF_V
c$$$      DOUBLE PRECISION LIMINF,LIMSUP
c$$$      INTEGER DIFF_NCHEBCOEFS
c$$$      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
c$$$      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
c$$$      DOUBLE PRECISION MFVAR_MAX
c$$$
c$$$      INTEGER DIFF_NCHEBCOEFS_TEMP
c$$$      DOUBLE PRECISION DIFF_MAXERR_TEMP
c$$$      DOUBLE PRECISION DIFF_PDFMIN
c$$$      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
c$$$      INTEGER DIFF_NCHEBCOEFS_EXTRA
c$$$      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
c$$$      DOUBLE PRECISION PMFPDFPASS
c$$$      COMMON/PMFPDFBLOK/PMFPDFPASS
c$$$      DOUBLE PRECISION PDF
c$$$
c$$$      MFVAR = PMFMFVARPASS
c$$$      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
c$$$      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
c$$$      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
c$$$      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
c$$$
c$$$C     DIIDV
c$$$      PMFMFMEANPASS = MFMEAN
c$$$      IF(PMFPDFPASS.LE.DIFF_PDFMIN) THEN
c$$$         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
c$$$      ELSE
c$$$         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
c$$$      ENDIF
c$$$      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIPMFPDF_V,
c$$$     $     MFVAR,1,IIPMFPDF_MV)
c$$$    
c$$$      RETURN
c$$$      END
c$$$
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION IIPMFPDF_VM(MFVAR)
c$$$C================================================================================== 
c$$$C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE MIXED 
c$$$C              DERIVATIVE OF II(ETA) WITH RESPECT TO THE [MIXTURE FRACTION 
c$$$C              MEAN AND THE MIXTURE FRACTION VARIANCE].
c$$$C              THIS FUNCTION IS INVOKED WHEN DIFFERENTIATION IS PERFORMED
c$$$C              USING RIDDERS' METHOD
c$$$C==================================================================================
c$$$C            VARIABLE  DESCRIPTION                    DATA TYPE
c$$$C            --------  -----------                    --------- 
c$$$C
c$$$C     INPUT:
c$$$C
c$$$C            MFMEAN    MIXTURE FRACTION MEAN          DOUBLE PRECISION 
c$$$C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION 
c$$$C==================================================================================
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION MFMEAN
c$$$      DOUBLE PRECISION MFVAR
c$$$      DOUBLE PRECISION PMFMFMEANPASS
c$$$      COMMON/MFMEANBLOK/PMFMFMEANPASS
c$$$      DOUBLE PRECISION PMFMFVARPASS
c$$$      COMMON/MFVARBLOK/PMFMFVARPASS
c$$$
c$$$      DOUBLE PRECISION H_MEAN,H_VAR,H_NEW
c$$$      DOUBLE PRECISION ERR,ERR_OLD
c$$$      DOUBLE PRECISION IIPMFPDF_M
c$$$      EXTERNAL IIPMFPDF_M
c$$$      DOUBLE PRECISION LIMINF,LIMSUP
c$$$      INTEGER DIFF_NCHEBCOEFS
c$$$      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
c$$$      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
c$$$      DOUBLE PRECISION MFVAR_MAX
c$$$
c$$$      DOUBLE PRECISION DIFF_NCHEBCOEFS_TEMP
c$$$      DOUBLE PRECISION DIFF_PDFMIN
c$$$      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
c$$$      INTEGER DIFF_NCHEBCOEFS_EXTRA
c$$$      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
c$$$      DOUBLE PRECISION PMFPDFPASS
c$$$      COMMON/PMFPDFBLOK/PMFPDFPASS
c$$$      DOUBLE PRECISION PDF
c$$$
c$$$      MFMEAN = PMFMFMEANPASS
c$$$      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
c$$$      LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
c$$$      LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
c$$$
c$$$C     DIIDM
c$$$      PMFMFVARPASS = MFVAR
c$$$      IF(PMFPDFPASS.LE.DIFF_PDFMIN) THEN
c$$$         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
c$$$      ELSE
c$$$         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
c$$$      ENDIF
c$$$      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIPMFPDF_M,
c$$$     $     MFMEAN,1,IIPMFPDF_VM)
c$$$
c$$$      RETURN
c$$$      END
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION IIEXTFUNC_GL(PHI)
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION PHI
c$$$      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
c$$$      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
c$$$      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2,INT_PHI
c$$$      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2,INT_PHI
c$$$      DOUBLE PRECISION R,X,TERM1,TERM2
c$$$      DOUBLE PRECISION ERF,NORMDIST
c$$$      EXTERNAL ERF,NORMDIST
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$
c$$$      TERM1 = PMFETAMPASS*ERF((PHI-INT_ALPHA1)
c$$$     $     /(D_TWO*DSQRT(INT_TAU)))
c$$$      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-INT_ALPHA2)
c$$$     $     /(D_TWO*DSQRT(INT_TAU)))
c$$$      X = D_HALF*(D_ONE+TERM1+TERM2)
c$$$      R = NORMDIST(PHI,D_ZERO,INT_SIGMA2)
c$$$
c$$$      IIEXTFUNC_GL = X*R*DEXP(PHI-INT_PHI)
c$$$
c$$$      RETURN
c$$$      END
c$$$
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION IIEXTFUNC_GH(PHI)
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION PHI
c$$$      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
c$$$      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
c$$$      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2,INT_PHI
c$$$      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2,INT_PHI
c$$$      DOUBLE PRECISION R,X,TERM1,TERM2
c$$$      DOUBLE PRECISION ERF,NORMDIST
c$$$      EXTERNAL ERF,NORMDIST
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$
c$$$      TERM1 = PMFETAMPASS*ERF((PHI-INT_ALPHA1)
c$$$     $     /(D_TWO*DSQRT(INT_TAU)))
c$$$      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-INT_ALPHA2)
c$$$     $     /(D_TWO*DSQRT(INT_TAU)))
c$$$      X = D_HALF*(D_ONE+TERM1+TERM2)
c$$$
c$$$      IIEXTFUNC_GH = X
c$$$
c$$$      RETURN
c$$$      END
c$$$
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C==================================================================================
c$$$C================================================================================== 
c$$$
c$$$      DOUBLE PRECISION FUNCTION IIEXTFUNC(PHI)
c$$$      IMPLICIT NONE
c$$$      DOUBLE PRECISION PHI
c$$$      DOUBLE PRECISION PMFETAPASS,PMFETAMPASS
c$$$      COMMON/ETABLOK/PMFETAPASS,PMFETAMPASS
c$$$      DOUBLE PRECISION INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2
c$$$      COMMON/IIINTEGBLOK/INT_TAU,INT_ALPHA1,INT_ALPHA2,INT_MFM,INT_MFV,
c$$$     $     INT_SIGMA2
c$$$      DOUBLE PRECISION R,X,TERM1,TERM2
c$$$      DOUBLE PRECISION ERF,NORMDIST
c$$$      EXTERNAL ERF,NORMDIST
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$
c$$$      TERM1 = PMFETAMPASS*ERF((PHI-INT_ALPHA1)
c$$$     $     /(D_TWO*DSQRT(INT_TAU)))
c$$$      TERM2 = (D_ONE-PMFETAMPASS)*ERF((PHI-INT_ALPHA2)
c$$$     $     /(D_TWO*DSQRT(INT_TAU)))
c$$$      X = D_HALF*(D_ONE+TERM1+TERM2)
c$$$      R = NORMDIST(PHI,D_ZERO,INT_SIGMA2)
c$$$
c$$$      IIEXTFUNC = X*R
c$$$
c$$$      RETURN
c$$$      END
