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
      DOUBLE PRECISION ETA_TEMP1(NETA-1),PDF_TEMP1(NETA-1)
      DOUBLE PRECISION ETA_TEMP2(NETA-2),PDF_TEMP2(NETA-2)
      DOUBLE PRECISION PDF_AREA,PDF_LA,PDF_RA,PDF_MA,PDF_SCALEFACT
      DOUBLE PRECISION BETA,BETAI,TRAP
      EXTERNAL BETA,BETAI,TRAP
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
    
      G = (MFMEAN*(D_ONE-MFMEAN)/MFVAR) - D_ONE
      V = MFMEAN*G
      W = (D_ONE-MFMEAN)*G
      B = BETA(V,W)

      DO IETA = 1,NETA
         PDF(IETA) = (ETA(IETA)**(V-D_ONE))
     $        *((D_ONE-ETA(IETA))**(W-D_ONE))/B
      ENDDO

C     DETERMINE SHAPE TO SET BOUNDARY CONDITIONS
C
      IF(((V.GT.D_ONE).AND.(W.GT.D_ONE))) THEN ! UNIMODAL:  POSITIVELY SKEWED  IF V < W
                                               !            NEGATIVELY SKEWED  IF V > W
                                               !            SYMMETRIC IF V = W
         PDF(1) = D_ZERO
         PDF(NETA) = D_ZERO
C
      ELSEIF(((V.LT.D_ONE).AND.(W.GE.D_ONE))   ! REVERSE J-SHAPED
     $     .OR.((V.EQ.D_ONE).AND.(W.GT.D_ONE))) THEN
C     USE THE INCOMPLETE BETA FUNCTION TO FIND THE AREA BETWEEN 
C     THE LEFT BOUNADRY AND THE FIRST INTERIOR GRID POINT
            PDF_LA = BETAI(V,W,ETA(2))
            PDF(1) = D_TWO*PDF_LA/(ETA(2)-ETA(1)) - PDF(2)        
C
      ELSEIF(((V.GE.D_ONE).AND.(W.LT.D_ONE))   ! J-SHAPED
     $        .OR.((V.GT.D_ONE).AND.(W.EQ.D_ONE))) THEN
C     USE THE INCOMPLETE BETA FUNCTION TO FIND THE AREA BETWEEN 
C     THE LEFT BOUNDARY AND THE LAST INTERIOR GRID POINT. 
            PDF_RA = D_ONE - BETAI(V,W,ETA(NETA-1))
            PDF(NETA) = D_TWO*PDF_RA/(ETA(NETA)-ETA(NETA-1))
     $           - PDF(NETA-1)
C
      ELSEIF((V.LT.D_ONE).AND.(W.LT.D_ONE)) THEN ! U-SHAPED
C     USE THE INCOMPLETE BETA FUNCTION TO FIND THE AREA BETWEEN 
C     THE LEFT BOUNADRY AND THE FIRST INTERIOR GRID POINT
         PDF_LA = BETAI(V,W,ETA(2))
         PDF(1) = D_TWO*PDF_LA/(ETA(2)-ETA(1)) - PDF(2)
C     USE THE INCOMPLETE BETA FUNCTION TO FIND THE AREA BETWEEN 
C     THE LEFT BOUNDARY AND THE LAST INTERIOR GRID POINT
         PDF_RA = D_ONE - BETAI(V,W,ETA(NETA-1))
         PDF(NETA) = D_TWO*PDF_RA/(ETA(NETA)-ETA(NETA-1))
     $        - PDF(NETA-1)
            
      ENDIF
      
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
      INTEGER I,IETA,NETA
      DOUBLE PRECISION MFMEAN,MFVAR,DT
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3),VEL(3)
      DOUBLE PRECISION ETA(NETA),CV(3,NETA)
      DOUBLE PRECISION G,V,W,B,PDF(NETA)
      DOUBLE PRECISION V_VAR,V_MEAN,W_VAR,W_MEAN,PDF_V,PDF_W,C1,C2
      DOUBLE PRECISION ETA_TEMP
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION BETA,DPSI
      EXTERNAL BETA,DPSI
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS

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
C         WRITE(*,*)'BETA CV ETA INDEX',IETA
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
C            ETA        MIXTURE FRACTION GRID           DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA AND CSDRH           INTEGER
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
      DOUBLE PRECISION DIIDV(NETA)
      DOUBLE PRECISION IM,H_NEW,ERR,ERR_OLD,H_MEAN,H_VAR    
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION IIBPDF_V
      EXTERNAL IIBPDF_V
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION ETAPASS
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINF,LIMSUP

      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH

      INTEGER DIFF_NCHEBCOEFS_TEMP
      
      CALL CHECKPARMS(MFMEAN,MFVAR)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)

C     DIIDV
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      MFMEANPASS = MFMEAN
      DO IETA = 2,NETA-1
C     
C     ADJUST NUMBER OF COEFS IF NECESSARY
         IF(PDF(IETA).LE.DIFF_PDFMIN) THEN
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
         ELSE
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
         ENDIF          
C         WRITE(*,*)'BETA CSDR-H',IETA,PDF(IETA),DIFF_NCHEBCOEFS
C     
         ETAPASS = ETA(IETA)
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIBPDF_V,
     $        MFVAR,1,DIIDV(IETA))
C          
      ENDDO
C       
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
C            ETA        MIXTURE FRACTION GRID           DOUBLE PRECISION (ARRAY,SIZE=NETA)
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
      DOUBLE PRECISION DIIDV(NETA)
      DOUBLE PRECISION DIIDM(NETA)
      DOUBLE PRECISION D2IIDV2(NETA)
      DOUBLE PRECISION D2IIDM2(NETA)
      DOUBLE PRECISION D2IIDMDV(NETA)
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION IIBPDF_M,IIBPDF_V,IIBPDF_MV
      EXTERNAL IIBPDF_M,IIBPDF_V,IIBPDF_MV

      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION ETAPASS
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION BPDFPASS
      COMMON/PPASS/BPDFPASS

      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINF,LIMSUP

      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
      INTEGER DIFF_NCHEBCOEFS_TEMP
      INTEGER FLAGOS,FLAGUS,IETAINTOS,IETAINTUS
      
      CALL CHECKPARMS(MFMEAN,MFVAR)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)

      DO IETA = 2,NETA-1
C  
C     ADJUST NUMBER OF COEFS IF NECESSARY
         IF(PDF(IETA).LE.DIFF_PDFMIN) THEN
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
         ELSE
            DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
         ENDIF  
C         WRITE(*,*)'BETA CSDR-I',IETA,PDF(IETA),DIFF_NCHEBCOEFS
C     
C     DIIDV AND D2IIDV2
         LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
         LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
         MFMEANPASS = MFMEAN
         ETAPASS = ETA(IETA)          
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIBPDF_V,
     $        MFVAR,1,DIIDV(IETA))
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIBPDF_V,
     $        MFVAR,2,D2IIDV2(IETA))
C     
C     D2IIDM2
         LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
         LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
         MFVARPASS = MFVAR
         ETAPASS = ETA(IETA)
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIBPDF_M,
     $        MFMEAN,2,D2IIDM2(IETA))
C     
C     D2IIDMDV
         LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
         LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
         MFVARPASS = MFVAR
         BPDFPASS = PDF(IETA)
         ETAPASS = ETA(IETA)
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIBPDF_MV,
     $        MFMEAN,1,D2IIDMDV(IETA))
C
      ENDDO
C     
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
     $        *(-DIIDV(IETA)*(-CHI + D_TWO*DT*PROD1)
     $        + DT*(PROD1*D2IIDM2(IETA) + PROD2*D2IIDV2(IETA)
     $        + D_TWO*PROD3*D2IIDMDV(IETA)))
                                         
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
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

      DOUBLE PRECISION FUNCTION IIBPDF_V(MFVAR)
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
      IIBPDF_V = (ETAPASS-MFMEANPASS)*BETAI(B1,B2,ETAPASS)+ MFVAR*P

      RETURN
      END

C===========================================================================
C===========================================================================
C===========================================================================
     
      DOUBLE PRECISION FUNCTION IIBPDF_M(MFMEAN)
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
      IIBPDF_M = (ETAPASS-MFMEAN)*BETAI(B1,B2,ETAPASS)+ MFVARPASS*P

      RETURN
      END


C===========================================================================
C===========================================================================
C===========================================================================

      DOUBLE PRECISION FUNCTION IIBPDF_MV(MFMEAN)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE MIXED 
C              DERIVATIVE OF II(ETA) WITH RESPECT TO THE [MIXTURE FRACTION 
C              MEAN AND THE MIXTURE FRACTION VARIANCE].
C              THIS FUNCTION IS INVOKED WHEN DIFFERENTIATION IS PERFORMED
C              USING CHEBYSHEV'S APPROXIMATION 
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
      DOUBLE PRECISION MFVAR
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION BPDFPASS
      COMMON/PPASS/BPDFPASS

      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION IIBPDF_V
      EXTERNAL IIBPDF_V
      DOUBLE PRECISION H_MEAN,H_VAR
          
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA

      DOUBLE PRECISION LIMINF,LIMSUP
      DOUBLE PRECISION H_NEW
      DOUBLE PRECISION ERR,ERR_OLD

      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      INTEGER DIFF_NCHEBCOEFS_TEMP

      MFVAR = MFVARPASS

      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
    
      MFMEANPASS = MFMEAN
C     ADJUST NUMBER OF COEFS IF NECESSARY
      IF(BPDFPASS.LE.DIFF_PDFMIN) THEN
         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS+DIFF_NCHEBCOEFS_EXTRA
      ELSE
         DIFF_NCHEBCOEFS_TEMP = DIFF_NCHEBCOEFS
      ENDIF 
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS_TEMP,IIBPDF_V,
     $     MFVAR,1,IIBPDF_MV)
      RETURN
      END
