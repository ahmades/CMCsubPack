C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CG_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      IMPLICIT NONE
      INCLUDE 'formats.h'
      INTEGER NETA,IETA
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
      DOUBLE PRECISION ETA(NETA),PDF(NETA)
      DOUBLE PRECISION G(NETA-2),GM1,GM2,ETA_TEMP(NETA-2)
      DOUBLE PRECISION BOUND,DUM
      INTEGER STATUS
      DOUBLE PRECISION AG,AREA,A1A2FLAG,A1,A2,FRAC
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION NORMDIST,CDFNORMDIST,TRAP
      EXTERNAL NORMDIST,CDFNORMDIST,TRAP    
      LOGICAL CHECK
      INTEGER NNLEQ,LDFJAC,INFO,LWA,IFLAG
      PARAMETER(NNLEQ = 2,LWA = (NNLEQ*(NNLEQ+13))/2+1,
     $     LDFJAC = NNLEQ)
      DOUBLE PRECISION FNORM
      DOUBLE PRECISION X(NNLEQ),FVEC(NNLEQ),FJAC(NNLEQ,NNLEQ),WA(LWA)
      DOUBLE PRECISION ENORM,DPMPAR
      EXTERNAL FCGPDF
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      INTEGER NLS_METH
      COMMON/NONLINEARSOLVERVARS1/NLS_METH
      DOUBLE PRECISION NLS_TOL
      COMMON/NONLINEARSOLVERVARS2/NLS_TOL

C     FIND THE MEAN AND VARIANCE OF THE CGPDF
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
    
C     G
      DO IETA = 2,NETA-1
         G(IETA-1) = NORMDIST(ETA(IETA),MFMEAN_G,MFVAR_G)
         ETA_TEMP(IETA-1) = ETA(IETA)
      ENDDO    
      
C     STRENGTH OF DELTA FUNCTION AT 0
C     GM1 = INT_-INF^0 [G(ETA)DETA]

      GM1 = CDFNORMDIST(D_ZERO,MFMEAN_G,MFVAR_G)
C      CALL CDFNOR(1,GM1,DUM,D_ZERO,MFMEAN_G,DSQRT(MFVAR_G),STATUS,BOUND)
      
C     STRENGTH OF DELTA FUNCTION AT 1
C     GM2 = INT_1^+INF [G(ETA)DETA]
C     USE THE FACT THAT INT_-INF^+INF[G(ETA)DETA] = 1
C     INT_-INF^1[G(ETA)DETA] + INT_1^+INF[G(ETA)DETA] = 1 =>
C     GM2 = INT_1^+INF[G(ETA)DETA] = 1 - INT_-INF^1[G(ETA)DETA]

      GM2 = D_ONE - CDFNORMDIST(D_ONE,MFMEAN_G,MFVAR_G)
C      CALL CDFNOR(1,GM2,DUM,D_ONE,MFMEAN_G,DSQRT(MFVAR_G),STATUS,BOUND)
C      GM2 = D_ONE-GM2

C     CONSTRUCT THE PDF        
C      IF(((GM1.EQ.D_ZERO).AND.(GM2.EQ.D_ZERO))
C     $     .OR.(GM1.LT.D_ZERO).OR.(GM2.LT.D_ZERO)) THEN 
      IF((GM1.EQ.D_ZERO).AND.(GM2.EQ.D_ZERO)) THEN
         A1 = D_ZERO
         A2 = D_ZERO
         A1A2FLAG = D_ZERO
      ELSE
         AG = TRAP(G,ETA_TEMP,NETA-2)
         AREA = D_ONE - AG        
         FRAC = GM1/(GM1+GM2)
         A1 = FRAC*AREA
         A2 = (D_ONE-FRAC)*AREA
         A1A2FLAG = D_ONE
      ENDIF
      
      DO IETA = 1,NETA
         IF(IETA.EQ.1) THEN
            PDF(IETA) = D_TWO*A1/(ETA(2)-ETA(1)) - A1A2FLAG*G(1)
         ELSEIF(IETA.EQ.NETA) THEN
            PDF(IETA) = D_TWO*A2/(ETA(NETA)-ETA(NETA-1)) 
     $           - A1A2FLAG*G(NETA-2)
         ELSE
            PDF(IETA) = G(IETA-1)
         ENDIF
      ENDDO 
           
      RETURN
      END


C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE CG_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,VEL,CV)
C========================================================================================== 
C     PURPOSE: COMPUTES THE CONDITIONAL VELOCITY USING THE PDF-GRADIENT 
C              MODEL BASED CG-PDF
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
      DOUBLE PRECISION ETA(NETA),CV(3,NETA),PDF(NETA)  
      DOUBLE PRECISION M,S,S2,MFMEAN_G,MFVAR_G
      DOUBLE PRECISION DPDM(NETA),DPDS2(NETA)
      DOUBLE PRECISION DMDMEAN,DMDVAR,ERR_DMDMEAN,ERR_DMDVAR
      DOUBLE PRECISION DS2DMEAN,DS2DVAR,ERR_DS2DMEAN,ERR_DS2DVAR
      DOUBLE PRECISION S2_V,S2_M,M_V,M_M
      EXTERNAL S2_V,S2_M,M_V,M_M
      DOUBLE PRECISION H_NEW,ERR,ERR_OLD,H_MEAN,H_VAR 
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION LIMINFM,LIMSUPM
      DOUBLE PRECISION LIMINFV,LIMSUPV
           
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION EPSL
      COMMON/EPSILON/EPSL

C     THESE BLOCKS ARE PASSED TO EXTERNAL FUNCTIONS
      DOUBLE PRECISION MFMEANPASS_CG
      COMMON/MFMEANBLOK_CG/MFMEANPASS_CG
      DOUBLE PRECISION MFVARPASS_CG
      COMMON/MFVARBLOK_CG/MFVARPASS_CG
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS

      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      

      CALL CHECKPARMS(MFMEAN,MFVAR)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL CG_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
      M  = MFMEAN_G
      S  = DSQRT(MFVAR_G)
      S2 = MFVAR_G

      DO IETA = 2,NETA-1
C
C     DPDM
         DPDM(IETA) = (D_ONE/(S*DSQRT(D_TWO*PI)))*((ETA(IETA)-M)/S2)
     $        * DEXP(-((ETA(IETA)-M)/(S*DSQRT(D_TWO)))**D_TWO)
C     
C     DPDS2
         DPDS2(IETA) = -(D_ONE/(D_TWO*DSQRT(D_TWO*PI)))
     $        * (D_ONE/(S**D_THREE))*(D_ONE-((ETA(IETA)-M)/S)**D_TWO)
     $        * DEXP(-((ETA(IETA)-M)/(S*DSQRT(D_TWO)))**D_TWO)
C     
      ENDDO

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
      LIMINFV = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUPV = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      LIMINFM = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUPM = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
C
C     DS2DVAR
      MFMEANPASS_CG = MFMEAN
      CALL CHEBDERIV(LIMINFV,LIMSUPV,DIFF_NCHEBCOEFS,S2_V,
     $     MFVAR,1,DS2DVAR)
C     
C     DS2DMEAN
      MFVARPASS_CG = MFVAR
      CALL CHEBDERIV(LIMINFM,LIMSUPM,DIFF_NCHEBCOEFS,S2_M,
     $     MFMEAN,1,DS2DMEAN)        
C     
C     DMDVAR
      MFMEANPASS_CG = MFMEAN
      CALL CHEBDERIV(LIMINFV,LIMSUPV,DIFF_NCHEBCOEFS,M_V,
     $     MFVAR,1,DMDVAR)
C     
C     DMDMEAN
      MFVARPASS_CG = MFVAR
      CALL CHEBDERIV(LIMINFM,LIMSUPM,DIFF_NCHEBCOEFS,M_M,
     $     MFMEAN,1,DMDMEAN)
C
C     COMPUTE THE CONDITIONAL VELOCITY
      DO IETA = 2,NETA-1
         DO I = 1,3
            CV(I,IETA) =  VEL(I) - (DT/ PDF(IETA))
     $           *((DPDS2(IETA)*DS2DMEAN + DPDM(IETA)*DMDMEAN)
     $           * MFMEANGRAD(I) 
     $           + (DPDS2(IETA)*DS2DVAR  + DPDM(IETA)*DMDVAR)
     $           * MFVARGRAD(I))
         ENDDO
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'CG: CONDITIONAL VELOCITY:'
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

      SUBROUTINE CG_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [HOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE CG-PDF
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
      DOUBLE PRECISION DIIDS2(NETA),DIIDM(NETA),
     $     ERR_DMDV,DMDV,ERR_DS2DV,DS2DV,
     $     MAXERR_DIIDV
      DOUBLE PRECISION H_NEW,ERR,ERR_OLD,H_MEAN,H_VAR    
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION S2_V,M_V,ERF
      EXTERNAL S2_V,M_V,ERF
      DOUBLE PRECISION MFMEANPASS_CG
      COMMON/MFMEANBLOK_CG/MFMEANPASS_CG
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINF,LIMSUP

      DOUBLE PRECISION MFMEAN_G,MFVAR_G,A1,A2,A3,A4,M,S,S2
      
      CALL CHECKPARMS(MFMEAN,MFVAR)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL CG_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      
C     DS2DV AND DMDV
      MFMEANPASS_CG = MFMEAN
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,S2_V,MFVAR,1,
     $     DS2DV)
      MFMEANPASS_CG = MFMEAN
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,M_V,MFVAR,1,
     $     DMDV)
C
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
      M  = MFMEAN_G
      S2 = MFVAR_G
      S  = DSQRT(MFVAR_G)
      DO IETA = 2,NETA-1
         DIIDS2(IETA) = D_ONE/(D_TWO*DSQRT(D_TWO*PI)*S)
     $        * (DEXP(-(((ETA(IETA)-M)/(S*DSQRT(D_TWO)))**D_TWO))
     $        - (D_ONE + ETA(IETA)*M/S2)
     $        * DEXP(-(M/(S*DSQRT(D_TWO)))**D_TWO))
         DIIDM(IETA) = -D_HALF*(ERF((ETA(IETA)-M)/(S*DSQRT(D_TWO)))
     $        + ERF(M/(S*DSQRT(D_TWO))))
     $        + ETA(IETA)/(DSQRT(D_TWO*PI)*S)
     $        * DEXP(-((M/(S*DSQRT(D_TWO)))**D_TWO))
      ENDDO
     
C     BOUNDARY VALUES ARE KNOWN
      CSDRH(1)    = D_ZERO
      CSDRH(NETA) = D_ZERO
C     COMPUTE CSDR AT INTERNAL GRID POINTS     
      DO IETA = 2,NETA-1
         CSDRH(IETA) = (D_TWO/PDF(IETA))*CHI
     $        *(DIIDS2(IETA)*DS2DV + DIIDM(IETA)*DMDV)
      ENDDO

      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE CG_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,
     $     MFVARGRAD,DT,CHI,CSDRI)
C==========================================================================================
C     PURPOSE: COMPUTES THE [INHOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
C              DISSIPATION RATE MODEL USING THE CG-PDF
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

      DOUBLE PRECISION IICGPDF_M,IICGPDF_V,IICGPDF_MV
      EXTERNAL IICGPDF_M,IICGPDF_V,IICGPDF_MV

      DOUBLE PRECISION MFMEANPASS_CG
      COMMON/MFMEANBLOK_CG/MFMEANPASS_CG
      DOUBLE PRECISION MFVARPASS_CG
      COMMON/MFVARBLOK_CG/MFVARPASS_CG
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINF,LIMSUP

      CALL CHECKPARMS(MFMEAN,MFVAR)

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL CG_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)    
C
C     DIIDV AND D2IIDV2
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      MFMEANPASS_CG = MFMEAN
      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,IICGPDF_V,
     $        MFVAR,1,DIIDV(IETA))
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,IICGPDF_V,
     $        MFVAR,2,D2IIDV2(IETA))
      ENDDO
C     
C     D2IIDM2
      LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
      MFVARPASS_CG = MFVAR
      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,IICGPDF_M,
     $        MFMEAN,2,D2IIDM2(IETA))
      ENDDO
      
      WRITE(*,*)'MORE WORK NEEDED FOR IN CG_CSDR_I'
      WRITE(*,*)'A STOP WAS INSERTED HERE'
      STOP

      
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
         WRITE(*,200)'INDEX','ETA','CSDRH'
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

      SUBROUTINE FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
C================================================================================== 
C     PURPOSE: COMPUTES THE MEAN AND VARIANCE OF THE CGPDF GIVEN THE MIXTURE
C              FRACTION MEAN AND VARIANCE
C==================================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            MMEAN     MIXTURE FRACTION MEAN          DOUBLE PRECISION 
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION 
C
C     OUTPUT:
C
C            MMEAN_G   MEAN OF THE CGPDF              DOUBLE PRECISION 
C            MFVAR_G   VARIANCE OF THE CGPDF          DOUBLE PRECISION 
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
C     NON-LINEAR SOLVER VARIABLES
      LOGICAL CHECK
      INTEGER NNLEQ,LDFJAC,INFO,LWA,IFLAG
      PARAMETER(NNLEQ = 2,LWA = (NNLEQ*(NNLEQ+13))/2+1,
     $     LDFJAC = NNLEQ)
      DOUBLE PRECISION FNORM
      DOUBLE PRECISION X(NNLEQ),FVEC(NNLEQ),FJAC(NNLEQ,NNLEQ),WA(LWA)
      DOUBLE PRECISION ENORM,DPMPAR
      EXTERNAL FCGPDF
      DOUBLE PRECISION NLS_TOL
      COMMON/NONLINEARSOLVERVARS/NLS_TOL

      MFMEANPASS = MFMEAN
      MFVARPASS  = MFVAR

C     INITIAL GUESS
C      X(1) = MFMEAN
C      X(2) = DSQRT(MFVAR)

      X(1) = MFMEAN
      X(2) = DSQRT(MFVAR)

C     COMPUTE MFMEAN_G AND MFVAR_G
C      IFLAG = 1
C      CALL FCGPDF(NNLEQ,X,FVEC,FJAC,LDFJAC,IFLAG)
      CALL HYBRJ1(FCGPDF,NNLEQ,X,FVEC,FJAC,LDFJAC,NLS_TOL,
     $     INFO,WA,LWA)
         
      MFMEAN_G = X(1)
      MFVAR_G  = X(2)*X(2)

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE FCGPDF(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C================================================================================== 
C     PURPOSE: COMPUTES THE THE FUNCTIONS AND THE JACOBIAN MATRIX REQUIRED BY THE
C              THE NON-LINERA EQUATION SOLVER HYBRJ1 FOR THE CALCULATION OF THE
C              MEAN AND VARIANCE OF THE CGPDF
C==================================================================================
C            VARIABLE  DESCRIPTION                 DATA TYPE
C            --------  -----------                 --------- 
C
C     INPUT:
C
C            N         NUMBER OF VARIABLES         INTEGER
C            X         X(1) = MEAN OF THE CGPDF    DOUBLE PRECISION (ARRAY,SIZE=2)
C                      X(2) = STANDARD DEVIATION
C                             OF THE CGPDF
C            FVEC      VECTOR CONTAINING THE       DOUBLE PRECISION (ARRAY,SIZE=2)
C                      SYSTEM OF NON-LINEAR
C                      EQUATIONS
C            FJAC      JACOBIAN MATRIX OF THE      DOUBLE PRECISION (ARRAY,SIZE=2:2)
C                      SYSTEM OF NON-LINEAR
C                      EQUATIONS
C            LDFJAC    LEADING DIMENSION OF FJAC   INTEGER (LDFJAC.GT.N)
C            IFLAG     CONTROL FLAG                INTEGER
C                      IFLAG = 1: CALCULATE THE 
C                        FUNCTIONS AT X AND RETURN 
C                        THE VECTOR IN FVEC.
C                      IFLAG = 2: CALCULATE THE 
C                        JACOBIAN AT X AND RETURN 
C                        THE MATRIX IN FJAC.
C================================================================================== 
      INTEGER N,LDFJAC,IFLAG
      DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
      IF (IFLAG.EQ.1) CALL CGPDFVEC(X,FVEC)
      IF (IFLAG.EQ.2) CALL CGPDFJAC(X,FJAC)
      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CGPDFVEC(X,FVEC)
C================================================================================== 
C     PURPOSE: SETS THE SYSTEM OF NON-LINEAR EQUATION REQUIRED FOR THE COMPUTATION
C              OF THE MEAN AND THE VARIANCE OF THE CGPDF
C==================================================================================
C            VARIABLE  DESCRIPTION                  DATA TYPE
C            --------  -----------                  --------- 
C
C     INPUT:
C
C            X         X(1) = MEAN OF THE CGPDF     DOUBLE PRECISION (ARRAY,SIZE=2)
C                      X(2) = STANDARD DEVIATION
C                             OF THE CGPDF
C
C     OUTPUT:
C
C            FVEC      VECTOR CONTAINING THE       DOUBLE PRECISION (ARRAY,SIZE=2)
C                      SYSTEM OF NON-LINEAR
C                      EQUATIONS
C================================================================================== 
      IMPLICIT NONE
      INTEGER N
      PARAMETER(N = 2)
      DOUBLE PRECISION X(N),FVEC(N)
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR

      DOUBLE PRECISION M,S,Z0,Z1,PZ0,PZ1

      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      
      DOUBLE PRECISION ERF
      EXTERNAL ERF

      M = X(1)
      S = X(2)

      Z0 = -M/S
      Z1 = (D_ONE-M)/S
      
      PZ0 = D_HALF * (D_ONE + ERF(Z0/DSQRT(D_TWO)))
      PZ1 = D_HALF * (D_ONE + ERF(Z1/DSQRT(D_TWO)))
      
      FVEC(1) = D_ONE + (M-D_ONE)*PZ1 - M*PZ0 + (S/(2.0D0*PI))
     $     * (DEXP(-(Z0**D_TWO)/2.0D0) - DEXP(-(Z1**D_TWO)/D_TWO)) 
     $     - MFMEANPASS
         
      FVEC(2) = D_ONE + ((S**D_TWO)+(M**D_TWO)-D_ONE)*PZ1
     $     - ((S**D_TWO)+(M**D_TWO))*PZ0 
     $     + ((S**D_TWO)/DSQRT(D_TWO*PI))
     $     * (Z0*DEXP(-(Z0**D_TWO)/D_TWO) 
     $     - Z1*DEXP(-(Z1**D_TWO)/D_TWO))
     $     + ((D_TWO*S*M)/DSQRT(D_TWO*PI))
     $     * (DEXP(-(Z0**D_TWO)/D_TWO) 
     $     - DEXP(-(Z1**D_TWO)/D_TWO)) 
     $     - (MFMEANPASS**D_TWO) - MFVARPASS
      
      RETURN
      END

C==========================================================================================
C=========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CGPDFJAC(X,FJAC)
C================================================================================== 
C     PURPOSE: COMPUTES THE (ANALYTICAL) JACOBIAN MATRIX OF THE SYSTEM OF NON-
C              LINEAR EQUATION REQUIRED FOR THE COMPUTATION OF THE MEAN AND THE 
C              VARIANCE OF THE CGPDF
C==================================================================================
C            VARIABLE  DESCRIPTION                 DATA TYPE
C            --------  -----------                 --------- 
C
C     INPUT:
C
C            X         X(1) = MEAN OF THE CGPDF    DOUBLE PRECISION (ARRAY,SIZE=2)
C                      X(2) = STANDARD DEVIATION
C                             OF THE CGPDF
C
C     OUTPUT:
C
C            FJAC      JACOBIAN MATRIX OF THE      DOUBLE PRECISION (ARRAY,SIZE=2:2)
C                      SYSTEM OF NON-LINEAR
C                      EQUATIONS
C================================================================================== 
      IMPLICIT NONE
      INTEGER N
      PARAMETER(N = 2)
      DOUBLE PRECISION X(N),FVEC(N),FJAC(N,N)    
      DOUBLE PRECISION M,S,Z0,Z1,PZ0,PZ1
      DOUBLE PRECISION PZ0_M,PZ1_M,PZ0_S,PZ1_S

      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
     
      DOUBLE PRECISION ERF
      EXTERNAL ERF

      M = X(1)
      S = X(2)
      
      Z0 = -M/S
      Z1 = (D_ONE-M)/S

C      PZ  = D_HALF * (D_ONE + ERF(Z/DSQRT(D_TWO)))
      PZ0 = D_HALF * (D_ONE + ERF(Z0/DSQRT(D_TWO)))
      PZ1 = D_HALF * (D_ONE + ERF(Z1/DSQRT(D_TWO)))

      PZ0_M = -(D_ONE/(S*DSQRT(D_TWO*PI)))*DEXP(-(Z0**D_TWO)/D_TWO)
      PZ1_M = -(D_ONE/(S*DSQRT(D_TWO*PI)))*DEXP(-(Z1**D_TWO)/D_TWO)
      PZ0_S = (M/((S**D_TWO)*DSQRT(D_TWO*PI)))
     $     *DEXP(-(Z0**D_TWO)/D_TWO)
      PZ1_S = -((D_ONE-M)/((S**D_TWO)*DSQRT(D_TWO*PI)))
     $     *DEXP(-(Z1**D_TWO)/D_TWO)
      
      FJAC(1,1) = (M-D_ONE)*PZ1_M + PZ1 - M*PZ0_M -PZ0
     $     - (D_ONE/(S*DSQRT(D_TWO*PI)))*
     $     (M*DEXP(-(Z0**D_TWO)/2.0D0) 
     $     + (D_ONE-M)*DEXP(-(Z1**D_TWO)/2.0D0))
      
      FJAC(1,2) = (M-D_ONE)*PZ1_S - M*PZ0_S
     $     + (D_ONE/((S**D_TWO)*DSQRT(D_TWO*PI)))
     $     *((M**D_TWO)*DEXP(-(Z0**D_TWO)/D_TWO)
     $     + (D_ONE-M)*DEXP(-(Z1**D_TWO)/D_TWO))
     $     + (D_ONE/DSQRT(D_TWO*PI))
     $     *(DEXP(-(Z0**D_TWO)/D_TWO)
     $     + DEXP(-(Z1**D_TWO)/D_TWO))
      
      
      FJAC(2,1) = ((S**D_TWO)+(M**D_TWO)-D_ONE)*PZ1_M + D_TWO*M*PZ1 
     $     - ((S**D_TWO)+(M**D_TWO))*PZ0_M - D_TWO*M*PZ0
     $     + (S/DSQRT(D_TWO*PI))
     $     *(-DEXP(-(Z0**D_TWO)/D_TWO)*(D_ONE+(Z0**2.0D0))
     $     + DEXP(-(Z1**D_TWO)/D_TWO)*(D_ONE-(Z1**2.0D0)))
     $     + (D_TWO*M/(S*DSQRT(D_TWO*PI)))
     $     * (M*DEXP(-(Z0**D_TWO)/D_TWO) 
     $     + (D_ONE-M)*DEXP(-(Z1**D_TWO)/D_TWO))
     $     + (D_TWO*S/(DSQRT(D_TWO*PI)))
     $     * (DEXP(-(Z0**D_TWO)/D_TWO)
     $     + DEXP(-(Z1**D_TWO)/D_TWO))
      
      FJAC(2,2) = ((S**D_TWO)+(M**D_TWO)-D_ONE)*PZ1_S + D_TWO*S*PZ1 
     $     - ((S**D_TWO)+(M**D_TWO))*PZ0_S - D_TWO*S*PZ0
     $     + (D_ONE/DSQRT(D_TWO*PI))
     $     * (M*(D_ONE-(Z0**D_TWO))*DEXP(-(Z0**D_TWO)/D_TWO)
     $     + (D_ONE-M)*(D_ONE-(Z1**D_TWO))*DEXP(-(Z1**D_TWO)/D_TWO))        
     $     + (D_TWO*S/DSQRT(D_TWO*PI))
     $     * (Z0*DEXP(-(Z0**D_TWO)/D_TWO) 
     $     + Z1*DEXP(-(Z1**D_TWO)/D_TWO))
     $     + (D_TWO*M/(DSQRT(D_TWO*PI)))
     $     * (DEXP(-(Z0**D_TWO)/D_TWO)
     $     + DEXP(-(Z1**D_TWO)/D_TWO))
     $     + (D_TWO*M/DSQRT(D_TWO*PI))
     $     * ((Z0**D_TWO)*DEXP(-(Z0**D_TWO)/D_TWO) 
     $     - (Z1**D_TWO)*DEXP(-(Z1**D_TWO)/D_TWO))
      
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION IICGPDF_M(MFMEAN)
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
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR
      DOUBLE PRECISION MFVARPASS_CG
      COMMON/MFVARBLOK_CG/MFVARPASS_CG
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION ERF
      EXTERNAL ERF
      DOUBLE PRECISION Q1,Q2
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
C     MFMEAN IS PROVIDED IN THE ARGUMENT OF FUNCTION IICGPDF_M
C     MFVAR IS PASSED FROM SUBROUTINE CG_CSDR (FIXED)
      MFVAR = MFVARPASS_CG  
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)      
      Q1 = (ETAPASS-MFMEAN_G)/DSQRT(D_TWO*MFVAR_G)
      Q2 = MFMEAN_G/DSQRT(D_TWO*MFVAR_G)
      IICGPDF_M = D_HALF*(ETAPASS-MFMEAN_G)*(ERF(Q1)+ERF(Q2)) 
     $   + DSQRT(MFVAR_G/(D_TWO*PI))
     $     *(DEXP(-(Q1**D_TWO))-DEXP(-(Q2**D_TWO)))
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION IICGPDF_V(MFVAR)
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
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR
      DOUBLE PRECISION MFMEANPASS_CG
      COMMON/MFMEANBLOK_CG/MFMEANPASS_CG
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION ERF
      EXTERNAL ERF
      DOUBLE PRECISION Q1,Q2
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
C     MFVAR IS PROVIDED IN THE ARGUMENT OF FUNCTION IICGPDF_V
C     MFMEAN IS PASSED FROM SUBROUTINE CG_CSDR (FIXED)
      MFMEAN = MFMEANPASS_CG  
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)  
      Q1 = (ETAPASS-MFMEAN_G)/DSQRT(D_TWO*MFVAR_G)
      Q2 = MFMEAN_G/DSQRT(D_TWO*MFVAR_G)
      IICGPDF_V = D_HALF*(ETAPASS-MFMEAN_G)*(ERF(Q1)+ERF(Q2)) 
     $     + DSQRT(MFVAR_G/(D_TWO*PI))
     $     *(DEXP(-(Q1**D_TWO))-DEXP(-(Q2**D_TWO)))

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION IICGPDF_MV(MFMEAN,MFVAR)
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
C            MFVAR     MIXTURE FRACTION VARIANCE      DOUBLE PRECISION
C================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION ERF
      EXTERNAL ERF
      DOUBLE PRECISION Q1,Q2
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
C     NOTHING TO BE TRANSFERRED HERE:
C     MFMEAN AND MFVAR ARE PROVIDED IN THE ARGUMENT OF IICGPDF_MV
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
      Q1 = (ETAPASS-MFMEAN_G)/DSQRT(D_TWO*MFVAR_G)
      Q2 = MFMEAN_G/DSQRT(D_TWO*MFVAR_G)
      IICGPDF_MV= D_HALF*(ETAPASS-MFMEAN_G)*(ERF(Q1)+ERF(Q2))
     $     + DSQRT(MFVAR_G/(D_TWO*PI))
     $     *(DEXP(-(Q1**D_TWO))-DEXP(-(Q2**D_TWO)))
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION IICGPDF_MV2(MFMEAN)
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
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR
      DOUBLE PRECISION MFVARPASS_CG
      COMMON/MFVARBLOK_CG/MFVARPASS_CG
      DOUBLE PRECISION MFMEANPASS_CG
      COMMON/MFMEANBLOK_CG/MFMEANPASS_CG
      DOUBLE PRECISION ETAPASS 
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION ERF,IICGPDF_V
      EXTERNAL ERF,IICGPDF_V
      DOUBLE PRECISION Q1,Q2
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION H_MEAN,H_VAR,LIMINF,LIMSUP
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
C     MFMEAN IS PROVIDED IN THE ARGUMENT OF FUNCTION IICGPDF_M
C     MFVAR IS PASSED FROM SUBROUTINE CG_CSDR (FIXED)
      MFVAR = MFVARPASS_CG  

      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
    
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)

      MFMEANPASS_CG = MFMEAN
      
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,IICGPDF_V,
     $           MFVAR,1,IICGPDF_MV2)

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION S2_V(MFVAR)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST-ORDER 
C              PARTIAL DERIVATIVE OF THE VARIANCE OF THE CG-PDF WITH RESPECT TO THE 
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
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR_X
      DOUBLE PRECISION MFMEANPASS_CG
      COMMON/MFMEANBLOK_CG/MFMEANPASS_CG
C     MFVAR IS PROVIDED IN THE ARGUMENT OF FUNCTION S2_V
C     MFMEAN IS PASSED FROM SUBROUTINE (FIXED)
      MFMEAN = MFMEANPASS_CG      
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
      S2_V = MFVAR_G
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION S2_M(MFMEAN)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST-ORDER 
C              PARTIAL DERIVATIVE OF THE VARIANCE OF THE CG-PDF WITH RESPECT TO THE 
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
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR_X
      DOUBLE PRECISION MFVARPASS_CG
      COMMON/MFVARBLOK_CG/MFVARPASS_CG
C     MFMEAN IS PROVIDED IN THE ARGUMENT OF FUNCTION S2_M
C     MFVAR IS PASSED FROM SUBROUTINE (FIXED)
      MFVAR = MFVARPASS_CG      
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
      S2_M = MFVAR_G
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION M_V(MFVAR)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST-ORDER 
C              PARTIAL DERIVATIVE OF THE MEAN OF THE CG-PDF WITH RESPECT TO THE 
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
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR_X
      DOUBLE PRECISION MFMEANPASS_CG
      COMMON/MFMEANBLOK_CG/MFMEANPASS_CG
C     MFVAR IS PROVIDED IN THE ARGUMENT OF FUNCTION M_V
C     MFMEAN IS PASSED FROM SUBROUTINE (FIXED)
      MFMEAN = MFMEANPASS_CG      
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
      M_V = MFMEAN_G
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION M_M(MFMEAN)
C================================================================================== 
C     PURPOSE: FUNCTION NECESSARY FOR THE COMPUTATION OF THE FIRST-ORDER 
C              PARTIAL DERIVATIVE OF THE MEAN OF THE CG-PDF WITH RESPECT TO THE 
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
      DOUBLE PRECISION MFMEAN,MFVAR,MFMEAN_G,MFVAR_G
C     THESE BLOCKS ARE PASSED FROM SUBROUTINE CG_CSDR_X
      DOUBLE PRECISION MFVARPASS_CG
      COMMON/MFVARBLOK_CG/MFVARPASS_CG
C     MFMEAN IS PROVIDED IN THE ARGUMENT OF FUNCTION M_M
C     MFVAR IS PASSED FROM SUBROUTINE (FIXED)
      MFVAR = MFVARPASS_CG      
      CALL FIND_CGPDF_PARMS(MFMEAN,MFVAR,MFMEAN_G,MFVAR_G)
      M_M = MFMEAN_G
      RETURN
      END
