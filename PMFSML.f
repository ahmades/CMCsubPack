C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE PMFSML_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
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

      DOUBLE PRECISION G,V,W
      DOUBLE PRECISION ETA_TEMP1(NETA-1),PDF_TEMP1(NETA-1)
      DOUBLE PRECISION PDF_AREA,PDF_LA,PDF_RA,PDF_SCALEFACT
      
      DOUBLE PRECISION ERFINV,TRAP,CDFNORMDIST
      EXTERNAL ERFINV,TRAP,CDFNORMDIST
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION TAUPASS
      COMMON/COMVARSTAU/TAUPASS

      CALL CHECKPARMS(MFMEAN,MFVAR)

      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)

      CALL FIND_TAU(MFMEAN,MFVAR,TAU)
     
      SIGMA2 = D_ONE - D_TWO*TAU

      DO IETA = 1,NETA
         E(IETA) = ERFINV(D_TWO*ETA(IETA)-D_ONE)
         PHI(IETA) = ALPHA + D_TWO*DSQRT(TAU)*E(IETA) 
         PDF(IETA) = DSQRT(D_TWO*TAU/SIGMA2)
     $        * DEXP(E(IETA)**D_TWO - (PHI(IETA)**D_TWO)/(D_TWO*SIGMA2))
      ENDDO

C     DETERMINE THE SHAPE OF THE PDF USING THE FIRST AND LAST TWO INTERIOR POINTS OF THE
C     PDF ARRAY IN ORDER TO SET BOUNDARY VALUES. THE PROCEDURE MAKES USE OF THE FACT
C     THAT R(PHI) IS THE CUMULATIVE DISTRIBUTION FUNCTION OF P(ETA).
C     R(PHI) IS THE CUMULATIVE NORMAL DISTRIBUTION WITH ZERO MEAN AND SIGMA2 VARIANCE
C     AND PHI IS RELATED TO ETA BY: PHI = ALPHA + 2*DSQRT(TAU)*ERFINV(2*ETA-1). 
C     R(PHI) IS USED TO COMPUTE THE AREAS UNDER THE PDF EXTENDING FROM ETA(1) TO ETA(2)
C     AND FROM ETA(NETA-1) TO ETA(NETA). THESE AREARS ARE THEN USED TO SET THE BOUNDARY 
C     VALUES OF THE PDF, PDF(1) AND PDF(NETA).
C
C     USE R(PHI) TO FIND THE AREA BETWEEN THE LEFT BOUNADRY AND THE FIRST INTERIOR GRID POINT
         PDF_LA = CDFNORMDIST(PHI(2),D_ZERO,SIGMA2)
         PDF(1) = D_TWO*PDF_LA/(ETA(2)-ETA(1)) - PDF(2)
         PDF(1) = DMAX1(PDF(1),D_ZERO)
C     USE R(PHI) TO FIND THE AREA BETWEEN THE LEFT BOUNADRY AND THE LAST INTERIOR GRID POINT
         PDF_RA = D_ONE - CDFNORMDIST(PHI(NETA-1),D_ZERO,SIGMA2)
         PDF(NETA) = D_TWO*PDF_RA/(ETA(NETA)-ETA(NETA-1)) - PDF(NETA-1)
         PDF(NETA) = DMAX1(PDF(NETA),D_ZERO)
C
      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF-B: PROBABILITY DENSITY FUNCTION:'
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

      SUBROUTINE PMFSML_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,VEL,CV)
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
      DOUBLE PRECISION MFMEAN,MFVAR,ALPHA,TAU,SIGMA2,DT
      DOUBLE PRECISION MFMEANGRAD(3),MFVARGRAD(3),VEL(3)
      DOUBLE PRECISION ETA(NETA),CV(3,NETA),E(NETA),PHI(NETA)
      DOUBLE PRECISION ALPHA_M,TAU_M,TAU_V
      DOUBLE PRECISION ERR1,ERR2
      DOUBLE PRECISION ERR_OLD, H_NEW,H_MEAN,H_VAR
      DOUBLE PRECISION ERFINV,TAUFUNC_M,TAUFUNC_V
      EXTERNAL ERFINV,TAUFUNC_M,TAUFUNC_V  
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION MFMEANPASS_DERIVS,MFVARPASS_DERIVS,ALPHAPASS
      COMMON/COMVARSMEAN/MFMEANPASS_DERIVS
      COMMON/COMVARSVAR/MFVARPASS_DERIVS
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINF,LIMSUP
   
      CALL CHECKPARMS(MFMEAN,MFVAR)

      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)

      ALPHA_M = -DSQRT(D_TWO*PI)*DEXP((ALPHA**D_TWO)/D_TWO)

      CALL FIND_TAU(MFMEAN,MFVAR,TAU)

      SIGMA2 = D_ONE - D_TWO*TAU


C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)     
C
C     TAU_M
      LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
      MFVARPASS_DERIVS = MFVAR
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,TAUFUNC_M,
     $     MFMEAN,1,TAU_M)
C     
C     TAU_V
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      MFMEANPASS_DERIVS = MFMEAN
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,TAUFUNC_V,
     $     MFVAR,1,TAU_V)
C     
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
         WRITE(*,50) 'PMF-B: CONDITIONAL VELOCITY:'
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

      SUBROUTINE PMFSML_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH)
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
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,ALPHA,TAU,SIGMA2
      DOUBLE PRECISION ETA(NETA),CSDRH(NETA),E(NETA),PHI(NETA)  
      DOUBLE PRECISION ERFINV
      EXTERNAL ERFINV
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
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
         CSDRH(IETA) =CHI * DSQRT((D_ONE-TAU)/TAU)
     $        * DEXP(-D_TWO*(E(IETA)**D_TWO)
     $        + (ALPHA**D_TWO)/(SIGMA2+D_ONE))
c         CSDRH(IETA) = DMAX1(CSDRH(IETA),D_ZERO)
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF-B: COND. SCAL. DISS. RATE (HOMOGENEOUS):'
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

      SUBROUTINE PMFSML_CSDR_I(ETA,NETA,MFMEAN,MFVAR,
     $     MFMEANGRAD,MFVARGRAD,DT,CHI,CSDRI)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [INHOMOGENEOUS] VERSION OF THE CONDITIONAL SCALAR 
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
      DOUBLE PRECISION PRODM,PRODV,PRODMV
      DOUBLE PRECISION T1,T2,T3,T4,T5,T6,FACT
      DOUBLE PRECISION ERR1,ERR2,ERR3,ERR4,ERR_TEMP,INT_FACT
      DOUBLE PRECISION ERR_OLD, H_NEW,H_MEAN,H_VAR
      DOUBLE PRECISION ERFINV,TAUFUNC_M,TAUFUNC_V
      EXTERNAL ERFINV,TAUFUNC_M,TAUFUNC_V  
      DOUBLE PRECISION MFMEANPASS_DERIVS,MFVARPASS_DERIVS,ALPHAPASS
      COMMON/COMVARSMEAN/MFMEANPASS_DERIVS
      COMMON/COMVARSVAR/MFVARPASS_DERIVS
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DOUBLE PRECISION MFVAR_MAX
      DOUBLE PRECISION LIMINF,LIMSUP

      CALL CHECKPARMS(MFMEAN,MFVAR)

C     ESTIMATION OF THE INITIAL STEPSIZE, H, IS ADOPTED FROM REF [2]
      CALL FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
      MFVAR_MAX = MFMEAN*(D_ONE-MFMEAN)
C
C     TAU_M
      LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
      MFVARPASS_DERIVS = MFVAR
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,TAUFUNC_M,
     $     MFMEAN,1,TAU_M)
C     
C     TAU_MM
      LIMINF = DMAX1(MFMEAN-H_MEAN,MFMEAN_MIN)
      LIMSUP = DMIN1(MFMEAN+H_MEAN,MFMEAN_MAX)
      MFVARPASS_DERIVS = MFVAR
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,TAUFUNC_M,
     $     MFMEAN,2,TAU_MM)
C     
C     TAU_V
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      MFMEANPASS_DERIVS = MFMEAN
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,TAUFUNC_V,
     $     MFVAR,1,TAU_V)
C     
C     TAU_VV
      LIMINF = DMAX1(MFVAR-H_VAR,MFVAR_MIN)
      LIMSUP = DMIN1(MFVAR+H_VAR,MFVAR_MAX)
      MFMEANPASS_DERIVS = MFMEAN
      CALL CHEBDERIV(LIMINF,LIMSUP,DIFF_NCHEBCOEFS,TAUFUNC_V,
     $     MFVAR,2,TAU_VV)
C
      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)

      ALPHA_M = -DSQRT(D_TWO*PI)*DEXP((ALPHA**D_TWO)/D_TWO)

      CALL FIND_TAU(MFMEAN,MFVAR,TAU)

      SIGMA2 = D_ONE - D_TWO*TAU
     
C     BOUNDARY VALUES ARE KNOWN
      CSDRI(1)    = D_ZERO
      CSDRI(NETA) = D_ZERO
      CSDRH(1)    = D_ZERO
      CSDRH(NETA) = D_ZERO
C     COMPUTE CSDR AT INTERNAL GRID POINTS
      PRODM  = D_ZERO    ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. MEAN
      PRODV  = D_ZERO    ! SCALAR PRODUCT OF THE GRADS OF MIX FRAC. VARIANCE
      PRODMV = D_ZERO    ! SCALAR PRODUCT OF MIX FRAC. MEAN AND MIX FRAC VARIANCE GRADS
      DO I = 1,3
         PRODM  = PRODM  + MFMEANGRAD(I)**D_TWO 
         PRODV  = PRODV  + MFVARGRAD(I)**D_TWO
         PRODMV = PRODMV + MFMEANGRAD(I)*MFVARGRAD(I)
      ENDDO
      DO IETA = 2,NETA-1
         E(IETA) =  ERFINV(D_TWO*ETA(IETA)-D_ONE)

         PHI(IETA) = ALPHA + D_TWO*DSQRT(TAU)*E(IETA) 

         A(IETA) = ALPHA/(D_ONE-TAU) - PHI(IETA)/SIGMA2

         B(IETA) = (ALPHA**D_TWO)/(D_TWO*((D_ONE-TAU)**D_TWO))
     $        - (D_ONE/SIGMA2)*(PHI(IETA)*E(IETA)/DSQRT(TAU) !DSQRT(D_TWO*TAU)
     $        + (PHI(IETA)**D_TWO)/SIGMA2
     $        - D_ONE/(D_ONE+SIGMA2))

         C(IETA) = DSQRT(D_TWO/PI) * (ALPHA_M
     $        + TAU_M*PHI(IETA)/(TAU*SIGMA2))
     $        * DEXP(-D_TWO*(E(IETA)**D_TWO) + (ALPHA**D_TWO)/D_TWO)

         CSDRH(IETA) = CHI * DSQRT((D_ONE-TAU)/TAU)
     $        * DEXP(-D_TWO*(E(IETA)**D_TWO)
     $        + (ALPHA**D_TWO)/(SIGMA2+D_ONE))

         T1 = PRODM*(D_TWO + ((TAU_M**D_TWO)*TAU_VV/(TAU_V**D_THREE))
     $        - TAU_MM/TAU_V)
         T2 = D_TWO*ALPHA_M*PRODMV*A(IETA)
         T3 = (PRODV*TAU_V + D_TWO*PRODMV*TAU_M 
     $        + PRODM*(TAU_M**D_TWO)/TAU_V)*B(IETA)
         T4 = D_TWO*DT*PRODM*C(IETA)

c$$$         CSDRI(IETA) = DMAX1((CHI - D_TWO*DT*(T1-T2-T3*B(IETA)))
c$$$     $        *(CSDRH(IETA)/CHI) - D_TWO*DT*T4, D_ZERO)

         CSDRI(IETA) = (CHI - D_TWO*DT*(T1-T2-T3))
     $        *(CSDRH(IETA)/CHI) - T4

         !CSDRI(IETA) = DMAX1(CSDRI(IETA), D_ZERO)
      ENDDO

      IF(VERBOSE) THEN
         WRITE(*,50)' '
         WRITE(*,50)'============================================='
         WRITE(*,50) 'PMF-B: COND. SCAL. DISS. RATE (INHOMOGENEOUS):'
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

      DOUBLE PRECISION FUNCTION PMF_PDF_FEXT(ETA)
      IMPLICIT NONE
      DOUBLE PRECISION MFMEAN,MFVAR,ALPHA,TAU,SIGMA2
      DOUBLE PRECISION ETA,PDF,E,PHI
      DOUBLE PRECISION ERFINV
      EXTERNAL ERFINV
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION TAUPASS
      COMMON/COMVARSTAU/TAUPASS

      MFMEAN = MFMEANPASS
      MFVAR  = MFVARPASS
      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)
      TAU = TAUPASS             !CALL FIND_TAU(MFMEAN,MFVAR,TAU)
      SIGMA2 = D_ONE - D_TWO*TAU
      E = ERFINV(D_TWO*ETA-D_ONE)
      PHI = ALPHA + D_TWO*DSQRT(TAU)*E
      PMF_PDF_FEXT = DSQRT(D_TWO*TAU/SIGMA2)
     $        * DEXP(E**D_TWO - (PHI**D_TWO)/(D_TWO*SIGMA2))
        
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION TAUFUNC_M(MFMEAN)
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
      DOUBLE PRECISION MFMEAN,MFVAR,TAU
      DOUBLE PRECISION MFVARPASS_DERIVS
      COMMON/COMVARSVAR/MFVARPASS_DERIVS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      MFVAR = MFVARPASS_DERIVS
      CALL FIND_TAU(MFMEAN,MFVAR,TAU)
      TAUFUNC_M = TAU

      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION TAUFUNC_V(MFVAR)
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
      DOUBLE PRECISION MFMEAN,MFVAR,TAU
      DOUBLE PRECISION MFMEANPASS_DERIVS
      COMMON/COMVARSMEAN/MFMEANPASS_DERIVS
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      MFMEAN = MFMEANPASS_DERIVS
      CALL FIND_TAU(MFMEAN,MFVAR,TAU)
      TAUFUNC_V = TAU

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
      DOUBLE PRECISION TAUFUNC,ERFINV,ZRIDDR,ZBRENT
      EXTERNAL TAUFUNC,ERFINV,ZRIDDR,ZBRENT
      DOUBLE PRECISION TAU,MFMEAN,MFVAR
      DOUBLE PRECISION TAU_MIN, TAU_MAX
      INTEGER IFLAG
      DOUBLE PRECISION MFMEANPASS,MFVARPASS
      COMMON/COMVARSMEANVAR/MFMEANPASS,MFVARPASS
      DOUBLE PRECISION ALPHAPASS
      COMMON/COMVARSALPHA/ALPHAPASS

      DOUBLE PRECISION ROOTF_RE,ROOTF_AE
      COMMON/ROOTFINDERVARS/ROOTF_RE,ROOTF_AE
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      MFMEANPASS = MFMEAN
      MFVARPASS  = MFVAR
      ALPHAPASS = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)
      
C     TAU VARIES BETWEEN 0 AND 0.5
      TAU_MIN = D_ZERO
      TAU_MAX = D_HALF
      CALL ZEROIN(TAUFUNC,TAU_MIN,TAU_MAX,ROOTF_RE,ROOTF_AE,IFLAG)
      TAU = TAU_MIN 
            
      RETURN 
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      DOUBLE PRECISION FUNCTION TAUFUNC(TAU)
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
    
      DOUBLE PRECISION TAUFUNC_INT
      EXTERNAL TAUFUNC_INT
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
      CALL DQAGI(TAUFUNC_INT,BOUND,INF,EPSABS,EPSREL,RESULT,
     $     INTEG_ABSERR,NEVAL,IER,INTEG_LIM,LENW,LAST,IWORK,WORK)
      INTEGRAL = RESULT
      TAUFUNC = MFMEANPASS**D_TWO + MFVARPASS - INTEGRAL

      RETURN 
      END

C==================================================================================
C==================================================================================
C==================================================================================
C==================================================================================  

      DOUBLE PRECISION FUNCTION TAUFUNC_INT(PHI)
C================================================================================== 
C     PURPOSE: COMPUTES THE TAUFUNC_INT OF THE SECOND TERM ON THE RIGHT-HAND SIDE 
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
      TAUFUNC_INT = (X**D_TWO)*R

      RETURN 
      END



