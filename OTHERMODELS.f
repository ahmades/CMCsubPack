C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE GIR_BETA_CSDR(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH)
C========================================================================================== 
C     PURPOSE: COMPUTES THE [HOMOGENEOUS] CONDITIONAL SCALAR  DISSIPATION RATE MODEL OF
C              GIRIMAJI [S.S. GIRIMAJI. ON THE MODELING OF SCALAR DIFFUSION IN ISOTROPIC 
C              FLOWS. PHYS. FLUIDS 4(11):2529–2537, 1992]. THE MODEL IS BASED ON THE 
C              DOUBLE INTEGRATION OF THE HOMOGENEOUS PDF TRANSPORT EQUATION AND PRESUMES
C              THE PDF USING THE BETA DISTRIBUTION.
C
C     NOTE:    THIS MODEL IS IDENTICAL TO THE HOMOGENEOUS VERSION OF THE BETA-PDF-BASED
C              MODEL DEVISED BY MORENSEN (SUBROUTINE BETA_CSDR_H). ANY POSSIBLE DIFFERENCES
C              BETWEEN THE OUTCOMES OF THE TWO MODELS SHOULD BE ATTRIBUTED TO NUMERICAL 
C              ERRORS.
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
      INTEGER NETA,IETA
      DOUBLE PRECISION MFMEAN,MFVAR,CHI
      DOUBLE PRECISION ETA(NETA),PDF(NETA),CSDRH(NETA),I(NETA)

      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      DOUBLE PRECISION INTI,I1,I2
      DOUBLE PRECISION IFUN,I1FUN,I2FUN
      EXTERNAL IFUN,I1FUN,I2FUN
      DOUBLE PRECISION LIMINF,LIMSUP
      INTEGER NEVAL,IER,LIMIT,LENW,LAST
      DOUBLE PRECISION EPSABS,EPSREL,RESULT
      INTEGER IWORK(INTEG_LIM)
      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR

      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION ETAPASS
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION I1PASS,I2PASS
      COMMON/IPASS/I1PASS,I2PASS

      MFMEANPASS = MFMEAN
      MFVARPASS  = MFVAR
  
      LIMIT   = INTEG_LIM
      LENW    = I_FOUR*INTEG_LIM
      EPSABS  = INTEG_EPSABS
      EPSREL  = INTEG_EPSREL

      LIMINF = ETA(1)
      LIMSUP = ETA(NETA)
      INTEG_ABSERR = D_ZERO
      CALL DQAGS(I1FUN,LIMINF,LIMSUP,EPSABS,EPSREL,I1,
     $     INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
      I1PASS = I1

C      WRITE(*,*)IER
C      PAUSE

      LIMINF = ETA(1)
      LIMSUP = ETA(NETA)
      INTEG_ABSERR = D_ZERO
      CALL DQAGS(I2FUN,LIMINF,LIMSUP,EPSABS,EPSREL,I2,
     $     INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK)
      I2PASS = I2
      
C      WRITE(*,*)IER
C      PAUSE

      DO IETA = 2,NETA-1
         ETAPASS = ETA(IETA)
         IF(ETA(IETA).LE.MFMEAN)THEN
            LIMINF = ETA(1)
            LIMSUP = ETA(IETA)
            INTEG_ABSERR = D_ZERO
            CALL DQAGS(IFUN,LIMINF,LIMSUP,EPSABS,EPSREL,INTI,
     $           INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK) 
            I(IETA) = INTI
C            WRITE(*,*)IER
C            PAUSE
         ELSE
            LIMINF = ETA(IETA)
            LIMSUP = ETA(NETA)
            INTEG_ABSERR = D_ZERO
            CALL DQAGS(IFUN,LIMINF,LIMSUP,EPSABS,EPSREL,INTI,
     $           INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK) 
            I(IETA) = -INTI
C            WRITE(*,*)IER
C            PAUSE
         ENDIF
      ENDDO

      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)

      CSDRH(1) = D_ZERO
      CSDRH(NETA) = D_ZERO
      DO IETA = 2,NETA-1
         CSDRH(IETA) = DMAX1(-D_TWO*CHI*(MFMEAN*(D_ONE-MFMEAN)
     $        /(MFVAR**D_TWO))*(I(IETA)/PDF(IETA)),D_ZERO) 
      ENDDO
          
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      DOUBLE PRECISION FUNCTION IFUN(ETA)
C========================================================================================== 
C     PURPOSE: FUNCTION NECESSARY TO COMPUTE THE INTEGRAND OF THE R.H.S. OF EQ. (30) IN 
C              [S.S. GIRIMAJI. ON THE MODELING OF SCALAR DIFFUSION IN ISOTROPIC FLOWS. PHYS. 
C              FLUIDS 4(11):2529–2537, 1992]
C==========================================================================================
C            VARIABLE   DESCRIPTION                        DATA TYPE
C            --------   -----------                        --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION VALUE             DOUBLE PRECISION
C
C     RETURNS:
C
C            IFUN      INTEGRAND OF THE R.H.S.OF EQ.(32)   DOUBLE PRECISION
C==========================================================================================
      IMPLICIT NONE
      DOUBLE PRECISION ETA
      DOUBLE PRECISION G,V,W,B,MFMEAN,MFVAR,PDF
      DOUBLE PRECISION BETA
      EXTERNAL BETA
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      DOUBLE PRECISION ETAPASS
      COMMON/EPASS/ETAPASS
      DOUBLE PRECISION I1PASS,I2PASS
      COMMON/IPASS/I1PASS,I2PASS
      MFMEAN = MFMEANPASS
      MFVAR  = MFVARPASS
      G = (MFMEAN*(D_ONE-MFMEAN)/MFVAR) - D_ONE
      V = MFMEAN*G
      W = (D_ONE-MFMEAN)*G
      B = BETA(V,W)
      PDF = (ETA**(V-D_ONE))*((D_ONE-ETA)**(W-D_ONE))/B
      IFUN = (MFMEAN*(DLOG(ETA)-I1PASS)
     $     + (1.0D0-MFMEAN)*(DLOG(1.0D0-ETA)-I2PASS))
     $     * PDF * (ETAPASS-ETA)
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      DOUBLE PRECISION FUNCTION I1FUN(ETA)    
C========================================================================================== 
C     PURPOSE: FUNCTION NECESSARY TO COMPUTE <ln(eta)> WHICH APPEARS IN THE INTEGRAND OF 
C              THE R.H.S. OF EQ. (30) IN [S.S. GIRIMAJI. ON THE MODELING OF SCALAR 
C              DIFFUSION IN ISOTROPIC FLOWS. PHYS. FLUIDS 4(11):2529–2537, 1992]
C==========================================================================================
C            VARIABLE   DESCRIPTION                        DATA TYPE
C            --------   -----------                        --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION VALUE             DOUBLE PRECISION
C
C     RETURNS:
C
C            I1FUN      ln(eta)P(eta)                      DOUBLE PRECISION
C==========================================================================================
      IMPLICIT NONE
      DOUBLE PRECISION ETA
      DOUBLE PRECISION G,V,W,B,MFMEAN,MFVAR,PDF
      DOUBLE PRECISION BETA
      EXTERNAL BETA
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      MFMEAN = MFMEANPASS
      MFVAR  = MFVARPASS
      G = (MFMEAN*(D_ONE-MFMEAN)/MFVAR) - D_ONE
      V = MFMEAN*G
      W = (D_ONE-MFMEAN)*G
      B = BETA(V,W)
      PDF = (ETA**(V-D_ONE))*((D_ONE-ETA)**(W-D_ONE))/B
      I1FUN = DLOG(ETA)*PDF
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      DOUBLE PRECISION FUNCTION I2FUN(ETA)
C========================================================================================== 
C     PURPOSE: FUNCTION NECESSARY TO COMPUTE <ln(1-eta)> WHICH APPEARS IN THE INTEGRAND OF 
C              THE R.H.S. OF EQ. (30) IN [S.S. GIRIMAJI. ON THE MODELING OF SCALAR 
C              DIFFUSION IN ISOTROPIC FLOWS. PHYS. FLUIDS 4(11):2529–2537, 1992]
C==========================================================================================
C            VARIABLE   DESCRIPTION                        DATA TYPE
C            --------   -----------                        --------- 
C
C     INPUT:
C
C            ETA        MIXTURE FRACTION VALUE             DOUBLE PRECISION
C
C     RETURNS:
C
C            I1FUN      ln(1-eta)P(eta)                    DOUBLE PRECISION
C==========================================================================================
      IMPLICIT NONE
      DOUBLE PRECISION ETA
      DOUBLE PRECISION G,V,W,B,MFMEAN,MFVAR,PDF
      DOUBLE PRECISION BETA
      EXTERNAL BETA
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      MFMEAN = MFMEANPASS
      MFVAR  = MFVARPASS
      G = (MFMEAN*(D_ONE-MFMEAN)/MFVAR) - D_ONE
      V = MFMEAN*G
      W = (D_ONE-MFMEAN)*G
      B = BETA(V,W)
      PDF = (ETA**(V-D_ONE))*((D_ONE-ETA)**(W-D_ONE))/B
      I2FUN = DLOG(1.0D0-ETA)*PDF
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE LINEAR_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,DT,VEL,CV)
C========================================================================================== 
C     PURPOSE: COMPUTES THE LINEAR CONDITIONAL VELOCITY MODEL [V. R. KUZNETSOV AND V. A. 
C              SABEL'NIKOV. TURBULENCE AND COMBUSTION. ENGLISH EDITION EDITOR: P. A. LIBBY. 
C              HEMISPHERE PUBLISHING CORPORATION, NEW YORK, U.S.A., 1990]
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
      DOUBLE PRECISION MFMEANGRAD(3),VEL(3)
      DOUBLE PRECISION ETA(NETA),CV(3,NETA)

      DO I = 1,3
         DO IETA = 1,NETA
            CV(I,IETA) = VEL(I) 
     $           - DT*MFMEANGRAD(I)*(ETA(IETA)-MFMEAN)/MFVAR
         ENDDO
      ENDDO

      RETURN
      END



C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  

      SUBROUTINE AMC_BETA_CSDR(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH)
C========================================================================================== 
C     PURPOSE: COMPUTES THE AMPLITUDE MAPPING CLOSURE [E. E. O'BRIAN AND T. L. JIANG. THE 
C              CONDITIONAL DISSIPATION RATE OF AN INITIALLY BINARY SCALAR IN HOMOGENEOUS 
C              TURBULENCE. PHYS. FLUIDS A, 3(12):3121–3123, 1991] WITH A PRESUMED [[BETA-PDF]]
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
      INTEGER NETA,IETA
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,CHI0,G
      DOUBLE PRECISION CSDRH(NETA),ETA(NETA),PDF(NETA)
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION LIMINF,LIMSUP
      INTEGER NEVAL,IER,LIMIT,LENW,LAST
      DOUBLE PRECISION EPSABS,EPSREL,RESULT
      INTEGER IWORK(INTEG_LIM)
      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION INTGP
      DOUBLE PRECISION ERFINV,AMCBETAFUN
      EXTERNAL ERFINV,AMCBETAFUN
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS

      MFMEANPASS = MFMEAN
      MFVARPASS  = MFVAR

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

      LIMINF  = D_ZERO
      LIMSUP  = D_ONE
      EPSREL  = INTEG_EPSREL
      EPSABS  = INTEG_EPSABS
      LIMIT   = INTEG_LIM
      LENW    = I_FOUR*INTEG_LIM
      CALL DQAGS(AMCBETAFUN,LIMINF,LIMSUP,EPSABS,EPSREL,INTGP,
     $     INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK) 

      CHI0 = CHI/INTGP

      CSDRH(1) = D_ZERO
      CSDRH(NETA) = D_ZERO
      DO IETA = 2,NETA-1
         G = DEXP(-D_TWO*((ERFINV(D_TWO*ETA(IETA)-D_ONE))**D_TWO))
         CSDRH(IETA) = DMAX1(CHI0*G,D_ZERO)
      ENDDO
        
      RETURN
      END

      DOUBLE PRECISION FUNCTION AMCBETAFUN(ETA)
      IMPLICIT NONE
      DOUBLE PRECISION ETA,MFMEAN,MFVAR
      DOUBLE PRECISION GAM,V,W,B,P,G
      DOUBLE PRECISION ERFINV,BETA
      EXTERNAL ERFINV,BETA
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      MFMEAN = MFMEANPASS
      MFVAR  = MFVARPASS
      GAM = (MFMEAN*(D_ONE-MFMEAN)/MFVAR) - D_ONE
      V = MFMEAN*GAM
      W = (D_ONE-MFMEAN)*GAM
      B = BETA(V,W)
      P = (ETA**(V-D_ONE))*((D_ONE-ETA)**(W-D_ONE))/B
      G = DEXP(-D_TWO*(ERFINV(D_TWO*ETA-D_ONE)**D_TWO))
      AMCBETAFUN = G*P
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================  
      SUBROUTINE AMC_PMFSML_CSDR(ETA,NETA,MFMEAN,MFVAR,CHI,CSDRH)
C========================================================================================== 
C     PURPOSE: COMPUTES THE AMPLITUDE MAPPING CLOSURE [E. E. O'BRIAN AND T. L. JIANG. THE 
C              CONDITIONAL DISSIPATION RATE OF AN INITIALLY BINARY SCALAR IN HOMOGENEOUS 
C              TURBULENCE. PHYS. FLUIDS A, 3(12):3121–3123, 1991] WITH THE [[PMF-PDF]]
C
C     NOTE:    THIS MODEL IS IDENTICAL TO THE HOMOGENEOUS VERSION OF THE PMF-PDF-BASED
C              MODEL DEVISED BY MORENSEN (SUBROUTINE PMF_CSDR_H). ANY POSSIBLE DIFFERENCES
C              BETWEEN THE OUTCOMES OF THE TWO MODELS SHOULD BE ATTRIBUTED TO NUMERICAL 
C              ERRORS.
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
      INTEGER NETA,IETA
      DOUBLE PRECISION MFMEAN,MFVAR,CHI,CHI0,G
      DOUBLE PRECISION CSDRH(NETA),ETA(NETA),PDF(NETA)
      LOGICAL VERBOSE_TEMP
      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION LIMINF,LIMSUP
      INTEGER NEVAL,IER,LIMIT,LENW,LAST
      DOUBLE PRECISION EPSABS,EPSREL,RESULT
      INTEGER IWORK(INTEG_LIM)
      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
      DOUBLE PRECISION INTEG_ABSERR
      COMMON/INTEGRATORABSERR/INTEG_ABSERR
      DOUBLE PRECISION INTGP
      DOUBLE PRECISION ERFINV,AMCPMFFUN
      EXTERNAL ERFINV,AMCPMFFUN
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS

      MFMEANPASS = MFMEAN
      MFVARPASS  = MFVAR

C     FIND THE PDF. TEMPORARILY DISABLE VERBOSITY IF ON.
      VERBOSE_TEMP = VERBOSE
      VERBOSE = .FALSE.
      CALL PMFSML_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
      VERBOSE = VERBOSE_TEMP

      LIMINF  = D_ZERO
      LIMSUP  = D_ONE
      EPSREL  = INTEG_EPSREL
      EPSABS  = INTEG_EPSABS
      LIMIT   = INTEG_LIM
      LENW    = I_FOUR*INTEG_LIM
      CALL DQAGS(AMCPMFFUN,LIMINF,LIMSUP,EPSABS,EPSREL,INTGP,
     $     INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,LAST,IWORK,WORK) 

      CHI0 = CHI/INTGP

      CSDRH(1) = D_ZERO
      CSDRH(NETA) = D_ZERO
      DO IETA = 2,NETA-1
         G = DEXP(-D_TWO*((ERFINV(D_TWO*ETA(IETA)-D_ONE))**D_TWO))
         CSDRH(IETA) = DMAX1(CHI0*G,D_ZERO)
      ENDDO
        
      RETURN
      END

      DOUBLE PRECISION FUNCTION AMCPMFFUN(ETA)
      IMPLICIT NONE
      DOUBLE PRECISION ETA,MFMEAN,MFVAR
      DOUBLE PRECISION ALPHA,TAU,SIGMA2,E,PHI,P,G
      DOUBLE PRECISION ERFINV
      EXTERNAL ERFINV
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      DOUBLE PRECISION MFMEANPASS
      COMMON/MFMEANBLOK/MFMEANPASS
      DOUBLE PRECISION MFVARPASS
      COMMON/MFVARBLOK/MFVARPASS
      MFMEAN = MFMEANPASS
      MFVAR  = MFVARPASS
      ALPHA = DSQRT(D_TWO)*ERFINV(D_ONE - D_TWO*MFMEAN)
      CALL FIND_TAU(MFMEAN,MFVAR,TAU)
      SIGMA2 = D_ONE - D_TWO*TAU
      E = ERFINV(D_TWO*ETA-D_ONE)
      PHI = ALPHA + D_TWO*DSQRT(TAU)*E
      P = DSQRT(D_TWO*TAU/SIGMA2)
     $     * DEXP(E**D_TWO - (PHI**D_TWO)/(D_TWO*SIGMA2))
      G = DEXP(-D_TWO*(ERFINV(D_TWO*ETA-D_ONE)**D_TWO))
      AMCPMFFUN = G*P
      RETURN
      END
      
c$$$      SUBROUTINE GIRIMAJI_CSDR(ETA,NETA,MFMEAN,MFVAR,CHI,CSDR)
c$$$      IMPLICIT NONE
c$$$      INTEGER NETA,IETA
c$$$      DOUBLE PRECISION MFMEAN,MFVAR,CHI
c$$$      DOUBLE PRECISION ETA(NETA),PDF(NETA),CSDR(NETA),J(NETA)
c$$$
c$$$      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
c$$$      INTEGER INTEG_LIM
c$$$      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
c$$$      COMMON/INTEGRATORVARS2/INTEG_LIM
c$$$      INTEGER I_TWO,I_FOUR
c$$$      COMMON/INTCONSTANTS/I_TWO,I_FOUR
c$$$      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
c$$$
c$$$      DOUBLE PRECISION INTJ,J1,J2
c$$$      DOUBLE PRECISION JFUN,J1FUN,J2FUN
c$$$      EXTERNAL JFUN,J1FUN,J2FUN
c$$$      DOUBLE PRECISION LIMINF,LIMSUP
c$$$      INTEGER NEVAL,IER,LIMIT,LENW,LAST
c$$$      DOUBLE PRECISION EPSABS,EPSREL,RESULT
c$$$      INTEGER IWORK(INTEG_LIM)
c$$$      DOUBLE PRECISION WORK(I_FOUR*INTEG_LIM)
c$$$      DOUBLE PRECISION INTEG_ABSERR
c$$$      COMMON/INTEGRATORABSERR/INTEG_ABSERR
c$$$      DOUBLE PRECISION ALF,BET
c$$$      INTEGER INTEGR 
c$$$
c$$$      DOUBLE PRECISION MFMEANPASS
c$$$      COMMON/MFMEANBLOK/MFMEANPASS
c$$$      DOUBLE PRECISION MFVARPASS
c$$$      COMMON/MFVARBLOK/MFVARPASS
c$$$      DOUBLE PRECISION ETAPASS
c$$$      COMMON/EPASS/ETAPASS
c$$$      DOUBLE PRECISION J1PASS,J2PASS
c$$$      COMMON/JPASS/J1PASS,J2PASS
c$$$
c$$$      MFMEANPASS = MFMEAN
c$$$      MFVARPASS  = MFVAR
c$$$  
c$$$      LIMIT   = INTEG_LIM
c$$$      LENW    = I_FOUR*INTEG_LIM
c$$$      EPSABS  = INTEG_EPSABS
c$$$      EPSREL  = INTEG_EPSREL
c$$$      INTEGR = 1
c$$$      ALF = 0.0D0
c$$$      BET = 0.0D0
c$$$
c$$$      LIMINF = ETA(1)
c$$$      LIMSUP = ETA(NETA)
c$$$      INTEG_ABSERR = D_ZERO
c$$$      CALL DQAWS(J1FUN,LIMINF,LIMSUP,ALF,BET,INTEGR,
c$$$     $        EPSABS,EPSREL,J1,INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,
c$$$     $        LAST,IWORK,WORK)
c$$$      J1PASS = J1
c$$$C      WRITE(*,*)IER,J1
c$$$C      PAUSE
c$$$
c$$$      LIMINF = ETA(1)
c$$$      LIMSUP = ETA(NETA)
c$$$      INTEG_ABSERR = D_ZERO
c$$$       CALL DQAWS(J2FUN,LIMINF,LIMSUP,ALF,BET,INTEGR,
c$$$     $        EPSABS,EPSREL,J2,INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,
c$$$     $        LAST,IWORK,WORK)
c$$$      J2PASS = J2
c$$$C      WRITE(*,*)IER,J2
c$$$C      PAUSE
c$$$
c$$$      DO IETA = 2,NETA-1
c$$$         ETAPASS = ETA(IETA)
c$$$         IF(ETA(IETA).LE.MFMEAN)THEN
c$$$            LIMINF = ETA(1)
c$$$            LIMSUP = ETA(IETA)
c$$$            INTEG_ABSERR = D_ZERO
c$$$            CALL DQAWS(JFUN,LIMINF,LIMSUP,ALF,BET,INTEGR,
c$$$     $           EPSABS,EPSREL,INTJ,INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,
c$$$     $           LAST,IWORK,WORK)
c$$$            J(IETA) = INTJ
c$$$C            WRITE(*,*)IER,J(IETA)
c$$$C            pause
c$$$         ELSE
c$$$             LIMINF = ETA(IETA)
c$$$             LIMSUP = ETA(NETA)
c$$$             INTEG_ABSERR = D_ZERO
c$$$             CALL DQAWS(JFUN,LIMINF,LIMSUP,ALF,BET,INTEGR,
c$$$     $           EPSABS,EPSREL,INTJ,INTEG_ABSERR,NEVAL,IER,LIMIT,LENW,
c$$$     $           LAST,IWORK,WORK)
c$$$             J(IETA) = -INTJ
c$$$C             WRITE(*,*)IER,J(IETA)
c$$$C             pause
c$$$          ENDIF
c$$$      ENDDO
c$$$
c$$$      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF)
c$$$
c$$$      CSDR(1) = D_ZERO
c$$$      CSDR(NETA) = D_ZERO
c$$$      DO IETA = 2,NETA-1
c$$$         CSDR(IETA) = -D_TWO*CHI*(MFMEAN*(D_ONE-MFMEAN)/(MFVAR**D_TWO))
c$$$     $        *(J(IETA)/PDF(IETA)) 
c$$$      ENDDO
c$$$          
c$$$      RETURN
c$$$      END
