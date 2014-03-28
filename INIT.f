C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      BLOCK DATA
C========================================================================================== 
C     SOME FREQUENTLY USED CONSTANTS
C     *****************************
C     ****    DO NOT CHANGE    ****
C     *****************************
C========================================================================================== 
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE,
     $     D_FOUR,D_EIGHT
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS2/D_FOUR,D_EIGHT
      INTEGER I_TWO,I_FOUR
      COMMON/INTCONSTANTS/I_TWO,I_FOUR
      DOUBLE PRECISION EPSL
      COMMON/EPSILON/EPSL
      DOUBLE PRECISION MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      COMMON/MFMINMAX/MFMEAN_MAX,MFMEAN_MIN,MFVAR_MIN
      DATA PI/3.141592653589793238462D0/
      DATA D_ZERO/0.0D0/
      DATA D_ONE/1.0D0/
      DATA D_HALF/0.5D0/
      DATA D_TWO/2.0D0/
      DATA D_THREE/3.0D0/
      DATA D_FOUR/4.0D0/
      DATA D_EIGHT/8.0D0/
      DATA I_TWO/2/
      DATA I_FOUR/4/
      DATA MFMEAN_MIN/1.0D-9/ !6/
      DATA MFMEAN_MAX/0.999999D0/
      DATA MFVAR_MIN/1.0D-6/ !4/
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE SETPARMS
C================================================================================== 
C     PURPOSE: SPECIFIES A SET OF OPTIONS AND/OR TOLERANCES FOR THE INTEGRATOR, 
C              DIFFERENTIATOR AND ROOT FINDER (MUST BE CALLED ONCE        
C              PRIOR TO CALLING ANY OF THE SUBMODEL SUBROUTINES AND BEFORE 
C              BEFORE LOOPING OVER ALL THE GRID POINTS OF THE PHYSICAL DOMAIN
C              IN ANY CFD CALCULATIONS EMPLOYING THE SUBROUTINES
C==================================================================================
      IMPLICIT NONE

      LOGICAL VERBOSE,PROGRESS
      COMMON/LOGICALVARS/VERBOSE,PROGRESS

      DOUBLE PRECISION ROOTF_RE,ROOTF_AE
      COMMON/ROOTFINDERVARS1/ROOTF_RE,ROOTF_AE

      DOUBLE PRECISION INTEG_EPSABS,INTEG_EPSREL
      INTEGER INTEG_LIM,INTEG_METH,INTEG_ORD
      COMMON/INTEGRATORVARS1/INTEG_EPSABS,INTEG_EPSREL
      COMMON/INTEGRATORVARS2/INTEG_LIM
      COMMON/INTEGRATORVARS3/INTEG_METH,INTEG_ORD

      INTEGER DIFF_NCHEBCOEFS
      COMMON/DIFFERENTIATORVARS1/DIFF_NCHEBCOEFS
      DOUBLE PRECISION DIFF_PDFMIN
      COMMON/DIFFERENTIATORVARS3/DIFF_PDFMIN
      INTEGER DIFF_NCHEBCOEFS_EXTRA
      COMMON/DIFFERENTIATORVARS4/DIFF_NCHEBCOEFS_EXTRA

      DOUBLE PRECISION NLS_TOL
      COMMON/NONLINEARSOLVERVARS/NLS_TOL

C===================================================================================================
C     VERBOSITY: [.FALSE.]OFF / [.TRUE.]ON
C===================================================================================================
C
      VERBOSE  = .FALSE.        ! SHOWS OUTPUT
      PROGRESS = .TRUE.         ! SHOWS PROGRESS (%) IN THE COMPUTATIONALLY DEMANDING SUBMODELS
C
C===================================================================================================
C     ROOT FINDER CONTROL 
C===================================================================================================
C
      ROOTF_RE = 1.0D-5         ! RELATIVE ERROR USED IN THE STOPPING CRITERION.
      ROOTF_AE = 1.0D-5         ! ABSOLUTE ERROR USED IN THE STOPPING CRITERION.
C
C===================================================================================================
C     INTEGRATOR CONTROL (REQUIRED FOR QUADPACK SUBROUTINES)
C===================================================================================================
C
      INTEG_METH = 2            ! INTEGRATION METHOD: [1] ADAPTIVE (QUADPACK)
                                !                     [2] NON-ADAPTIVE (IPACK)
                                ! NOTE: 
                                ! A) SPECIFYING [2] ONLY AFFECTS THE COMPUTATION OF THE PARAMETER 
                                !    TAU IN THE **TRINARY PMF CLOSURES**. ALL OTHER INTEGRATIONS  
                                !    IN THE CODE USE ADAPTIVE QUADRATURE METHOD (QUADPACK).
                                ! B) OPTION [2] EMPLOYS A GAUSS-HERMITE QUADRATURE RULE USING 
                                !    INTEG_ORD POINTS (SPECIFIED BELOW)
                                ! C) OPTION [1] IS MORE ACCURATE THAN OPTION [2] BUT LESS EFFICIENT,
                                !    SPECIALLY INTEG_LIM (SEE BELOW) IS LARGE. 
                                ! D) THE ACCURACY OF INTEGRATION USING OPTION [2] DEPENDS ON THE 
                                !    THE NUMBER OF POINT, TEG_ORD. LARGE TEG_ORD IMPROVES ACCURACY
                                !    BUT INCREASES THE COMPUTATIONAL COST
C
C===================================================================================================
C     PARAMETERS AFFECTING THE ADAPTIVE INTEGRATION METHOD
C===================================================================================================
C
      INTEG_LIM = 1000 !5       ! MAXIMUM NUMBER OF SUBINTERVALS IN THE PARTITION OF THE GIVEN 
                                ! INTEGRATION INTERVAL.
C
      INTEG_EPSABS = 1.0D-5     ! REQUESTED ABSOLUTE ACCURACY.
C
      INTEG_EPSREL = 1.0D-5     ! REQUESTED RELATIVE ACCURACY.
C
C===================================================================================================
C     PARAMETERS AFFECTING THE NON-ADAPTIVE INTEGRATION METHOD
C===================================================================================================
C
      INTEG_ORD = 16            ! NUMBER OF POINTS USED TO GENERATE THE WEIGHTS OF THE GAUSS-HERMITE
                                ! QUADRATURE RULE
C
C===================================================================================================
C     DIFFERENTIATOR CONTROL
C===================================================================================================
C
C     PARAMETERS FOR CHEBYSHEV'S APPROXIMATION [1]
C
      DIFF_NCHEBCOEFS = 6       ! NUMBER OF CHEBYSHEV COEFFICIENTS.
C
      DIFF_PDFMIN = 1.0D-5      ! MINIMUM PDF VALUE BELOW WHICH DIFF_NCHEBCOEFS IS INCREASED BY
                                ! DIFF_NCHEBCOEFS_EXTRA (SEE BELOW). THIS ENSURES ACCURACY WHEN
                                ! THERE IS DIVISION BY SMALL PDF(IETA) VALUES BY MEANS OF REDUCTION
                                ! OF TRUNCATION ERROR.
                                ! NOTE 1: ONLY AFFECT THE PARTIAL DIFFERENTIATION OF [PDF(IETA)]
                                !         AND [II(ETA)] WITH RESPECT TO THE MIXTURE FRACTION MEAN 
                                !         AND VARIANCE. DIFF_NCHEBCOEFS IS USED FOR THE EVALUATION
                                !         OF THE PARTIAL DERIVATIVES OF [TAU] WITH RESPECT TO THE  
                                !         MIXTURE FRACTION MEAN AND VARIANCE.
                                ! NOTE 2: SETTING LARGE DIFF_PDFMIN CAN SLOW DOWN PERFORMANCE 
                                !         SIGNIFICANTLY, DEPENDING ON THE SHAPE OF THE PDF 
                                !         (NUMBER OF ETA GRID PTS WHERE PDF <= DIFF_PDFMIN BECOMES 
                                !         LARGER).
C
      DIFF_NCHEBCOEFS_EXTRA = 8 ! NUMBER OF CHEBYSHEV COEFFICIENTS ADDED TO DIFF_NCHEBCOEFS
                                ! WHEN PDF(IETA) <= DIFF_PDFMIN. SET DIFF_NCHEBCOEFS_EXTRA 
                                ! CONSERVATIVELY. SETTING LARGE DIFF_NCHEBCOEFS_EXTRA CAN SLOW DOWN 
                                ! PERFORMANCE SIGNIFICANTLY, DEPENDING ON THE SHAPE OF THE PDF.
C
C===================================================================================================
C     NON-LINEAR EQUATION SOLVER CONTROL (REQUIRED FOR MINPACK SUBROUTINES)
C===================================================================================================
C
      NLS_TOL  = 1.0D-5        ! NON-LINEAR EQUATION SOLVER TOLERANCE.
C
C===================================================================================================

      RETURN
      END


