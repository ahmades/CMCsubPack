
      PROGRAM EXAMPLE
C================================================================================== 
C     PURPOSE: DRIVER PROGRAM
C==================================================================================
      IMPLICIT NONE
      INCLUDE 'formats.h'
      INCLUDE 'COM'
      
      INTEGER IETA        ! MIXTURE FRACTION SPACE CONTER
C
      DOUBLE PRECISION
     $     MFMEAN,        ! MIXTURE FRACTION MEAN
     $     MFVAR,         ! MIXTURE FRACTION VARIANCE
     $     MFMEANGRAD(3), ! GRADIENT OF THE MIXTURE FRACTION MEAN (COMPONENTS: 1=X, 2=Y, 3=Z)
     $     MFVARGRAD(3),  ! GRADIENT OF THE MIXTURE FRACTION VARIANCE (COMPONENTS: 1=X, 2=Y, 3=Z)
     $     DT,            ! TURBULENT DIFFUSIVITY
     $     VEL(3),        ! VELOCITY VECTOR (COMPONENTS: 1=X, 2=Y, 3=Z)
     $     CHIMEAN,       ! MEAN SCALAR DISSIPATION RATE
     $     PILOTMF,       ! PILOT MIXTURE FRACTION
     $     C1MEAN,C2MEAN,C3MEAN,CMEAN  ! RELATIVE MASS INJECTIONS OF AIR, PILOT AND FUEL STREAM
C
      DOUBLE PRECISION ETA(NETA)! MIXTURE FRACTION

C     PMF SML ARRAYS
      DOUBLE PRECISION
     $     PDF_PMFSML(NETA),     ! PROBABILITY DENSITY FUNCTION
     $     CV_PMFSML(3,NETA),    ! CONDITIONAL VELOCITY
     $     CSDRH_PMFSML(NETA),   ! CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS VERSION]]
     $     CSDRI_PMFSML(NETA)    ! CONDITIONAL SCALAR DISSIPATION RATE [[INHOMOGENEOUS VERSION]]

C     PMF DSML ARRAYS
      DOUBLE PRECISION
     $     PDF_PMFDSML(NETA),     ! PROBABILITY DENSITY FUNCTION
     $     CV_PMFDSML(3,NETA),    ! CONDITIONAL VELOCITY
     $     CSDRH_PMFDSML(NETA),   ! CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS VERSION]]
     $     CSDRI_PMFDSML(NETA)    ! CONDITIONAL SCALAR DISSIPATION RATE [[INHOMOGENEOUS VERSION]]

C     BETA ARRAYS
      DOUBLE PRECISION 
     $     PDF_BETA(NETA),    ! PROBABILITY DENSITY FUNCTION
     $     CV_BETA(3,NETA),   ! CONDITIONAL VELOCITY
     $     CSDRH_BETA(NETA),  ! CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS VERSION]]
     $     CSDRI_BETA(NETA)   ! CONDITIONAL SCALAR DISSIPATION RATE [[INHOMOGENEOUS VERSION]]

C     CLIPPED GAUSSIAN ARRAYS
      DOUBLE PRECISION 
     $     PDF_CG(NETA),      ! CLIPPED GAUSSIAN PDF
     $     CV_CG(3,NETA),      ! CONDITIONAL VELOCITY CG-PDF
     $     CSDRH_CG(NETA),    ! CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS VERSION]]
     $     CSDRI_CG(NETA)     ! CONDITIONAL SCALAR DISSIPATION RATE [[INHOMOGENEOUS VERSION]]
     
C     UNIFORM WITH DELTA-PDF
      DOUBLE PRECISION 
     $     PDF_DELTA(NETA),   ! DIRAC DELTA PDF
     $     CV_UNFRM(3,NETA),  ! UNIFRM CONDITIONAL VELOCITY
     $     CSDR_UNFRM(NETA)   ! UNIFROM CONDITIONAL SCALAR DISSIPATION RATE

C     OTHER MODELS ARRAYS
      DOUBLE PRECISION
     $     CV_LIN(3,NETA),       ! LINEAR CONDITIONAL VELOCITY MODEL
     $     CSDRH_GIR_BETA(NETA), ! GIRIMAJI'S CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS]]
     $     CSDRH_AMC_BETA(NETA), ! BETA-AMC CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS]]
     $     CSDRH_AMC_PMFSML(NETA)   ! PMF-AMC CONDITIONAL SCALAR DISSIPATION RATE [[HOMOGENEOUS]]
     

      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE

      DOUBLE PRECISION INTEGRAND(NETA),INTEGRAL,RELERR,ERR ! A UTILITY ARRAY AND VARIABLES
      DOUBLE PRECISION TSTART,TFINISH,CPUT_PMFSML,CPUT_PMFDSML,
     $     CPUT_BETA,CPUT_CG
      DOUBLE PRECISION AREA_PMFSML,AREA_PMFDSML,AREA_BETA,AREA_CG,
     $     ERR_PMFSML,ERR_PMFDSML,ERR_BETA
      DOUBLE PRECISION TRAP,RELERRP
      EXTERNAL TRAP,RELERRP

      DOUBLE PRECISION C1L,C1U,C2L,C2U,C3L,C3U


C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                           *****INPUT*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      PILOTMF       = 0.27D0
      MFMEAN        = 0.4440126121D0
      MFVAR         = 0.2559307963D-1
      CHIMEAN       = 0.2305485382D+3
      MFMEANGRAD(1) = 0.2829907608D2
      MFMEANGRAD(2) = -0.5249074097D3
      MFMEANGRAD(3) = 0.0D0
      MFVARGRAD(1)  = 0.2471093655D1
      MFVARGRAD(2)  = -0.5286500931D2
      MFVARGRAD(3)  = 0.0D0
      C1MEAN        = 0.6990208476D-5
      VEL(1)        = 0.1966718483D2
      VEL(2)        = -0.5240319967D0
      VEL(3)        = 0.0D0
      DT            = 0.1144416248D-3

      C2MEAN        = (MFMEAN+C1MEAN-1.0D0)/(PILOTMF-1.0D0)
      C3MEAN        = 1.0D0 - (C1MEAN+C2MEAN)
      C1L           = 0.0D0
      C1U           = 1.0D0-MFMEAN 
      C2L           = 0.0D0 
      C2U           = (1.0D0-MFMEAN)/(1.0D0-PILOTMF)
      C3L           = (MFMEAN-PILOTMF)/(1.0D0-PILOTMF)
      C3U           = MFMEAN


      WRITE(*,'(A)')' '
      WRITE(*,'(A)')' =========================='
      WRITE(*,'(A)')'|      INPUT PARAMETERS    |'
      WRITE(*,'(A)')' =========================='
      WRITE(*,201)'C1(O) =',C1MEAN,' RANGE: [',C1L,',',C1U,']'
      WRITE(*,201)'C2(P) =',C2MEAN,' RANGE: [',C2L,',',C2U,']'
      WRITE(*,201)'C3(F) =',C3MEAN,' RANGE: [',C3L,',',C3U,']'
      WRITE(*,301)'MEAN  =',MFMEAN
      WRITE(*,301)'VAR   =',MFVAR
      WRITE(*,301)'CHI   =',CHIMEAN
      WRITE(*,301)'DMDX  =',MFMEANGRAD(1)
      WRITE(*,301)'DMDY  =',MFMEANGRAD(2)
      WRITE(*,301)'DVDX  =',MFVARGRAD(1)
      WRITE(*,301)'DVDY  =',MFVARGRAD(2)
      WRITE(*,301)'DT    =',DT
      WRITE(*,301)'U    =',VEL(1)
      WRITE(*,301)'V    =',VEL(2)
      WRITE(*,'(A)')' ========================='
      WRITE(*,'(A)')'| END OF INPUT PARAMETERS |'
      WRITE(*,'(A)')' ========================='
      WRITE(*,'(A)')' '

 101  FORMAT(A,E18.10E3,A,E18.10E3)
 201  FORMAT(A,E18.10E3,A,E18.10E3,A,E18.10E3,A)
 301  FORMAT(A,E18.10E3)


    
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                           *****SETUP*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      WRITE(*,50)'============================================='

C     GENERATE UNIFORM MIXTURE FRACTION GRID.
C     OTHER NON-UNIFORM GRID GENERATION TECHNIQUES ARE AVAILABLE IN GRID.f
      WRITE(*,*)'CREATING MIXTURE FRACTION GRID'
      CALL ETAGRID_UNIFORM(NETA,ETA)

C     INITIALISE
      WRITE(*,*)'SETTING PARAMETERS'
      CALL SETPARMS

C     CHECK THE VALIDITY OF THE INPUT PARAMETERS (MIX. FRAC. MEAN NAD VARIANCE)
      CALL CHECKPARMS(MFMEAN,MFVAR)

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                     *****BINARY PMF APPROACH*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      WRITE(*,50)'============================================='
      WRITE(*,50) 'PRESUMED MAPPING FUNCTION APPROACH - BINARY'
      WRITE(*,50)'============================================='

      CALL CPU_TIME(TSTART)

C     COMPUTE THE PMF-PDF
      WRITE(*,50) '*PDF'
      CALL PMFSML_PDF(ETA,NETA,MFMEAN,MFVAR,PDF_PMFSML)

C     COMPUTE THE PMF-CV
      WRITE(*,50) '*CV'
      CALL PMFSML_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,DT,
     $     VEL,CV_PMFSML)

C     COMPUTE THE [[HOMOGENEOUS]] PMF-CSDR
      WRITE(*,50) '*CSDR-HOMOGENEOUS'
      CALL PMFSML_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHIMEAN,CSDRH_PMFSML)

C     COMPUTE THE [[INHOMOGENEOUS]] PMF-CSDR
      WRITE(*,50) '*CSDR-INHOMOGENEOUS'
      CALL PMFSML_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,CHIMEAN,CSDRI_PMFSML)

      CALL CPU_TIME(TFINISH)
    
      CPUT_PMFSML = TFINISH - TSTART

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                     *****TRINARY PMF APPROACH*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      WRITE(*,50)'============================================='
      WRITE(*,50) 'PRESUMED MAPPING FUNCTION APPROACH - TRINARY'
      WRITE(*,50)'============================================='

      CALL CPU_TIME(TSTART)

C     COMPUTE THE PMF-PDF
      WRITE(*,50) '*PDF'
      CALL PMFDSML_PDF(ETA,NETA,PILOTMF,C1MEAN,MFMEAN,MFVAR,
     $     PDF_PMFDSML)

C     COMPUTE THE PMF-CV
      WRITE(*,50) '*CV'
      CALL PMFDSML_CV(ETA,NETA,PILOTMF,C1MEAN,MFMEAN,MFVAR,
     $     MFMEANGRAD,MFVARGRAD,DT,VEL,CV_PMFDSML)

C     COMPUTE THE [[HOMOGENEOUS]] PMF-CSDR
      WRITE(*,50) '*CSDR-HOMOGENEOUS'
      CALL PMFDSML_CSDR_H(ETA,NETA,PILOTMF,C1MEAN,MFMEAN,MFVAR,
     $     CHIMEAN,CSDRH_PMFDSML)

C     COMPUTE THE [[INHOMOGENEOUS]] PMF-CSDR
      WRITE(*,50) '*CSDR-INHOMOGENEOUS'
      CALL PMFDSML_CSDR_I(ETA,NETA,PILOTMF,C1MEAN,MFMEAN,MFVAR,
     $     MFMEANGRAD,MFVARGRAD,DT,CHIMEAN,CSDRI_PMFDSML)

      CALL CPU_TIME(TFINISH)
      
      CPUT_PMFDSML = TFINISH - TSTART

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                       *****BETA APPROACH*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      WRITE(*,50)'============================================='
      WRITE(*,50) 'BETA PROBABILITY DENSITY FUNCTION APPROACH'
      WRITE(*,50)'============================================='

      CALL CPU_TIME(TSTART)

C     COMPUTE THE BETA-PDF
      WRITE(*,50) '*PDF'
      CALL BETA_PDF(ETA,NETA,MFMEAN,MFVAR,PDF_BETA)

C     COMPUTE THE BETA-CV
      WRITE(*,50) '*CV'
      CALL BETA_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,DT,
     $     VEL,CV_BETA)

C     COMPUTE THE [[HOMOGENEOUS]] BETA-CSDR
      WRITE(*,50) '*CSDR-HOMOGENEOUS'
      CALL BETA_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHIMEAN,CSDRH_BETA)

C     COMPUTE THE [[INHOMOGENEOUS]] BETA-CSDR
      WRITE(*,50) '*CSDR-INHOMOGENEOUS'
      CALL BETA_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,CHIMEAN,CSDRI_BETA)

      CALL CPU_TIME(TFINISH) 

      CPUT_BETA = TFINISH -TSTART

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                   *****UNIFORM WITH DELTA PDF*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      WRITE(*,50)'============================================='
      WRITE(*,50) 'DELTA PROBABILITY DENSITY FUNCTION          '
      WRITE(*,50) 'UNIFORM DISTRIBUTIONS                       '
      WRITE(*,50)'============================================='

C     COMPUTE THE DELTA-PDF
      WRITE(*,50) '*DELTA PDF'
      CALL DELTA_PDF(ETA,NETA,MFMEAN,PDF_DELTA)

C     COMPUTE THE UNIFORM CV
      WRITE(*,50) '*UNIFORM CV'
      CALL UNFRM_CV(ETA,NETA,VEL,CV_UNFRM)

C     COMPUTE THE UNIFORM CSDR
      WRITE(*,50) '*UNIFORM CSDR'
      CALL UNFRM_CSDR(ETA,NETA,CHIMEAN,CSDR_UNFRM)

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                     *****CLIPPED GAUSSIAN*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      WRITE(*,50)'============================================='
      WRITE(*,50) 'C.G. PROBABILITY DENSITY FUNCTION APPROACH  '
      WRITE(*,50)'============================================='

      CALL CPU_TIME(TSTART)

C     COMPUTE THE CG-PDF
      WRITE(*,50) '*PDF'
      CALL CG_PDF(ETA,NETA,MFMEAN,MFVAR,PDF_CG)

C     COMPUTE THE CG-CV
      WRITE(*,50) '*CV'
      CALL CG_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,MFVARGRAD,
     $     DT,VEL,CV_CG)

C     COMPUTE THE [[HOMOGENEOUS]] CG-CSDR
      WRITE(*,50) '*CSDR-HOMOGENEOUS'
      CALL CG_CSDR_H(ETA,NETA,MFMEAN,MFVAR,CHIMEAN,CSDRH_CG)

C     COMPUTE THE [[INHOMOGENEOUS]] CG-CSDR
      WRITE(*,50) '*CSDR-INHOMOGENEOUS'
c      CALL CG_CSDR_I(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,
c     $     MFVARGRAD,DT,CHIMEAN,CSDRI_CG)

      CALL CPU_TIME(TFINISH)

      CPUT_CG = TFINISH -TSTART

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                       *****OTHER MODELS*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================

      WRITE(*,50)'============================================='
      WRITE(*,50) 'OTHER SUBMODELS                             '
      WRITE(*,50)'============================================='

      WRITE(*,50) '*CV-LINEAR MODEL'
      CALL LINEAR_CV(ETA,NETA,MFMEAN,MFVAR,MFMEANGRAD,DT,VEL,CV_LIN)

      WRITE(*,50) '*CSDR-MODEL OF GIRIMAJI USING BETA PDF'
      CALL GIR_BETA_CSDR(ETA,NETA,MFMEAN,MFVAR,CHIMEAN,CSDRH_GIR_BETA)

      WRITE(*,50) '*CSDR-AMC USUNG BETA PDF'
      CALL AMC_BETA_CSDR(ETA,NETA,MFMEAN,MFVAR,CHIMEAN,CSDRH_AMC_BETA)

      WRITE(*,50) '*CSDR-AMC USUNG PMF PDF'
      CALL AMC_PMFSML_CSDR(ETA,NETA,MFMEAN,MFVAR,CHIMEAN,
     $     CSDRH_AMC_PMFSML)

C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C                       *****VALIDATION*****       
C==========================================================================
C==========================================================================
C==========================================================================
C==========================================================================
C
C     COMPARISON OF THE MEAN VALUE OF CHIMEAN WITH THE INTEGRAL
C     OF THE PDF-WEIGHTED CSDR OVER MIXTURE FRACTION SPACE
C
      WRITE(*,*)' '
      WRITE(*,50)'============================================='
      WRITE(*,*)' '
      WRITE(*,*)'==========='
      WRITE(*,*)'VALIDATION |'
      WRITE(*,*)'==========='
      WRITE(*,*)' '
      WRITE(*,*)'NOTE: THE ERRORS CALCULATED BELOW CAN BE LARGE'
      WRITE(*,*)'      WHEN A COARSE ETA GRID IS EMPLOYED'
      WRITE(*,*)'      (SIMPLE TRAPEZOIDAL INTEGRATION IS USED)'
      WRITE(*,*)' '
      WRITE(*,*)'===================================='
      WRITE(*,*)'CONDITIONAL SCALAR DISSIPATION RATE |'
      WRITE(*,*)'===================================='
      WRITE(*,*)' '
      WRITE(*,*)'             1'
      WRITE(*,*)'CHECKING IF |  <CHI|ETA>P(ETA)dETA = CHIMEAN'
      WRITE(*,*)'             0'
      WRITE(*,*)' '
      WRITE(*,'(A)')'                    +----------------
     $-----------------------+'
      WRITE(*,'(A)')'                    |      VALUE     
     $  |     REL ERR (%)    |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
      WRITE(*,'(A,F17.12,A)')'| CHIMEAN MEAN      |',CHIMEAN,
     $     ' |          -         |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
C
      CALL UNCONDAVG(CSDRH_PMFDSML,PDF_PMFDSML,ETA,NETA,INTEGRAL)
      WRITE(*,'(A,F17.12,A,F17.12,A)')'| PMF-PDF CSDR-T (H)|', INTEGRAL,
     $     ' | ',RELERRP(CHIMEAN,INTEGRAL),'  |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
C      
      CALL UNCONDAVG(CSDRI_PMFDSML,PDF_PMFDSML,ETA,NETA,INTEGRAL)
      WRITE(*,'(A,F17.12,A,F17.12,A)')'| PMF-PDF CSDR-T (I)|', INTEGRAL,
     $     ' | ',RELERRP(CHIMEAN,INTEGRAL) ,'  |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
C
      CALL UNCONDAVG(CSDRH_PMFSML,PDF_PMFSML,ETA,NETA,INTEGRAL)
      WRITE(*,'(A,F17.12,A,F17.12,A)')'| PMF-PDF CSDR-B (H)|', INTEGRAL,
     $     ' | ',RELERRP(CHIMEAN,INTEGRAL),'  |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
C      
      CALL UNCONDAVG(CSDRI_PMFSML,PDF_PMFSML,ETA,NETA,INTEGRAL)
      WRITE(*,'(A,F17.12,A,F17.12,A)')'| PMF-PDF CSDR-B (I)|', INTEGRAL,
     $     ' | ',RELERRP(CHIMEAN,INTEGRAL) ,'  |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
C
      CALL UNCONDAVG(CSDRH_BETA,PDF_BETA,ETA,NETA,INTEGRAL)
      WRITE(*,'(A,F17.12,A,F17.12,A)')'| BETA-PDF CSDR (H) |', INTEGRAL,
     $     ' | ',RELERRP(CHIMEAN,INTEGRAL),'  |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
C
      CALL UNCONDAVG(CSDRI_BETA,PDF_BETA,ETA,NETA,INTEGRAL)
      WRITE(*,'(A,F17.12,A,F17.12,A)')'| BETA-PDF CSDR (I) |', INTEGRAL,
     $     ' | ',RELERRP(CHIMEAN,INTEGRAL),'  |'
      WRITE(*,'(A)')'+-------------------+----------------
     $--+--------------------+'
C     
      WRITE(*,'(A)')' B = BINARY        T = TRINARY'
      WRITE(*,'(A)')' H = HOMOGENEOUS   I = INHOMOGENEOUS'
      
C
C     CHECK IF THE PDFS INTEGRATE TO UNITY
C
      WRITE(*,*)' '
      WRITE(*,*)' '
      WRITE(*,*)'============================='
      WRITE(*,*)'PROBABILITY DENISTY FUNCTION |'
      WRITE(*,*)'============================='
      WRITE(*,*)' '
      WRITE(*,*)'             1'
      WRITE(*,*)'CHECKING IF |  P(ETA)dETA = 1'
      WRITE(*,*)'             0'
      WRITE(*,*)' '
      CALL PDFAREA(PDF_PMFDSML,ETA,NETA,AREA_PMFDSML)
      CALL PDFAREA(PDF_PMFSML,ETA,NETA,AREA_PMFSML)
      CALL PDFAREA(PDF_BETA,ETA,NETA,AREA_BETA)
      CALL PDFAREA(PDF_CG,ETA,NETA,AREA_CG)
      WRITE(*,'(A)')'+--------------+-------------------------------+'
      WRITE(*,'(A)')'|   APPROACH   |     INTEGRAL    | REL ERR (%) |'
      WRITE(*,'(A)')'+--------------+-------------------------------+'
      WRITE(*,'(A,F12.8,A,F10.7,A)'),'|     PMF-T    |  ', 
     $     AREA_PMFDSML ,'   | ',RELERRP(1.0D0,AREA_PMFDSML),'  |'
      WRITE(*,'(A)')'+--------------+-------------------------------+'
      WRITE(*,'(A,F12.8,A,F10.7,A)'),'|     PMF-B    |  ', 
     $     AREA_PMFSML ,'   | ',RELERRP(1.0D0,AREA_PMFSML),'  |'
      WRITE(*,'(A)')'+--------------+-------------------------------+'
      WRITE(*,'(A,F12.8,A,F10.7,A)'),'|     BETA     |  ', 
     $     AREA_BETA,'   | ',RELERRP(1.0D0,AREA_BETA),'  |'
      WRITE(*,'(A)')'+--------------+-------------------------------+'
      WRITE(*,'(A,F12.8,A,F10.7,A)'),'|      CG      |  ', 
     $     AREA_CG,'   | ',RELERRP(1.0D0,AREA_CG),'  |'
      WRITE(*,'(A)')'+--------------+-------------------------------+'
C
C     CPU TIME
C
      WRITE(*,*)' '
      WRITE(*,*)' '
      WRITE(*,'(A)')'+--------------+----------------+'
      WRITE(*,'(A)')'|        CPU STATISTICS         |'
      WRITE(*,'(A)')'+--------------+----------------+'
      WRITE(*,'(A)')'|   APPROACH   | CPU TIME (SEC) |'
      WRITE(*,'(A)')'+--------------+----------------+'
      WRITE(*,'(A,F8.4,A)')'|     PMF-T    |  ', CPUT_PMFDSML ,'      |'
      WRITE(*,'(A)')'+--------------+----------------+'
      WRITE(*,'(A,F8.4,A)')'|     PMF-B    |  ', CPUT_PMFSML ,'      |'
      WRITE(*,'(A)')'+--------------+----------------+'
      WRITE(*,'(A,F8.4,A)')'|     BETA     |  ', CPUT_BETA,'      |'
      WRITE(*,'(A)')'+--------------+----------------+'
      WRITE(*,*)' '
      WRITE(*,'(1X,A,F8.4)')'PMF-B TO BETA CPU TIME RATIO  =',
     $     CPUT_PMFSML/CPUT_BETA
      WRITE(*,'(1X,A,F8.4)')'PMF-T TO PMF-B CPU TIME RATIO =  ',
     $     CPUT_PMFDSML/CPUT_PMFSML
      WRITE(*,*)' '
      WRITE(*,*)' '
C
C     WRITE RESULTS TO FILES
C
      OPEN(1,FILE='res/eta.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(2,FILE='res/pmf1PDF.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(3,FILE='res/pmf1CV.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(4,FILE='res/pmf1CSDRH.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(5,FILE='res/pmf1CSDRI.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(6,FILE='res/betaPDF.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(7,FILE='res/betaCV.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(8,FILE='res/betaCSDRH.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(9,FILE='res/betaCSDRI.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(10,FILE='res/cgPDF.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(11,FILE='res/cgCV.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
       OPEN(111,FILE='res/cgCSDRH.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(14,FILE='res/deltaPDF.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(15,FILE='res/unfrmCV.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(16,FILE='res/unfrmCSDR.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(17,FILE='res/girCSDRH.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(18,FILE='res/amcBetaCSDRH.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(19,FILE='res/amcPmfCSDRH.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(20,FILE='res/linCV.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(21,FILE='res/pmf2PDF.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(22,FILE='res/pmf2CV.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(23,FILE='res/pmf2CSDRH.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      OPEN(24,FILE='res/pmf2CSDRI.dat',FORM='FORMATTED',
     $     STATUS='UNKNOWN')
      DO IETA = 1,NETA
         WRITE(1,110)  ETA(IETA)
         WRITE(2,110)  PDF_PMFSML(IETA)
         WRITE(3,120)  CV_PMFSML(1:3,IETA)
         WRITE(4,110)  CSDRH_PMFSML(IETA)
         WRITE(5,110)  CSDRI_PMFSML(IETA)
         WRITE(6,110)  PDF_BETA(IETA)
         WRITE(7,120)  CV_BETA(1:3,IETA)
         WRITE(8,110)  CSDRH_BETA(IETA)
         WRITE(9,110)  CSDRI_BETA(IETA)
         WRITE(10,110) PDF_CG(IETA)
         WRITE(11,120) CV_CG(1:3,IETA)
         WRITE(111,120) CSDRH_CG(IETA)
         WRITE(14,110) PDF_DELTA(IETA)
         WRITE(15,120) CV_UNFRM(1:3,IETA)
         WRITE(16,110) CSDR_UNFRM(IETA)
         WRITE(17,110) CSDRH_GIR_BETA(IETA)
         WRITE(18,110) CSDRH_AMC_BETA(IETA)
         WRITE(19,110) CSDRH_AMC_PMFSML(IETA)
         WRITE(20,120) CV_LIN(1:3,IETA)
         WRITE(21,110) PDF_PMFDSML(IETA)
         WRITE(22,120) CV_PMFDSML(1:3,IETA)
         WRITE(23,110) CSDRH_PMFDSML(IETA)
         WRITE(24,110) CSDRI_PMFDSML(IETA)
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
      CLOSE(10)
      CLOSE(11)
      CLOSE(111)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
      CLOSE(17)
      CLOSE(18)
      CLOSE(19)
      CLOSE(20)
      CLOSE(21)
      CLOSE(22)
      CLOSE(23)
      CLOSE(24)
C
      WRITE(*,*)' '
      WRITE(*,'(A)') 'THE RESULTS ARE WRITTEN TO ../res/' 
      WRITE(*,*)' '

      STOP
      END

      
