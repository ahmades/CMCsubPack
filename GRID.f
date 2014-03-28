C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE ETAGRID_BOUNDARYCLUSTER(BETA,NETA,ETA)
C================================================================================== 
C     PURPOSE: CREATES A MIXTURE FRACTION GRID ETA(NETA) CLUSTERED NEAR THE 
C              BOUNDARIES OF THE MIXTURE FRACTION SPACE GIVEN THE STRETCHING 
C              PARAMETER BETA.
C==================================================================================
C            VARIABLE  DESCRIPTION                              DATA TYPE
C            --------  -----------                              --------- 
C
C     INPUT:
C
C            BETA      STRETCHING PARAMETER (BETA > 1)          DOUBLE PRECISION
C                      GRADUALLY INCREASING BETA DECREASES
C                      THE LEVEL OF REFINEMENT NEAR THE 
C                      BOUNDARIES UP TO THE POINT WHERE
C                      THE GRID BECOMES UNIFORM. INCREASING 
C                      BETA FURHTER WILL NOT HAVE ANY EFFECT
C            NETA      NUMBER OF GRID POINTS                    INTEGER
C
C     OUTPUT:
C
C            ETA       MIXTURE FRACTION GRID CLUSTERED          DOUBLE PRECISION 
C                      NEAR THE BOUNDARIES                      (ARRAY,SIZE=NETA)
C================================================================================== 
C     SOURCE: TRANSFORMATION 2 IN [TANNEHILL, J.C., ANDERSON, D.A., PLETCHER, 
C             R.H., 1997. COMPUTATIONAL FLUID MECHANICS AND HEAT TRANSFER. IN: 
C             MINKOWYCZ, W.J., SPARROW, E.M. (EDS.), COMPUTATIONAL AND PHYSICAL 
C             PROCESSES IN MECHANICS AND THERMAL SCIENCES. TAYLOR & FRANCIS, 
C             LONDON, UK, P. 336].
C================================================================================== 
      IMPLICIT NONE
      INTEGER NETA,IETA
      DOUBLE PRECISION ETA(NETA),ETAU(NETA),ETAMIN,ETAMAX,
     $     ALPHA,BETA,DELTA,H,TERM
      PARAMETER(ALPHA = 0.5D0)
      IF(BETA.LE.1.0D0) THEN
         WRITE(*,*)'ERROR IN SUBROUTINE ETAGRID_BOUNDARYCLUSTER'
         WRITE(*,*)'BETA MUST BE STRICTLY GREATER THAN 1.0'
         WRITE(*,*)'COMPUTATIONS ABORTED'
         STOP
      ENDIF
      ETAMIN = 0.0D0
      ETAMAX = 1.0D0
      H = ETAMAX-ETAMIN
      CALL ETAGRID_UNIFORM(NETA,ETAU)
      DO IETA = 2,NETA-1
         TERM = ((BETA+1.0D0)/(BETA-1.0D0))**((ETAU(IETA)-ALPHA)
     $        /(1.0D0-ALPHA))
         ETA(IETA) = H*((BETA+2.0D0*ALPHA)*TERM - BETA + 2.0D0*ALPHA)
     $        /((2.0D0*ALPHA+1.0D0)*(1.0D0+TERM))
      ENDDO
      ETA(1)  = ETAMIN
      ETA(NETA) = ETAMAX
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE ETAGRID_INTERIORCLUSTER(ETAC,TAU,NETA,ETA)
C================================================================================== 
C     PURPOSE: CREATES A MIXTURE FRACTION GRID ETA(NETA) CLUSTERED AROUND ETAC 
C              GIVEN THE STRETCHING PARAMETER TAU.
C==================================================================================
C            VARIABLE  DESCRIPTION                   DATA TYPE
C            --------  -----------                   --------- 
C
C     INPUT:
C
C            ETAC      POINT WHERE THE GRID IS       DOUBLE PRECISION
C                      TO BE CLUSTERED
C            TAU       STRETCHING PARAMETER          DOUBLE PRECISION 
C                      TAU = 0 STRETCHING IS OFF      
C                      TAU > 0 STRETCHING IS ON
C                              LARGER TAU INCREASES 
C                              THE GRID DENSITY 
C                              ABOUT ETAC        
C            NETA      NUMBER OF GRID POINTS         INTEGER
C
C     OUTPUT:
C
C            ETA       MIXTURE FRACTION GRID         DOUBLE PRECISION 
C                      CLUSTERED ABOUT ETAC          (ARRAY,SIZE=NETA)
C================================================================================== 
C     SOURCE: TRANSFORMATION 3 IN [TANNEHILL, J.C., ANDERSON, D.A., PLETCHER, 
C             R.H., 1997. COMPUTATIONAL FLUID MECHANICS AND HEAT TRANSFER. IN: 
C             MINKOWYCZ, W.J., SPARROW, E.M. (EDS.), COMPUTATIONAL AND PHYSICAL 
C             PROCESSES IN MECHANICS AND THERMAL SCIENCES. TAYLOR & FRANCIS, 
C             LONDON, UK. P. 337].
C================================================================================== 
      IMPLICIT NONE
      INTEGER NETA,IETA
      DOUBLE PRECISION ETA(NETA),ETAU(NETA),ETAC,ETAMIN,ETAMAX,
     $     TAU,B,H
      IF((ETAC.LE.0.0D0).OR.(ETAC.GE.1.0D0)) THEN
         WRITE(*,*)'ERROR IN SUBROUTINE ETAGRID_INTERIORCLUSTER'
         WRITE(*,*)'ETAC MUST BE IN ]0,1['
         WRITE(*,*)'COMPUTATIONS ABORTED'
         STOP
      ENDIF
      ETAMIN = 0.0D0
      ETAMAX = 1.0D0
      H = ETAMAX-ETAMIN
      CALL ETAGRID_UNIFORM(NETA,ETAU)
      B = (1.0D0/(2.0D0*TAU))
     $     *DLOG((1.0D0 + (DEXP(TAU)-1.0D0)*(ETAC/H))
     $     /(1.0D0 + (DEXP(-TAU)-1.0D0)*(ETAC/H)))
      DO IETA = 2,NETA-1
         ETA(IETA) = ETAC*(1.0D0 + DSINH(TAU*(ETAU(IETA)-B))
     $        /DSINH(TAU*B))
      ENDDO
      ETA(1)  = ETAMIN
      ETA(NETA) = ETAMAX
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE ETAGRID_UNIFORM(NETA,ETA)
C================================================================================== 
C     PURPOSE: CREATES A UNIFORMLY DISTRIBUTED MIXTURE FRACTION GRID ETA(NETA)
C==================================================================================
C            VARIABLE  DESCRIPTION                DATA TYPE
C            --------  -----------                --------- 
C
C     INPUT:
C
C            NETA        NUMBER OF GRID POINTS    INTEGER
C
C     OUTPUT:
C
C            ETA         UNIFORMLY DISTRIBUTED    DOUBLE PRECISION 
C                        MIXTURE FRACTION GRID    (ARRAY,SIZE=NETA)
C================================================================================== 
      IMPLICIT NONE
      INTEGER IETA,NETA
      DOUBLE PRECISION ETA(NETA),ETAMIN,ETAMAX,DETA
      ETAMIN = 0.0D0
      ETAMAX = 1.0D0
      DETA = (ETAMAX-ETAMIN)/DBLE(NETA-1)
      ETA(1) = ETAMIN
      ETA(NETA) = ETAMAX
      DO IETA = 2,NETA-1
         ETA(IETA) = ETA(IETA-1)+DETA
      ENDDO
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 

      SUBROUTINE ETAGRID_EXPONENTIAL(MS,NETA,ETA)
C================================================================================== 
C     PURPOSE: CREATES AN EXPONENTIALLY DISTRIBUTED MIXTURE FRACTION GRID ETA(NETA)  
C              GIVEN A STRETCHING PARAMETER MS.A UNIFORM GRID IS RETURNED IF MS IS 
C              ZERO. 
C==================================================================================
C            VARIABLE  DESCRIPTION                             DATA TYPE
C            --------  -----------                             --------- 
C
C     INPUT:
C
C            MS        STRETCHING PARAMETER (MAPPING           DOUBLE PRECISION
C                      STRENGTH)
C                      MS = 0 GRID IS UNIFORM
C                      MS > 0 GRID INCREASES EXPONENTIALLY
C                             LARGER MS YIELDS MORE POINTS 
C                             NEAR ETAMIN
C                      MS < 0 GRID DECREASES EXPONENTIALLY
C                             SMALLER MS YIELDS MORE POINTS 
C                             NEAR ETAMAX
C            NETA     NUMBER OF GRID POINT                    INTEGER
C
C     OUTPUT:
C
C            ETA      EXPONENTIALLY DISTRIBUTED               DOUBLE PRECISION 
C                     MIXTURE FRACTION GRID                   (ARRAY,SIZE=NETA)
C==================================================================================
      IMPLICIT NONE
      INTEGER NETA,IETA
      DOUBLE PRECISION ETA(NETA),MS,ETAMIN,ETAMAX,LINSPACE(NETA-1),
     $     DETA,H
      ETAMIN = 0.0D0
      ETAMAX = 1.0D0
      H = ETAMAX-ETAMIN
      DETA  = (ETAMAX-ETAMIN)/DBLE(NETA-1)
      DO IETA = 1,NETA-1
         IF(MS.EQ.0) THEN
            ETA(IETA) = ETAMIN + DBLE(IETA-1)*DETA
         ELSE
            LINSPACE(IETA) = ETAMIN+DBLE(IETA-1)*DETA
            ETA(IETA) = ETAMIN + H*(DEXP(MS*LINSPACE(IETA))-1.0D0)
     $           /(DEXP(MS)-1.0D0)
         ENDIF
      ENDDO
      ETA(1)  = ETAMIN
      ETA(NETA) = ETAMAX
      RETURN 
      END
