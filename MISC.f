C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE FINDDIFFDELTA(MFMEAN,MFVAR,H_MEAN,H_VAR)
C================================================================================== 
C     PURPOSE: NUMERICAL DIFFERENTIATION WITH RESPECT TO THE MIXTURE FRACTION MEAN
C              (VARIANCE) REQUIRES THE SPECIFICATION OF A RANGE ABOUT THE MEAN 
C              (VARIANCE). 
C              GIVEN THE MIXTURE FRACTION MEAN (MFMEAN) AND VARIANCE (MFVAR),
C              THIS SUBROUTINE COMPUTES H_MEAN AND H_VAR WHICH ARE NECESSARY TO 
C              DEFINE THE RANGE REQUIRED FOR DIFFERENTIATION.              
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
      DOUBLE PRECISION MFMEAN,MFVAR,H_MEAN,H_VAR
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE,
     $     D_FOUR,D_EIGHT
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS2/D_FOUR,D_EIGHT

C      H_MEAN = (1.0D0 - DSQRT(D_ONE - 4.0D0*MFVAR))/2.0D0
C      H_VAR  = MFVAR/2.0D0

      H_MEAN = (D_ONE - DSQRT(D_ONE - D_FOUR*MFVAR))/D_TWO
      H_VAR  = MFVAR/D_TWO

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CHECKPARMS(MFMEAN,MFVAR)
C================================================================================== 
C     PURPOSE: CHECKS THE VALIDITY OF THE VALUES OF THE MIXTURE FRACTION MEAN,
C              THE MIXTURE FRACTION AND VARIANCE, AND THE INTENSITY OF SEGREGATION
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
      INCLUDE 'formats.h'
      DOUBLE PRECISION MFMEAN,MFVAR,ISEG
      DOUBLE PRECISION PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      COMMON/DBLECONSTANTS/PI,D_ZERO,D_ONE,D_HALF,D_TWO,D_THREE
      LOGICAL VERBOSE
      COMMON/LOGICALVARS/VERBOSE

      ISEG = MFVAR/(MFMEAN*(D_ONE-MFMEAN))

      IF(VERBOSE) THEN
         WRITE(*,50)'============================================='
         WRITE(*,50)'PARAMETERS CHECK:'
         WRITE(*,50)'============================================='
         WRITE(*,100)'MIXTURE FRACTION MEAN     =', MFMEAN
         WRITE(*,100)'MIXTURE FRACTION VARIANCE =', MFVAR
         WRITE(*,100)'INTENSITY OF SEGREGATION  =', ISEG
         WRITE(*,50)' '
      ENDIF
      IF((MFMEAN.LT.D_ZERO).OR.(MFMEAN.GT.D_ONE)) THEN
         WRITE(*,50)'============================================='
         WRITE(*,50)'WRONG INPUT:'
         WRITE(*,50)'MIXTURE FRACTION MEAN IS NOT BETWEEN 0 AND 1.'
         WRITE(*,50)'COMPUTATIONS ABORTED.'
         WRITE(*,50)'============================================='
         STOP
      ELSE
         IF(VERBOSE)
     $        WRITE(*,50)'MIXTUR FRACTION MEAN      -> VALUE IS VALID.'
      ENDIF

      IF(MFVAR.LT.D_ZERO) THEN
         WRITE(*,50)'============================================='
         WRITE(*,50)'WRONG INPUT:'
         WRITE(*,50)'MIXTURE FRACTION VARIANCE IS LESS THAN 0.'
         WRITE(*,50)'COMPUTATIONS ABORTED.'
         WRITE(*,50)'============================================='
         STOP
      ELSE
         IF(VERBOSE)
     $        WRITE(*,50)'MIXTUR FRACTION VARAINCE  -> VALUE IS VALID.'   
      ENDIF

      IF((ISEG.LT.D_ZERO).OR.(ISEG.GT.D_ONE)) THEN
         WRITE(*,50)'WRONG INPUT:'
         WRITE(*,50)'INTENSITY OF SEGREGATION IS NOT BETWEEN 0 AND 1.'
         WRITE(*,50)'COMPUTATIONS ABORTED.'
         WRITE(*,50)'============================================='
         STOP
      ELSE
         IF(VERBOSE)
     $        WRITE(*,50)'INTENSITY OF SEGREGATION  -> VALUE IS VALID.'
      ENDIF
      IF(VERBOSE)
     $     WRITE(*,50)'============================================='

      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      SUBROUTINE PDFAREA(PDF,ETA,NETA,AREA)
C========================================================================================== 
C     PURPOSE: COMPUTES THE AREA UNDER THE PDF BY INTEGRATING PDF(ETA) OVER ETA SPACE. 
C
C     NOTE:    TRAPEZOIDAL INTEGRATION IS USED. LARGE ERRORS MAY ARISE
C              WHEN COARSE MIXTURE FRACTION GRIDS ARE EMPLOYED.
C==========================================================================================
C            VARIABLE   DESCRIPTION                     DATA TYPE
C            --------   -----------                     --------- 
C
C     INPUT:
C
C            PDF        PROBABILITY DENSITY FUNCTION    DOUBLE PRECISION (ARRAY,SIZE=NETA)      
C            ETA        MIXTURE FRACTION GRID           DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA AND CSDRI           INTEGER
C            
C
C     OUTPUT:
C
C            AREA       AREA UNDER PDF                  DOUBLE PRECISION
C==========================================================================================       
      IMPLICIT NONE
      INTEGER IETA,NETA
      DOUBLE PRECISION PDF(NETA),ETA(NETA)
      DOUBLE PRECISION AREA
      DOUBLE PRECISION TRAP
      EXTERNAL TRAP
      AREA = TRAP(PDF,ETA,NETA)
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      SUBROUTINE UNCONDAVG(X_COND,PDF,ETA,NETA,X_UNCOND)
C========================================================================================== 
C     PURPOSE: COMPUTES THE UNCONDITIONAL (FAVRE) AVERAGE OF QUANTITY X 
C              BY INTEGRATING <X|ETA>*PDF(ETA) OVER ETA SPACE. 
C
C     NOTE:    TRAPEZOIDAL INTEGRATION IS USED. LARGE ERRORS MAY ARISE
C              WHEN COARSE MIXTURE FRACTION GRIDS ARE EMPLOYED.
C==========================================================================================
C            VARIABLE   DESCRIPTION                         DATA TYPE
C            --------   -----------                         --------- 
C
C     INPUT:
C    
C
C            X_COND     CONDITIONAL AVERAGE OF X  .         DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            PDF        PROBABILITY DENSITY FUNCTION        DOUBLE PRECISION (ARRAY,SIZE=NETA)      
C            ETA        MIXTURE FRACTION GRID               DOUBLE PRECISION (ARRAY,SIZE=NETA)
C            NETA       SIZE OF ETA AND CSDRI               INTEGER
C            
C
C     OUTPUT:
C
C            X_UNCOND   UNCONDITIONAL AVERAGE OF X          DOUBLE PRECISION
C==========================================================================================       
      IMPLICIT NONE
      INTEGER IETA,NETA
      DOUBLE PRECISION X_COND(NETA),PDF(NETA),ETA(NETA),X_UNCOND 
      DOUBLE PRECISION INTEGRAND(NETA)
      DOUBLE PRECISION TRAP
      EXTERNAL TRAP
      DO IETA = 1,NETA
         INTEGRAND(IETA) = X_COND(IETA)*PDF(IETA)
      ENDDO
      X_UNCOND = TRAP(INTEGRAND,ETA,NETA)
      RETURN
      END

C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION TRAP(F,X,N)
C========================================================================================== 
C     PURPOSE: COMPUTES THE INTEGRAL OF A DISCRETE FUNCTION USING TRAPEZOIDAL INTEGRATION
C
C     NOTE:    LARGE ERRORS MAY ARISE WHEN COARSE MIXTURE FRACTION GRIDS ARE EMPLOYED.
C==========================================================================================
C            VARIABLE   DESCRIPTION                DATA TYPE
C            --------   -----------                --------- 
C
C     INPUT:
C    
C
C            F          DISCRETE FUNCTION OF X     DOUBLE PRECISION (ARRAY,SIZE=N)      
C            X          X VALUES                   DOUBLE PRECISION (ARRAY,SIZE=N)
C            N          SIZE OF F AND X            INTEGER
C            
C
C     OUTPUT:
C
C            TRAP       INTEGRAL OF F(X)           DOUBLE PRECISION
C==========================================================================================
      IMPLICIT NONE
      INTEGER N,I
      DOUBLE PRECISION F(N),X(N)	
      TRAP = 0.0D0
      DO I = 1,N-1
         TRAP = TRAP + (X(I+1)-X(I)) * (F(I)+F(I+1))
      ENDDO
      TRAP = TRAP * 0.50D0

      RETURN
      END
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      DOUBLE PRECISION FUNCTION RELERRP(X1,X2)
C========================================================================================== 
C     PURPOSE: COMPUTES THE RELATIVE ERROR (%) OF X2 WITH RESPECT TO X1
C
C==========================================================================================
C            VARIABLE   DESCRIPTION                DATA TYPE
C            --------   -----------                --------- 
C
C     INPUT:
C    
C
C            X1          ACTUAL VALUE OF X        DOUBLE PRECISION
C            X2          ESTIMATED VALUE OF X     DOUBLE PRECISION
C
C     OUTPUT:
C
C            RELERRP     RELATIVE ERROR * 100     DOUBLE PRECISION
C==========================================================================================
      IMPLICIT NONE
      DOUBLE PRECISION X1,X2	
      RELERRP = 100.0D0*(X2-X1)/X1
      RETURN
      END
    
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 

      LOGICAL FUNCTION ISNAN(X)
C========================================================================================== 
C     PURPOSE: DETECTS THE OCCURRENCE OF NaN
C==========================================================================================
C            VARIABLE   DESCRIPTION                DATA TYPE
C            --------   -----------                --------- 
C
C     INPUT:
C    
C
C            X          VARIABLE TO BE TESTED      DOUBLE PRECISION
C
C     OUTPUT:
C
C            ISNAN      TEST RESULT                LOGICAL
C==========================================================================================
      IMPLICIT NONE
      DOUBLE PRECISION X
      IF(DABS(X).LE.HUGE(X))THEN 
         ISNAN = .FALSE.
      ELSE
         ISNAN = .TRUE.
      END IF
      RETURN
      END
     
C==================================================================================
C==================================================================================
C==================================================================================
C================================================================================== 
 
      LOGICAL FUNCTION ISINF(X)
C========================================================================================== 
C     PURPOSE: DETECTS THE OCCURRENCE OF INF
C==========================================================================================
C            VARIABLE   DESCRIPTION                DATA TYPE
C            --------   -----------                --------- 
C
C     INPUT:
C    
C
C            X          VARIABLE TO BE TESTED      DOUBLE PRECISION
C
C     OUTPUT:
C
C            IINF       TEST RESULT                LOGICAL
C========================================================================================== 
      IMPLICIT NONE
      DOUBLE PRECISION X
      IF ((X+1.0D0).EQ.X) THEN
         ISINF = .TRUE.
      ELSE
         ISINF = .FALSE.
      END IF
      RETURN
      END


C==========================================================================================
C==========================================================================================
C==========================================================================================
C========================================================================================== 
      SUBROUTINE LOCATE(X,NX,XX,INDX)
C========================================================================================== 
C     PURPOSE: GIVEN A MONOTONICALLY INCREASING REAL ARRAY X(1:NX) AND A REAL VALUE XX,
C              THIS SUBROUTINE RETURNS INDX SUCH THAT XX IS CLOSETS TO X(INDX).
C==========================================================================================
C            VARIABLE   DESCRIPTION                       DATA TYPE
C            --------   -----------                       --------- 
C
C     INPUT:
C
C            X          ARRAY TO SEARCH                   DOUBLE PRECISION (ARRAY,SIZE=NX)
C            NX         SIZE OF X                         INTEGER
C            XX         VALUE TO LOCATE IN ARRAY          DOUBLE PRECISION
C
C     OUTPUT:
C
C            INDX       INDEX IN ARRAY X SUCH THAT        DOUBLE PRECISION (ARRAY,SIZE=N)
C                       XX IS CLOSEST TO X(INDX)  
C========================================================================================== 
      INTEGER NX,INDX,IX
      DOUBLE PRECISION X(NX),XX
      IF((XX.LT.X(1)).OR.(XX.GT.X(NX))) THEN
         WRITE(*,*)'ERROR IN SUBROUTINE SEARCH: XX IS OUT OF BOUNDS'
         WRITE(*,*)'COMPUTATIONS ABORTED'
         STOP
      ENDIF
      IF(XX.EQ.X(NX)) THEN
         INDX = NX
      ELSE
         DO IX = 1,NX-1
            IF((XX.GE.X(IX)).AND.(XX.LT.X(IX+1))) INDX = IX                      
         ENDDO
         IF(DABS(XX-X(INDX+1)).LT.DABS(XX-X(INDX))) INDX = INDX+1
      ENDIF

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================
      SUBROUTINE MAXIMUM(ARRAY,NSIZE,MAXVAL,MAXVALIND)
      IMPLICIT NONE
      INTEGER J,NSIZE,MAXVALIND
      DOUBLE PRECISION ARRAY(NSIZE),MAXVAL
      MAXVAL = ARRAY(1)
      MAXVALIND = 1
      DO J = 2,NSIZE
         IF(MAXVAL.LT.ARRAY(J)) THEN
            MAXVAL = ARRAY(J)
            MAXVALIND = J
         ENDIF
      ENDDO
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================
      SUBROUTINE MINIMUM(ARRAY,NSIZE,MINVAL,MINVALIND)
      IMPLICIT NONE
      INTEGER J,NSIZE,MINVALIND
      DOUBLE PRECISION ARRAY(NSIZE),MINVAL
      MINVAL = ARRAY(1)
      DO J = 2,NSIZE
         IF(MINVAL.GT.ARRAY(J)) THEN 
            MINVAL = ARRAY(J)
            MINVALIND = J
         ENDIF
      ENDDO
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================
      SUBROUTINE PROD1D(X,NX,Y,NY,Z)
C     RETURNS THE ELEMENTWISE PRODUCT OF 1D VECTORS X(NX) AND Y(NY) IN Z(NX).
C     AN ERROR IS PRODUCED IF NX.NE.NY
      IMPLICIT NONE
      INTEGER I,NX,NY
      DOUBLE PRECISION X(NX),Y(NY),Z(NX)
      IF(NX.NE.NY) THEN
         WRITE(*,*)'ERROR IN SUBROUTINE PROD1D'
         WRITE(*,*)'THE ARRAYS X AND Y ARE NOT OF THE SAME SIZE'
         WRITE(*,*)'COMPUTAIONS ABORTED'
         STOP
      ENDIF
      DO I = 1,NX
         Z(I) = X(I)*Y(I)
      ENDDO
      
      RETURN
      END


C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================
      SUBROUTINE PROGRESSPC(J,N)
C========================================================================================== 
C     PURPOSE: SHOWS THE PROGRESS (PERCENTAGE) AT INDEX J OF A DO LOOP THAT EXECUTES N TIMES
C==========================================================================================
      IMPLICIT NONE
      INTEGER I,J,K,N
      CHARACTER(LEN=17) BAR
C     SET UP BAR
      BAR = " PROGRESS: ???% "
      DO I = 1,N
         BAR = BAR // " "
      ENDDO
      BAR(17:17) = " "
C     FILL OUT BAR
      WRITE(UNIT=BAR(12:14),FMT="(I3)") 100*J/N
C     PRINT BAR
      WRITE(UNIT=6,FMT="(A1,A17)",ADVANCE="NO") CHAR(13),BAR
      IF (J.NE.N) THEN
         FLUSH(UNIT = 6)
      ELSE
         WRITE(UNIT = 6,FMT=*)
      ENDIF
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================
      SUBROUTINE PROGRESSBAR(J,N)
C========================================================================================== 
C     PURPOSE: SHOWS THE PROGRESS (PERCENTAGE AND BAR) AT INDEX J OF A DO LOOP THAT 
C              EXECUTES N TIMES
C==========================================================================================
      IMPLICIT NONE
      INTEGER I,J,K,N
      CHARACTER(LEN=N+18) BAR
      CHARACTER(LEN=1000) FRM
C     SET UP BAR
      BAR = "PROGRESS: ???% |"
      DO I = 1,N
         BAR = BAR // " "
      ENDDO
      BAR(N+17:N+17) = "|"
      BAR(N+18:N+18) = " "
C     SET UP FORMAT
      WRITE(FRM, '(I3)')N+18
      FRM = ADJUSTL(FRM)
      FRM = "(A1,A" // TRIM(FRM) // ")"
C     FILL OUT BAR
      WRITE(UNIT=BAR(11:13),FMT="(I3)") 100*J/N
      DO K = 1, J
         BAR(16+K:16+K) = "*"
      ENDDO
C     PRINT BAR
      WRITE(UNIT=6,FMT=FRM(1:LEN_TRIM(FRM)),ADVANCE="NO") CHAR(13),BAR
      IF (J.NE.N) THEN
         FLUSH(UNIT = 6)
      ELSE
         WRITE(UNIT = 6,FMT=*)
      ENDIF
      RETURN
      END
