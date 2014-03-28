C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CHEBDERIV(A,B,N,F,X,ORD,DFDX)
      IMPLICIT NONE
      INTEGER N,M,ORD
      DOUBLE PRECISION A,B,F,X,DFDX,C(N),CD(N)
      EXTERNAL F
      M = N
      CALL CHEBYSHEV_COEFFICIENTS(A,B,N,F,C)
      CALL CHEBYSHEV_DERIVATIVE(A,B,N,C,ORD,CD)    
      CALL CHEBYSHEV_EVALUATION(A,B,CD,N,M,X,DFDX)
      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CHEBYSHEV_COEFFICIENTS(A,B,N,F,C)
C==========================================================================================
C     CHEBYSHEV_COEFFICIENTS DETERMINES CHEBYSHEV INTERPOLATION COEFFICIENTS.
C     PARAMETERS:
C     INPUT,DOUBLE PRECISION A,B,THE DOMAIN OF DEFINITION.
C     INPUT,INTEGER N,THE ORDER OF THE INTERPOLANT.
C     INPUT,DOUBLE PRECISION,EXTERNAL :: F (X),AN EXTERNAL FUNCTION.
C     OUTPUT,DOUBLE PRECISION C(N),THE CHEBYSHEV COEFFICIENTS.
C==========================================================================================
      IMPLICIT NONE
      INTEGER N,I,J
      DOUBLE PRECISION A,B,X,ANGLE,C(N),F,FX(N)
      EXTERNAL F
      DOUBLE PRECISION PI
      PARAMETER (PI = 3.141592653589793238462D0)

      DO I = 1,N
        FX(I) = F(0.5D0*(A+B) 
     $        + DCOS(DBLE(2*I-1)*PI/DBLE(2*N))*0.5D0*(B-A))
      ENDDO

      DO I = 1,N
        C(I) = 0.0D0
        DO J = 1,N
          C(I) = C(I) + FX(J)*DCOS(DBLE((I-1)*(2*J-1))*PI/DBLE(2*N))
        ENDDO
        C(I) = 2.0D0 * C(I)/DBLE(N)
      ENDDO

      RETURN
      END

C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

c$$$      SUBROUTINE CHEBYSHEV_COEFFICIENTS(A,B,N,F,C)
c$$$C==========================================================================================
c$$$C     CHEBYSHEV_COEFFICIENTS DETERMINES CHEBYSHEV INTERPOLATION COEFFICIENTS.
c$$$C     PARAMETERS:
c$$$C     INPUT,DOUBLE PRECISION A,B,THE DOMAIN OF DEFINITION.
c$$$C     INPUT,INTEGER N,THE ORDER OF THE INTERPOLANT.
c$$$C     INPUT,DOUBLE PRECISION,EXTERNAL :: F (X),AN EXTERNAL FUNCTION.
c$$$C     OUTPUT,DOUBLE PRECISION C(N),THE CHEBYSHEV COEFFICIENTS.
c$$$C==========================================================================================
c$$$      IMPLICIT NONE
c$$$      INCLUDE 'omp_lib.h'
c$$$      INTEGER N,I,J
c$$$      DOUBLE PRECISION A,B,X,ANGLE,C(N),F,FX(N)
c$$$      EXTERNAL F
c$$$      DOUBLE PRECISION PI
c$$$      PARAMETER (PI = 3.141592653589793238462D0)
c$$$
c$$$      CALL OMP_SET_NUM_THREADS(8)
c$$$
c$$$C$OMP PARALLEL
c$$$C$OMP& DEFAULT(SHARED) 
c$$$C$OMP& PRIVATE (I)
c$$$C$OMP DO   
c$$$      DO I = 1,N
c$$$        FX(I) = F(0.5D0*(A+B) 
c$$$     $        + DCOS(DBLE(2*I-1)*PI/DBLE(2*N))*0.5D0*(B-A))
c$$$      ENDDO
c$$$C$OMP END DO
c$$$C$OMP END PARALLEL
c$$$
c$$$      DO I = 1,N
c$$$        C(I) = 0.0D0
c$$$        DO J = 1,N
c$$$          C(I) = C(I) + FX(J)*DCOS(DBLE((I-1)*(2*J-1))*PI/DBLE(2*N))
c$$$        ENDDO
c$$$        C(I) = 2.0D0 * C(I)/DBLE(N)
c$$$      ENDDO
c$$$
c$$$      RETURN
c$$$      END

c$$$C=============================================================================
c$$$C=============================================================================
c$$$C=============================================================================
c$$$C=============================================================================
c$$$
c$$$      SUBROUTINE CHEBYSHEV_COEFFICIENTS(A,B,N,F,C)
c$$$C==========================================================================================
c$$$C     CHEBYSHEV_COEFFICIENTS DETERMINES CHEBYSHEV INTERPOLATION COEFFICIENTS.
c$$$C     PARAMETERS:
c$$$C     INPUT,DOUBLE PRECISION A,B,THE DOMAIN OF DEFINITION.
c$$$C     INPUT,INTEGER N,THE ORDER OF THE INTERPOLANT.
c$$$C     INPUT,DOUBLE PRECISION,EXTERNAL :: F (X),AN EXTERNAL FUNCTION.
c$$$C     OUTPUT,DOUBLE PRECISION C(N),THE CHEBYSHEV COEFFICIENTS.
c$$$C==========================================================================================
c$$$      !Use omp_lib
c$$$      IMPLICIT NONE
c$$$      INCLUDE 'omp_lib.h'
c$$$      INTEGER N,I,J
c$$$      DOUBLE PRECISION A,B,X,ANGLE,C(N),F,FX(N),CTEMP
c$$$      EXTERNAL F
c$$$      DOUBLE PRECISION PI
c$$$      PARAMETER (PI = 3.141592653589793238462D0)
c$$$      
c$$$      CALL OMP_SET_NUM_THREADS(8)
c$$$
c$$$c$$$C$OMP PARALLEL
c$$$c$$$C$OMP& DEFAULT(SHARED) 
c$$$c$$$C$OMP& PRIVATE (I)
c$$$c$$$C$OMP DO ORDERED      
c$$$c$$$      DO I = 1,N
c$$$c$$$!$OMP ORDERED
c$$$c$$$         FX(I) = F(0.5D0*(A+B) 
c$$$c$$$     $        + DCOS(DBLE(2*I-1)*PI/DBLE(2*N))*0.5D0*(B-A))
c$$$c$$$!$OMP END ORDERED
c$$$c$$$      ENDDO
c$$$c$$$C$OMP END DO
c$$$c$$$C$OMP END PARALLEL
c$$$
c$$$C$OMP PARALLEL
c$$$C$OMP& DEFAULT(SHARED) 
c$$$C$OMP& PRIVATE (I)
c$$$C$OMP DO   
c$$$      DO I = 1,N
c$$$         FX(I) = F(0.5D0*(A+B) 
c$$$     $        + DCOS(DBLE(2*I-1)*PI/DBLE(2*N))*0.5D0*(B-A))
c$$$      ENDDO
c$$$C$OMP END DO
c$$$C$OMP END PARALLEL
c$$$
c$$$      
c$$$      DO I = 1,N
c$$$         CTEMP = 0.0D0
c$$$         DO J = 1,N
c$$$            CTEMP = CTEMP + FX(J)*DCOS(DBLE((I-1)*(2*J-1))*PI/DBLE(2*N))
c$$$         ENDDO
c$$$         C(I) = 2.0D0 * CTEMP/DBLE(N)
c$$$      ENDDO
c$$$C     
c$$$      RETURN
c$$$      END

C=============================================================================
C=============================================================================
C=============================================================================
C=============================================================================

      SUBROUTINE CHEBYSHEV_DERIVATIVE(A,B,N,C,ORD,CD) 
C=============================================================================
C     PURPOSE: CALCULATES THE CHEBYSHEV COEFFICIENTS OF THE DERIVATIVE OF 
C              ORDER ORD OF THE FUNCTION WHOSE COEFFICIENTS ARE GIVEN BY C
C     NOTES: 
C           [1] THIS SUBROUTINE IS BASED ON SUBROUTINE CHDER [W.H. PRESS, B.P. 
C               FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING, NUMERICAL RECIPES 
C               IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING, CAMBRIDGE 
C               UNIVERSITY PRESS, 1992].
C           [2] THIS SUBROUTINE IS ACCURATE FOR THE COMPUTATION OF THE FIRST
C               FEW DERIVATIVES OF A GIVEN FUNCTION. LARGE ERRORS MAY RESULT 
C               IN THE COMPUTATION OF HIGHER-ORDER DERIVATIVES DUE TO THE 
C               MAGNIFICATION OF ROUNDOFF ERROR BY MEANS OF RECUSRSION. THIS 
C               EFFECT IS INSIGNIFICANT IN THE FORTRAN IMPLEMENTATION OF 
C               PMFPACK SINCE ONLY THE FIRST AND SECOND DERIVATIVES OF SOME 
C               RELEVANT FUNCTIONS ARE NEEDED.
C=============================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            A         LOWER BOUND OF THE DOMAIN     DOUBLE PRECISION
C                      OF EVALUATION
C            B         UPPER BOUND OF THE DOMAIN     DOUBLE PRECISION
C                      OF EVALUATION
C            N         DEGREE OF APPROXIMATION       INTEGER
C            C         CHEBYSHEV COEFFICIENTS        DOUBLE PRECISION
C                      OF THE FUNCTION TO BE         (ARRAY,SIZE=N)
C                      DIFFERENTIATED
C            ORD       ORDER OF THE DERIVATIVE        INTEGER
C
C     OUTPUT:
C
C            CD       CHEBYSHEV COEFFICIENTS        DOUBLE PRECISION
C                     OF THE DERIVATIVE OF          (ARRAY,SIZE=N)
C                     ORDER ORD OF THE FUNCTION
C============================================================================== 
      IMPLICIT NONE
      INTEGER N,I,ORD,IORD
      DOUBLE PRECISION A,B,C(N),CD(N),CD1(N),CD2(N)

      DO I = 1,N
         CD1(I) = C(I)
      ENDDO

      DO IORD = 1,ORD
         CD2(N) = 0.0D0             
         CD2(N-1) = 2.0D0*DBLE(N-1)*CD1(N)
         DO I = N-2,1,-1
            CD2(I) = CD2(I+2) + 2.0D0*DBLE(I)*CD1(I+1)
         ENDDO
         DO I = 1,N                 
            CD2(I) = (2.0D0/(B-A))*CD2(I)
            CD1(I)  = CD2(I)
         ENDDO
      ENDDO

      DO I = 1,N
         CD(I)  = CD2(I)
      ENDDO
            
      RETURN
      END      

C=============================================================================
C=============================================================================
C=============================================================================
C=============================================================================
    
      SUBROUTINE CHEBYSHEV_EVALUATION(A,B,C,N,M,X,FX)
C=============================================================================
C     PURPOSE: EVALUATES A FUNCTION WHOSE CHEBYSHEV COEFFICIENTS ARE GIVEN BY C
C 
C     NOTE:    THIS SUBROUTINE IS BASED ON FUNCTION CHEBEV [W.H. PRESS, B.P. 
C              FLANNERY, S.A. TEUKOLSKY, W.T. VETTERLING, NUMERICAL RECIPES 
C              IN FORTRAN: THE ART OF SCIENTIFIC COMPUTING, CAMBRIDGE 
C              UNIVERSITY PRESS, 1992].
C=============================================================================
C            VARIABLE  DESCRIPTION                    DATA TYPE
C            --------  -----------                    --------- 
C
C     INPUT:
C
C            A         LOWER BOUND OF THE DOMAIN     DOUBLE PRECISION
C                      OF EVALUATION
C            B         UPPER BOUND OF THE DOMAIN     DOUBLE PRECISION
C                      OF EVALUATION
C            N         DEGREE OF APPROXIMATION       INTEGER
C            C         CHEBYSHEV COEFFICIENTS        DOUBLE PRECISION
C                      OF THE FUNCTION TO BE         (ARRAY,SIZE=N)
C                      DIFFERENTIATED
C            X         ORDER OF THE DERIVATIVE       DOUBLE PRECISION
C
C     OUTPUT:
C
C            FX       APPROXIMATION OF THE           DOUBLE PRECISION
C                     FUNCTION AT X 
C============================================================================== 
      IMPLICIT NONE
      INTEGER N,M,I
      DOUBLE PRECISION A,B,C(N),X,FX
      DOUBLE PRECISION R,D1,D2,S,XX

      R = (X-A)*(X-B)
      IF (R.GT.0.0D0) THEN
         WRITE(*,*)'ERROR IN CHEBYSHEV_EVALUATION: INVALID X VALUE'
         WRITE(*,*)'COMPUTATIONS ABORTED'
         STOP
      ENDIF
 
      D1 = 0.0D0
      D2 = 0.0D0
      XX = (2.0D0*X-A-B)/(B-A)   
     
      DO I = M,2,-1            
         S = D1
         D1 = 2.0D0*XX*D1 - D2+C(I)
         D2 = S
      ENDDO

      FX = C(1)/2.0D0 + D1*XX - D2

      RETURN
      END

C=============================================================================
C=============================================================================
C=============================================================================
C=============================================================================  



C==========================================================================================
C==========================================================================================
C==========================================================================================
C==========================================================================================

      SUBROUTINE CHEBDERIVADAPT(A,B,N,F,X,ORD,DFDX)
      IMPLICIT NONE
      INTEGER N,ORD,I,NTEMP,COUNTER,NMAX
      DOUBLE PRECISION A,B,X,DFDX,DFDXOLD,ERR,ERROLD,MP
      DOUBLE PRECISION F,D1MACH
      EXTERNAL F,D1MACH
      DOUBLE PRECISION , ALLOCATABLE :: C(:)
      DOUBLE PRECISION , ALLOCATABLE :: CD(:)

      INTEGER NCHEB2DER
      COMMON/NCHEB2DERCOM/NCHEB2DER
      
C     STORE NCOEFS
      NTEMP = N

      COUNTER = 0
      ERR = HUGE(1.0D0)
      MP = 100.0D0*D1MACH(4) !1.0D-13
      NMAX = 60 !50 !80
      
      DO WHILE(ERR.GT.MP)

         NCHEB2DER = N

         IF(COUNTER.EQ.0) THEN
            DFDXOLD = HUGE(1.0D0)
            ERROLD = HUGE(1.0D0)
         ELSE
            DFDXOLD = DFDX
            ERROLD = ERR
         ENDIF

         ALLOCATE (C(N))
         ALLOCATE (CD(N))

         CALL CHEBYSHEV_COEFFICIENTS(A,B,N,F,C)
         CALL CHEBYSHEV_DERIVATIVE(A,B,N,C,ORD,CD)    
         CALL CHEBYSHEV_EVALUATION(A,B,CD,N,N,X,DFDX)
                  
         ERR = DABS(DFDX-DFDXOLD)
         
         IF(COUNTER.GT.0)WRITE(*,*)'COUNT =',COUNTER,'N =',N,'ERR =',ERR
                  
         N = N+10 !5 !1 !2 !10 !5

         COUNTER = COUNTER+1

         DEALLOCATE (C)
         DEALLOCATE (CD)

C     IF ERROR IS INCREASING CHOOSE BEST (PREVIOUS) ANSWER
         IF((COUNTER.GT.0).AND.(ERR.GT.ERROLD)) THEN
            DFDX = DFDXOLD
            GOTO 100 
         ENDIF

         IF(N.GT.NMAX) GOTO 100 !EXIT

      ENDDO

 100  CONTINUE

C     RESET NCOEFS
      N = NTEMP
      
      RETURN
      END
