!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the position of the nodes
!The formula is y = sqrt(x)*(1-x)
!  NELEM    Figure of elements (odd) (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  T        Thickness                (IN)
!  C        Chord                    (IN)
!It use an iterative method to make the elements of the same lenght
SUBROUTINE GEOMWING(NELEM, XNODE, C, T)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL(KIND=8), INTENT(IN) :: T, C
      REAL(KIND=8), DIMENSION(NELEM,2), INTENT(OUT) :: XNODE
      REAL(KIND=8), DIMENSION(NELEM/2) :: DS, DX, PHI
      REAL(KIND=8), PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I, J, NITER = 15
      REAL(KIND=8) :: X, Y, SSQ, SIGMA, MEAN !, THETA
      !First and mid nodes
      XNODE(1,1) = C
      XNODE(1,2) = 0.
      XNODE(NELEM/2+1,1) = 0.
      XNODE(NELEM/2+1,2) = 0.
      XNODE(NELEM/2+2,1) = 0.
      XNODE(NELEM/2+2,2) = 0.
      !Symmetric profile
      DX(:) = C/REAL(NELEM/2,8)
      DO I = 1,NITER
       DO J = 1,NELEM/2
        IF (J .LT. NELEM/2) THEN
         X = C-SUM(DX(:J))
         XNODE(J+1,1) = X
         Y = T / DSQRT(16.D0/27.D0) * DSQRT(X/C)*(1.-X/C)
         XNODE(J+1,2) = Y
         XNODE(NELEM-J+1,1) = X
         XNODE(NELEM-J+1,2) = -Y
        END IF
        DS(J) = DSQRT( (XNODE(J+1,1)- XNODE(J,1))**2 + (XNODE(J+1,2)-XNODE(J,2))**2 )
        PHI(J) = DATAN2( XNODE(J,2) - XNODE(J+1,2) , XNODE(J,1) - XNODE(J+1,1) )
       END DO
       !Mean and stdev
       MEAN = SUM(DS)/DFLOAT(NELEM/2)
       SSQ = 0.
       DO J = 1,NELEM/2
        SSQ = SSQ + (DS(J)-MEAN)**2
       END DO
       SIGMA = DSQRT( 1./DFLOAT(NELEM/2-1) * SSQ )
       !Iteration
       DO J = 0,NELEM/2-2
        DX(NELEM/2-J) = MEAN*DCOS(PHI(NELEM/2-J))
       END DO
       DX(1) = C - SUM(DX(2:NELEM/2))
      END DO
      XNODE(NELEM/2+1,2) = 0.01*MEAN/2.
      XNODE(NELEM/2+2,2) = -0.01*MEAN/2.
END SUBROUTINE

