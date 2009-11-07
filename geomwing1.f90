!Copyright (c) 2009 Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the position of the nodes
!  NELEM    Figure of elements       (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  R        Radius                   (IN)
SUBROUTINE GEOMWING1(NELEM, XNODE, C, T)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL(KIND=8), INTENT(IN) :: C,T
      REAL(KIND=8), DIMENSION(NELEM,3), INTENT(OUT) :: XNODE
      REAL(KIND=8), DIMENSION(NELEM/2) :: LENGTH
      REAL(KIND=8), DIMENSION(NELEM) :: DIF
      REAL(KIND=8), DIMENSION(NELEM/2) :: DX,PHI
      REAL(KIND=8), PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: J, I, NITER = 10
      REAL(KIND=8) :: DS, THETA, MEAN, SIGMA, X
      THETA = (PI)/((NELEM/2.)+1.)

      XNODE(1,1) = 0.
      XNODE(1,2) = 0.
      XNODE(NELEM/2+1,1) = C
      XNODE(NELEM/2+1,2) = 0.

      DX(:) = 1.*C/REAL(NELEM/2)
      DO I = 1,NITER
       DO J = 1,NELEM/2
        X = SUM(DX(:J))
        XNODE(J+1,1) = X
        XNODE(J+1,2) = -(T/(2.*SQRT(4./27.)))*SQRT(X/C)*(1-X/C)
        PHI(J) = ATAN2(XNODE(J+1,2)-XNODE(J,2),XNODE(J+1,1)-XNODE(J,1))
        XNODE(NELEM-J+1,1) = XNODE(J+1,1)
        XNODE(NELEM-J+1,2) = -XNODE(J+1,2)
        DS = SQRT((XNODE(J+1,1)-XNODE(J,1))**2+(XNODE(J+1,2)-XNODE(J,2))**2)
        LENGTH(J) = DS
       END DO

       MEAN = SUM(LENGTH)/REAL(NELEM/2,8)
       DO J = 1,NELEM/2
        DIF(J) = (LENGTH(J)-MEAN)**2
       END DO
       SIGMA = DSQRT( 1./REAL(NELEM/2-1,8) * SUM(DIF) )

       DO J = 1,NELEM/2-1
        DX(J) = MEAN*COS(PHI(J))
       END DO
       DX(NELEM/2) = C - SUM(DX(:NELEM/2-1))

      END DO
END SUBROUTINE

