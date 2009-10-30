!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the position of the nodes
!  NELEM    Figure of elements       (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  R        Radius                   (IN)
SUBROUTINE CIRCLEGEOM(NELEM, XNODE, R)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL, INTENT(IN) :: R
      REAL, DIMENSION(NELEM,3), INTENT(OUT) :: XNODE
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I
      REAL :: DTH, THETA
      DTH = 2.*PI/REAL(NELEM)
      DO I = 1,NELEM
       THETA = DTH*REAL(I-1)
       XNODE(I,1) = R*COS(THETA)
       XNODE(I,2) = R*SIN(THETA)
       XNODE(I,3) = THETA+DTH*0.5
      END DO
END SUBROUTINE

