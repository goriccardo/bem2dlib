!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the position of the nodes
!The formula is y = sqrt(x)*(1-x)
!  NELEM    Figure of elements       (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  T        Thickness                (IN)
!  C        Chord                    (IN)
SUBROUTINE GEOMWING(NELEM, XNODE, C, T)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL(KIND=8), INTENT(IN) :: T, C
      REAL(KIND=8), DIMENSION(NELEM,3), INTENT(OUT) :: XNODE
      REAL(KIND=8), PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I
      REAL(KIND=8) :: DX !, THETA
      DX = 2.*C/REAL(NELEM,8)
      !First and last nodes
      XNODE(1,1) = C
      XNODE(1,2) = 0.
      XNODE(NELEM/2+1,1) = C
      XNODE(NELEM/2+1,2) = 0.
      !Symmetric profile
      DO I = 1,NELEM/2
       XNODE(I+1,1) = C-DX*I
       XNODE(I+1,2) = SQRT(XNODE(I+1,1))*(1-XNODE(I+1,1))
       XNODE(NELEM-I+1,1) = XNODE(I+1,1)
       XNODE(NELEM-I+1,2) = -XNODE(I+1,2)
       !XNODE(I,3) = ATAN2(XNODE(I,2),XNODE(I,1))+PI/2.
      END DO
END SUBROUTINE

