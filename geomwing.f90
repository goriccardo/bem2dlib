!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the position of the nodes
!  NELEM    Figure of elements       (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  R        Radius                   (IN)
SUBROUTINE GEOMWING(NELEM, XNODE, C)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL, INTENT(IN) :: C
      REAL, DIMENSION(NELEM,3), INTENT(OUT) :: XNODE
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I
      REAL :: DX !, THETA
      DX = 1./REAL(NELEM)
      XNODE(1,1) = 0.
      XNODE(1,2) = 0.
      XNODE(NELEM/2+1,1) = C
      XNODE(NELEM/2+1,2) = 0.
      DO I = 1,NELEM/2
       XNODE(I+1,1) = DX*I
       XNODE(I+1,2) = SQRT(DX*I)*(1-DX*I)
       XNODE(NELEM-I+1,1) = DX*I
       XNODE(NELEM-I+1,2) = -XNODE(I+1,2)
       !XNODE(I,3) = ATAN2(XNODE(I,2),XNODE(I,1))+PI/2.
      END DO
END SUBROUTINE


!  R        Circle radius
!  NHX      Half of points on X figure
!  DELTAX   Size on X
!  NHY      Half of points on Y figure
!  DELTAY   Size on Y
SUBROUTINE CIRCLEFLD(R, NHX, DELTAX, NHY, DELTAY, XFIELD)
      IMPLICIT NONE
      REAL, INTENT(IN) :: R, DELTAX, DELTAY
      INTEGER, INTENT(IN) :: NHX, NHY
      REAL, DIMENSION(NHX*NHY*4,2), INTENT(OUT) :: XFIELD
      REAL, PARAMETER :: PI = 4.*ATAN(1.)
      REAL :: X, Y, DX, DY
      INTEGER :: IY, IX, NX, NY
      NX = NHX*2
      NY = NHY*2
      X = -DELTAX/2.
      DX = DELTAX/REAL(NX-1)
      DO IY = 1, NHY
       X = -DELTAX/2.
       DO IX = 1, NHX
        IF (ABS(X) .GE. R) THEN
          DY = (DELTAY)/REAL(NY-1)
        ELSE
          DY = (DELTAY - 2.*SIN(ACOS(ABS(X)/R)))/REAL(NY-1)
        END IF
        Y = -DELTAY/2. + (IY-1)*DY
!       Bottom left
        XFIELD((IY-1)*NX + IX,:) = (/X, Y/)
!       Bottom right
        XFIELD((IY)*NX - IX + 1,:) = (/-X, Y/)
!       Top left
        XFIELD(NHX*NHY*4 - (IY)*NX + IX,:) = (/X, -Y/)
!       Top right
        XFIELD(NHX*NHY*4 - (IY-1)*NX - IX + 1,:) = (/-X, -Y/)
        X = X + DX
       END DO
      END DO
END SUBROUTINE
