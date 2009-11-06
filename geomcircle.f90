!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the position of the nodes
!  NELEM    Figure of elements       (IN)
!  XNODE    Nodes and angles vector  (OUT)
!  R        Radius                   (IN)
SUBROUTINE CIRCLEGEOM(NELEM, XNODE, R)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NELEM
      REAL(KIND=8), INTENT(IN) :: R
      REAL(KIND=8), DIMENSION(NELEM,3), INTENT(OUT) :: XNODE
      REAL(KIND=8), PARAMETER :: PI = 4.*ATAN(1.)
      INTEGER :: I
      REAL(KIND=8) :: DTH, THETA
      DTH = 2.*PI/REAL(NELEM,8)
      DO I = 1,NELEM
       THETA = DTH*REAL(I-1,8)
       XNODE(I,1) = R*COS(THETA)
       XNODE(I,2) = R*SIN(THETA)
       XNODE(I,3) = THETA+DTH*0.5
      END DO
END SUBROUTINE


!  R        Circle radius
!  NHX      Half of points on X figure
!  DELTAX   Size on X
!  NHY      Half of points on Y figure
!  DELTAY   Size on Y
SUBROUTINE CIRCLEFLD(R, NHX, DELTAX, NHY, DELTAY, XFIELD)
      IMPLICIT NONE
      REAL(KIND=8), INTENT(IN) :: R, DELTAX, DELTAY
      INTEGER, INTENT(IN) :: NHX, NHY
      REAL(KIND=8), DIMENSION(NHX*NHY*4,2), INTENT(OUT) :: XFIELD
      REAL(KIND=8), PARAMETER :: PI = 4.*ATAN(1.)
      REAL(KIND=8) :: X, Y, DX, DY
      INTEGER :: IY, IX, NX, NY
      NX = NHX*2
      NY = NHY*2
      X = -DELTAX/2.
      DX = DELTAX/REAL(NX-1,8)
      DO IY = 1, NHY
       X = -DELTAX/2.
       DO IX = 1, NHX
        IF (ABS(X) .GE. R) THEN
          DY = (DELTAY)/REAL(NY-1,8)
        ELSE
          DY = (DELTAY - 2.*SIN(ACOS(ABS(X)/R)))/REAL(NY-1,8)
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
