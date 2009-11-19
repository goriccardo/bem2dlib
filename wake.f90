!Copyright (c) 2009    Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE

SUBROUTINE blabla(NELEM, TIMESTEP, PHIT) !TEMPORARY RANDOM PHI MATRIX DEPENDENT ON TIME TO SEE IF PROGRAM WORKS
      INTEGER, INTENT(IN) :: NELEM, TIMESTEP
      REAL(KIND=8), DIMENSION(NELEM,TIMESTEP), INTENT(OUT) :: PHIT
      INTEGER :: I
      PHIT(:,:) = 0.
      DO I = 1,TIMESTEP
       PHIT(1,I) = I*1
       PHIT(NELEM,I) = I*2
      END DO
END SUBROUTINE


SUBROUTINE WAKE(NELEM, TIMESTEP, PHIT, DPHIW)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, timestep
      !real(kind=8), dimension(timestep,2), intent(IN) :: XWNODE
      real(kind=8), dimension(Nelem,timestep), intent(IN) :: PHIT
      !REAL(KIND=8), DIMENSION(NELEM,NELEM) :: B, C
      !real(kind=8), dimension(NTstep,2), intent(IN) :: USCALAR
      real(kind=8), dimension(TIMESTEP,TIMESTEP), intent(OUT) :: DPHIW
      integer :: I, J
      CALL BLABLA(NELEM, TIMESTEP, PHIT)
      !CALL WAKEGRID(NELEM, TIMESTEP, XNODE, USCALAR, XWNODE)
      !CALL BCONDVEL(NELEM, XNODE, U, CHI)
      !CALL SOLVEPHI(NELEM, B, C, PHI, CHI)
      !CALC THE DPHI FOR EACH WAKE (TIME) PART
      DPHIW(:,:) = 0.
      DO I = 2,TIMESTEP
        DPHIW(1,I) = PHIT(1,I-1) - PHIT(NELEM,I-1)
        DO J = 2,TIMESTEP
          DPHIW(J,I) = DPHIW(J-1,I-1)
        END DO
      END DO
END SUBROUTINE


SUBROUTINE WAKEGRID(NELEM, TIMESTEP, XNODE, USCALAR, XWNODE)
       IMPLICIT NONE
       REAL(KIND=8), INTENT(IN) :: USCALAR
       INTEGER, INTENT(IN) :: TIMESTEP, NELEM
       REAL(KIND=8), DIMENSION(NELEM,2), INTENT(IN) :: XNODE
       REAL(KIND=8), DIMENSION(TIMESTEP,2), INTENT(OUT) :: XWNODE
       REAL(kind=8), DIMENSION(NELEM,2) :: Cpoint
       REAL(KIND=8), DIMENSION(2) :: XHALF, CVERSOR
       REAL(KIND=8) :: DXW, DIST
       INTEGER :: I
       !Determine direction of wake by cord direction
       CALL COLLOCATION(NELEM, XNODE, CPOINT)
       XHALF = CPOINT(NELEM/2+1,:)
       CVERSOR = (XNODE(1,:) - XHALF)/ DIST(XNODE(1,:), XHALF)
       DXW = USCALAR / DBLE(TIMESTEP)
       XWNODE(1,:) = XNODE(1,:)
       DO I = 1,TIMESTEP
        IF (I .GT. 1) THEN
          XWNODE(I,:) = XNODE(1,:) + (I-1)*DXW*CVERSOR
        END IF
       END DO
END SUBROUTINE
