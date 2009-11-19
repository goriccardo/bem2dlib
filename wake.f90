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

subroutine WAKE(NELEM, TIMESTEP, PHIT, DPHIW)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, timestep
      !real(kind=8), dimension(timestep,2), intent(IN) :: XWNODE
      real(kind=8), dimension(Nelem,timestep), intent(IN) :: PHIT
      !REAL(KIND=8), DIMENSION(NELEM,NELEM) :: B, C
      !real(kind=8), dimension(NTstep,2), intent(IN) :: USCALAR
      real(kind=8), dimension(TIMESTEP,TIMESTEP), intent(OUT) :: DPHIW
      integer :: I, J
      CALL BLABLA(NELEM, TIMESTEP, PHIT)
      !CALL WAKEGRID(TIMESTEP, USCALAR, XNODE, NELEM, XWNODE)
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
end subroutine
