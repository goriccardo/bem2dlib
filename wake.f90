!Copyright (c) 2009 Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE

SUBROUTINE blabla(NELEM, TIMESTEP, PHIT) !TEMPORARY RANDOM PHI MATRIX DEPENDENT ON TIME TO SEE IF PROGRAM WORKS
      INTEGER, intent(IN) :: NELEM, TIMESTEP
      real(kind=8), dimension(NELEM,TIMESTEP), intent(OUT) :: PHIT
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
      !real(kind=8), dimension(NELEM,NELEM) :: B, C
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


subroutine WAKEGRID(Nelem, TimeStep, Xnode, Uscalar, XWnode)
      IMPLICIT NONE
      real(kind=8), intent(IN) :: Uscalar
      integer, intent(IN) :: TIMESTEP, Nelem
      real(kind=8), dimension(NELEM,2), intent(IN) :: Xnode
      real(kind=8), dimension(TIMESTEP,2), intent(OUT) :: XWnode
      real(kind=8), dimension(NELEM,2) :: Cpoint
      real(kind=8), dimension(2) :: Xhalf, Cversor
      real(kind=8) :: DXW, DIST
      integer :: i
      !Determine direction of wake by cord direction
      call collocation(Nelem, Xnode, Cpoint)
      Xhalf = Cpoint(Nelem/2+1,:)
      Cversor = (Xnode(1,:) - Xhalf)/ dist(Xnode(1,:), Xhalf)
      DXW = Uscalar / DBLE(TimeStep)
      XWnode(1,:) = Xnode(1,:)
      do i = 1,TimeStep
       if (i .gt. 1) then
         XWnode(i,:) = Xnode(1,:) + (i-1)*DXW*Cversor
       end if
      end do
end subroutine
