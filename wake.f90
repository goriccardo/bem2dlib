!Copyright (c) 2009 Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


subroutine WAKEGRID(Nelem, Xnode, Uscalar, DT, NTime, XWnode)
      IMPLICIT NONE
      real(kind=8), intent(IN) :: Uscalar
      real(kind=8), intent(IN) :: DT
      integer, intent(IN) :: Nelem, NTime
      real(kind=8), dimension(NELEM,2), intent(IN) :: Xnode
      real(kind=8), dimension(NTIME,2), intent(OUT) :: XWnode
      real(kind=8), dimension(NELEM,2) :: Cpoint
      real(kind=8), dimension(2) :: Xhalf, Cversor
      real(kind=8) :: DXW, DIST
      integer :: i
      !Determine direction of wake by cord direction
      call collocation(Nelem, Xnode, Cpoint)
      Xhalf = Cpoint(Nelem/2+1,:)
      Cversor = (Xnode(1,:) - Xhalf)/dist(Xnode(1,:), Xhalf)
      DXW = Uscalar / DT
      do i = 1,NTime
       XWnode(i,:) = Xnode(1,:) + (i)*DXW*Cversor
      end do
end subroutine


SUBROUTINE WAKE(NELEM, Ntime, ITime, PHITime, DPHIW)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, Ntime, ITime      
      real(kind=8), dimension(Nelem,NTime), intent(IN) :: PHITime
      real(kind=8), dimension(NTime) :: DPhiW
      DPhiW(2:ITime) = DPhiW(:ITime-1)
      DPhiW(1) = PhiTime(1,ITime) - PhiTime(Nelem,ITime)
!      write(*,*) DPhiW
END SUBROUTINE
