!Copyright (c) 2009 Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


subroutine WakeGrid(Nelem, Xnode, Uscalar, DT, NWake, XWnode)
      IMPLICIT NONE
      real(kind=8), intent(IN) :: Uscalar
      real(kind=8), intent(IN) :: DT
      integer, intent(IN) :: Nelem, NWake
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(NWake,2), intent(OUT) :: XWnode
      real(kind=8), dimension(Nelem,2) :: Cpoint
      real(kind=8), dimension(2) :: Xhalf, Cversor
      real(kind=8) :: DXW, Dist
      integer :: i
      !Determine direction of wake by cord direction
      call collocation(Nelem, Xnode, Cpoint)
      Xhalf = Cpoint(Nelem/2+1,:)
      Cversor = (Xnode(1,:) - Xhalf)/dist(Xnode(1,:), Xhalf)
      DXW = Uscalar * DT
      do i = 1,NWake
       XWnode(i,:) = Xnode(1,:) + (i)*DXW*Cversor
      end do
end subroutine


SUBROUTINE Wake(Nelem, NWake, NTime, ITime, PhiTime, DPhiW)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NWake, NTime, ITime
      real(kind=8), dimension(Nelem,NTime), intent(IN) :: PHITime
      real(kind=8), dimension(NWake,NTime) :: DPhiW
      DPhiW(2:,ITime+1) = DPhiW(:NWake-1,ITime)
      DPhiW(1,ITime+1) = PhiTime(1,ITime) - PhiTime(Nelem,ITime)
END SUBROUTINE


!Super cool matrix D*R*S. D is aka F
!R is the retard (delays) matrix
!S is the stupid matrix
!DRS is a 'temple' matrix (name by Robin & Martin)
subroutine MatDRS(Nelem, NWake, D, DT, s, DRS)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NWake
      real(kind=8), dimension(Nelem,NWake), intent(IN) :: D
      real(kind=8), intent(IN) :: DT
      complex, intent(IN) :: s
      complex, dimension(NWake) :: Rvec
      complex, dimension(Nelem,Nelem), intent(OUT) :: DRS
      integer :: i
      DRS(:,:) = 0.
      do i = 1,NWake
       Rvec(i) = cdexp(-s*DT*i)
      end do
      do i = 1,Nelem
       DRS(i,1) = dot_product(D(i,:),Rvec)
       DRS(i,Nelem) = -DRS(i,1)
      end do
end subroutine
