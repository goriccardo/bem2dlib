!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE


!Surface velocity in the Air Frame of Reference
subroutine calcsrfvel(Nelem, Xnode, TEat1, phi, chi, srfvel)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem), intent(IN) :: phi, chi
      real(kind=8), dimension(Nelem,2), intent(OUT) :: srfvel
      real(kind=8), dimension(Nelem) :: dphids
      real(kind=8), dimension(2) :: t, n
      real(kind=8) :: dist
      integer :: k, kp1, NOE
      call CalcDPhiDs(Nelem, Xnode, TEat1, phi, dphids)
      do k = 1,Nelem
       kp1 = NOE(Nelem, k, 1)
       t = (Xnode(kp1,:) - Xnode(k,:))/dist(Xnode(k,:), Xnode(kp1,:))
       n = (/t(2), -t(1)/)
       srfvel(k,:) = chi(k)*n + dphids(k)*t
      end do
end subroutine


!Surface velocity in the Air Frame of Reference
subroutine calcSrfVelLap(Nelem, Xnode, TEat1, Nfreq, phiLap, chiLap, srfvelLap)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, Nfreq
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      complex(kind=8), dimension(Nelem,Nfreq), intent(IN) :: phiLap, chiLap
      complex(kind=8), dimension(Nelem,Nfreq,2), intent(OUT) :: srfvelLap
      complex(kind=8), dimension(Nelem) :: dphids
      real(kind=8), dimension(2) :: t, n
      real(kind=8) :: dist
      integer :: i, k, kp1, NOE
      do k = 1,Nelem
       kp1 = NOE(Nelem, k, 1)
       t = (Xnode(kp1,:) - Xnode(k,:))/dist(Xnode(k,:), Xnode(kp1,:))
       n = (/t(2), -t(1)/)
       do i = 1,Nfreq
        call CalcDPhiDsLap(Nelem, Xnode, TEat1, phiLap(:,i), dphids)
        srfvelLap(k,i,:) = chiLap(k,i)*n + dphids(k)*t
       end do
      end do
end subroutine


!Finite difference derivative of the potential on the surface
subroutine calcdphids(Nelem, Xnode, TEat1, phi, dphids)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem), intent(IN) :: phi
      real(kind=8), dimension(Nelem), intent(OUT) :: dphids
      real(kind=8), dimension(Nelem) :: DS
      real(kind=8) :: dist, h
      integer :: NOE, k, km1, kp1, stk = 1
      !If there is a trailing edge we do a special derivative for first and last
      call DeltaS(Nelem, Xnode, DS)
      if (TEat1) then
       h = dist( (Xnode(3,:) + Xnode(2,:))/2., (Xnode(2,:) + Xnode(1,:))/2. )
       dphids(1) = (phi(2) - phi(1)) / h
       h = dist( (Xnode(1,:) + Xnode(Nelem,:))/2., (Xnode(Nelem,:) + Xnode(Nelem-1,:))/2. )
       dphids(Nelem) = (phi(Nelem) - phi(Nelem-1)) / h
       stk = 2
      end if
      do k = stk,Nelem+1-stk
       km1 = NOE(Nelem, k,-1)
       kp1 = NOE(Nelem, k, 1)
       h = (DS(km1)+DS(kp1))/dble(2)+DS(k)
       dphids(k) = ( phi(kp1) - phi(km1) ) / h
      end do
end subroutine


!Finite difference derivative of the potential on the surface
!in the Laplace domai
subroutine calcdphidsLap(Nelem, Xnode, TEat1, PhiLap, dphids)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      complex(kind=8), dimension(Nelem), intent(IN) :: PhiLap
      complex(kind=8), dimension(Nelem), intent(OUT) :: dphids
      real(kind=8), dimension(Nelem) :: DS
      real(kind=8) :: dist, h
      integer :: NOE, k, km1, kp1, stk = 1
      !If there is a trailing edge we do a special derivative for first and last
      call DeltaS(Nelem, Xnode, DS)
      if (TEat1) then
       h = dist( (Xnode(3,:) + Xnode(2,:))/2., (Xnode(2,:) + Xnode(1,:))/2. )
       dphids(1) = (phiLap(2) - phiLap(1)) / h
       h = dist( (Xnode(1,:) + Xnode(Nelem,:))/2., (Xnode(Nelem,:) + Xnode(Nelem-1,:))/2. )
       dphids(Nelem) = (phiLap(Nelem) - phiLap(Nelem-1)) / h
       stk = 2
      end if
      do k = stk,Nelem+1-stk
       km1 = NOE(Nelem, k,-1)
       kp1 = NOE(Nelem, k, 1)
       h = (DS(km1)+DS(kp1))/dble(2)+DS(k)
       dphids(k) = ( phiLap(kp1) - phiLap(km1) ) / h
      end do
end subroutine


!Finite difference x component of the potential gradient
subroutine calcdphidx(Nelem, Xnode, TEat1, phi, dphidx)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem), intent(IN) :: phi
      real(kind=8), dimension(Nelem), intent(OUT) :: dphidx
      real(kind=8) :: dx
      integer :: NOE, k, km1, kp1, kp2, stk = 1
      !If there is a trailing edge we do a special derivative for first and last
      if (TEat1) then
       dx = (Xnode(1,1) - Xnode(3,1))/2.
       dphidx(1) = (phi(1) - phi(2)) / dx
       dx = (Xnode(Nelem-1,1) - Xnode(1,1))/2.
       dphidx(Nelem) = (phi(Nelem-1) - phi(Nelem)) / dx
       stk = 2
      end if
      do k = stk,Nelem+1-stk
       km1 = NOE(Nelem, k,-1)
       kp1 = NOE(Nelem, k, 1)
       kp2 = NOE(Nelem, k, 2)
       dx = (Xnode(k,1) + Xnode(km1,1) - Xnode(kp1,1) - Xnode(kp2,1))/2.
       dphidx(k) = ( phi(km1) - phi(kp1) ) / dx
      end do
end subroutine
