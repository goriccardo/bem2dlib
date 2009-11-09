!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE


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
      CALL calcdphids(Nelem, Xnode, TEat1, phi, dphids)
      do k = 1,Nelem
       kp1 = NOE(Nelem, k, 1)
       t = (Xnode(k,:) - Xnode(kp1,:))/dist(Xnode(k,:), Xnode(kp1,:))
       n = (/-t(2), t(1)/)
       srfvel(k,:) = chi(k)*n + dphids(k)*t
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
      real(kind=8) :: dist, h1, h2
      integer :: NOE, k, km1, kp1, kp2, stk = 1
      !If there is a trailing edge we do a special derivative for first and last
      if (TEat1) then
       h1 = dist( (Xnode(3,:) + Xnode(2,:))/2., (Xnode(2,:) + Xnode(1,:))/2. )
       dphids(1) = (phi(1) - phi(2)) / h1
       h1 = dist( (Xnode(1,:) + Xnode(Nelem,:))/2., (Xnode(Nelem,:) + Xnode(Nelem-1,:))/2. )
       dphids(Nelem) = (phi(Nelem-1) - phi(Nelem)) / h1
       stk = 2
      end if
      do k = stk,Nelem+1-stk
       km1 = NOE(Nelem, k,-1)
       kp1 = NOE(Nelem, k, 1)
       kp2 = NOE(Nelem, k, 2)
       h1 = dist( (Xnode(kp1,:) + Xnode(k,:))/2., (Xnode(k,:) + Xnode(km1,:))/2. )
       h2 = dist( (Xnode(kp1,:) + Xnode(k,:))/2., (Xnode(kp2,:) + Xnode(kp1,:))/2. )
       dphids(k) = ( phi(km1) - phi(kp1) ) / (h1 + h2)
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
