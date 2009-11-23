!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Collocation points matrix
subroutine collocation(N, Xnode, Cpoint)
      integer, intent(IN) :: N
      real(kind=8), dimension(N,2), intent(IN) :: Xnode
      real(kind=8), dimension(N,2), intent(OUT) :: Cpoint
      integer :: k, kp1, NOE
      do k = 1,N
       kp1 = NOE(N, k, 1)
       Cpoint(k,:) = (Xnode(k,:) + Xnode(kp1,:)) / 2.
      end do
end subroutine

!Computes the distance from vector X to Y
double precision function dist(X, Y)
      IMPLICIT NONE
      real(kind=8), dimension(2) :: X, Y
      dist = dsqrt( (X(1)-Y(1))**2 + (X(2)-Y(2))**2 )
      return
end function

!Calculate the timestep for the wake
double precision function CalcDT(Nelem, Xnode, UHoriz)
      IMPLICIT NONE
      INTEGER :: Nelem
      real(kind=8) :: UHoriz
      real(KIND=8), dimension(Nelem,2) :: XNODE
      real(kind=8) :: dist
      CalcDT = dist(Xnode(1,:),Xnode(2,:)) / UHoriz
      return
end function

!Node of element
!  N    # of elements
!  el   index of element
!  idx  index of node respect to the element
!e.g
! I need the first node of element k
!  NOE(N, k, 0)
integer function NOE(N, el, idx)
      IMPLICIT NONE
      integer :: N, el, idx
      NOE = el + idx
      if ((el + idx) .lt. 1) then
        NOE = N + el + idx
        return
      else if ((el + idx) .gt. N) then
        NOE = el + idx - N
        return
      end if
      return
end function


!Surface normals
subroutine normals(Nelem, Xnode, n)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(Nelem,2), intent(OUT) :: n
      real(kind=8), dimension(2) :: t
      real(kind=8) :: dist
      integer :: k, kp1, NOE
      do k = 1,Nelem
       kp1 = NOE(Nelem, k, 1)
       t = (Xnode(k,:) - Xnode(kp1,:))/dist(Xnode(k,:), Xnode(kp1,:))
       n(k,:) = (/-t(2), t(1)/)
      end do
end subroutine


!Rotation of the body with angle alpha (being positive counterclockwise) BFR
subroutine bodyrotation(uscalar, alpha, u)
       IMPLICIT NONE
       real(kind=8), intent(IN) :: uscalar , alpha
       real(kind=8), dimension(2), intent(OUT) :: u
       REAL(KIND=8), PARAMETER :: PI = 4.D0*datan(1.D0)
       u(1) = -uscalar*cos(alpha*pi/dble(180))
       u(2) = -uscalar*sin(alpha*pi/dble(180))
end subroutine


subroutine DeltaS(Nelem, Xnode, ds)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(Nelem), intent(OUT) :: ds
      real(kind=8) :: dist
      integer :: k, kp1, NOE
      do k = 1,Nelem
       kp1 = NOE(Nelem, k, 1)
       ds(k) = dist(Xnode(k,:), Xnode(kp1,:))
      end do
end subroutine

