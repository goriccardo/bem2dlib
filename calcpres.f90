!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the pressure in on the surface
!  Nelem    # of elements
!  Xnode    Coordinates of xnode
!  TEat1    Node 1 is a trailing edge
!  NTstep   # of timesteps
!  U        Vector speed
!  phi      Potential
!  chi      Normal wash
!  pres     Pressure (out)
subroutine CALCPRES(Nelem, Xnode, TEat1, NTstep, U, phi, pres)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTstep
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem,NTstep), intent(IN) :: phi
      real(kind=8), dimension(Nelem,NTstep), intent(OUT) :: pres
      real(kind=8), dimension(NTstep,2), intent(IN) :: U
      real(kind=8), dimension(Nelem,2) :: dphidx
      real(kind=8), parameter :: PI = 4.*atan(1.)
      integer :: K, N
      !For each timestep
      do n = 2,NTstep
       !Foreach element
       CALL calcdphidx(Nelem, Xnode, TEat1, phi(:,n), dphidx(:,2))
       do k = 1,Nelem/2

       end do
      end do
end subroutine
