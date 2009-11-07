!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE

!Calculate the pressure in on the surface
!  Nelem    # of elements
!  Xnode    Coordinates of xnode
!  TEnode   Node index of the trailing edge
!  NTstep   # of timesteps
!  U        Vector speed
!  phi      Potential
!  chi      Normal wash
!  pres     Pressure (out)
subroutine CALCPRES(Nelem, Xnode, TEnode, NTstep, U, phi, chi, pres)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTstep
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: TEnode
      real(kind=8), dimension(Nelem,NTstep), intent(IN) :: phi, chi
      real(kind=8), dimension(Nelem,NTstep), intent(OUT) :: pres
      real(kind=8), dimension(NTstep,2), intent(IN) :: U
      real(kind=8), parameter :: PI = 4.*atan(1.)
      real(kind=8) :: UINF, DX, DT, VB, TIME
      integer :: K, N
      !For each timestep
      do N=2, NTstep
       UINF = SQRT(U(N,1)**2+U(N,2)**2)
       DX = SQRT((Xnode(1,1)-Xnode(2,1))**2+(Xnode(1,2)-Xnode(2,2))**2)
       pres(1,N) = 0.5*(((phi(1,N)-phi(1,N-1))/DT)+((phi(Nelem-1,N)-phi(Nelem-1,N-1))/DT))+UINF*(phi(1,N)-phi(Nelem,N))/DX 
       pres(Nelem,N) = pres(1,N)
       !Foreach element
       do K=2, Nelem/2
        DX = SQRT((Xnode(K-1,1)-Xnode(K,1))**2+(Xnode(K-1,2)-Xnode(K,2))**2)
        pres(K,N) = 0.5*(((phi(K,N)-phi(K,N-1))/DT)+((phi(K-1,N)-phi(K-1,N-1))/DT))+UINF*(phi(1,N)-phi(Nelem,N))/DX
        pres(Nelem-K+1,N) = pres(K,N)
       end do
      end do
end subroutine

