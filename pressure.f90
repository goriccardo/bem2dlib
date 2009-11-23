!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!                   Jelle Reichert <jellereichert@gmail.com>
!Released under BSD license, see LICENSE


!Calculate the pressure in on the surface
!  Nelem    # of elements
!  Xnode    Coordinates of xnode
!  TEat1    Node 1 is a trailing edge
!  NTstep   # of timesteps
!  phi      Potential
!  Uinf     Speed
!  pinfr    Pinf/rho
!  presr     Pressure/rho (out)
subroutine calcpres(Nelem, Xnode, TEat1, NTstep, dt, U, phi, chi, presr)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTstep
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem,NTstep), intent(IN) :: phi, chi
      real(kind=8), dimension(NTstep,2), intent(IN) :: U
      real(kind=8), dimension(Nelem,NTstep), intent(OUT) :: presr
      real(kind=8), intent(IN) :: dt
      real(kind=8), dimension(Nelem) :: dphidt
      real(kind=8), dimension(Nelem,2) :: srfvel
      real(kind=8) :: v2, timeder
      integer :: k, n
      !For each timestep
      dphidt(:) = 0.
      presr(:,1) = 0.
      do n = 1,NTstep
       !Foreach element
       if (n .GT. 1) then
        dphidt = (phi(:,n) - phi(:,n-1)) / dt
       end if
       if (NTstep .GT. 1) then
        CALL CalcSrfVel(Nelem, Xnode, TEat1, phi(:,n), chi(:,n), srfvel)
        do k = 1,Nelem
         v2 = dot_product(srfvel(k,:), srfvel(k,:))
         timeder = dphidt(k)
         presr(k,n) = dot_product(U(n,:),srfvel(k,:)) - v2/dble(2) - timeder
        end do
       end if
      end do
end subroutine


subroutine calccp(Nelem, Xnode, TEat1, NTstep, dt, U, phi, chi, cp)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTstep
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem,NTstep), intent(IN) :: phi, chi
      real(kind=8), dimension(NTstep,2), intent(IN) :: U
      real(kind=8), dimension(Nelem,NTstep) :: presr
      real(kind=8), dimension(Nelem,NTstep), intent(OUT) :: cp
      real(kind=8), intent(IN) :: dt
      real(kind=8) :: ru2
      real(kind=8), parameter :: EPS = 1.D-14
      integer :: n
      CALL calcpres(Nelem, Xnode, TEat1, NTstep, dt, U, phi, chi, presr)
      do n = 1,NTstep
       ru2 = 0.5D0*dot_product(U(n,:),U(n,:))
       if (ru2 .LT. EPS) then
        cp(:,n) = 0.D0
       else
        cp(:,n) = presr(:,n)/ru2
       end if
      end do
end subroutine


!Lift on the upper side of a symmetric profile
subroutine calcHalfLift(Nelem, Xnode, TEat1, NTstep, dt, U, phi, chi, lift)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTstep
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem,NTstep), intent(IN) :: phi, chi
      real(kind=8), dimension(NTstep,2), intent(IN) :: U
      real(kind=8), dimension(Nelem,NTstep) :: presr
      real(kind=8), dimension(NTstep), intent(OUT) :: lift
      real(kind=8), intent(IN) :: dt
      real(kind=8), dimension(Nelem) :: ds
      real(kind=8), dimension(2):: F
      real(kind=8), dimension(Nelem,2) :: nrm
      integer :: k, n
      CALL calcpres(Nelem, Xnode, TEat1, NTstep, dt, U, phi, chi, presr)
      CALL normals(Nelem, Xnode, nrm)
      CALL deltas(Nelem, Xnode, ds)
      do n = 1,NTstep
       F = 0.
       do k = 1,Nelem/2
        F = F - nrm(k,:)*presr(k,n)*ds(k)
       end do
       lift(n) = F(2)
      end do
end subroutine


!Lift on the upper side of a symmetric profile
subroutine calcHalfCl(Nelem, Xnode, TEat1, NTstep, dt, U, phi, chi, cl)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTstep
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem,NTstep), intent(IN) :: phi, chi
      real(kind=8), dimension(NTstep,2), intent(IN) :: U
      real(kind=8), dimension(Nelem,NTstep) :: cp
      real(kind=8), dimension(NTstep), intent(OUT) :: cl
      real(kind=8), intent(IN) :: dt
      real(kind=8), dimension(Nelem) :: ds
      real(kind=8), dimension(2):: F
      real(kind=8), dimension(Nelem,2) :: nrm
      integer :: k, n
      CALL calccp(Nelem, Xnode, TEat1, NTstep, dt, U, phi, chi, cp)
      CALL normals(Nelem, Xnode, nrm)
      CALL deltas(Nelem, Xnode, ds)
      do n = 1,NTstep
       F = 0.
       do k = 1,Nelem/2
        F = F - nrm(k,:)*cp(k,n)*ds(k)
       end do
       cl(n) = F(2)
      end do
end subroutine



!Lift on the upper side of a symmetric profile
subroutine CalcCl(Nelem, Xnode, TEat1, NTime, dt, U, phi, chi, cl)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTime
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      logical, intent(IN) :: TEat1
      real(kind=8), dimension(Nelem,NTime), intent(IN) :: phi, chi
      real(kind=8), dimension(NTime,2), intent(IN) :: U
      real(kind=8), dimension(Nelem,NTime) :: cp
      real(kind=8), dimension(NTime), intent(OUT) :: cl
      real(kind=8), intent(IN) :: dt
      real(kind=8), dimension(Nelem) :: ds
      real(kind=8), dimension(2):: F
      real(kind=8), dimension(Nelem,2) :: nrm
      integer :: k, n
      CALL calccp(Nelem, Xnode, TEat1, NTime, dt, U, phi, chi, cp)
      CALL normals(Nelem, Xnode, nrm)
      CALL deltas(Nelem, Xnode, ds)
      do n = 1,NTime
       F = 0.
       do k = 1,Nelem
        F = F - nrm(k,:)*cp(k,n)*ds(k)
       end do
       cl(n) = F(2)
      end do
end subroutine


subroutine CalcDPhiBody(Nelem, NTime, PhiT, DPhiBody)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NTime
      real(kind=8), dimension(Nelem,NTime), intent(IN) :: PhiT
      real(kind=8), dimension(Nelem/2,NTime), intent(OUT) :: DPhiBody
      integer :: i, j
      do i = 1,NTime
       do j = 1,Nelem/2
        DPhiBody(Nelem-j+1,i) = PhiT(j,i) - PhiT(Nelem-j,i)
       end do
      end do
end subroutine

