!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Impose boundary conditions for a moving body with U velocity
!in the air frame of reference
subroutine BCondVel(Nelem, Xnode, U, Chi)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(2), intent(IN) :: U
      real(kind=8), dimension(Nelem), intent(OUT) :: Chi
      real(kind=8), dimension(Nelem,2) :: n
      integer :: i
      call normals(Nelem,XNode,n)
      do i = 1, Nelem
       Chi(i) = dot_product(U,n(i,:))
      end do
end subroutine


subroutine BCondStraight(Nelem, Xnode, Uscalar, alpha, Ntime, Ut, ChiTime)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, Ntime
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: Uscalar, alpha
      real(kind=8), dimension(Nelem,Ntime), intent(OUT) :: ChiTime
      real(kind=8), dimension(Ntime,2), intent(OUT) :: Ut
      real(kind=8), dimension(2) :: U
      integer :: i
      call bodyrotation(Uscalar, alpha, U)
      Ut(:,1) = U(1)
      Ut(:,2) = U(2)
      do i = 1,Ntime
        call BCondVel(Nelem, Xnode, U, ChiTime(:,i))
      end do
end subroutine


!Impose boundary conditions for a roto-translating body
!in the air frame of reference
! N         # of elements
! Xnode     Node vector
! Xo        Rotation center
! Uo        Body Speed Vector
! alpha     Angle of incidence
! alphaAmpl Amplitude
! Freq      Frequency [rad/step]
! ChiTime   Boundary Conditions Matrix
subroutine BCondRot(Nelem, Xnode, Xo, UScalar, alpha, alphaAmpl, Freq, DT, Ntime, Ut, ChiTime)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, Ntime
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: alpha, alphaAmpl, freq, DT, UScalar
      real(kind=8), dimension(2), intent(IN) :: Xo
      real(kind=8), dimension(Nelem,Ntime), intent(OUT) :: ChiTime
      real(kind=8), dimension(Ntime,2), intent(OUT) :: Ut
      real(kind=8), dimension(2) :: U, Ur, R, T
      real(kind=8), dimension(Nelem,2) :: n, Cpoint
      real(kind=8) :: alphat, ws, wt
      real(KIND=8), parameter :: PI = 4.D0*datan(1.D0)
      integer :: i, j
      call normals(Nelem, Xnode, n)
      call collocation(Nelem, Xnode, Cpoint)
      ws = dble(2)*PI*freq*DT ![rad/step]
      wt = dble(2)*PI*freq    ![rad/s]
      Ut(:,1) = U(1)
      Ut(:,2) = U(2)
      do i = 1,NTime
       alphat = alpha + alphaAmpl*dsin(ws*i)
       call bodyrotation(uscalar, alphat, U)
       do j = 1,Nelem
!       Constant component
        R = (Cpoint(j,:) - Xo)
        call rotateVec90(1,R,T)
        Ur = U + alphaAmpl*PI/dble(180)*T*wt*dcos(ws*i)
        ChiTime(j,i) = dot_product(Ur, n(j,:))
       end do
      end do
end subroutine


subroutine BCondRotLap(Nelem, Xnode, Xo, Uscalar, alpha, alphaAmpl, DT, NFreq, s, ChiLap)
      IMPLICIT NONE
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(2), intent(IN) :: Xo
      real(kind=8), intent(IN) :: alphaAmpl, Uscalar, DT, alpha
      real(kind=8) :: alpharad, wt, vxampl, vyampl, R, dist
      integer, intent(IN) :: NFreq
      complex(kind=8), dimension(Nfreq), intent(OUT) :: s
      complex(kind=8), dimension(Nelem,Nfreq), intent(OUT) :: ChiLap
      real(kind=8), dimension(Nelem,2) :: n, CPoint
      integer :: I
      call normals(Nelem, Xnode, n)
      call collocation(Nelem, Xnode, Cpoint)
!     Zeroth freq
      
      
      wt = PI !*freq    ![rad/s]
      do I = 1, Nelem
       R = dist(Cpoint(i,:), Xo)
       vxampl = Uscalar*dcos(PI/dble(180)*alphaAmpl)
       vyampl = Uscalar*dsin(PI/dble(180)*alphaAmpl)
       ChiLap(i,1) = dcmplx(wt*R + vxampl*n(i,1) + vyampl*n(i,2), 0.)
      end do
end subroutine


!Impose boundary conditions for vertical oscillating body
!in the air frame of reference
! N         # of elements
! Uscalar   Module of Body Velocity [m/s]
! alpha     Angle of incidence [deg°]
! VelAmpl   Vertical Velocity Amplitude [m/s]
! Freq      Frequency [Hz]
! Dt        Time step size [s]
! NTime     # time steps
! ChiTime   Boundary Conditions Matrix
subroutine BCondOscil(Nelem, Xnode, Uscalar, alpha, VelAmpl, Freq, DT, Ntime, Ut, ChiTime)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, Ntime
      real(kind=8), intent(IN) :: alpha, VelAmpl, Freq, Uscalar, DT
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(Ntime,2), intent(OUT) :: Ut
      real(kind=8), dimension(2) :: U
      real(kind=8), dimension(Nelem,Ntime), intent(OUT) :: ChiTime
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
!     w <- Angular Velocity [rad/step]
      real(kind=8) :: w, alpharad, Utmp
      integer :: i
      alpharad = alpha*PI/dble(180)
      call bodyrotation(uscalar, alpha, U)
      w = dble(2)*PI*freq*DT ![rad/step]
      do i = 1,NTime
       Utmp = -VelAmpl*dsin(w*dble(i))
       Ut(i,1) = U(1) + Utmp*dsin(alpharad)
       Ut(i,2) = U(2) - Utmp*dcos(alpharad)
       call BCondVel(Nelem, Xnode, Ut(i,:), ChiTime(:,i))
      end do
end subroutine


subroutine BCondOscilLap(Nelem, Xnode, Uscalar, alpha, VelAmpl, Freq, NFreq, s, ChiLap)
      IMPLICIT NONE
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: VelAmpl, alpha, Freq, Uscalar
      integer, intent(IN) :: Nfreq
      complex(kind=8), dimension(Nfreq), intent(OUT) :: s
      complex(kind=8), dimension(Nelem,Nfreq), intent(OUT) :: ChiLap
      real(kind=8), dimension(Nelem,2) :: n
      real(kind=8) :: alpharad
      integer :: I
      if (Nfreq .ne. 2) then
        write(*,*) "WARNING, BCondOscilLap needs Nfreq == 2"
        s(:) = dcmplx(0)
        ChiLap(:,:) = 0.
        return
      end if
      s(1) = dcmplx(0)
      s(2) = dcmplx(0,Freq)
      call normals(Nelem, Xnode, n)
      alpharad = alpha*PI/dble(180)
      do I = 1, Nelem
       ChiLap(i,1) = dcmplx(-UScalar*(n(i,1)*dcos(alpharad) + n(i,2)*dsin(alpharad)))
       ChiLap(i,2) = dcmplx(-VelAmpl*(-n(i,2)*dcos(alpharad) + n(i,1)*dsin(alpharad)))
      end do
end subroutine