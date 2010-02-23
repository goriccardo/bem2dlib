!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Solve the (NxN) linear system (0.5*I - C)*phi = C*chi
subroutine SolvePhi(N, B, C, PHI, CHI)
      IMPLICIT NONE
      integer, intent(IN) :: N
      real(kind=8), dimension(N,N), intent(IN) :: B, C
      real(kind=8), dimension(N), intent(OUT) :: PHI
      real(kind=8), dimension(N), intent(IN) :: CHI
      integer :: I, INFO
      real(kind=8), dimension(N,N) :: A
      real(kind=8), dimension(N) :: RHS, IPIV
      A = -C
      do I = 1, N
       A(I,I) = A(I,I) + 0.5
      end do
      RHS = matmul(B,CHI)
      call DGESV(N, 1, A, N, IPIV, RHS, N, INFO)
      if (INFO .NE. 0) then
        write(*,*) "ERROR IN LINEAR SYSTEM"
      end if
      PHI = RHS
end subroutine


subroutine SolvePhiTime(N, B, C, NWake, D, NTime, ChiTime, PhiTime, DPhiW)
      IMPLICIT NONE
      integer, intent(IN) :: N, NWake, NTime
      real(kind=8), dimension(N,N), intent(IN) :: B, C
      real(kind=8), dimension(N,NWake), intent(IN) :: D
      real(kind=8), dimension(N,NTime), intent(IN) :: ChiTime
      real(kind=8), dimension(NWake,NTime), intent(OUT) :: DPhiW
      real(kind=8), dimension(N,NTime), intent(OUT) :: PhiTime
      integer :: I, INFO
      real(kind=8), dimension(N,N) :: A
      real(kind=8), dimension(N) :: RHS, IPIV
      DPhiW(:,:) = 0.
      A = -C
      do i = 1, N
       A(i,i) = A(i,i) + DBLE(0.5)
      end do
      call DGETRF(N, N, A, N, IPIV, INFO)
      do i = 1,NTime
       RHS = matmul(B,ChiTime(:,i)) + matmul(D, DPhiW(:,i))
       call DGETRS('N', N, 1, A, N, IPIV, RHS, N, INFO)
       if (INFO .NE. 0) then
         write(*,*) "ERROR IN LINEAR SYSTEM"
       end if
       PhiTime(:,i) = RHS
       call Wake(N, NWake, NTime, i, PhiTime, DPhiW)
      end do
end subroutine


subroutine SolvePhiLap(N, B, C, Nwake, D, Nfreq, s, DT, ChiLap, PhiLap)
      IMPLICIT NONE
      integer, intent(IN) :: N, NWake, Nfreq
      real(kind=8), dimension(N,N), intent(IN) :: B, C
      real(kind=8), dimension(N,NWake), intent(IN) :: D
      complex(kind=8), dimension(N,N) :: DRS
      complex(kind=8), dimension(N,Nfreq), intent(IN) :: ChiLap
      complex(kind=8), dimension(N,Nfreq), intent(OUT) :: PhiLap
      complex(kind=8), dimension(Nfreq), intent(IN) :: s
      complex(kind=8) :: ss
      real(kind=8), intent(IN) :: DT
      integer :: i, j, INFO
      complex(kind=8), dimension(N,N) :: A
      complex(kind=8), dimension(N) :: RHS, IPIV
      real(KIND=8), parameter :: PI = 4.D0*datan(1.D0)
      do i = 1, Nfreq
       ss = 2.D0*PI*s(i)
       call MatDRS(N, NWake, D, ss, DT, DRS)
       A = - dcmplx(C) - DRS
       do j = 1, N
        A(j,j) = A(j,j) + dcmplx(0.5,0.)
       end do
       RHS = matmul(dcmplx(B),ChiLap(:,i))
       call ZGESV(N, 1, A, N, IPIV, RHS, N, INFO)
       if (INFO .NE. 0) then
         write(*,*) "ERROR IN LINEAR SYSTEM"
       end if
       PhiLap(:,i) = RHS
      end do
end subroutine


!Calculate the potential in the field
!  PHISRF   On surface                          (in)
!  CHISRF   On surface                          (in)
!  B, C     In the field                        (in)
!  PHifLD   In the field                        (out)
subroutine CalcPhiFld(N, PhiSrf, ChiSrf, NX, B, C, PHiFld)
      IMPLICIT NONE
      integer, intent(IN) :: N, NX
      real(kind=8), dimension(N), intent(IN) :: PhiSrf, ChiSrf
      real(kind=8), dimension(NX,N) :: B, C
      real(kind=8), dimension(NX), intent(OUT) :: PhiFld
      PhiFld = matmul(B, ChiSrf) + matmul(C, PhiSrf)
end subroutine
