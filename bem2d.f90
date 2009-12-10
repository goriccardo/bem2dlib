!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Impose boundary conditions for a moving circle with u velocity
!in the air frame of reference
subroutine BCondVel(N, XNODE, U, CHI)
      IMPLICIT NONE
      integer, intent(IN) :: N
      real(kind=8), dimension(N,2), intent(IN) :: XNODE
      real(kind=8), dimension(2), intent(IN) :: U
      real(kind=8), dimension(N), intent(OUT) :: CHI
      real(kind=8), dimension(2) :: T
      integer :: I, NOE
      do I = 1, N
       T = XNODE(NOE(N,I,1),:2) - XNODE(NOE(N,I,0),:2)
       T = T/SQRT(T(1)**2 + T(2)**2)
       CHI(I) = U(1)*T(2) - U(2)*T(1)
      end do
end subroutine


!We assume that U Horizontal is constant
subroutine BCondLap(Nelem, Xnode, Ampl, ChiLap)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: Ampl
      complex(kind=8), dimension(Nelem), intent(OUT) :: ChiLap
      real(kind=8), dimension(Nelem,2) :: n
      integer :: I
      call normals(Nelem, Xnode, n)
      do I = 1, Nelem
       ChiLap(i) = Ampl*n(i,2)
      end do
end subroutine


subroutine BCondRotLap(Nelem, Xnode, alphaAmpl, ChiLap)
      IMPLICIT NONE
      real(KIND=8), parameter :: PI = 4.D0*datan(1.D0)
      integer, intent(IN) :: Nelem
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: alphaAmpl
      real(kind=8) :: alpharad
      complex(kind=8), dimension(Nelem), intent(OUT) :: ChiLap
      real(kind=8), dimension(Nelem,2) :: n
      integer :: I
      call normals(Nelem, Xnode, n)
      alpharad = dble(2)*PI*alphaAmpl/dble(360)
      do I = 1, Nelem
       ChiLap(i) = dcmplx(alpharad*n(i,2),alpharad*n(i,1))
      end do
end subroutine


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
      real(kind=8), dimension(N,N) :: A, ATEMP
      real(kind=8), dimension(N) :: RHS, IPIV
      DPhiW(:,:) = 0.
      A = -C
      do i = 1, N
       A(i,i) = A(i,i) + DBLE(0.5)
      end do
      do i = 1,NTime
       ATemp = A
       RHS = matmul(B,ChiTime(:,i)) + matmul(D, DPhiW(:,i))
       call DGESV(N, 1, ATEMP, N, IPIV, RHS, N, INFO)
       if (INFO .NE. 0) then
         write(*,*) "ERROR IN LINEAR SYSTEM"
       end if
       PhiTime(:,i) = RHS
       call Wake(N, NWake, NTime, i, PhiTime, DPhiW)
      end do
end subroutine


subroutine SolvePhiLap(N, B, C, DRS, ChiLap, PhiLap)
      IMPLICIT NONE
      integer, intent(IN) :: N
      real(kind=8), dimension(N,N), intent(IN) :: B, C
      complex(kind=8), dimension(N,N), intent(IN) :: DRS
      complex(kind=8), dimension(N), intent(IN) :: ChiLap
      complex(kind=8), dimension(N), intent(OUT) :: PhiLap
      integer :: I, INFO
      complex(kind=8), dimension(N,N) :: A
      complex(kind=8), dimension(N) :: RHS, IPIV
      A = - dcmplx(C) - DRS
      do i = 1, N
       A(i,i) = A(i,i) + dcmplx(0.5,0.)
      end do
      RHS = matmul(dcmplx(B),ChiLap)
      call ZGESV(N, 1, A, N, IPIV, RHS, N, INFO)
      if (INFO .NE. 0) then
        write(*,*) "ERROR IN LINEAR SYSTEM"
      end if
      PhiLap = RHS
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
