!Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Return the adimensional matrix E_BC for the wing -> E(p)
! p : Complex variable (IN)
subroutine EMatrixWing(Nelem, Xnode, NWake, p, E)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NWake
      complex(kind=8), intent(IN) :: p
      complex(kind=8), dimension(2,2), intent(OUT) :: E
      complex(kind=8) :: EBC(Nelem,2), EIE(Nelem,Nelem), EBT(Nelem,Nelem), EGF(2,Nelem)
      real(kind=8), parameter :: chord = 1.D0, Uinf = 1.D0, alpha = 0.D0
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(NWake,2) :: XWnode
      real(kind=8), dimension(Nelem,2) :: n0, Cpoint
      real(kind=8), dimension(Nelem) :: ds
      real(kind=8), dimension(2) :: nk
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(Nelem, NWake) :: D
      complex(kind=8), dimension(Nelem,Nelem) :: DRS
      complex(kind=8), dimension(Nelem,Nelem) :: A
      complex(kind=8), dimension(Nelem) :: IPIV
      real(kind=8) :: DT, CalcDT
      integer :: i, INFO

      ! Geometry
      DT = CalcDT(Nelem, Xnode, Uinf)
      call WakeGrid(Nelem, Xnode, Uinf, DT, NWake, XWnode)

      call collocation(Nelem, Xnode, Cpoint)
      call normals(Nelem, Xnode, n0)

      ! Create E_BC(1)
      do i = 1, Nelem
       call rotateVec90(1, n0(i,:), nk)
       EBC(i,1) = p*n0(i,2)
       EBC(i,2) = -nk(1)+p*n0(i,2)*CPoint(i,1)
      end do

      ! Create E_IE
      call SrfMatBCD(Nelem, Xnode, NWake, XWnode, B, C, D)
      call MatDRS(Nelem, NWake, D, p, DT, DRS)
      A = - dcmplx(C) - DRS
      EIE = dcmplx(B)
      do i = 1, Nelem
       A(i,i) = A(i,i) + dcmplx(0.5D0,0.)
      end do
      call ZGESV(Nelem, Nelem, A, Nelem, IPIV, EIE, Nelem, INFO)
      if (INFO .NE. 0) then
        write(*,*) "ERROR IN LINEAR SYSTEM"
      end if

      ! Create E_BT
      call EmatrixBT(Nelem, Xnode, p, EBT)

      ! Create EGF
      EGF(:,:) = dcmplx(0)
      do i = 1, Nelem
       EGF(1,i) = n0(i,2)*ds(i)
       EGF(2,i) = n0(i,2)*ds(i)*Cpoint(i,1)
      end do

      E = matmul(EGF,matmul(matmul(EBT,EIE),EBC))
end subroutine


subroutine EmatrixBT(Nelem, Xnode, p, EBT)
      implicit none
      integer, intent(IN) :: Nelem
      complex(kind=8), intent(IN) :: p
      complex(kind=8), intent(OUT) :: EBT(Nelem,Nelem)
      real(kind=8), parameter :: Uinf = 1.D0, alpha = 0.D0
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(Nelem,2) :: Cpoint
      real(kind=8), dimension(Nelem) :: ds
      real(kind=8) :: h
      integer :: i, k, km1, kp1, NOE

      call DeltaS(Nelem, Xnode, ds)
      EBT(:,:) = dcmplx(0)
      h = Cpoint(2,1)-Cpoint(1,1)
      EBT(1,1) = -1/h
      EBT(1,2) = 1/h
      h = Cpoint(Nelem,1)-Cpoint(Nelem-1,1)
      EBT(Nelem,Nelem-1) = -1/h
      EBT(Nelem,Nelem) = 1/h
      do k = 2, Nelem-1
       km1 = NOE(Nelem, k,-1)
       kp1 = NOE(Nelem, k, 1)
       h = Cpoint(kp1,1)-Cpoint(km1,1)
       if (dabs(h) .gt. 1.D-10) then
         EBT(k,k-1) = -1/h
         EBT(k,k+1) =  1/h
       end if
      end do
      EBT = -2*EBT
      do i = 1, Nelem
       EBT(i,i) = EBT(i,i)-2.D0*p
      end do
end subroutine


!subroutine EmatrixStripTheory(Nelem, Xnode, alpha, Nwake, Nlen, L, Nmodes, F, p, E)
!      implicit none
!      integer, intent(IN) :: Nelem, Nwake, Nlen, Nmodes
!      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
!      real(kind=8), dimension(Nmodes,Nlen), intent(IN) :: F
!      real(kind=8), intent(IN) :: L, alpha
!      complex(kind=8), intent(IN) :: p
!      real(kind=8), dimension(Nmodes, Nmodes), intent(OUT) :: E
!      real(kind=8), parameter :: uscalar = 1.D0
!      real(kind=8) :: UHoriz
!      integer, dimension(2), parameter :: Nfr = (/2,3/)
!      complex(kind=8) :: s1(2), s2(3)
!      integer :: i, j
!      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
!      UHoriz = Uscalar*dcos(alpha*PI/dble(180))
!      DT = CalcDT(Nelem, Xnode, UHoriz)

!     Oscillation
!      call BCondOscilLap(Nelem, Xnode, uscalar, alpha, VelAmpl, freq, 2, s1, Us, ChiLap1)
!     Rotation
!      call BCondRotLap(Nelem, Xnode, Xo, Uscalar, alpha, alphaAmpl, freq, 3, s2, Us, ChiLap2)

!      do i = 1,Nmodes
!       do j = 1,Nmodes
!        E(i,j) = 1.
!       end do
!      end do
!end subroutine

