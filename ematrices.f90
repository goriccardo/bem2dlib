!Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Return the adimensional matrix E_BC for the wing -> E(p)
! xe: Rotation center
! p : Complex variable (IN)
subroutine EMatrixWingPres(Nelem, Xnode, alpha, NWake, xe, p, E)
      IMPLICIT NONE
      integer, intent(IN) :: Nelem, NWake
      complex(kind=8), intent(IN) :: p
      real(kind=8), intent(IN) :: xe, alpha
      complex(kind=8), dimension(Nelem,2), intent(OUT) :: E
      complex(kind=8) :: EBC(Nelem,2), EIE(Nelem,Nelem), EBT(Nelem,Nelem)
      real(kind=8), parameter :: chord = 1.D0, Uinf = 1.D0
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), dimension(NWake,2) :: XWnode
      real(kind=8), dimension(Nelem,2) :: n0, Cpoint
      real(kind=8), dimension(2) :: nk, Uvec
      real(kind=8), dimension(Nelem,Nelem) :: B, C
      real(kind=8), dimension(Nelem, NWake) :: D
      complex(kind=8), dimension(Nelem,Nelem) :: DRS
      complex(kind=8), dimension(Nelem,Nelem) :: A
      complex(kind=8), dimension(Nelem) :: IPIV
      real(kind=8) :: DT, CalcDT, h
      real(KIND=8), parameter :: PI = 4.D0*datan(1.D0)
      integer :: i, k, km1, kp1, NOE, INFO

      ! Geometry
      DT = CalcDT(Nelem, Xnode, Uinf)
      call WakeGrid(Nelem, Xnode, Uinf, DT, NWake, XWnode)

      call collocation(Nelem, Xnode, Cpoint)
      call normals(Nelem, Xnode, n0)

      Uvec = (/dcos(alpha)*PI/dble(180), -dsin(alpha)*PI/dble(180)/)

      ! Create E_BC(1)
      do i = 1, Nelem
       call rotateVec90(1, n0(i,:), nk)
       EBC(i,1) = p*n0(i,2)
       EBC(i,2) = sum(Uvec*nk) + p*n0(i,2)*(CPoint(i,1)-xe)
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

      E = matmul(matmul(EBT,EIE),EBC)
end subroutine


subroutine EMatrixWing2D(Nelem, Xnode, alpha, NWake, p, E)
      implicit none
      integer, intent(IN) :: Nelem, NWake
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: alpha
      complex(kind=8), intent(IN) :: p
      complex(kind=8), dimension(2,2), intent(OUT) :: E
      complex(kind=8), dimension(Nelem,2) :: EE
      complex(kind=8), dimension(2,Nelem) :: EGF
      real(kind=8), dimension(Nelem) :: ds
      real(kind=8), dimension(Nelem,2) :: n0, Cpoint
      integer :: i
      call EMatrixWingPres(Nelem, Xnode, alpha, NWake, 0.D0, p, EE)
      call DeltaS(Nelem, Xnode, ds)
      call normals(Nelem, Xnode, n0)
      call collocation(Nelem, Xnode, Cpoint)
      ! Create EGF
      EGF(:,:) = dcmplx(0)
      do i = 1, Nelem
       EGF(1,i) = -n0(i,2)*ds(i)
       EGF(2,i) = n0(i,2)*ds(i)*Cpoint(i,1)
      end do
      do i = 1, Nelem
       write(*,*) EE(i,:), Cpoint(i,1)
      end do
      E = matmul(EGF,EE)
end subroutine


subroutine EMatrixWing3DA(Nelem, Xnode, alpha, Nlen, L, Nwake, p, E)
      implicit none
      integer, intent(IN) :: Nelem, NWake, Nlen
      real(kind=8), dimension(Nelem,2), intent(IN) :: Xnode
      real(kind=8), intent(IN) :: alpha, L
      complex(kind=8), intent(IN) :: p
      complex(kind=8), dimension(2,2), intent(OUT) :: E
      complex(kind=8), dimension(Nelem,2) :: EE
      real(kind=8) :: dL, q1, q2
      real(kind=8), dimension(Nelem) :: ds
      real(kind=8), dimension(Nelem,2) :: n0, Cpoint
      integer :: i
      call EMatrixWingPres(Nelem, Xnode, alpha, NWake, 5.D-1, p, EE)
      call DeltaS(Nelem, Xnode, ds)
      call normals(Nelem, Xnode, n0)
      call collocation(Nelem, Xnode, Cpoint)
      dL = L/dble(Nlen)
      E(:,:) = 0.
      do i = 1, Nlen
       q1 = (i/dble(Nlen))**2
       q2 = i/dble(Nlen)
       E(1,1) = E(1,1) - sum(n0(:,2)*ds(:)*dL*EE(:,1)*q1)
       E(1,2) = E(1,2) - sum(n0(:,2)*ds(:)*dL*EE(:,2)*q2)
       E(2,1) = E(2,1) + sum(n0(:,2)*ds(:)*dL*EE(:,1)*q1*(Cpoint(:,1) - 5.D-1))
       E(2,2) = E(2,2) + sum(n0(:,2)*ds(:)*dL*EE(:,2)*q2*(Cpoint(:,1) - 5.D-1))
      end do
end subroutine

