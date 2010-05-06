!For a naca profile
! Nup:      Nuber of elements of a 2D section
! Tr:       Thickness ratio
! Lr:       Lenght rato
! alpha:    Angle of incidence
! mu:       2*m/(rho*c**3)
! M:        Mass Matrix
! Kmat:     Stiffness matrix
! Vf:       Flutter speed
subroutine flutterstriptheory(Nup, Tr, Lr, alpha, mu, M, Kmat, Vf)
      implicit none
      integer, intent(IN) :: Nup
      integer :: Nelem, NWake, Nlen
      real(kind=8), dimension(2,2), intent(IN) :: M, Kmat
      complex(kind=8), dimension(2,2) :: A, A1, B, E
      complex(kind=8), dimension(2) :: l1, l2
      real(kind=8), dimension(2,2) :: ss
      complex(kind=8) :: vl, vr
      complex(kind=8), dimension(4) :: work, rrwork
      real(kind=8), dimension(16) :: rwork
      real(kind=8), intent(IN) :: Tr, Lr, alpha, mu
      real(kind=8), intent(OUT) :: Vf
      real(kind=8), allocatable :: Xnode(:,:)
      real(kind=8), parameter :: kmin = -1.D-2, kmax = 1.D-1
      real(kind=8) :: k, dk
      integer, parameter :: Nk = 100
      integer, dimension(2) :: IPIV
      integer :: i, j, INFO
      Nelem = Nup*2+1
      NWake = Nup*10
      Nlen = int(Lr*dble(Nup))

      allocate(Xnode(Nelem, 2))

      call GeomNACA00xx(Nup, 1.D0, Tr, Xnode)
      dk = (kmax - kmin)/dble(Nk)
      ss(:,:) = 0.D0
      do i = 1,Nk
       k = kmin + i*dk
       call EMatrixWing3DA(Nelem, Xnode, alpha, Nlen, Lr, Nwake, dcmplx(0.,k), E)
       A = - M - E/(mu*k**2)
       B = Kmat
       !call ZGESV( 2, 2, B, 2, IPIV, A, 2, INFO )
       !call ZGEEV( 'N', 'N', 2, A, 2, w, VL, 1, VR, 1, WORK, 4, RRWORK, INFO )
       !A = - M + E/(mu*k**2)
       !B = Kmat
       call ZGGEV('N', 'N', 2, A, 2, B, 2, l1, l2, VL, 1, VR, 1, WORK, 4, RWORK, INFO)
       l1 = l1/l2
       ss(:,2) = sign(1.D0,imagpart(l1))
       if (any((ss(:,1)+ss(:,2)) .eq. 0.D0)) then
         write(*,*) "Flutter"
         return
       end if
       write(*,*) k
       do j = 1,2
        write(*,*) "L:", l1(j), ss(j,2)
        Vf = dsqrt(1.D0/(k**2*realpart(-l1(j))))
        write(*,*) "U:", Vf
       end do
       write(*,*)
      end do


end subroutine

