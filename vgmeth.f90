!For a naca profile
subroutine flutterstriptheory(Nup, C, T, L, alpha, M, Kmat, Vf)
      implicit none
      integer, intent(IN) :: Nup
      integer :: Nelem, NWake, Nlen
      real(kind=8), dimension(2,2), intent(IN) :: M, Kmat
      complex(kind=8), dimension(2,2) :: A, B, E
      complex(kind=8), dimension(2) :: l1, l2
      complex(kind=8) :: vl, vr
      complex(kind=8), dimension(4) :: work
      real(kind=8), dimension(16) :: rwork
      real(kind=8), intent(IN) :: C, T, L, alpha
      real(kind=8), intent(OUT) :: Vf
      real(kind=8), allocatable :: Xnode(:,:)
      real(kind=8), parameter :: kmin = 0.1D0, kmax = 2.D0
      real(kind=8) :: k, dk
      integer, parameter :: Nk = 1000
      integer :: i, INFO
      Nelem = Nup*2+1
      NWake = Nup*10
      Nlen = int(L*dble(Nup))

      allocate(Xnode(Nelem, 2))

      call GeomNACA00xx(Nup, C, T, Xnode)
      dk = (kmax - kmin)/dble(Nk)
      do i = 1,Nk
       k = kmin + i*dk
       call EMatrixWing3DA(Nelem, Xnode, alpha, Nlen, L, Nwake, dcmplx(0.,k), E)
       A = M + E/(a*k**2)
       B = -K
       call ZGGEV('N', 'N', A, 2, B, 2, l1, l2, VL, 1, VR, 1, WORK, 4, RWORK, INFO)

      end do


end subroutine

