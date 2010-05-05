!Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Test program
PROGRAM etest
      IMPLICIT NONE
      integer, parameter :: Nup = 100
      integer, parameter :: Nelem = Nup*2+1
      integer, parameter :: NWake = Nup*10
      integer, parameter :: Nlen = 100
      real(kind=8), dimension(Nelem, 2) :: Xnode
      real(kind=8), parameter :: rt = 6.D-2, chord = 1.D0, alpha = 0.D0, L = 1.0D0
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      real(kind=8) :: thick = chord*rt
      complex(kind=8) :: p = 0*2.D0*PI*dcmplx(0,0.05D0)
      complex(kind=8), dimension(2) :: q
      complex(kind=8), dimension(2) :: f
      complex(kind=8), dimension(2,2) :: E
      integer :: i
      q(1) = 0!-dcmplx(0.1D0)/abs(p)
      q(2) = dcmplx(4.D0)
      E(:,:) = 0
      call GeomNACA00xx(Nup, Chord, Thick, Xnode)
      !call EMatrixWing2d(Nelem, Xnode, alpha, NWake, p, E)
      call EMatrixWing3DA(Nelem, Xnode, alpha, Nlen, L, Nwake, p, E)
      write(*,*) "BEM:"
      do i = 1,2
       write(*,1001) E(i,:)
      end do
      write(*,*)
      f = matmul(E,q)
      do i = 1,2
       write(*,1001) f(i)
      end do
1001 FORMAT('',100(F15.8))
END PROGRAM

