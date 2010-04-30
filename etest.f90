!Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Test program
PROGRAM etest
      IMPLICIT NONE
      integer, parameter :: Nup = 100
      integer, parameter :: Nelem = Nup*2+1
      integer, parameter :: NWake = Nup*10
      real(kind=8), dimension(Nelem, 2) :: Xnode
      real(kind=8), parameter :: rt = 0.12D0, chord = 1.
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      real(kind=8) :: thick = chord*rt
      complex(kind=8) :: p = 2.D0*PI*dcmplx(0,0.05D0)
      complex(kind=8), dimension(2) :: f, q
      complex(kind=8), dimension(2,2) :: E
      integer :: i
      q(1) = -dcmplx(0.1D0)/abs(p)
      q(2) = dcmplx(0)
      E(:,:) = 0
      call GeomNACA00xx(Nup, Chord, Thick, Xnode)
      call EMatrixWing(Xnode, Nelem, NWake, p, E)
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
