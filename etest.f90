!Copyright (c) 2010 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Test program
PROGRAM etest
      IMPLICIT NONE
      integer, parameter :: Nup = 100
      integer, parameter :: Nelem = Nup*2+1
      integer, parameter :: NWake = Nup*10
      real(kind=8), dimension(Nelem, 2) :: Xnode
      real(kind=8), parameter :: rt = 6.D-2, chord = 1.D0, alpha = 0.D0, L = 1.0D0
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      real(kind=8), dimension(2,2) :: M, K
      real(kind=8) :: thick = chord*rt, vf
      complex(kind=8) :: p = 0*2.D0*PI*dcmplx(0,0.05D0)
      complex(kind=8), dimension(2) :: q
      complex(kind=8), dimension(2) :: f
      complex(kind=8), dimension(2,2) :: E, Et
      integer :: i
      q(1) = 0!-dcmplx(0.1D0)/abs(p)
      q(2) = dcmplx(pi*4.D0/180.,0.)
      E(:,:) = 0
      call GeomNACA00xx(Nup, Chord, Thick, Xnode)
      !call EMatrixWing2d(Nelem, Xnode, alpha, NWake, p, E)
      call EMatrixWing3DA(Nelem, Xnode, alpha, L, Nwake, p, E)
      call theodorsen(p, L, Et)
      write(*,*) "BEM:"
      do i = 1,2
       write(*,1001) E(i,:)
      end do
      write(*,*)
      f = matmul(E,q)
      do i = 1,2
       write(*,1001) f(i)
      end do
      write(*,*) "Thoedorsen:"
      do i = 1,2
       write(*,1001) Et(i,:)
      end do
      write(*,*)
      f = matmul(Et,q)
      do i = 1,2
       write(*,1001) f(i)
      end do
      return
      M(1,:) = (/1.D0, 0.D0/)
      M(2,:) = (/0.D0, .7D0/)
      K(1,:) = (/.0312D0, 0.D0/)
      K(2,:) = (/0.D0, .7D0/)
      call flutterstriptheory(Nup, rt, L, alpha, 5.D0, M, K, Vf)
1001 FORMAT('',100(F15.8))
END PROGRAM

