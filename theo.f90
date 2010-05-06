subroutine theodorsen(pbem, L, E)
      implicit none
      complex(kind=8), intent(in) :: pbem
      complex(kind=8), dimension(2,2), intent(out) :: E
      real(kind=8) :: k, L
      complex(kind=8) :: c, p
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      complex(kind=8), parameter :: j = dcmplx(0.D0,1.D0)
      p = pbem/2.D0
      k = imagpart(p)
      E(1,1) = 0.5D0*(k**2/5.D0 - 0.4D0*j*k*c(k))
      E(1,2) = (-k*(k/8.D0 - j/4.D0) + c(k)*0.5D0*(1.D0+j*k))
      E(2,1) = 0.5D0*(-c(k)*j*k/4.D0)
      E(2,2) = (k*(k/24.D0-j/6.D0) + c(k)/3.D0*(1.D0+j*k))
      E = 2.D0*PI*L*E
end subroutine


subroutine theodorsenST(pbem, L, E)
      implicit none
      complex(kind=8), intent(in) :: pbem
      complex(kind=8), dimension(2,2), intent(out) :: E
      real(kind=8) :: k, L
      complex(kind=8) :: c, p
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      complex(kind=8), parameter :: j = dcmplx(0.D0,1.D0)
      p = pbem/2.D0
      k = imagpart(p)
      E(1,1) = 0.5D0*(pi*k**2/32.D0 - pi*j*k*c(k)/16.D0)
      E(1,2) = ( 2.D0/15.D0*(-k*(k/2.D0 - j)) + c(k)*4.D0/15.D0*(1.D0+j*k) )
      E(2,1) = 0.5D0*(-c(k)*j*k*2.D0/15.D0)
      E(2,2) = (pi/32.D0*k*(k/4.D0 - j) + c(k)*pi/16.D0*(1.D0+j*k))
      E = 2.D0*PI*L*E
end subroutine


function c(k)
      implicit none
      complex(kind=8) :: c
      real(kind=8) :: k
      complex(kind=8) :: jk
      jk = dcmplx(0.D0,k)
      c = 0.5D0*(jk+0.135D0)*(jk+0.651D0)/((jk+0.0956D0)*(jk+0.4555D0))
      return
end function

