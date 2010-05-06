
subroutine theodorsen(pbem, E)
      implicit none
      complex(kind=8), intent(in) :: pbem
      complex(kind=8), dimension(2,2), intent(out) :: E
      real(kind=8) :: k
      complex(kind=8) :: c, p
      real(kind=8), parameter :: PI = 4.D0*datan(1.D0)
      complex(kind=8), parameter :: j = dcmplx(0.D0,1.D0)
      p = pbem/2.D0
      k = cdabs(p)
      E(1,1) = 2.D0*pi*(k**2/5.D0 - 0.4D0*j*k*c(k))
      !E(1,1) = PI*2.D0*(p*(0.4D0*c(k) + p*0.2D0)) !Moltiplicata per 2 PI
      E(1,2) = pi*(-k*(k/8.D0 - j/4.D0) + c(k)*0.5D0*(1.D0+j*k))
      !E(1,2) = -0.5D0*PI*c(k) + p*(-0.25D0-c(k)/2.D0 - 0.125D0*p) !Moltiplicata per PI
      E(2,1) = pi*(-c(k)*j*k/4.D0)
      !E(2,1) = PI*(p*0.25D0*c(k)) !Moltiplicata per PI
      E(2,2) = pi/2.D0*(k*(k/24.D0-j/6.D0) + c(k)/3.D0*(1.D0+j*k))
      !E(2,2) = -0.5D0*PI/3.D0*c(k) + p*(1.D0/6.D0 - c(k) + 1.D0/24.D0*p) !Moltiplicata per PI/2
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

