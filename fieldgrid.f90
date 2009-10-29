!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Create a vector representing a grid for field contour plotting
subroutine fieldgrid(xmin, xmax, nx, ymin, ymax, ny, xyvec)
  integer, intent(in) :: nx, ny
  real, intent(in) :: xmin, xmax, ymin, ymax
  real, dimension(nx*ny,2), intent(out) :: xyvec
  integer :: i,j
  real :: dx, dy

  dx = (xmax-xmin)/real(nx-1)
  dy = (ymax-ymin)/real(ny-1)
  do i = 1,nx
    xyvec((i-1)*ny:(i)*ny,1) = xmin + dx*(i-1)
    do j = 1,ny
      xyvec((i-1)*ny+j,2) = ymin + dy*(j-1)
    end do
  end do
end subroutine
