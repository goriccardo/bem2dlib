!Copyright (c) 2009 Riccardo Gori <goriccardo@gmail.com>
!Released under BSD license, see LICENSE

!Create a vector representing a grid for field contour plotting
subroutine fieldgrid(xmin, xmax, nx, ymin, ymax, ny, xyvec)
  integer, intent(in) :: nx, ny
  real(kind=8), intent(in) :: xmin, xmax, ymin, ymax
  real(kind=8), dimension(nx*ny,2), intent(out) :: xyvec
  integer :: i,j
  real(kind=8) :: dx, dy
  dx = (xmax-xmin)/real(nx-1,8)
  dy = (ymax-ymin)/real(ny-1,8)
  do i = 1,ny
   xyvec((i-1)*nx+1:i*nx,2) = ymin + dy*(i-1)
   do j = 1,nx
    xyvec((i-1)*nx+j,1) = xmin + dx*(j-1)
   end do
  end do
end subroutine

