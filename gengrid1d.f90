subroutine gengrid1d(iv,fv,ngrid,grid)
implicit none
real(8), intent(in) :: iv, fv
integer, intent(in) :: ngrid
real(8), intent(out) :: grid(ngrid)
integer i
do i=1,ngrid
  grid(i)=iv+(i-1)*(fv-iv)/(ngrid-1)
enddo
end subroutine
