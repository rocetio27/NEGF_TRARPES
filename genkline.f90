subroutine genkline(hspts1,kvecs1,kline1)
use modmain
implicit none
real(8), intent(in) :: hspts1(nhspts,3)
real(8), intent(out) :: kvecs1(nkline,3)
real(8), intent(out) :: kline1(nkline)
integer :: i1,i2,i3
real(8) :: temp_v3(3)

  i3=1
  kvecs1(i3,:)=hspts1(1,:)
  kline1(i3)=0
  do i1=1,(nhspts-1)
    temp_v3(:)=(hspts1(i1+1,:)-hspts1(i1,:))/n_p_to_p_grid(i1)
    do i2=1,n_p_to_p_grid(i1)
      i3=i3+1
      kvecs1(i3,:)=kvecs1(i3-1,:)+temp_v3(:)
      kline1(i3)=kline1(i3-1)+sqrt(temp_v3(1)**2+temp_v3(2)**2+temp_v3(3)**2)
    enddo
  enddo
  
endsubroutine


