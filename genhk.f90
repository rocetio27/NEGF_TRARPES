subroutine genhk(k_vec,h)
use modmain
implicit none
complex(8), intent(in) :: k_vec(3)
complex(8), intent(inout) :: h(norb,norb)
integer :: i1,i2,ir1,ir2
real(8) :: d_vec(3)
real(8) :: lr_vec(3)
complex(8) :: cphase
h(:,:)=0.d0
do i1=1,norb
  do i2=1,norb
    if (i1.ne.i2)then
      do ir1=-nrmax,nrmax
        do ir2=-nrmax,nrmax
          lr_vec=ir1*a_1+ir2*a_2
          d_vec=lr_vec+atom_pos(i2,:)-atom_pos(i1,:)
          if(norm2(d_vec).gt.(1.1d0*norm2(tau_1_vec))) cycle
          cphase=k_vec(1)*d_vec(1)+k_vec(2)*d_vec(2)
           h(i1,i2)=h(i1,i2)-hpp(d_vec)*cmplx(cos(dble(cphase)),sin(dble(cphase)),8)*(dexp(-aimag(cphase)))
        enddo
      enddo
    endif
  enddo
enddo
do i1=1,norb
  h(i1,i1)=h(i1,i1)-e_fermi-w_func
enddo
contains
include 'hpp.f90'

end subroutine genhk
