subroutine u_transform(ns,theta_arg_vec,h_k_pw,h_k_pw_ut)
use modconstants
use modmain
use modpulse
implicit none
integer   , intent(in)  :: ns
real(8), intent(in)  :: theta_arg_vec(ns)
complex(8), intent(in)  :: h_k_pw(ns,ns)
complex(8), intent(inout) :: h_k_pw_ut(ns,ns)
integer :: i,j
complex(8) :: exp_i_theta_vec(ns)
complex(8) :: conjg_exp_i_theta_vec(ns)
exp_i_theta_vec=cdexp(ci*theta_arg_vec)
conjg_exp_i_theta_vec=conjg(exp_i_theta_vec)

do i=1,ns
do j=1,ns
  h_k_pw_ut(i,j)=exp_i_theta_vec(i)*conjg_exp_i_theta_vec(j)*h_k_pw(i,j)
enddo
enddo
end subroutine u_transform
