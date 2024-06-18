subroutine init_theta_arg_vec(ns,pp_vec,gpz_vecs,theta_arg_vec)
use modmain
use modconstants
use modpulse
implicit none
integer, intent(in) :: ns
real(8), intent(in) :: pp_vec(3)
real(8), intent(in) :: gpz_vecs(ns,3)
complex(8), intent(inout) :: theta_arg_vec(ns)
integer :: igpz
real(8) :: p_vec(3)
complex(8) :: norm2_p_vec_cmplx_sq, p_vec_cmplx(3)
complex(8) :: integrand_it, integrand_it_p_1
!when itex equal to 0 then initialize the value of theta_arg_vec as the value at t=dt/2d0
do igpz=1,ns
  p_vec               =pp_vec+gpz_vecs(igpz,:)

  !define integrand at itex=1
  integrand_it=c0

  !define integrand at itex=2
  p_vec_cmplx         =cmplx(p_vec,0d0,8)-pump(2,:)/sol
  norm2_p_vec_cmplx_sq=p_vec_cmplx(1)**2+p_vec_cmplx(2)**2+p_vec_cmplx(3)**2
  integrand_it_p_1=(norm2_p_vec_cmplx_sq/2d0-norm2(p_vec)**2/2d0)

  !sum of trapzoidal area to theta_arg_vec
  theta_arg_vec(igpz)=dt/2d0*(integrand_it+integrand_it_p_1)
enddo
end subroutine
