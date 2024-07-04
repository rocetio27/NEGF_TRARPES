subroutine init_theta_arg_vec(ns,pp_vec,gpz_vecs,theta_arg_vec)
use modmain
use modconstants
use modpulse
implicit none
integer, intent(in) :: ns
real(8), intent(in) :: pp_vec(3)
real(8), intent(in) :: gpz_vecs(ns,3)
real(8), intent(inout) :: theta_arg_vec(ns)
integer :: igpz
real(8) :: p_vec(3)
real(8) :: p_vec_a_field(3)
real(8) :: integrand_it, integrand_it_p_1
do igpz=1,ns
  p_vec               =pp_vec+gpz_vecs(igpz,:)

  !define integrand at itex=1
  integrand_it=c0

  !define integrand at itex=2
  p_vec_a_field         =p_vec-pump(2,:)/sol
  integrand_it_p_1=(sum(p_vec_a_field**2)-sum(p_vec**2))/2d0
  
  !sum of trapzoidal area to theta_arg_vec
   theta_arg_vec(igpz)=(dt/2d0)/2d0*(integrand_it+integrand_it_p_1)
enddo
end subroutine
