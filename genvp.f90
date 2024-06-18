subroutine genvp(p_vec,iextp,v_mat)
use modconstants
use modmain
use modpulse
real(8), intent(in) :: p_vec(3)
integer, intent(in) :: iextp
complex(8), intent(inout) :: v_mat(ns,ns)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: ig, iqz, jg, jqz, i, j
real(8) :: qi_vec(3), qj_vec(3)
real(8) ::    e_pqi   , e_pqj
complex(8) :: arg_pqi , arg_pqj
complex(8) ::   phase_pqi, phase_pqj_conjg
v_mat=c0
do ig=1,ngpt
do iqz=1,nqzpt
i=nqzpt*(ig-1)+iqz
qi_vec=g_vecs(ig,:)+(/0d0,0d0,pz_grid(iqz)/)
e_pqi=norm2(p_vec+qi_vec)**2/2d0+w_func
arg_pqi=cmplx(e_pqi*(dt/2d0)*(iextp-1),0d0,8)+&
        ((p_vec(1)+qi_vec(1))*pump_integral_vec(iextp,1)+&
        (p_vec(2)+qi_vec(2))*pump_integral_vec(iextp,2)+&
        (p_vec(3)+qi_vec(3))*pump_integral_vec(iextp,3))/sol
phase_pqi=cmplx(cos(dble(arg_pqi)),sin(dble(arg_pqi)),8)*(dexp(-aimag(arg_pqi)))
do jg=1,ngpt
do jqz=1,nqzpt
j=nqzpt*(jg-1)+jqz
qj_vec=g_vecs(jg,:)+(/0d0,0d0,pz_grid(jqz)/)
e_pqj=norm2(p_vec+qj_vec)**2/2d0+w_func
arg_pqj=cmplx(e_pqj*(dt/2d0)*(iextp-1),0d0,8)+&
        ((p_vec(1)+qj_vec(1))*pump_integral_vec(iextp,1)+&
        (p_vec(2)+qj_vec(2))*pump_integral_vec(iextp,2)+&
        (p_vec(3)+qj_vec(3))*pump_integral_vec(iextp,3))/sol
phase_pqj_conjg=cmplx(cos(dble(arg_pqj)),-sin(dble(arg_pqj)),8)*(dexp(aimag(arg_pqj)))

!v_mat(i,j)=conjg(lvq(-qi_vec+qj_vec))*phase_pqi*phase_pqj_conjg
enddo
enddo
enddo
enddo
v_mat=v_mat/cell_area/h_box
end subroutine
