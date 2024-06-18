subroutine genhk_pw_basis_and_u_trans(ns,theta_arg_vec,pp_vec,gpz_vecs,h_k_tb,h_k_pw)
use modconstants
use modmain
use modpulse
implicit none
integer, intent(in) :: ns
complex(8), intent(in) :: theta_arg_vec(ns)
real(8), intent(in) :: pp_vec(3)
real(8), intent(in) :: gpz_vecs(ns,3)
complex(8), intent(in) :: h_k_tb(norb,norb)
complex(8), intent(out) :: h_k_pw(ns,ns)
integer :: iatom,imu,jnu,igpz, jgpz
real(8) :: gpz_vec_i(3), gpz_vec_j(3)
real(8) :: arg_imu, arg_jnu, arg
complex(8) :: arg_igpz, arg_jgpz, theta_igpz, theta_jgpz
complex(8) :: kinetic_mat_el, pot_mat_el
real(8) :: q_vec(3)
real(8) :: phase
real(8) :: h_box, dpz
h_k_pw=c0
dpz=gpz_vecs(2,3)-gpz_vecs(1,3)
h_box=2d0*pi/dpz
if (.false.) then !
do igpz=1,ns
do jgpz=1,ns
gpz_vec_i=gpz_vecs(igpz,:)
gpz_vec_j=gpz_vecs(jgpz,:)
arg_igpz=theta_arg_vec(igpz)
arg_jgpz=theta_arg_vec(jgpz)
theta_igpz=(dcos(dble(arg_igpz))+ci*dsin(dble(arg_igpz)))*dexp(-aimag(arg_igpz))
theta_jgpz=(dcos(dble(arg_jgpz))+ci*dsin(dble(arg_jgpz)))*dexp(-aimag(arg_jgpz))
do imu=1,nband
do jnu=1,nband
  arg_imu=gpz_vec_i(1)*atom_pos(imu,1)+gpz_vec_i(2)*atom_pos(imu,2)-gpz_vec_i(3)*atom_pos(imu,3)
  arg_jnu=gpz_vec_j(1)*atom_pos(jnu,1)+gpz_vec_j(2)*atom_pos(jnu,2)-gpz_vec_j(3)*atom_pos(jnu,3)
  arg=arg_imu-arg_jnu
  h_k_pw(igpz,jgpz)=h_k_pw(igpz,jgpz)+&
                    h_k_tb(imu,jnu)*(dcos(arg)+ci*dsin(arg))

enddo
enddo
 h_k_pw(igpz,jgpz)=h_k_pw(igpz,jgpz)*theta_igpz*conjg(theta_jgpz)*&
                        phi_p(pp_vec+gpz_vec_i)*conjg(phi_p(pp_vec+gpz_vec_j))/cell_area/h_box
enddo
enddo
endif !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!free electron potential
if (.false.) then
do igpz=1,ns
do jgpz=1,ns
  if (igpz.ne.jgpz) cycle
  gpz_vec_i=gpz_vecs(igpz,:)
  h_k_pw(igpz,jgpz)=((pp_vec(1)+gpz_vec_i(1))**2+(pp_vec(2)+gpz_vec_i(2))**2+(pp_vec(3)+gpz_vec_i(3))**2)/2d0
enddo
enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!phenomenological potential
if (.true.) then
do igpz=1,ns
do jgpz=1,ns
  gpz_vec_i     =gpz_vecs(igpz,:)
  gpz_vec_j     =gpz_vecs(jgpz,:)

  arg_igpz=theta_arg_vec(igpz)
  arg_jgpz=theta_arg_vec(jgpz)

  theta_igpz=(dcos(dble(arg_igpz))+ci*dsin(dble(arg_igpz)))*dexp(-aimag(arg_igpz))
  theta_jgpz=(dcos(dble(arg_jgpz))+ci*dsin(dble(arg_jgpz)))*dexp(-aimag(arg_jgpz))

  kinetic_mat_el=c0
  pot_mat_el    =c0
  if(igpz.eq.jgpz) then
    kinetic_mat_el=((pp_vec(1)+gpz_vec_i(1))**2+(pp_vec(2)+gpz_vec_i(2))**2+(pp_vec(3)+gpz_vec_i(3))**2)/2d0
  endif
  q_vec=gpz_vec_i-gpz_vec_j
  do iatom=1,natom
    phase=-(q_vec(1)*atom_pos(iatom,1)+q_vec(2)*atom_pos(iatom,2)+q_vec(3)*atom_pos(iatom,3))
    pot_mat_el=pot_mat_el-4d0*pi*atomic_nums(iatom)/(norm2(q_vec)**2+mu(iatom)**2)*(dcos(phase)+ci*dsin(phase))
  enddo
  !h_k_pw(igpz,jgpz)=kinetic_mat_el+pot_mat_el/cell_area
  h_k_pw(igpz,jgpz)=kinetic_mat_el+theta_igpz*conjg(theta_jgpz)*pot_mat_el/cell_area/h_box
enddo
enddo
endif

contains
include 'phi_p.f90'
end subroutine genhk_pw_basis_and_u_trans
