function pw_coeff_bloch(k_vec,g_vec,pz,eigvec_nk)
use modmain
use modconstants
!this function calculate < k+G+pz*z\hat | nk >
implicit none
complex(8) :: pw_coeff_bloch
real(8) :: k_vec(3)
real(8) ::g_vec(3)
real(8) :: pz
complex(8) :: eigvec_nk(nband), phi_p
integer :: j
real(8) :: phase
pw_coeff_bloch=c0
do j=1,nband
 phase=-(g_vec(1)*atom_pos(j,1)+g_vec(2)*atom_pos(j,2)+pz*atom_pos(j,3))
 pw_coeff_bloch=pw_coeff_bloch+eigvec_nk(j)*(dcos(phase)+ci*dsin(phase))*phi_p(k_vec+g_vec+pz*(/0d0,0d0,1d0/))
enddo
pw_coeff_bloch=pw_coeff_bloch/dsqrt(cell_area)
end function

