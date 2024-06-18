subroutine get_atom_pos(tmp_atom_pos)
use modconstants
use modmain
implicit none
real(8), intent(inout) :: tmp_atom_pos(natom,3)
!CAUTION, mind atomital species.
tmp_atom_pos(1,:)=tau_1_vec
tmp_atom_pos(2,:)=0d0
!tmp_atom_pos(3,:)=tau_1_vec-(/0d0,0d0,d_0/)
!tmp_atom_pos(4,:)=-(/0d0,0d0,d_0/)


end subroutine
