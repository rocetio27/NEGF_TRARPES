module modmain
use modconstants

integer, parameter :: nrmax=3
real(8), parameter :: a=2.46d0*astr
real(8), parameter :: a_1(3)=a*(/1d0,0d0,0d0/)
real(8), parameter :: a_2(3)=a*(/1d0/2d0,dsqrt(3d0)/2d0,0d0/)
real(8), parameter :: b_1(3)=(2*pi/a)*(/1d0,-1d0/dsqrt(3d0),0d0/)
real(8), parameter :: b_2(3)=(2*pi/a)*(/0d0, 2d0/dsqrt(3d0),0d0/)
real(8), parameter :: tau_1_vec(3)=(2d0*a_2-a_1)/3d0
real(8), parameter :: nnv(3,3)=transpose(reshape(a*(/0.d0     ,  1.d0/dsqrt(3.d0)      ,0.d0&
                                                    ,1.d0/2.d0, -1.d0/2.d0/dsqrt(3.d0) ,0.d0&
                                                   ,-1.d0/2.d0, -1.d0/2.d0/dsqrt(3.d0) ,0.d0/),(/3,3/)))
real(8), parameter :: v_pppi_0=-2.7d0*eV
real(8), parameter :: v_ppsgm_0=0.48d0*eV
real(8), parameter :: dlt_0=0.184d0*a
real(8), parameter :: a_0=a/dsqrt(3d0)
real(8), parameter :: d_0=3.35d0*astr

real(8), parameter :: cell_area=&
                      norm2((/a_1(2)*a_2(3)-a_1(3)*a_2(2),-a_1(1)*a_2(3)+a_1(3)*a_2(1),a_1(1)*a_2(2)-a_1(2)*a_2(1)/))
! LINEBANDPLOT
integer :: nkline, nhspts
integer, allocatable :: n_p_to_p_grid(:)
real(8), allocatable :: hspts(:,:), kline(:)

! LAPACK
integer :: lwork,lainfo
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)

! PARAMETERS FOR Hamiltonian
integer, parameter   :: natom=2
real(8), parameter   :: atomic_nums(natom)=6d0
real(8), parameter   :: mu(natom)=2.1d0
real(8)              :: atom_pos(natom,3)
integer              :: norb_per_atom(0:natom)
integer              :: norb ! initialized in main.f90
!USE FOLLOWING ITERATION FOR ORBITALS.
!do i=1,natom
!do imu=1,norb_per_atom(i)
!   index=norb_per_atom(i-1)*(i-1)+imu
!enddo
!enddo
integer              :: nband
! PARAMETERS FOR PES CAL
integer, parameter :: iorbstart=1
logical :: lineplotmode=.true.
integer, parameter :: nppx=41, nppy=41, nenrg=31

!G-vector parameters
integer, parameter :: ngmax=1
integer, parameter :: ngtotax=2*ngmax+1
integer, parameter :: ngpt=ngtotax**2
integer, parameter :: igc=(ngpt-1)/2+1 !center index of g-vecs
real(8) :: g_vecs(ngpt,3)

!pz-grid paramters
integer, parameter :: m_for_pz_max=3     ! From this, pz_max will be determined in the main code as (pz_max=m/2*pz)
integer, parameter :: n_pz_half_grid=100 ! npzpt = 2*m*k + 1
integer, parameter :: npzpt=2*n_pz_half_grid*m_for_pz_max+1

!integer, parameter :: ipzc=(npzpt-1)/2+1 !center index of pz-vecs

! ARRAYS
real(8) :: pxg(nppx),pyg(nppy),enrg(nenrg)
integer :: npppt, nppt
real(8), allocatable :: pp_vecs(:,:)
real(8), allocatable :: p_vecs(:,:)
end module
