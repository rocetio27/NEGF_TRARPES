module modpulse
use modconstants
implicit none
real(8) :: t_pump_start
real(8) :: t_probe_start
real(8) :: t_pump_end
real(8) :: t_error
integer :: exitmax
integer :: itmax
real(8) :: dt=0.1d0
real(8), parameter :: theta=0d0/180d0*pi !0: normal incidence with arb.inplane pol 
real(8), parameter :: phi=0d0/180d0*pi    !0: p_pol=x_hat s_pol=y_hat
!incidence vector of incident light
complex(8), parameter :: i_hat(3)=-c1*(/cos(phi)*sin(theta),sin(phi)*sin(theta),cos(theta)/)
!useable polarization vectors of incident light
complex(8), parameter :: s_pol(3)=c1*(/-sin(phi),cos(phi),0d0/)
complex(8), parameter :: p_pol(3)=c1*(/-cos(phi)*cos(theta),-sin(phi)*cos(theta),sin(theta)/)
complex(8), parameter :: lcp_pol(3)=(s_pol-ci*p_pol)
complex(8), parameter :: rcp_pol(3)=(s_pol+ci*p_pol)

real(8), parameter :: w_pump=0.4*eV
real(8), parameter :: e_to_v_pump=eV/astr*137/w_pump ! E-field(eV/astr) to A-field(a.u.) unit conversion factor
real(8), parameter :: field_intensity_pump=0.005d0*e_to_v_pump
complex(8), parameter :: pol_vec_pump(3)=p_pol
real(8), parameter :: phase_pump=0d0
real(8), parameter :: hd_pump=20*fs ! half-duration

!integer, parameter :: min_time_step_pump=20
!real(8), parameter :: t_period_pump=2d0*pi/w_pump
!real(8), parameter :: dt_pump=t_period_pump/min_time_step_pump

real(8), parameter :: w_probe=30*eV
real(8), parameter :: e_to_v_probe=eV/astr*137/w_probe ! E-field to A-field unit conversion factor
real(8), parameter :: field_intensity_probe=0.00005d0*e_to_v_probe
complex(8), parameter :: pol_vec_probe(3)=p_pol
real(8), parameter :: phase_probe=0d0
real(8), parameter :: hd_probe=20*fs ! half-duration

real(8) :: delta_t=0*fs ! delta_t=probe pusle center - pump pulse center
real(8) :: tshift=0d0

real(8), allocatable :: pump(:,:)      ,probe(:,:)       !timeindex,vectorindex
complex(8), allocatable :: pump_profile(:),probe_profile(:) !timeindex,vectorindex
real(8)   , allocatable :: e_pump(:,:)    ,e_probe(:,:)     !timeindex,vectorindex
end module
