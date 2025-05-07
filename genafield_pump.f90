subroutine genafield_pump(pumptmp,profile)
use modmain
use modconstants
use modpulse
implicit none
integer it
real(8)    :: t,t1,t3, t3p, phase,is
complex(8) :: t2, t2p
complex(8) :: shape_func, shape_funcp
real(8), intent(inout) :: pumptmp(exitmax,3)
complex(8), intent(inout) :: profile(exitmax)
is=2d0/3d0
phase=0d0
do it=1,exitmax
  !generate pulse shape
  t=(it-1)*(dt/2d0)
  t1=t+tshift
  t2 =cos(w_pump*t1)+ci*sin(-w_pump*t1) !exp(-iwt)
  !t2p=is*(cos(2d0*w_pump*t1+phase)+ci*sin(-2d0*w_pump*t1+phase)) !intensity_scale*exp(-2iwt)
  t3 =cos(0.5d0*pi*t1/hd_pump)
  shape_func = t2*(t3**2)
  !shape_funcp=t2p*(t3**2)
  if ((abs(t1).ge.hd_pump)) shape_func=0.d0
  if (abs(shape_func).lt.1.d-20) shape_func=0.d0
  !if ((abs(t1).ge.hd_pump)) shape_funcp=0.d0
  !if (abs(shape_funcp).lt.1.d-20) shape_funcp=0.d0
  !define vector field
  !pumptmp(it,:)=dble(pol_vec_pump*field_intensity_pump*shape_func+conjg(pol_vec_pump)*field_intensity_pump*shape_funcp)
  pumptmp(it,:)=dble(pol_vec_pump*field_intensity_pump*shape_func)
  profile(it)=field_intensity_pump*shape_func
enddo
end subroutine
