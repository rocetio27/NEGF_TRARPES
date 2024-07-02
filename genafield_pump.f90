subroutine genafield_pump(pumptmp,profile)
use modmain
use modconstants
use modpulse
implicit none
integer it
real(8)    :: t,t1,t3
complex(8) :: t2
complex(8) :: shape_func
complex(8), intent(inout) :: pumptmp(exitmax,3)
complex(8), intent(inout) :: profile(exitmax)
do it=1,exitmax
  !generate pulse shape
  t=(it-1)*(dt/2d0)
  t1=t+tshift
  t2=cos(w_pump*t1)+ci*sin(-w_pump*t1) !exp(-iwt)
  t3=cos(0.5d0*pi*t1/hd_pump)
  shape_func=t2*(t3**2)
  if ((abs(t1).ge.hd_pump)) shape_func=0.d0
  if (abs(shape_func).lt.1.d-20) shape_func=0.d0
  !define vector field
  pumptmp(it,:)=cmplx(dble(pol_vec_pump*field_intensity_pump*shape_func),0d0,8)
  profile(it)=field_intensity_pump*shape_func
enddo
end subroutine
