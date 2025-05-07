subroutine genafield_probe(probetmp,probetmp_conjg,profile)
use modmain
use modconstants
use modpulse
implicit none
integer it
real(8)    :: t,t1,t3
complex(8) :: t2
complex(8) :: shape_func
real(8), intent(inout) :: probetmp(exitmax,3)
real(8), intent(inout) :: probetmp_conjg(exitmax,3)
complex(8), intent(inout) :: profile(exitmax)

if(w_probe.lt.0d0) then
  write(*,*) "(error): w_probe is negative value"
  stop
endif

do it=1,exitmax
  !generate pulse shape
  t=(it-1)*(dt/2d0)
  t1=t-delta_t+tshift
  t2=cos(-w_probe*t1)+ci*sin(-w_probe*t1) !exp(-iwt)
  t3=cos(0.5d0*pi*t1/hd_probe)
  shape_func=t2*(t3**2) 
  if ((abs(t1).ge.hd_probe)) shape_func=0d0
  if (abs(shape_func).lt.1.d-20) shape_func=0d0
  !define vector field
  probetmp(it,:)=dble(pol_vec_probe*field_intensity_probe*shape_func)
  probetmp_conjg(it,:)=dble(conjg(pol_vec_probe)*field_intensity_probe*shape_func)
  profile(it)=field_intensity_probe*shape_func
enddo
end subroutine
