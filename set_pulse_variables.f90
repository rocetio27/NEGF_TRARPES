subroutine set_pulse_variables(itmax1,t_pump_start1,t_pump_end1,t_probe_start1,&
                               it_observe1,t_error1,exitmax1)
use modpulse
implicit none
integer, intent(out) :: itmax1
real(8), intent(out) :: t_pump_start1,t_pump_end1,t_probe_start1
integer, intent(out) :: it_observe1
real(8), intent(out) :: t_error1
integer, intent(out) :: exitmax1
if((-1d0*hd_pump).le.(delta_t-hd_probe)) then
  tshift=-1d0*hd_pump !coordinate shift to fit the left-side first
else
  tshift=delta_t-hd_probe !coordinate shift to fit the left-side first
endif
itmax1=ceiling((max(-tshift-hd_pump,-tshift+hd_pump,delta_t-tshift-hd_probe,delta_t-tshift+hd_probe)&
              -min(-tshift-hd_pump,-tshift+hd_pump,delta_t-tshift-hd_probe,delta_t-tshift+hd_probe))/dt)+1
t_pump_start1=-tshift-hd_pump
t_pump_end1=-tshift+hd_pump
t_probe_start1=(t_pump_start1+t_pump_end1)/2d0+delta_t-hd_probe
it_observe1=floor(t_observe/dt)
t_error1=ceiling(t_pump_start1/dt)*dt-t_pump_start1
exitmax1=2*itmax1-1
if(t_probe_start1.gt.t_observe) then
  write(*,*) "observing time is ealier than the probe-start time"
  stop
endif
if(it_observe1.gt.itmax1) then
  write(*,*) "observing time exceeds itmax1"
  stop
endif
end subroutine
