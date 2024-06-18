subroutine genefield_probe(probetmp)
use modmain
use modconstants
use modpulse
implicit none
integer it
real(8) gs,t,t1,t2,t3
real(8) gsp,tp,t1p,t2p,t3p
real(8), intent(inout) :: probetmp(exitmax,3)
do it=1,exitmax
  t=(it-1)*(dt/2d0)
  tp=it*(dt/2d0)

  t1=t-delta_t+tshift
  t1p=tp-delta_t+tshift
  
  t2=cos(w_probe*t1)
  t2p=cos(w_probe*t1p)

  t3=cos(0.5*pi*t1/hd_probe)
  t3p=cos(0.5*pi*t1p/hd_probe)
  
  if(circular_probe) then
    gs=t3**2
    gsp=t3p**2
  else
    gs=t2*(t3**2)
    gsp=t2p*(t3p**2)
  endif
  
  if ((abs(t1).ge.hd_probe)) gs=0.d0
  if ((abs(t1p).ge.hd_probe)) gsp=0.d0
  if (abs(gs).lt.1.d-20) gs=0.d0
  if (abs(gsp).lt.1.d-20) gsp=0.d0
  
  if(circular_probe) then
    probetmp(it,:)=-afield_intensity_cprobe*((/cos(w_cprobe*t1p),sin(w_cprobe*t1p),0.d0/)*gsp&
                                             -(/cos(w_cprobe*t1),sin(w_cprobe*t1),0.d0/)*gs)/137d0/(dt/2d0)
  else
    probetmp(it,:)=-afield_vec_probe(:)*(gsp-gs)/137d0/(dt/2d0)
  endif
enddo
end subroutine
