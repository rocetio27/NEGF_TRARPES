function hpp(d_vec)
use modmain
use modconstants
implicit none
real(8) d_vec(3), d, dz
real(8) hpp
d=norm2(d_vec)
dz=d_vec(3)
hpp=-v_pppi_0 *dexp(-(d-a_0)/dlt_0)*(1-(dz/d)**2)&
    -v_ppsgm_0*dexp(-(d-d_0)/dlt_0)*   (dz/d)**2
end function
