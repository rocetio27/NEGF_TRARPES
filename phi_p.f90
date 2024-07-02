function phi_p(p_vec)
   use modconstants
   implicit none
   real(8) :: p_vec(3)
   real(8) :: a=6d0/2d0
   complex(8) :: phi_p
   complex(8) :: coeff
   coeff=-ci*32d0*dsqrt(pi)*dsqrt(a)**7
   phi_p=coeff*cmplx(p_vec(3)/(a**2+p_vec(1)**2+p_vec(2)**2+p_vec(3)**2)**3,0d0,8)
end function
