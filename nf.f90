function nf(enrg)
    use modconstants
    implicit none
    real(8) :: enrg
    real(8) :: nf
    nf=1d0/(exp(beta*enrg)+1d0)
end function
