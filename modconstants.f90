module modconstants
implicit none
real(8), parameter :: astr=1.8897259886d0
real(8), parameter :: eV=0.036749323858378d0
real(8), parameter :: pi=4.d0*ATAN(1.d0)
real(8), parameter :: fs= 41.341374575751d0
real(8), parameter :: sol=137d0
real(8), parameter :: w_func=4.30d0*eV
real(8), parameter :: e_fermi=0.0d0*eV
real(8), parameter :: temperature=100d0
real(8), parameter :: kb=8.6173d0*1e-5*eV
real(8), parameter :: beta=1d0/kb/temperature
real(8), parameter :: t_observe=40d0*fs
real(8), parameter :: radius=100d0
real(8), parameter :: delta=1e-6
complex(8), parameter :: ci=cmplx(0d0,1d0,8)
complex(8), parameter :: c0=cmplx(0d0,0d0,8)
complex(8), parameter :: c1=cmplx(1d0,0d0,8)

contains
include 'r3cross.f90'
end module
