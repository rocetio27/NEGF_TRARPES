program main
use modmain
use modconstants
use modpulse
use modmpi
use modomp
implicit none
integer :: ie, itex, ipp, ip, ipz, in, inp, ig, ig1, ig2, ipx, ipy, i, j, imy
integer :: it, it_observe,itex_half, igpz, igpz1
integer :: igc2
integer, allocatable :: set_ig(:)
real(8) :: p_vec(3)=0d0, p1_vec(3)=0d0
real(8) :: pp_vec(3)=0d0
real(8) :: g_vec(3)=0d0
real(8) :: origin_vec(3)=0d0
real(8) :: ec, er
real(8) :: k_vec(3)=0d0
!pzgrid variables-----------------------------------------
real(8) :: delta_e=0d0, pzmax=0.1d0
integer :: ns=0, npzpt=0, ipzc, igpzc, ngpt2
real(8), allocatable :: pz_grid(:)
!pzgrid variables-----------------------------------------
complex(8) :: t_integral_vec1(3)=c0
complex(8) :: t_integral_vec2(3)=c0
complex(8), allocatable :: pes(:)
real(8), allocatable :: pes_real(:), pes_imag(:)
complex(8), allocatable :: eigvec(:)
real(8)   , allocatable :: eigvals(:)
complex(8), allocatable :: eigvec_mat(:,:)
complex(8), allocatable :: temp_hamilt(:,:)
real(8) :: height, width
integer :: nmpi, ompthnum, nthd
real(8) :: normalizer, temp
!!!!!!!!!!!!!!!!!!!!!!!!!!VERSION2 VARIABLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8) :: pz=0d0
complex(8) :: theta_p_t, theta_arg_p_t
complex(8), allocatable :: hk_t_lesser(:,:)
complex(8), allocatable :: hk_t_greater(:,:)
complex(8), allocatable :: hk_pw(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!UNIRATY MATRIX RELATED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8), allocatable :: gpz_vecs(:,:)
complex(8), allocatable :: temp_col_vec_lesser(:)
complex(8), allocatable :: pes_t_integrand(:)
complex(8), allocatable :: pes_t_integrand_t_p_dt(:)
complex(8), allocatable :: pes_t_integral_val1(:)
complex(8) :: dipole_mat_el, greater_mat_el, lesser_mat_el, greater_dagger_mat_el
complex(8), allocatable :: theta_arg_vec(:)
complex(8), allocatable :: u_mat_greater(:,:), u_mat_greater_dagger(:,:)
complex(8), allocatable :: u_mat_lesser(:,:)
complex(8), allocatable :: u_mat_l_greater(:,:), u_mat_r_greater(:,:)
complex(8), allocatable :: u_mat_l_lesser(:,:), u_mat_r_lesser(:,:)
complex(8), allocatable :: u_mat_evolver_greater(:,:)
complex(8), allocatable :: u_mat_evolver_lesser(:,:)
complex(8), allocatable :: identity_lesser(:,:)
complex(8), allocatable :: identity_greater(:,:)
complex(8), allocatable :: work_lesser(:)
complex(8), allocatable :: work_greater(:)
integer, allocatable :: ipiv_lesser(:), ipiv_greater(:)
integer :: info
complex(8), allocatable :: temp_cmat_lesser(:,:)
complex(8), allocatable :: temp_cmat_greater(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP INITIALIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call mpi_init(ierror)
call mpi_comm_dup(mpi_comm_world,mpicom,ierror)
call mpi_comm_size(mpicom,np_mpi,ierror)
call mpi_comm_rank(mpicom,lp_mpi,ierror)
if (lp_mpi.eq.0) then
  mp_mpi=.true.
else
  mp_mpi=.false.
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP INITIALIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*---------------------------------------------------------------------------------------------------------
!   define in-plane p-vectors
!*---------------------------------------------------------------------------------------------------------
if(lineplotmode.eqv.(.true.)) then
    nhspts=2 !!
    allocate(hspts(nhspts,3)) !3 means px,py,pz
    allocate(n_p_to_p_grid(nhspts-1))

    hspts(1,:)=0.66666666d0*b_1+0.33333333d0*b_2-(/0.0d0, 0.12d0, 0.0d0/) !!
    hspts(2,:)=0.66666666d0*b_1+0.33333333d0*b_2+(/0.0d0, 0.12d0, 0.0d0/) !!
    write(*,*) "(info) origin_vec:", (hspts(1,:)+hspts(2,:))/2d0
    !hspts(1,:)=-(/1.1d0, 0d0, 0.0d0/) !!
    !hspts(2,:)= (/1.1d0, 0d0, 0.0d0/) !!
    n_p_to_p_grid=(/48/) !!
    nkline=sum(n_p_to_p_grid)+1
    npppt=nkline

    allocate(pp_vecs(npppt,3))
    allocate(kline(npppt))
    call genkline(hspts,pp_vecs,kline)
elseif(lineplotmode.eqv.(.false.)) then
    npppt=nppx*nppy
    allocate(pp_vecs(npppt,3))

    origin_vec=0.66666666d0*b_1+0.33333333d0*b_2 !!
    write(*,*) "(info) origin_vec:", origin_vec
    width=0.06d0                               !!
    height=0.06d0                              !!
    call gengrid1d(origin_vec(1)-width/2d0 ,origin_vec(1)+width/2d0 ,nppx,pxg)
    call gengrid1d(origin_vec(2)-height/2d0,origin_vec(2)+height/2d0,nppy,pyg)
    do ipy=1,nppy
    do ipx=1,nppx
      pp_vecs(nppx*(ipy-1)+ipx,:)=(/pxg(ipx),pyg(ipy),0.d0/)
    enddo
    enddo
endif
!*---------------------------------------------------------------------------------------------------------
!   define energy grid
!*---------------------------------------------------------------------------------------------------------
if(lineplotmode.eqv.(.false.)) then
  enrg=w_probe-w_func !here
else
  er=1.2d0*eV !energy range
  ec=w_probe-w_func !energy center
  call gengrid1d(ec-er/2d0,ec+er/2d0,nenrg,enrg)
endif
!*---------------------------------------------------------------------------------------------------------
!   define p-vectors
!*---------------------------------------------------------------------------------------------------------
if(lineplotmode.eqv..false.) then
nppt=npppt*1
else
nppt=npppt*nenrg
endif

allocate(p_vecs(nppt,3))
ip=1
do ipp=1,npppt
do ie=1,nenrg
  if((2*enrg(ie)-pp_vecs(ipp,1)**2-pp_vecs(ipp,2)**2).lt.0d0) then
  write(*,*) "(warning)"
  p_vecs(ip,:)=999d0
  ip=ip+1
  else
  p_vecs(ip,:)=(/pp_vecs(ipp,1),pp_vecs(ipp,2),dsqrt(2*enrg(ie)-pp_vecs(ipp,1)**2-pp_vecs(ipp,2)**2)/)
  ip=ip+1
  endif
enddo
enddo
write(*,*) "nkline:", nkline
write(*,*) "nenrg:", nenrg
write(*,*) "nppt:", nppt
!*---------------------------------------------------------------------------------------------------------
!   define G-vectros
!*---------------------------------------------------------------------------------------------------------
do ig1=-ngmax,ngmax
do ig2=-ngmax,ngmax
  g_vecs((ig1+ngmax)*ngtotax+(ig2+ngmax+1),:)=ig1*b_1+ig2*b_2
enddo
enddo
!*---------------------------------------------------------------------------------------------------------
!   define atom positions and norb
!*---------------------------------------------------------------------------------------------------------
norb_per_atom(0)=0
norb_per_atom(1)=1
norb_per_atom(2)=1
norb=sum(norb_per_atom)
call get_atom_pos(atom_pos)
!*---------------------------------------------------------------------------------------------------------
!   pump - probe generation
!*---------------------------------------------------------------------------------------------------------
call set_pulse_variables(itmax,t_pump_start,t_pump_end,t_probe_start,it_observe,t_error,exitmax)
allocate(pump(exitmax,3),probe(exitmax,3))
allocate(pump_profile(exitmax),probe_profile(exitmax))
call genafield_pump(pump,pump_profile)
call genafield_probe(probe,probe_profile)
!*---------------------------------------------------------------------------------------------------------
!   generation of pump probe raw data
!*---------------------------------------------------------------------------------------------------------
if(mp_mpi) then
 open(50,file='pump.txt',form='FORMATTED')
  do itex=1,exitmax
  write(50,'(5G18.10)'),itex, itex*(dt/2.d0), real(pump(itex,1)), real(pump(itex,2)), real(pump(itex,3))
  enddo
 close(50)

 open(50,file='probe.txt',form='FORMATTED')
  do itex=1,exitmax
  write(50,'(5G18.10)'),itex, itex*(dt/2.d0), real(probe(itex,1)), real(probe(itex,2)), real(probe(itex,3))
  enddo
 close(50)
endif
!*---------------------------------------------------------------------------------------------------------
!   lesser variables                                              
!*---------------------------------------------------------------------------------------------------------
nband=norb
allocate(temp_col_vec_lesser(nband))
allocate(ipiv_lesser(nband),work_lesser(nband))
allocate(hk_t_lesser(nband,nband))
allocate(identity_lesser(nband,nband))
allocate(u_mat_lesser(nband,nband))
allocate(u_mat_l_lesser(nband,nband))
allocate(u_mat_r_lesser(nband,nband))
allocate(u_mat_evolver_lesser(nband,nband))
allocate(temp_cmat_lesser(nband,nband))
allocate(eigvec(nband),eigvals(nband),eigvec_mat(nband,nband),temp_hamilt(nband,nband))
if (allocated(work)) deallocate(work)
if (allocated(rwork)) deallocate(rwork)
lwork = 2*nband-1
allocate(work(lwork),rwork(3*nband-2))
eigvec(:)=c0
eigvals(:)=0d0
eigvec_mat(:,:)=c0
temp_hamilt(:,:)=c0
identity_lesser=c0
do i=1,nband
  identity_lesser(i,i)=c1
enddo
!*---------------------------------------------------------------------------------------------------------
!   PES variables               
!*---------------------------------------------------------------------------------------------------------
allocate(pes(nppt))
allocate(pes_real(nppt),pes_imag(nppt))
pes=c0
pes_real=c0
pes_imag=c0
call mpi_barrier(mpicom,ierror)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip,it,ig,ig2,igc2,in,inp,i,j,ipz,igpz,igpz1,imy,itex_half,set_ig)&
!$OMP PRIVATE(p_vec, p1_vec, pp_vec, pz, g_vec, k_vec, gpz_vecs)&
!$OMP PRIVATE(delta_e,ngpt2,npzpt,ns,identity_greater)&
!$OMP PRIVATE(temp_hamilt,eigvec_mat,eigvals,normalizer,temp)&
!$OMP PRIVATE(t_integral_vec1,t_integral_vec2)&
!$OMP PRIVATE(lwork,work,rwork,lainfo,info)&
!$OMP PRIVATE(work_lesser,ipiv_lesser,work_greater,ipiv_greater)&
!$OMP PRIVATE(theta_arg_vec,theta_p_t,theta_arg_p_t)&
!$OMP PRIVATE(pes_t_integral_val1,pes_t_integrand,pes_t_integrand_t_p_dt)&
!$OMP PRIVATE(dipole_mat_el, greater_mat_el, greater_dagger_mat_el, lesser_mat_el)&
!$OMP PRIVATE(u_mat_lesser,u_mat_greater,u_mat_greater_dagger)&
!$OMP PRIVATE(u_mat_evolver_lesser,u_mat_evolver_greater)&
!$OMP PRIVATE(temp_cmat_lesser,temp_cmat_greater,temp_col_vec_lesser)&
!$OMP PRIVATE(hk_t_lesser,hk_t_greater,hk_pw)&
!$OMP PRIVATE(u_mat_l_lesser,u_mat_r_lesser,u_mat_l_greater,u_mat_r_greater)&
!$OMP NUM_THREADS(24)
ompthnum=omp_get_thread_num()
nthd=omp_get_num_threads()
info=0
lwork=2*nband-1
rwork=0d0
work_lesser=c0
ipiv_lesser=0
temp_hamilt=c0
eigvec_mat=c0
eigvals=0d0
u_mat_evolver_lesser=c0
theta_arg_p_t=c0



!$OMP DO
do ip=1,nppt
  if (mod(ip-1,np_mpi).ne.lp_mpi) cycle
  write(*,*) "(info) ip:", ip, ompthnum, np_mpi
  p_vec=p_vecs(ip,:)
  if(p_vec(1).eq.999d0) cycle
  pp_vec=(/p_vec(1),p_vec(2),0d0/)
  pz=p_vec(3)
  !*-------------------------------------------------------
  !    Define ngpt2, set_ig and igc2
  !*-------------------------------------------------------
  ngpt2=0
  do ig=1,ngpt
  !do ig=igc,igc
    g_vec=g_vecs(ig,:)
    if (sum(p_vec**2)-sum((pp_vec+g_vec)**2).ge.0d0) ngpt2=ngpt2+1
  enddo
  write(*,*) "(info) ngpt2:", ngpt2
  allocate(set_ig(ngpt2))
  i=1
  do ig=1,ngpt
  !do ig=igc,igc
    if (sum(p_vec**2)-sum((pp_vec+g_vec)**2).ge.0d0) then
      set_ig(i)=ig
      if (ig.eq.igc) igc2=i
      i=i+1
    endif
  enddo
  !*-------------------------------------------------------
  npzpt=201
  ns=ngpt2*npzpt
  !*-------------------------------------------------------
  ! allocate memories for variables depends on ns
  !*-------------------------------------------------------
  allocate(gpz_vecs(ns,3))
  allocate(ipiv_greater(ns),work_greater(ns))
  allocate(theta_arg_vec(ns))
  allocate(hk_t_greater(ns,ns))
  allocate(hk_pw(ns,ns))
  allocate(identity_greater(ns,ns))
  allocate(u_mat_greater(ns,ns))
  allocate(u_mat_greater_dagger(ns,ns))
  allocate(u_mat_l_greater(ns,ns))
  allocate(u_mat_r_greater(ns,ns))
  allocate(u_mat_evolver_greater(ns,ns))
  allocate(temp_cmat_greater(ns,ns))
  allocate(pes_t_integrand(nband))
  allocate(pes_t_integrand_t_p_dt(nband))
  allocate(pes_t_integral_val1(nband))
  identity_greater=c0
  hk_pw=c0
  work_greater=c0
  ipiv_greater=0
  u_mat_evolver_greater=c0
  do igpz=1,ns
    identity_greater(igpz,igpz)=c1
  enddo
  !*-------------------------------------------------------------------------------
  !   Define gpz_vecs
  !*-------------------------------------------------------------------------------
  do ig2=1,ngpt2
    g_vec=g_vecs(set_ig(ig2),:)
    do ipz=1,npzpt
       gpz_vecs(npzpt*(ig2-1)+ipz,:)=g_vec
    enddo
    call gengrid1d(-4d0*pz,4d0*pz,npzpt,gpz_vecs(npzpt*(ig2-1)+1:npzpt*ig2,3)) !here
  enddo
  igpzc=npzpt*(igc2-1)+(5*npzpt+3)/8
  !*-------------------------------------------------------------------------------
  !   Ground state calculation
  !*-------------------------------------------------------------------------------
  call genhk(cmplx(pp_vec,0d0,8),temp_hamilt)
  call genhk_pw_basis(ns,pp_vec,gpz_vecs,temp_hamilt,hk_pw)
  eigvec_mat=temp_hamilt
  call zheev('V','U',nband,eigvec_mat,nband,eigvals,work,lwork,rwork,lainfo)
  !*-------------------------------------------------------------------------------
  !   normalization and gauge setting of groundstates
  !*-------------------------------------------------------------------------------
    do i=1,nband
    normalizer=0d0
    do j=1,nband
    normalizer=normalizer+dble(eigvec_mat(j,i))**2+aimag(eigvec_mat(j,i))**2
    enddo
    normalizer=dsqrt(normalizer)
    eigvec_mat(:,i)=eigvec_mat(:,i)/normalizer
    eigvec_mat(:,i)=eigvec_mat(:,i)/eigvec_mat(1,i)*dsqrt(dble(eigvec_mat(1,i))**2+aimag(eigvec_mat(1,i))**2)
    enddo
  !*-----INITIALIZATION OF U(t,t0)--------------------------
  u_mat_lesser=identity_lesser
  u_mat_greater=identity_greater
  !*-----INITIALIZATION OF pes_t_integrand_t_p_dt------------
  pes_t_integrand       =c0
  pes_t_integrand_t_p_dt=c0
  pes_t_integral_val1   =c0
  !*---------------------------------------------------------
  ! codes below generate Ug(t,t0)
  !*---------------------------------------------------------
  call init_theta_arg_vec(ns,pp_vec,gpz_vecs,theta_arg_vec)
  do it=1,it_observe-1
  call u_transform(ns,theta_arg_vec,hk_pw,hk_t_greater) !!!
  u_mat_r_greater=identity_greater-ci*(dt/2d0)*hk_t_greater
  u_mat_l_greater=identity_greater+ci*(dt/2d0)*hk_t_greater
  call zgetrf(ns,ns,u_mat_l_greater,ns,ipiv_greater,info)
  call zgetri(ns,u_mat_l_greater,ns,ipiv_greater,work_greater,ns,info)
  call zgemm('n','n',ns,ns,ns,c1,u_mat_l_greater,ns,u_mat_r_greater,ns,c0,u_mat_evolver_greater,ns)
  call zgemm('n','n',ns,ns,ns,c1,u_mat_evolver_greater,ns,u_mat_greater,ns,c0,temp_cmat_greater,ns)
  u_mat_greater=temp_cmat_greater
  call evolve_theta_arg_vec(ns,2*it,pp_vec,gpz_vecs,theta_arg_vec)
  call evolve_theta_arg_vec(ns,2*it+1,pp_vec,gpz_vecs,theta_arg_vec)
  !endif
  enddo
  !*------------------------------------------------------------------------
  ! codes below do time integration for the time vairable t1.
  !   for each t1,
  !   Ug(t,t1) is generated by evolving Ug(t,t0) step by step 
  !   and Ul(t1,t0) from the Ul(t0,t0)
  !*------------------------------------------------------------------------
  call init_theta_arg_vec(ns,pp_vec,gpz_vecs,theta_arg_vec)
  do it=1,it_observe-1
  call genhk(cmplx(pp_vec,0d0,8)-pump(2*it,:)/sol,hk_t_lesser)
  u_mat_r_lesser=identity_lesser-ci*(dt/2d0)*hk_t_lesser
  u_mat_l_lesser=identity_lesser+ci*(dt/2d0)*hk_t_lesser
  call zgetrf(nband,nband,u_mat_l_lesser,nband,ipiv_lesser,info)
  call zgetri(nband,u_mat_l_lesser,nband,ipiv_lesser,work_lesser,nband,info)
  call zgemm('n','n',nband,nband,nband,c1,u_mat_l_lesser,nband,u_mat_r_lesser,nband,c0,u_mat_evolver_lesser,nband)
  call zgemm('n','n',nband,nband,nband,c1,u_mat_evolver_lesser,nband,u_mat_lesser,nband,c0,temp_cmat_lesser,nband)
  u_mat_lesser=temp_cmat_lesser

  call u_transform(ns,theta_arg_vec,hk_pw,hk_t_greater) !!!
  u_mat_r_greater=identity_greater-ci*(dt/2d0)*hk_t_greater
  u_mat_l_greater=identity_greater+ci*(dt/2d0)*hk_t_greater
  call zgetrf(ns,ns,u_mat_l_greater,ns,ipiv_greater,info)
  call zgetri(ns,u_mat_l_greater,ns,ipiv_greater,work_greater,ns,info)
  call zgemm('n','n',ns,ns,ns,c1,u_mat_l_greater,ns,u_mat_r_greater,ns,c0,u_mat_evolver_greater,ns)
  u_mat_evolver_greater=transpose(conjg(u_mat_evolver_greater))
  call zgemm('n','n',ns,ns,ns,c1,u_mat_greater,ns,u_mat_evolver_greater,ns,c0,temp_cmat_greater,ns)
  u_mat_greater=temp_cmat_greater
  !*------------------------------------------------------------------------------------------------------------------------
  ! code below define integrand function for t1 corresponding to (it+1) - pes_t_integrand_t_p_dt
  ! and do trapzoidal summation with integrand at it
  !*------------------------------------------------------------------------------------------------------------------------
  call evolve_theta_arg_vec(ns,2*it,pp_vec,gpz_vecs,theta_arg_vec)
  !*------------------------------------------------------------------------------------------------------------------------
  do in=1,nband
  do inp=1,nband
     lesser_mat_el=c0
     call zgemv('N', nband, nband, c1, u_mat_lesser, nband, eigvec_mat(:,in), 1, c0, temp_col_vec_lesser, 1)
     do i=1,nband
     lesser_mat_el=lesser_mat_el+conjg(eigvec_mat(i,inp))*temp_col_vec_lesser(i)
     enddo
  do igpz1=1,ns
     p1_vec=pp_vec+gpz_vecs(igpz1,:)
     theta_p_t=(dcos(dble(theta_arg_vec(igpz1)))+ci*dsin(dble(theta_arg_vec(igpz1))))*dexp(-aimag(theta_arg_vec(igpz1)))
     dipole_mat_el=pw_coeff_bloch(pp_vec,(/gpz_vecs(igpz1,1),gpz_vecs(igpz1,2),0d0/),gpz_vecs(igpz1,3),eigvec_mat(:,inp))*&
     (probe(2*it+1,1)*p1_vec(1)+probe(2*it+1,2)*p1_vec(2)+probe(2*it+1,3)*p1_vec(3))
     
     pes_t_integrand_t_p_dt(in)=pes_t_integrand_t_p_dt(in)&
                               +theta_p_t*u_mat_greater(igpzc,igpz1)*dipole_mat_el*lesser_mat_el
  enddo !do for p1
  enddo !do for band index n'
     pes_t_integral_val1(in)=pes_t_integral_val1(in)+dt/2d0*(pes_t_integrand(in)+pes_t_integrand_t_p_dt(in))
     pes_t_integrand(in)=pes_t_integrand_t_p_dt(in)
  enddo !do for band index n
  !*------------------------------------------------------------------------------------------------------------------------
  call evolve_theta_arg_vec(ns,2*it+1,pp_vec,gpz_vecs,theta_arg_vec)
  !*------------------------------------------------------------------------------------------------------------------------
  enddo !do for t-integral
  !*------------------------------------------------------------------------------------------------------------------------
  do in=1,nband
     pes_real(ip)=pes_real(ip)+nf(eigvals(in)+w_func)*(dble(pes_t_integral_val1(in))**2+aimag(pes_t_integral_val1(in))**2)
  enddo
  !*------------------------------------------------------------------------------------------------------------------------
deallocate(gpz_vecs)
deallocate(ipiv_greater,work_greater)
deallocate(hk_t_greater)
deallocate(hk_pw)
deallocate(u_mat_greater)
deallocate(u_mat_greater_dagger)
deallocate(u_mat_l_greater)
deallocate(u_mat_r_greater)
deallocate(u_mat_evolver_greater)
deallocate(identity_greater)
deallocate(temp_cmat_greater)
deallocate(pes_t_integral_val1)
deallocate(pes_t_integrand)
deallocate(pes_t_integrand_t_p_dt)
deallocate(theta_arg_vec)
deallocate(set_ig)
enddo !do for p_vec
!$OMP END DO
!$OMP END PARALLEL
call mpi_barrier(mpicom,ierror)
if (np_mpi.gt.1) then
  nmpi=nppt
  call mpi_allreduce(mpi_in_place,pes_real,nmpi,mpi_double_precision,mpi_sum,mpicom,ierror)
endif
!write pes
if (mp_mpi) then
if(lineplotmode.eqv..false.) then
    open(50,file='pes.txt',form='FORMATTED')
      ip=1
      do ipp=1,npppt
      do ie=1,nenrg
        write(50,'(5G18.10)') pp_vecs(ipp,1),pp_vecs(ipp,2),enrg(ie),pes_real(ip), pes_imag(ip)
        ip=ip+1
      enddo
      enddo
    close(50)
else
    open(50,file='pes.txt',form='FORMATTED')
      ip=1
      do ipp=1,npppt
      do ie=1,nenrg
        write(50,'(5G18.10)') kline(ipp), enrg(ie), pes_real(ip)
        ip=ip+1
      enddo
      enddo
    close(50)
endif

write(*,*) "calculation finished"

endif

if(lineplotmode.eqv.(.true.)) then
deallocate(hspts)
deallocate(n_p_to_p_grid)
deallocate(kline)
endif
deallocate(work,rwork)
deallocate(pump,probe,pump_profile,probe_profile)
deallocate(pp_vecs,p_vecs)
deallocate(pes)
deallocate(pes_real,pes_imag)
deallocate(eigvec,eigvals,eigvec_mat,temp_hamilt)
deallocate(ipiv_lesser,work_lesser)
deallocate(hk_t_lesser)
deallocate(u_mat_lesser)
deallocate(u_mat_l_lesser)
deallocate(u_mat_r_lesser)
deallocate(u_mat_evolver_lesser)
deallocate(identity_lesser)
deallocate(temp_cmat_lesser)
deallocate(temp_col_vec_lesser)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! terminate MPI execution environment
call mpi_finalize(ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
stop
contains
include 'nf.f90'
include 'phi_p.f90'
include 'pw_coeff_bloch.f90'
end program
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
