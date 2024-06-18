program main
use modmain
use modconstants
use modpulse
use modmpi
use modomp
implicit none
integer :: ie, itex, ipp, ip, ipz, in, inp, ig, ig1, ig2, ipx, ipy, i, j, imy
integer :: it, it_pb_start, it_observe,igpz,igpzp, igpz1
real(8) :: p_vec(3)=0d0, p1_vec(3)=0d0
real(8) :: pp_vec(3)=0d0
real(8) :: g_vec(3)=0d0
real(8) :: origin_vec(3)=0d0
real(8) :: ec, er
real(8) :: k_vec(3)=0d0
!pzgrid variables-----------------------------------------
real(8) :: delta_e=0d0, pzmax=0.1d0
integer :: pzcount=0, ns=0, npzpt=0, ipzc, igpzc, num_positive_delta_e
real(8), allocatable :: pz_grid(:)
!pzgrid variables-----------------------------------------
complex(8) :: t_integral_vec1(3)=c0
complex(8) :: t_integral_vec2(3)=c0
complex(8), allocatable :: pes(:)
real(8), allocatable :: pes_real(:), pes_imag(:)
complex(8), allocatable :: eigvec(:)
real(8)   , allocatable :: eigvals(:)
complex(8), allocatable :: eigvec_mat(:,:)
complex(8), allocatable :: eigvect_mat(:,:)
complex(8), allocatable :: temp_hamilt(:,:)
real(8) :: height, width
integer :: nmpi, ompthnum, nthd
real(8) :: normalizer, temp
!!!!!!!!!!!!!!!!!!!!!!!!!!VERSION2 VARIABLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8) :: pz=0d0
complex(8) :: theta_p_t, theta_arg_p_t
complex(8), allocatable :: hk_t_lesser(:,:)
complex(8), allocatable :: hk_t_greater(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!UNIRATY MATRIX RELATED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(8), allocatable :: gpz_vecs(:,:)
complex(8), allocatable :: temp_col_vec_lesser(:)
complex(8), allocatable :: pes_t_integrand(:,:)
complex(8), allocatable :: pes_t_integrand_t_p_dt(:,:)
complex(8), allocatable :: pes_t_integral_val1(:,:)
complex(8), allocatable :: pes_t_integral_val2(:)
complex(8) :: dipole_mat_el, greater_mat_el, lesser_mat_el, greater_dagger_mat_el
complex(8), allocatable :: theta_arg_vec(:)
complex(8), allocatable :: u_mat_greater(:,:), u_mat_greater_dagger(:,:), u_mat_greater_ptime(:,:)
complex(8), allocatable :: u_mat_lesser(:,:), u_mat_lesser_ptime(:,:)
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
write(*,*) "pol_vec_probe", pol_vec_probe
!define p-parallel_vecs
!----------------------------------------------------------
if(lineplotmode.eqv.(.true.)) then
!----------------------------------------------------------
    nhspts=2 !!
    allocate(hspts(nhspts,3)) !3 means px,py,pz
    allocate(n_p_to_p_grid(nhspts-1))

    hspts(1,:)=0.66666666d0*b_1+0.33333333d0*b_2-(/0.0d0, 0.12d0, 0.0d0/) !!
    hspts(2,:)=0.66666666d0*b_1+0.33333333d0*b_2+(/0.0d0, 0.12d0, 0.0d0/) !!
    write(*,*) "(info) origin_vec:", (hspts(1,:)+hspts(2,:))/2d0
    !hspts(1,:)=-(/1.1d0, 0d0, 0.0d0/) !!
    !hspts(2,:)= (/1.1d0, 0d0, 0.0d0/) !!
    n_p_to_p_grid=(/96/) !!
    nkline=sum(n_p_to_p_grid)+1
    npppt=nkline

    allocate(pp_vecs(npppt,3))
    allocate(kline(npppt))
    call genkline(hspts,pp_vecs,kline)

!--------------------------------------------------------
elseif(lineplotmode.eqv.(.false.)) then
!--------------------------------------------------------
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
!----------------------------------------------------------
endif
!define energy grid
!CONSTANT ENERGY----------------------------------------------
if(lineplotmode.eqv.(.false.)) then
  enrg=w_probe-w_func !here
!LINE ENERGY----------------------------------------------
elseif(lineplotmode.eqv.(.true.)) then
  er=1.2d0*eV                       !!energy range
  ec=w_probe-w_func                    !!energy center
  call gengrid1d(ec-er/2d0,ec+er/2d0,nenrg,enrg)
endif

!gernerate p_vecs from pp_vecs and enrg
if(lineplotmode.eqv..false.) then
nppt=npppt*1
else
nppt=npppt*nenrg
endif

allocate(pes(nppt))
allocate(pes_real(nppt),pes_imag(nppt))
pes=c0
pes_real=c0
pes_imag=c0
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

!define g_vecs
do ig1=-ngmax,ngmax
do ig2=-ngmax,ngmax
  g_vecs((ig1+ngmax)*ngtotax+(ig2+ngmax+1),:)=ig1*b_1+ig2*b_2
enddo
enddo
!define pz_grid
!call gengrid1d(-pzmax,pzmax,npzpt,pz_grid)
!dpz=pz_grid(2)-pz_grid(1)
!h_box=2d0*pi/dpz

!define g+pz vecs
!do ig=1,ngpt
! do ipz=1,npzpt
!   imy=npzpt*(ig-1)+ipz
!   gpz_vecs(imy,:)=g_vecs(ig,:)+(/0d0,0d0,pz_grid(ipz)/)
! enddo
!enddo

!define orb_pos
norb_per_atom(0)=0
norb_per_atom(1)=1
norb_per_atom(2)=1
norb=sum(norb_per_atom)
nband=norb
write(*,*) "norb:",norb
allocate(temp_col_vec_lesser(nband))
allocate(ipiv_lesser(nband),work_lesser(nband))
allocate(hk_t_lesser(nband,nband))
allocate(identity_lesser(nband,nband))
allocate(u_mat_lesser(nband,nband))
allocate(u_mat_lesser_ptime(nband,nband))
allocate(u_mat_l_lesser(nband,nband))
allocate(u_mat_r_lesser(nband,nband))
allocate(u_mat_evolver_lesser(nband,nband))
allocate(temp_cmat_lesser(nband,nband))

!allocate(ipiv_greater(ns),work_greater(ns))
!allocate(theta_arg_vec(ns))
!allocate(hk_t_greater(ns,ns))
!allocate(identity_greater(ns,ns))
!allocate(u_mat_greater(ns,ns))
!allocate(u_mat_greater_dagger(ns,ns))
!allocate(u_mat_greater_ptime(ns,ns))
!allocate(u_mat_l_greater(ns,ns))
!allocate(u_mat_r_greater(ns,ns))
!allocate(u_mat_evolver_greater(ns,ns))
!allocate(temp_cmat_greater(ns,ns))
!allocate(pes_t_integrand(nband,ns))
!allocate(pes_t_integrand_t_p_dt(nband,ns))
!allocate(pes_t_integral_val1(nband,ns))
allocate(pes_t_integral_val2(nband))

allocate(eigvec(nband),eigvals(nband),eigvec_mat(nband,nband),eigvect_mat(nband,nband),temp_hamilt(nband,nband))
eigvec(:)=c0
eigvals(:)=0d0
eigvec_mat(:,:)=c0
eigvect_mat(:,:)=c0
temp_hamilt(:,:)=c0
call get_atom_pos(atom_pos)

!pump-probe generation
  if((-1d0*hd_pump).le.(delta_t-hd_probe)) then
    tshift=-1d0*hd_pump !coordinate shift to fit the left-side first
  else
    tshift=delta_t-hd_probe !coordinate shift to fit the left-side first
  endif
  itmax=ceiling((max(-tshift-hd_pump,-tshift+hd_pump,delta_t-tshift-hd_probe,delta_t-tshift+hd_probe)&
       -min(-tshift-hd_pump,-tshift+hd_pump,delta_t-tshift-hd_probe,delta_t-tshift+hd_probe))/dt)+1
  t_pump_start=-tshift-hd_pump
  t_pump_end=-tshift+hd_pump
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t_probe_start=(t_pump_start+t_pump_end)/2d0+delta_t-hd_probe
  write(*,*) "t_probe_start", t_probe_start/fs
  write(*,*) "t_probe_start", (delta_t-tshift-hd_probe)/fs
  if(t_probe_start.gt.t_observe) then
  write(*,*) "observing time is ealier than the probe-start time"
  stop
  endif
  it_pb_start=floor(t_probe_start/dt)
  it_observe=floor(t_observe/dt)
  
  if(it_observe.gt.itmax) then
  write(*,*) "observing time exceeds itmax"
  stop
  endif
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  t_error=ceiling(t_pump_start/dt)*dt-t_pump_start

  exitmax=2*itmax-1
  allocate(pump(exitmax,3),probe(exitmax,3))
  allocate(pump_profile(exitmax),probe_profile(exitmax))
  !allocate(e_probe(exitmax,3))
  call genafield_pump(pump,pump_profile)
  call genafield_probe(probe,probe_profile)
  !call genefield_probe(e_probe)

  !allocation for the zheev
  if (allocated(work)) deallocate(work)
  if (allocated(rwork)) deallocate(rwork)
  lwork = 2*nband-1
  allocate(work(lwork),rwork(3*nband-2))
!MESSAGES AND PLOTS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if(mp_mpi) then
     write(*,27) "(Info) pump start(fs):", t_pump_start/fs
     27 format (a40,g18.10,g18.10)
     write(*,27) "(Info) probe start(fs):", (delta_t-tshift-hd_probe)/fs
     write(*,27) "(Info) pump half duration(fs):",hd_pump/fs
     write(*,27) "(Info) probe half duration(fs):",hd_probe/fs
     write(*,27) "(Info) pol_vec_pump:",  dble(pol_vec_pump), aimag(pol_vec_probe)
     write(*,27) "(Info) pol_vec_probe:", dble(pol_vec_probe), aimag(pol_vec_probe)
     write(*,27) "(Info) nwpt:", nwpt
     write(*,27) "(Info) delay(fs):",delta_t/fs
     write(*,27) "(Info) itmax:", itmax
     write(*,27) "(Info) exitmax:",exitmax
     write(*,27) "(Info) lineplotmode:", lineplotmode
     write(*,27) "(Info) nenrg:", nenrg
     write(*,27) "(Info) omega_half_range:", omega_half_range/eV
    if(lineplotmode) then
     write(*,27) "(Info) nkline :",nkline
     else
     write(*,27) "(Info) npx,npy:",nppx,nppy
     endif
     write(*,27) "(Info) the number of in-plane p points:", npppt
     write(*,27) "(Info) ngmax, ngpt:",ngmax, ngpt
     write(*,27) "(Info) t_observe(fs), radius(astr):", t_observe/fs, radius/astr
     !-------------------------------------------------------
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

     !open(50,file='eprobe.txt',form='FORMATTED')
     !  do itex=1,exitmax
     !  write(50,'(5G18.10)'),itex, itex*(dt/2.d0), e_probe(itex,1), e_probe(itex,2), e_probe(itex,3)
     !  enddo
     !close(50)
     if(lineplotmode) then
     !open(50,file='eigenvalues.txt',form='FORMATTED')
     !  do i=1,norb
     !  do ip=1,nppt
     !      write(50,'(5G18.10)'), kline(ip), eigvals(ip,i)
     !  enddo
     !  enddo
     !close(50)
     endif
     !-------------------------------------------------------
   endif
!==========================================================
identity_lesser=c0
do i=1,nband
  identity_lesser(i,i)=c1
enddo
!photocurrent calculation
call mpi_barrier(mpicom,ierror)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ip,it,ig,in,inp,i,j,ipz,igpz,igpz1,igpzp,imy)&
!$OMP PRIVATE(p_vec, p1_vec, pp_vec, pz, g_vec, k_vec, gpz_vecs)&
!$OMP PRIVATE(delta_e,num_positive_delta_e,pzcount,npzpt,ns,identity_greater)&
!$OMP PRIVATE(temp_hamilt,eigvec_mat,eigvect_mat,eigvals,normalizer,temp)&
!$OMP PRIVATE(t_integral_vec1,t_integral_vec2)&
!$OMP PRIVATE(lwork,work,rwork,lainfo,info)&
!$OMP PRIVATE(work_lesser,ipiv_lesser,work_greater,ipiv_greater)&
!$OMP PRIVATE(theta_arg_vec,theta_p_t,theta_arg_p_t)&
!$OMP PRIVATE(pes_t_integral_val1,pes_t_integral_val2,pes_t_integrand,pes_t_integrand_t_p_dt)&
!$OMP PRIVATE(dipole_mat_el, greater_mat_el, greater_dagger_mat_el, lesser_mat_el)&
!$OMP PRIVATE(u_mat_lesser,u_mat_lesser_ptime,u_mat_greater,u_mat_greater_dagger,u_mat_greater_ptime)&
!$OMP PRIVATE(u_mat_evolver_lesser,u_mat_evolver_greater)&
!$OMP PRIVATE(temp_cmat_lesser,temp_cmat_greater,temp_col_vec_lesser)&
!$OMP PRIVATE(hk_t_lesser,hk_t_greater)&
!$OMP PRIVATE(u_mat_l_lesser,u_mat_r_lesser,u_mat_l_greater,u_mat_r_greater)&
!$OMP NUM_THREADS(20)
ompthnum=omp_get_thread_num()
nthd=omp_get_num_threads()
info=0
lwork=2*nband-1
rwork=0d0
work_lesser=c0
ipiv_lesser=0

temp_hamilt=c0
eigvec_mat=c0
eigvect_mat=c0
eigvals=0d0

u_mat_evolver_lesser=c0
theta_arg_p_t=c0
!$OMP DO
do ip=1,nppt
if (mod(ip-1,np_mpi).ne.lp_mpi) cycle
  write(*,*) "(info) ip:", ip, ompthnum, np_mpi
  p_vec=p_vecs(ip,:)
  pp_vec=(/p_vec(1),p_vec(2),0d0/)
  pz=p_vec(3)
  !count npzpt define ns
  num_positive_delta_e=0
  do ig=1,ngpt
    g_vec=g_vecs(ig,:)
    delta_e=norm2(p_vec)**2-norm2(pp_vec+g_vec)**2
    if (delta_e.ge.0d0) then
    num_positive_delta_e=num_positive_delta_e+1
    endif
  enddo
  !ns=num_positive_delta_e*2 !here
   ns=193                   !here
! allocate memories for variables depends on ns
allocate(gpz_vecs(ns,3))
allocate(ipiv_greater(ns),work_greater(ns))
allocate(theta_arg_vec(ns))
allocate(hk_t_greater(ns,ns))
allocate(identity_greater(ns,ns))
allocate(u_mat_greater(ns,ns))
allocate(u_mat_greater_dagger(ns,ns))
allocate(u_mat_greater_ptime(ns,ns))
allocate(u_mat_l_greater(ns,ns))
allocate(u_mat_r_greater(ns,ns))
allocate(u_mat_evolver_greater(ns,ns))
allocate(temp_cmat_greater(ns,ns))
allocate(pes_t_integrand(nband,ns))
allocate(pes_t_integrand_t_p_dt(nband,ns))
allocate(pes_t_integral_val1(nband,ns))
identity_greater=c0
do igpz=1,ns
  identity_greater(igpz,igpz)=c1
enddo
work_greater=c0
ipiv_greater=0
u_mat_evolver_greater=c0
  !define g+pz vecs
  pzcount=0
  !do ig=1,ngpt
   do ig=igc,igc
   g_vec=g_vecs(ig,:) 
   delta_e=norm2(p_vec)**2-norm2(pp_vec+g_vec)**2
   if (delta_e.ge.0d0) then
   !gpz_vecs(pzcount+1,:)=g_vec+(/0d0,0d0,dsqrt(delta_e)/)
   !gpz_vecs(pzcount+2,:)=g_vec+(/0d0,0d0,-dsqrt(delta_e)/)
   !if(ig.eq.igc) igpzc=pzcount+1
   !pzcount=pzcount+2
   do igpz=1,ns
     gpz_vecs(igpz,:)=g_vec
   enddo
   call gengrid1d(-4d0*pz,4d0*pz,ns,gpz_vecs(:,3)) !here
   igpzc=(5*ns+3)/8
   endif
  enddo
  !do ig=1,ngpt
  ! do ipz=1,npzpt
  !   imy=npzpt*(ig-1)+ipz
  !   gpz_vecs(imy,:)=g_vecs(ig,:)+(/0d0,0d0,pz+pz_grid(ipz)/)
  ! enddo
  !enddo

  if(p_vec(1).eq.999d0) cycle
  !generate ground state hamiltonian to temp_hamilt
  call genhk(cmplx(pp_vec,0d0,8),temp_hamilt)
  eigvec_mat=temp_hamilt
  !generate ground eigen states of temp_hamilt
  call zheev('V','U',nband,eigvec_mat,nband,eigvals,work,lwork,rwork,lainfo)
  !normalization and gauge setting of groundstates
  !*------------------------------------------------
  do i=1,nband
    normalizer=0d0
    do j=1,nband
    normalizer=normalizer+dble(eigvec_mat(j,i))**2+aimag(eigvec_mat(j,i))**2
    enddo
    normalizer=dsqrt(normalizer)
    eigvec_mat(:,i)=eigvec_mat(:,i)/normalizer
    eigvec_mat(:,i)=eigvec_mat(:,i)/eigvec_mat(1,i)*dsqrt(dble(eigvec_mat(1,i))**2+aimag(eigvec_mat(1,i))**2)
  enddo
  eigvect_mat=eigvec_mat
  !*-----INITIALIZATION OF Up(t,t0)--------------------------
  u_mat_lesser=identity_lesser
  u_mat_greater=identity_greater
  !*-----INITIALIZATION OF pes_t_integrand_t_p_dt------------
  pes_t_integrand       =c0
  pes_t_integrand_t_p_dt=c0
  pes_t_integral_val1   =c0
  pes_t_integral_val2   =c0
  !*-----TIME EVOLUTION OF Up(t,t0)--------------------------
  call init_theta_arg_vec(ns,pp_vec,gpz_vecs,theta_arg_vec)
  do it=1,it_observe-1
  ! in the do it=1,~ iteration u_mat_lesser_ptime corresponds to it, and
  ! u_mat_laesser to it+1
  u_mat_lesser_ptime=u_mat_lesser
  call genhk(cmplx(pp_vec,0d0,8)-pump(2*it,:)/sol,hk_t_lesser)
  u_mat_r_lesser=identity_lesser-ci*(dt/2d0)*hk_t_lesser
  u_mat_l_lesser=identity_lesser+ci*(dt/2d0)*hk_t_lesser
  call zgetrf(nband,nband,u_mat_l_lesser,nband,ipiv_lesser,info)
  call zgetri(nband,u_mat_l_lesser,nband,ipiv_lesser,work_lesser,nband,info)
  call zgemm('n','n',nband,nband,nband,c1,u_mat_l_lesser,nband,u_mat_r_lesser,nband,c0,u_mat_evolver_lesser,nband)
  call zgemm('n','n',nband,nband,nband,c1,u_mat_evolver_lesser,nband,u_mat_lesser,nband,c0,temp_cmat_lesser,nband)
  u_mat_lesser=temp_cmat_lesser

  u_mat_greater_ptime=u_mat_greater !save the matrix of it-1 time
  call genhk_pw_basis_and_u_trans(ns,theta_arg_vec,pp_vec,gpz_vecs,temp_hamilt,hk_t_greater)
  u_mat_r_greater=identity_greater-ci*(dt/2d0)*hk_t_greater
  u_mat_l_greater=identity_greater+ci*(dt/2d0)*hk_t_greater
  call zgetrf(ns,ns,u_mat_l_greater,ns,ipiv_greater,info)
  call zgetri(ns,u_mat_l_greater,ns,ipiv_greater,work_greater,ns,info)
  call zgemm('n','n',ns,ns,ns,c1,u_mat_l_greater,ns,u_mat_r_greater,ns,c0,u_mat_evolver_greater,ns)
  call zgemm('n','n',ns,ns,ns,c1,u_mat_evolver_greater,ns,u_mat_greater,ns,c0,temp_cmat_greater,ns)
  u_mat_greater=temp_cmat_greater
  u_mat_greater_dagger=transpose(conjg(u_mat_greater))
  call evolve_theta_arg_vec(ns,2*it,pp_vec,gpz_vecs,theta_arg_vec)
  !*--------------------------------------------------------
  do in=1,nband
  do igpzp=1,ns
  do inp=1,nband
     lesser_mat_el=c0
     call zgemv('N', nband, nband, c1, u_mat_lesser, nband, eigvec_mat(:,in), 1, c0, temp_col_vec_lesser, 1)
     do i=1,nband
     lesser_mat_el=lesser_mat_el+conjg(eigvec_mat(i,inp))*temp_col_vec_lesser(i)
     enddo
  do igpz1=1,ns
     !generation pes_t_integrand at at it+1, the value at it is initialized c0
     p1_vec=pp_vec+gpz_vecs(igpz1,:)
     theta_p_t=(dcos(dble(theta_arg_vec(igpz1)))+ci*dsin(dble(theta_arg_vec(igpz1))))*dexp(-aimag(theta_arg_vec(igpz1)))
     dipole_mat_el=pw_coeff_bloch(pp_vec,(/gpz_vecs(igpz1,1),gpz_vecs(igpz1,2),0d0/),gpz_vecs(igpz1,3),eigvec_mat(:,inp))*&
     (probe(2*it+1,1)*p1_vec(1)+probe(2*it+1,2)*p1_vec(2)+probe(2*it+1,3)*p1_vec(3))
     greater_dagger_mat_el=u_mat_greater_dagger(igpzp,igpz1)
     
     pes_t_integrand_t_p_dt(in,igpzp)=&
     pes_t_integrand_t_p_dt(in,igpzp)+theta_p_t*greater_dagger_mat_el*dipole_mat_el*lesser_mat_el
  enddo !do for p1
  enddo !do for band index n'
     pes_t_integral_val1(in,igpzp)=pes_t_integral_val1(in,igpzp)+dt/2d0*(pes_t_integrand(in,igpzp)+pes_t_integrand_t_p_dt(in,igpzp))
     pes_t_integrand(in,igpzp)=pes_t_integrand_t_p_dt(in,igpzp)
  enddo !do for p'
  enddo !do for band index n
  !*---------------------------------------------------------
  call evolve_theta_arg_vec(ns,2*it+1,pp_vec,gpz_vecs,theta_arg_vec)
  enddo !do for t-integral
  !*---------------------------------------------------------
  do in=1,nband
  do igpzp=1,ns
     pes_t_integral_val2(in)=pes_t_integral_val2(in)+u_mat_greater(igpzc,igpzp)*pes_t_integral_val1(in,igpzp)
  enddo
     pes_real(ip)=pes_real(ip)+nf(eigvals(in)+w_func)*(dble(pes_t_integral_val2(in))**2+aimag(pes_t_integral_val2(in))**2)
  enddo
deallocate(gpz_vecs)
deallocate(ipiv_greater,work_greater)
deallocate(hk_t_greater)
deallocate(u_mat_greater)
deallocate(u_mat_greater_dagger)
deallocate(u_mat_greater_ptime)
deallocate(u_mat_l_greater)
deallocate(u_mat_r_greater)
deallocate(u_mat_evolver_greater)
deallocate(identity_greater)
deallocate(temp_cmat_greater)
deallocate(pes_t_integral_val1)
deallocate(pes_t_integrand)
deallocate(pes_t_integrand_t_p_dt)
deallocate(theta_arg_vec)
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
!deallocate(pump,probe,e_probe)
deallocate(pump,probe,pump_profile,probe_profile)
deallocate(pp_vecs,p_vecs)
deallocate(pes)
deallocate(pes_real,pes_imag)
deallocate(eigvec,eigvals,eigvec_mat,eigvect_mat,temp_hamilt)
deallocate(ipiv_lesser,work_lesser)
deallocate(hk_t_lesser)
deallocate(u_mat_lesser)
deallocate(u_mat_lesser_ptime)
deallocate(u_mat_l_lesser)
deallocate(u_mat_r_lesser)
deallocate(u_mat_evolver_lesser)
deallocate(identity_lesser)
deallocate(temp_cmat_lesser)
deallocate(pes_t_integral_val2)
deallocate(temp_col_vec_lesser)

!deallocate(ipiv_greater,work_greater)
!deallocate(hk_t_greater)
!deallocate(u_mat_greater)
!deallocate(u_mat_greater_dagger)
!deallocate(u_mat_greater_ptime)
!deallocate(u_mat_l_greater)
!deallocate(u_mat_r_greater)
!deallocate(u_mat_evolver_greater)
!deallocate(identity_greater)
!deallocate(temp_cmat_greater)
!deallocate(pes_t_integral_val1)
!deallocate(pes_t_integrand)
!deallocate(pes_t_integrand_t_p_dt)
!deallocate(theta_arg_vec)

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
