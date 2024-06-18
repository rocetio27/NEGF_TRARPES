
! Copyright (C) 2002-2011 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! main routine for the Elk code
program elk
use modomp
use modmpi
implicit none
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP INITIALIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialise MPI execution environment
call mpi_init(ierror)
! duplicate mpi_comm_world
call mpi_comm_dup(mpi_comm_world,mpicom,ierror)
! determine the number of MPI processes
call mpi_comm_size(mpicom,np_mpi,ierror)
! determine the local MPI process number
call mpi_comm_rank(mpicom,lp_mpi,ierror)
! determine if the local process is the master
if (lp_mpi.eq.0) then
  mp_mpi=.true.
  if (np_mpi.gt.1) then
    write(*,*)
    write(*,'("Using MPI, number of processes : ",I8)') np_mpi
  end if
else
  mp_mpi=.false.
end if

! initialise OpenMP variables
call omp_init
if (mp_mpi) then
  write(*,*)
  write(*,'("Number of CPUs per MPI process : ",I4)') ncpu
  write(*,'("Number of OpenMP threads per MPI process : ",I4)') maxthd
  write(*,'("Nested OpenMP : ",L1)') omp_get_nested()
  write(*,'("Maximum active nested levels : ",I10)') omp_get_max_active_levels()
  write(*,'("OpenMP dynamic thread adjustment : ",L1)') omp_get_dynamic()
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP INITIALIZATION!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
! reset the OpenMP thread variables
call omp_reset
! terminate MPI execution environment
call mpi_finalize(ierror)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!MPI AND OMP STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
stop
end program
