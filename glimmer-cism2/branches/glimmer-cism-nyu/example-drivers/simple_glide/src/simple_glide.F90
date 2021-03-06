! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  simple_glide.f90 - part of the Glimmer-CISM ice model    + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

program simple_glide
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases
  use glimmer_global, only:rk
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init

  use glide_diagnostics

  implicit none

#ifdef GPTL
#include <gptl.inc>
#endif
#ifdef PAPI
#include <f90papi.h>
#endif

#ifdef GLIMMER_MPI
#include <mpif.h>
  integer nproc,ierr,irank
#endif

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret

  integer :: tstep_count

  logical,parameter :: do_restore_inflow_thk = .true.

  real(kind=dp),dimension(:,:),allocatable :: prev_thk
  logical :: is_steady = .false.

  ! start gptl
#ifdef GPTL
  ret = gptlsetoption (gptlprint_method,gptlfull_tree)
  ret = gptlsetoption (PAPI_FP_OPS, 1)
  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('total')
#endif

#ifdef GLIMMER_MPI
 call MPI_Init(ierr)
! call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
! call MPI_Comm_rank(MPI_COMM_WORLD,irank,ierr)
#endif
  call glimmer_GetCommandline()

  ! start logging
  call open_log(unit=50, fname=logname(commandline_configname))
  
  ! setup paths
  call filenames_init(commandline_configname)

  ! read configuration
  print *, commandline_configname
  call ConfigRead(trim(commandline_configname),config)

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! initialise GLIDE
  call glide_config(model,config)
  call simple_initialise(climate,config)
  call glide_initialise(model)

  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart

  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)

  tstep_count = 0

  allocate(prev_thk(model%general%ewn,model%general%nsn))
  prev_thk = 0.d0

  do while(time.le.model%numerics%tend .and. .not.(is_steady))

     !have to specify temperature of inflowing ice for ice_stream case
     model%temper%temp(:,:,1:3) = climate%airt(1)
     ! reset western boundary temp
     model%temper%temp(:,1:2,:) = climate%airt(1)
     !  and also eastern boundary temp
     model%temper%temp(:,(model%general%ewn-2):(model%general%ewn-1),:) = climate%airt(1)

     call glide_tstep_p1(model,time)
   
     call glide_tstep_p2(model,no_write=.true.)

     ! zero out the efvs in the inflow edge
     model%stress%efvs(:,:,1:2) = 0.d0

     call glide_tstep_p3(model,no_write=.false.)

     ! override masking stuff for now
     tstep_count = tstep_count + 1
     if (mod(tstep_count, model%numerics%ndiag) == 0) then
        call glide_write_diag(model, time, model%numerics%idiag, &
                                           model%numerics%jdiag )
     endif

     time = time + model%numerics%tinc
     
     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)     

     call check_for_steady(model,is_steady)

  end do

  ! finalise GLIDE
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log

  deallocate(prev_thk)

#ifdef GLIMMER_MPI
 call MPI_Finalize(ierr)
#endif

  ! stop gptl
#ifdef GPTL
  ret = gptlstop ('total')
  ret = gptlpr (0)
  ret = gptlfinalize ()
#endif


contains  

subroutine check_for_steady(model, is_steady)

 type(glide_global_type),intent(in) :: model
 logical,intent(out)                :: is_steady

 real(kind=dp),dimension(model%general%ewn,model%general%nsn) :: rel_thk_change
 real(kind=dp) :: max_rel_thk_change
 real(kind=dp),parameter :: steadiness_tol = 5.0d-4

 where (model%geometry%thck > 0.d0)
     rel_thk_change = abs(model%geometry%thck - prev_thk)/model%geometry%thck
 else where
     rel_thk_change = 0.d0
 end where

 max_rel_thk_change = maxval(rel_thk_change)

 prev_thk = model%geometry%thck

 print *, 'max_rel_thk_change', max_rel_thk_change

 if (max_rel_thk_change < steadiness_tol) then
    is_steady = .true.
 else
    is_steady = .false. 
 end if

 end subroutine check_for_steady


end program simple_glide