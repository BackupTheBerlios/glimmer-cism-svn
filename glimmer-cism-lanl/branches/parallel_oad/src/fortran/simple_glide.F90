! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  simple_glide.f90 - part of the GLIMMER ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif


subroutine numericalCore(model, climate)
  use glimmer_global, only:rk, dp
  use glide
  use simple_forcing
  use glide_diagnostics

  implicit none

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  real(kind=rk) time
  integer :: tstep_count

  time = model%numerics%tstart

  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)

  tstep_count = 0

  do while(time.le.model%numerics%tend)
     call glide_tstep_p1(model,time)
     call glide_tstep_p2(model)
     call glide_tstep_p3(model)
     ! override masking stuff for now

     tstep_count = tstep_count + 1
     if (mod(tstep_count, model%numerics%ndiag) == 0) then
        call glide_write_diag(model, time, model%numerics%idiag, &
                                           model%numerics%jdiag )
     endif

     time = time + model%numerics%tinc

     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)     

  end do

  call glide_finalise(model)

end subroutine

subroutine set_ctrl(inittemp,inCtrl)
  use glimmer_global, only:dp
  REAL(dp), DIMENSION(:, :, :) :: inittemp
  REAL(dp) inCtrl
  ! stub for automatic differentiation
  inittemp(1,1,1)=inittemp(1,1,1)+inCtrl
end subroutine 

subroutine numericalCoreWithModelInit(model, climate, config, inCtrl, outCtrl)
  use glide
  use simple_forcing
  use glimmer_config
  use glide_diagnostics 

  implicit none

  type(glide_global_type)      :: model 
  type(ConfigSection), pointer :: config
  type(simple_climate)         :: climate 
  real(dp) inCtrl, outCtrl

  interface 
     subroutine numericalCore(model, climate)
       use glide
       use simple_forcing
       type(glide_global_type) :: model 
       type(simple_climate) :: climate 
     end subroutine numericalCore
  end interface
  interface
     subroutine set_ctrl(inittemp,inCtrl)
       use glimmer_global, only:dp
       REAL(dp), DIMENSION(:, :, :) :: inittemp
       REAL(dp) inCtrl
     end subroutine set_ctrl
  end interface

  call glide_config(model,config)
  call simple_initialise(climate,config)
  call glide_initialise(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  call set_ctrl(model%tempwk%inittemp,inCtrl)
  call numericalCore(model, climate)
  outCtrl=tot_energy
end subroutine

program simple_glide
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases
  use glimmer_global, only:rk, dp
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats_module

  use glide_diagnostics

  implicit none
 
  interface 
     subroutine numericalCoreWithModelInit(model, climate, config, inCtrl, outCtrl)
       use glide
       use simple_forcing
       use glimmer_config
       type(glide_global_type)      :: model 
       type(ConfigSection), pointer :: config
       type(simple_climate)         :: climate 
       real(dp) inCtrl, outCtrl
     end subroutine
  end interface

#ifdef GPTL
#include <gptl.inc>
#endif
#ifdef PAPI
#include <f90papi.h>
#endif

#ifdef HAVE_MPI
#include <mpif.h>

  integer npes,ierr,iam,lengthofconfigname
#endif

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret
  real(dp) inCtrl, outCtrl

  ! start gptl
#ifdef GPTL
  ret = gptlsetoption (gptlprint_method,gptlfull_tree)
  ret = gptlsetoption (PAPI_FP_OPS, 1)
  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('total')
#endif

#ifdef HAVE_MPI
 call MPI_Init(ierr)
 call MPI_Comm_size(MPI_COMM_WORLD,npes,ierr)
 call MPI_Comm_rank(MPI_COMM_WORLD,iam,ierr)

!! RN_20100115: An attempt to do MPI in Fortran
!if (iam .eq. 0) then
!  write(*,*), 'i am in here'
!  call glimmer_GetCommandline()
!  lengthofconfigname = len_trim(commandline_configname)
!  write(*,*), 'configname: ', commandline_configname
!  call MPI_Bcast(lengthofconfigname, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!  ! send commandline_configname
!  call MPI_Bcast(commandline_configname, lengthofconfigname, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
!else
!  write(*,*) "hello from ", iam
!! receive commandline_configname
!  call MPI_Bcast(lengthofconfigname, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!  call MPI_Bcast(commandline_configname, lengthofconfigname, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
!  write(*,*), 'configname: ', commandline_configname  
!endif  
!#else
!  call glimmer_GetCommandline()

#endif
  call glimmer_GetCommandline()



  ! start logging
  call open_log(unit=50, fname=logname(commandline_configname))
  
  ! read configuration
  call ConfigRead(commandline_configname,config)

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! initialise GLIDE
  inCtrl=0.0D0
  call numericalCoreWithModelInit(model, climate,config,inCtrl,outCtrl)
  print *,"DEPENDENT output", outCtrl
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_writestats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log

#ifdef HAVE_MPI
 call MPI_Finalize(ierr)
#endif

  ! stop gptl
#ifdef GPTL
  ret = gptlstop ('total')
  ret = gptlpr (0)
  ret = gptlfinalize ()
#endif
  
end program simple_glide
