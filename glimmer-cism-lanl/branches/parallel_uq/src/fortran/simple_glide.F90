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

program simple_glide
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases
  use glimmer_global, only:rk
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats_module

  use glide_diagnostics

  implicit none

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
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret

  integer :: tstep_count

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
  call glide_config(model,config)
  call simple_initialise(climate,config)
  call glide_initialise(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart

  tstep_count = 0

  do while(time.le.model%numerics%tend)
     call glide_tstep_p1(model,time)
     call glide_tstep_p2(model)
     call glide_tstep_p3(model)
     ! override masking stuff for now

!     tstep_count = tstep_count + 1
!     if (mod(tstep_count, model%numerics%ndiag) == 0) then
!        call glide_write_diag(model, time, model%numerics%idiag, &
!                                           model%numerics%jdiag )
!     endif

     time = time + model%numerics%tinc

  end do

  ! finalise GLIDE
  call glide_finalise(model)
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
