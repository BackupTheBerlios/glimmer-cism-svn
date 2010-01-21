! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  shelf_driver.F90 - part of the GLIMMER ice model         + 
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

program shelf_driver
  !*FD This is a simple GLIDE driver. It can be used to run
  !*FD a Petermann shelf simulation
  use glimmer_global, only:rk
  use glide
  use shelf_forcing
 ! use glimmer_log
 ! use glimmer_config
  use glimmer_commandline
  use glimmer_writestats_module
  use plume
  use plume_io
  use plume_global, only:bpos

  use glimmer_paramets, only: thk0

  implicit none

  type(glide_global_type) :: model        ! model instance
  type(shelf_climate) :: climate_cfg      ! climate configuration info
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) :: time
  real(kind=dp) :: t1,t2
  integer :: clock,clock_rate
  character(len=512) :: plume_nl
  logical :: plume_suppress_output

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

  ! read [plume] section of glimmer config file
  call plume_read_print_config(model,config,plume_nl,plume_suppress_output)
!  model%options%which_bmlt = 1 !use plume
  model%options%which_bmlt = 0 !use plume
	 	
  ! config the [Peterman shelf] section of the config file
  call shelf_config_initialise(climate_cfg,config)

  !initialise glide
  call glide_initialise(model)

  ! this is a hack to get the lower surface into glide
  ! without changing the code that reads in initial data
  model%geometry%lsrf(:,:) = -1.0*920.0/1030.0 * model%geometry%thck(:,:) / thk0 
  model%geometry%topg = -2000.0/thk0

  call CheckSections(config)

  ! fill dimension variables
  call glide_nc_fillall(model)

  time = model%numerics%tstart

  call spinup_lithot(model)

!  call plume_io_initialize('/home/gladish/computation/glimmer-cism/test_run/', &
!                           '999',.false.)

!  call plume_initialise(plume_nl,model%geometry%lsrf * thk0)

 !run the plume to steady state with given shelf draft
 !call plume_iterate(time, &
 !                   model%numerics%tinc*1.d0, &
 !                   model%geometry%lsrf(:,:) * thk0, &
 !                   model%temper%temp(:,:,model%general%upn-1), &
 !                   (model%numerics%sigma(model%general%upn) - &
 !                   model%numerics%sigma(model%general%upn -1)) * &
!	             model%geometry%thck(:,:), &
!                    model%temper%bmlt(:,:), &
!                    model%temper%temp(:,:,model%general%upn), &            
!                    .true., &
!                    plume_stopping_tol = 1.d-10)

 write(*,*) model%temper%bmlt(10,:)


 do while(time.le.model%numerics%tend)

     write(*,*)'beginning of timestep at time =', time

     model%climate%acab(:,:) = climate_cfg%accumulation_rate
     model%climate%artm(:,:) = climate_cfg%artm

     write(*,*)'ice step'
     call glide_tstep_p1(model,time) ! temp evolution, basal traction
     call glide_tstep_p2(model)      ! velocities, thickness advection, basal stress
     call glide_tstep_p3(model)      ! isostasy, upper/lower surfaces

     time = time + model%numerics%tinc
     
!     write(*,*)'plume_iterate'
!     call plume_iterate(time, &
!                       model%numerics%tinc*1.d0, &
!                       model%geometry%lsrf(:,:) * thk0, &
!                       model%temper%temp(:,:,model%general%upn-1), &
!                      (model%numerics%sigma(model%general%upn) - &
!                       model%numerics%sigma(model%general%upn -1)) * &
!	               model%geometry%thck, &
!                       model%temper%bmlt, &
!                       model%temper%temp(:,:,model%general%upn), &                    
!                       run_plume_to_steady = .false., &
!                       plume_stopping_tol = 1.0d-10)

  end do

  ! finalise GLIDE
  call glide_finalise(model)
 ! call plume_finalise()	
 ! call plume_io_finalize()
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_writestats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log

contains

subroutine plume_read_print_config(model,config,plume_nl_file,plume_supp_output)

    !*FD read configuration file for plume-related things

    use glide_types
    use glimmer_config

    implicit none

    type(glide_global_type) :: model       
    type(ConfigSection), pointer :: config 
                                           
    character(len=*),intent(out) :: plume_nl_file
    logical,intent(out) :: plume_supp_output

    ! local variables
    type(ConfigSection), pointer :: section
    character(len=512) :: message

    call GetSection(config,section,'plume')
    if (associated(section)) then
        call GetValue(section, 'plume_nl_file', plume_nl_file)
	call GetValue(section, 'suppress_output', plume_supp_output)	
    else
	call write_log('for shelf driver there must be a [plume] &
                        section in config file',  &
	          GM_FATAL)
    endif 	          		

    call write_log('Plume config')
    write(message,*) 'plume namelist file:', trim(plume_nl_file)
    call write_log(message)
    write(message,*) 'suppressing plume output:',plume_supp_output
 	
end subroutine

end program shelf_driver
