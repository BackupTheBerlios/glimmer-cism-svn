#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

program shelf_driver

  use glide
  use glide_mask
  use glimmer_commandline
  use glimmer_writestats_module
  use plume
  use plume_io

  use glimmer_paramets, only: thk0
  use glide_temp, only: timeevoltemp
  use glam_strs2, only: use_shelf_bc_1
  use glimmer_scales, only: scale2d_f1

  implicit none

  type shelf_climate
     ! holds parameters for the shelf climate
     real(kind=sp) :: artm = 0.0
     real(kind=sp) :: accumulation_rate = 0.0
     real(kind=sp) :: eus = 0.0
  end type shelf_climate

  integer,parameter :: USE_PLUME = 1
  integer,parameter :: fake_landw = 2

  !local variables
  type(glide_global_type) :: model        ! model instance
  type(shelf_climate) :: climate_cfg      ! climate configuration info
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) :: time
  real(kind=dp) :: t1,t2
  integer :: clock,clock_rate
  logical :: is_steady = .false.

  logical :: check_for_steady = .false.
  real(kind=dp) :: thk_steady_tol = 1.0e-5

  real(kind=rk),dimension(:),allocatable :: upstream_thck
  logical,      dimension(:,:),allocatable :: plume_land_mask
  real(kind=dp),dimension(:,:),allocatable :: plume_lsrf_ext,plume_ice_dz,plume_t_interior
  real(kind=dp),dimension(:,:),allocatable :: plume_bmelt_out, plume_btemp_out
  real(kind=dp),dimension(:,:),allocatable :: prev_ice_thk

  !PLUME configuration values
  character(len=512) :: plume_nl,plume_output_nc_file,plume_output_prefix,plume_ascii_output_dir
  logical :: plume_suppress_ascii_output,plume_suppress_logging
  logical :: plume_write_all_states= .false.
  real(kind=dp) :: plume_min_subcycle_time,plume_min_spinup_time,plume_steadiness_tol
  integer :: plume_imin,plume_imax,plume_kmin,plume_kmax
  logical :: plume_const_bmlt
  real(kind=dp) :: plume_const_bmlt_rate

  use_shelf_bc_1 = .true.

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

  if (model%options%use_plume == USE_PLUME) then
     ! read [plume] section of glimmer config file
     call plume_read_print_config(model,config) 

  end if

  ! config the [Peterman shelf] section of the config file
  call shelf_config_initialise(climate_cfg,config)

  call CheckSections(config)

  !initialise sea level here so that ice mask get set correctly (ie floating)
  model%climate%eus = climate_cfg%eus 

  !initialise glide
  call glide_initialise(model)

  !initialise climate in model
  !these have to be done after since acab and artm are allocated space
  !inside glide_initialise
  model%climate%acab(:,:) = climate_cfg%accumulation_rate
  model%climate%artm(:,:) = climate_cfg%artm

  allocate(upstream_thck(model%general%ewn))
  upstream_thck = model%geometry%thck(:,model%general%nsn)

  ! fill dimension variables
  call glide_nc_fillall(model)

  time = model%numerics%tstart

  if ((model%options%use_plume == USE_PLUME) .and. .not. plume_const_bmlt) then
     ! we are using the plume

     allocate(plume_land_mask(  model%general%ewn+2*fake_landw, &
          model%general%nsn+1*fake_landw))
     allocate(plume_lsrf_ext(   model%general%ewn+2*fake_landw, &
          model%general%nsn+1*fake_landw))
     allocate(plume_ice_dz(     model%general%ewn+2*fake_landw, &
          model%general%nsn+1*fake_landw))
     allocate(plume_t_interior( model%general%ewn+2*fake_landw, &
          model%general%nsn+1*fake_landw)) 
     allocate(plume_bmelt_out(  model%general%ewn+2*fake_landw, &
          model%general%nsn+1*fake_landw))
     allocate(plume_btemp_out(  model%general%ewn+2*fake_landw, &
          model%general%nsn+1*fake_landw))

     call write_logical_plume_array(GLIDE_IS_LAND(model%geometry%thkmask) .or. &
          GLIDE_IS_GROUND(model%geometry%thkmask), &       
          plume_land_mask, .true., &
          model%general%ewn, model%general%nsn, fake_landw)
     call write_real_plume_array(model%geometry%lsrf *thk0, plume_lsrf_ext, 0.0, &
          model%general%ewn, model%general%nsn, fake_landw)
     call write_real_plume_array( model%temper%temp(model%general%upn-1,:,:), &
          plume_t_interior, 0.0, &
          model%general%ewn, model%general%nsn, fake_landw)
     call write_real_plume_array((model%numerics%sigma(model%general%upn) - &
          model%numerics%sigma(model%general%upn -1)) * &
          model%geometry%thck * thk0, &
          plume_ice_dz, 0.0, &
          model%general%ewn, model%general%nsn, fake_landw)   

     call plume_logging_initialize(trim(plume_ascii_output_dir), &
          trim(plume_output_prefix), &
          plume_suppress_logging)

     call plume_initialise(trim(plume_nl), &
          plume_suppress_ascii_output, &
          plume_suppress_logging,&
          plume_lsrf_ext,&
          plume_land_mask, &
          plume_imin,plume_imax,plume_kmin,plume_kmax)

     call plume_io_initialize(plume_suppress_ascii_output, &
          trim(plume_output_nc_file))

     !run the plume to steady state with given shelf draft
     call plume_iterate(time, &
          model%numerics%tinc*1.d0, &
          plume_lsrf_ext, &
          plume_land_mask, &
          plume_t_interior, &
          plume_ice_dz, &
          plume_bmelt_out, &
          plume_btemp_out, &
          plume_steadiness_tol, &
          plume_min_spinup_time, &
          .true.,&
          plume_write_all_states)


     call write_real_ice_array(plume_bmelt_out / scale2d_f1,model%temper%bmlt, &
          model%general%ewn, model%general%nsn, fake_landw)
     call write_real_ice_array(plume_btemp_out,model%temper%temp(model%general%upn,:,:), &
          model%general%ewn, model%general%nsn, fake_landw)

  end if

  if (plume_const_bmlt) then
     ! NB: a positive plume_const_bmlt_rate should indicate
     !     mass being lost from the ice shelf
     model%temper%bmlt = plume_const_bmlt_rate / scale2d_f1
  end if

  if (plume_write_all_states) then
     write(*,*) 'Stopping because plume_write_all_states is true'
     call plume_finalise()
     call plume_io_finalize()
     call plume_logging_finalize()

     stop
  end if

  allocate(prev_ice_thk(model%general%ewn,model%general%nsn))
  prev_ice_thk = 0.0

  do while(time .le. model%numerics%tend  .and. &
       .not. (check_for_steady .and. is_steady))

     model%climate%acab(:,:) = climate_cfg%accumulation_rate
     model%climate%artm(:,:) = climate_cfg%artm
     model%climate%eus = climate_cfg%eus

     call glide_tstep_p1(model,time) ! temp evolution
     call glide_tstep_p2(model)      ! velocities, thickness advection
     ! adjust the heights in the upstream row
     model%geometry%thck(:,model%general%nsn) = upstream_thck
     model%geometry%thck(:,model%general%nsn-1) = upstream_thck

     ! impose dh/dx = 0 along the lateral side walls
     model%geometry%thck(model%general%ewn,:) = model%geometry%thck(model%general%ewn-1,:)
     model%geometry%thck(                1,:) = model%geometry%thck(                  2,:)

     call glide_tstep_p3(model)      ! isostasy, upper/lower surfaces

     is_steady =  check_for_steady .and. &
          all(  abs(prev_ice_thk - model%geometry%thck) .le. &
          abs(model%geometry%thck)*thk_steady_tol  )
     prev_ice_thk = model%geometry%thck

     time = time + model%numerics%tinc

     if ((model%options%use_plume == USE_PLUME) .and.  (.not. plume_const_bmlt)) then	

        ! We would normally expect the plume to reach a steady state
        ! with respect to the current ice after only a few days.
        ! If necessary we run the plume over the entire ice timestep
        ! but if the basal melt rate becomes constant to the specified
        ! tolerance than we assume the plume melts the ice at the
        ! steady rate for the rest of the ice timestep and we stop
        ! timestepping the plume.
        call write_logical_plume_array(GLIDE_IS_LAND(model%geometry%thkmask) .or. &
             GLIDE_IS_GROUND(model%geometry%thkmask), &       
             plume_land_mask, .true.,&
             model%general%ewn, model%general%nsn, fake_landw)
        call write_real_plume_array(model%geometry%lsrf *thk0, &
             plume_lsrf_ext, 0.0, &                              
             model%general%ewn, model%general%nsn, fake_landw)
        call write_real_plume_array( model%temper%temp(model%general%upn-1,:,:), &
             plume_t_interior, 0.0, &
             model%general%ewn, model%general%nsn, fake_landw)

        call write_real_plume_array((model%numerics%sigma(model%general%upn) - &
             model%numerics%sigma(model%general%upn -1)) * &
             model%geometry%thck * thk0, &
             plume_ice_dz, 0.0, &
             model%general%ewn, model%general%nsn, fake_landw)

        call plume_iterate(time, &
             model%numerics%tinc*1.d0, &
             plume_lsrf_ext, &
             plume_land_mask, &                      
             plume_t_interior, &
             plume_ice_dz, &
             plume_bmelt_out, &
             plume_btemp_out, &            
             plume_steadiness_tol, &
             plume_min_subcycle_time, &
             .false., &
             plume_write_all_states)

        call write_real_ice_array(plume_bmelt_out / scale2d_f1,model%temper%bmlt, &
             model%general%ewn, model%general%nsn, fake_landw)
        call write_real_ice_array(plume_btemp_out,model%temper%temp(model%general%upn,:,:), &
             model%general%ewn, model%general%nsn, fake_landw)
     end if

     if (plume_const_bmlt) then
        model%temper%bmlt = plume_const_bmlt_rate / scale2d_f1
     end if

  end do

  if (is_steady) then
     call write_log('Stopped time-stepping after detecting steady ice profile')
  end if

  if ((model%options%use_plume .eq. USE_PLUME) .and. (.not. plume_const_bmlt)) then 

     call plume_finalise()	
     call plume_io_finalize()
     call plume_logging_finalize()

     deallocate(plume_land_mask,plume_lsrf_ext,plume_ice_dz, &
          plume_t_interior, plume_bmelt_out, plume_btemp_out)

  end if

  ! finalise GLIDE
  deallocate(upstream_thck, prev_ice_thk)
  call glide_finalise(model)

  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_writestats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log()

contains

  subroutine write_real_plume_array(a_in,a_out,padding_val, m,n,padw)

    real(kind=dp),dimension(:,:),intent(in) :: a_in
    real(kind=dp),dimension(:,:),intent(out):: a_out
    real,                        intent(in) :: padding_val
    integer,                     intent(in) :: m,n,padw

    a_out = padding_val
    a_out( (1+padw):(m+padw), 1:n ) = a_in

  end subroutine write_real_plume_array

  subroutine write_logical_plume_array(a_in,a_out,padding_val, m,n,padw)

    logical,dimension(:,:),intent(in) :: a_in
    logical,dimension(:,:),intent(out):: a_out
    logical,               intent(in) :: padding_val
    integer,               intent(in) :: m,n,padw

    a_out = padding_val
    a_out( (1+padw) : (m+padw), 1:n ) = a_in

  end subroutine write_logical_plume_array

  subroutine write_real_ice_array(a_in,a_out, m,n,padw)

    real(kind=dp),dimension(:,:),intent(in) :: a_in
    real(kind=dp),dimension(:,:),intent(out):: a_out
    integer,                     intent(in) :: m,n,padw

    a_out = a_in( (1+padw):(m+padw), : )

  end subroutine write_real_ice_array

  subroutine shelf_config_initialise(climate_cfg,config)
    ! initialise shelf climate model

    use glimmer_paramets, only: thk0, acc0, scyr
    use glimmer_config
    use glimmer_log

    type(shelf_climate) :: climate_cfg         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   

    ! local variables
    type(ConfigSection), pointer :: section
    character(len=100) :: message

    call GetSection(config,section,'Petermann shelf')
    if (associated(section)) then
       call GetValue(section,'air_temperature',climate_cfg%artm)
       call GetValue(section,'accumulation_rate',climate_cfg%accumulation_rate)   
       call GetValue(section,'eustatic_sea_level',climate_cfg%eus)
       call GetValue(section,'check_for_steady',check_for_steady)
       call GetValue(section,'thk_steady_tol',thk_steady_tol)
    else
       !log error
       call write_log('No Petermann section',GM_FATAL)
    end if

    call write_log_div
    call write_log('Petermann shelf configuration')
    call write_log('------------------------------------')
    write(message,*) 'air temperature  : ',climate_cfg%artm
    call write_log(message)
    write(message,*) 'accumulation rate: ',climate_cfg%accumulation_rate
    call write_log(message)
    write(message,*) 'eustatic sea level:',climate_cfg%eus
    call write_log(message)

  end subroutine shelf_config_initialise

  subroutine plume_read_print_config(model,config) 

    !*FD read configuration file for plume-related things

    use glide_types
    use glimmer_config

    type(glide_global_type) :: model       
    type(ConfigSection), pointer :: config                                            

    ! local variables
    type(ConfigSection), pointer :: section
    character(len=512) :: message

    call GetSection(config,section,'plume')
    if (.not. associated(section)) then
       call write_log('for shelf driver there must be a [plume] section &
            & in config file',  &
            GM_FATAL)
    end if
    call GetValue(section, 'plume_nl_file', plume_nl)
    call GetValue(section, 'plume_output_dir', plume_ascii_output_dir)
    call GetValue(section, 'plume_output_file', plume_output_nc_file)
    call GetValue(section, 'suppress_ascii_output', plume_suppress_ascii_output)
    call GetValue(section, 'suppress_logging', plume_suppress_logging)
    call GetValue(section, 'plume_output_prefix',plume_output_prefix)
    call GetValue(section, 'plume_write_all_states',plume_write_all_states)
    call GetValue(section, 'plume_min_spinup_time',plume_min_spinup_time)
    call GetValue(section, 'plume_min_subcycle_time',plume_min_subcycle_time)
    call GetValue(section, 'plume_steadiness_tol', plume_steadiness_tol)
    call GetValue(section, 'plume_imin',plume_imin)
    call GetValue(section, 'plume_imax',plume_imax)
    call GetValue(section, 'plume_kmin',plume_kmin)
    call GetValue(section, 'plume_kmax',plume_kmax)
    call GetValue(section, 'plume_const_bmlt', plume_const_bmlt)
    call GetValue(section, 'plume_const_bmlt_rate', plume_const_bmlt_rate)

    call write_log('Plume config')

    write(message,*) 'plume namelist file:', trim(plume_nl)
    call write_log(message)
    write(message,*) 'suppressing plume output:',plume_suppress_ascii_output
    call write_log(message)
    write(message,*) 'suppressing plume screen/file logging:',plume_suppress_logging
    call write_log(message)
    write(message,*) 'plume output written to netcdf file:', trim(plume_output_nc_file)
    call write_log(message)
    write(message,*) 'plume old-style output written to directory:', trim(plume_ascii_output_dir)
    call write_log(message)
    write(message,*) 'plume output file prefix (jobid):', trim(plume_output_prefix)
    call write_log(message)
    write(message,*) 'plume_write_all_states:', plume_write_all_states
    call write_log(message)
    write(message,*) 'plume_min_spinup_time',plume_min_spinup_time
    call write_log(message)
    write(message,*) 'plume_min_subcycle_time', plume_min_subcycle_time
    call write_log(message)
    write(message,*) 'imin',plume_imin,'imax',plume_imax,'kmin',plume_kmin,'kmax',plume_kmax
    call write_log(message)
    write(message,*) 'plume_const_bmlt', plume_const_bmlt, 'plume_const_bmlt_rate', plume_const_bmlt_rate
    call write_log(message)

  end subroutine plume_read_print_config


end program shelf_driver
