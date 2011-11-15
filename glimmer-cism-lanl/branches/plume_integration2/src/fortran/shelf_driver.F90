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
  integer,parameter :: thk_zero_margin = 3
  logical,parameter :: use_thk_zero_margin = .true.

  real(kind=kind(1.d0)),parameter :: min_melt_depth = -50.d0
!  real(kind=kind(1.d0)),parameter :: min_melt_depth = 0.d0

  !local variables
  type(glide_global_type) :: model        ! model instance
  type(shelf_climate) :: climate_cfg      ! climate configuration info
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) :: time
  real(kind=dp) :: t1,t2
  integer :: clock,clock_rate
  logical :: is_steady = .false.
  logical :: do_plume_coupling = .true.
  logical :: plume_reached_steady

  logical :: check_for_steady = .false.
  logical :: doCrossShelfAvg = .false.
  logical :: hide_shelf_inflow_row = .true.

  real(kind=dp) :: thk_steady_tol = 1.0e-2
  real(kind=dp) :: mean_thk_steady_tol = 1.0e-6
  real(kind=dp) :: mean_rel_thk_change = 0.d0
  real(kind=dp) :: max_rel_thk_change = 0.d0	

  real(kind=dp) :: last_bmlt_time = 0.d0
  real(kind=dp),parameter :: big_tinc = 0.0d0

  ! smooth_diff is the diffusivity of thickness in x direction, in (m^2/a)
  real(kind=dp)  :: smooth_diff = 0.d0 
  real(kind=dp)  :: non_dim_smooth_diff

  integer :: i_smooth_col 
  integer :: ice_ntimestep = 0	

  real(kind=rk),dimension(:,:),allocatable :: upstream_thck
  real(kind=rk),dimension(:,:),allocatable :: smoothed_thck
  logical,      dimension(:,:),allocatable :: plume_land_mask,no_plume
  real(kind=dp),dimension(:,:),allocatable :: plume_lsrf_ext,plume_ice_dz,plume_t_interior
  real(kind=dp),dimension(:,:),allocatable :: plume_bmelt_out, plume_btemp_out
  real(kind=dp),dimension(:,:),allocatable :: prev_ice_thk

  !PLUME configuration values
  character(len=1048) :: plume_nl,plume_output_nc_file,plume_output_prefix,plume_ascii_output_dir
  logical :: plume_suppress_ascii_output,plume_suppress_logging
  logical :: plume_write_all_states= .false.
  integer :: plume_write_every_n = 1
  real(kind=dp) :: plume_output_frequency = 0.d0

  real(kind=dp) :: plume_min_subcycle_time,plume_max_subcycle_time
  real(kind=dp) :: plume_min_spinup_time,plume_max_spinup_time
  real(kind=dp) :: plume_steadiness_tol,plume_speed_steadiness_tol

  integer :: plume_imin,plume_imax,plume_kmin,plume_kmax
  logical :: plume_const_bmlt,plume_initial_bmlt, plume_delayed_coupling
  real(kind=dp) :: plume_homotopy_frac = 1.d0
  real(kind=dp) :: plume_homotopy_ramp = 0.001d0

  real(kind=dp) :: plume_const_bmlt_rate

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

  model%climate%acab(:,:) = climate_cfg%accumulation_rate / scale2d_f1
  model%climate%artm(:,:) = climate_cfg%artm

  !initialise ice temperature to air temperature
  call timeevoltemp(model,0)

  allocate(upstream_thck(model%general%ewn-2*thk_zero_margin,1))
  allocate(smoothed_thck(model%general%ewn,model%general%nsn))

  upstream_thck(:,1) = &
       model%geometry%thck(thk_zero_margin+1:model%general%ewn-thk_zero_margin, &
       	                   1)  !south edge is inflow
  smoothed_thck = 0.d0

  ! fill dimension variables
  call glide_nc_fillall(model)

  time = model%numerics%tstart

  if ((model%options%use_plume == USE_PLUME) .and. .not. plume_const_bmlt) then
     ! we are using the plume

     ! these are plume grid fields
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

     ! no_plume is an ice-grid variable
     allocate(no_plume(model%general%ewn, &
                       model%general%nsn))

     no_plume = GLIDE_IS_LAND(model%geometry%thkmask) .or. &
                GLIDE_IS_GROUND(model%geometry%thkmask)

     if (use_thk_zero_margin) then
     	no_plume(1:thk_zero_margin,:) = .true.
	no_plume(model%general%ewn-(thk_zero_margin-1):model%general%ewn,:) = .true.
     end if

     if (hide_shelf_inflow_row) then
        no_plume(:,1) = .true.    !hide southern row of ice
     end if

     call write_logical_plume_array(no_plume,plume_land_mask, .true., &
          model%general%ewn, model%general%nsn, fake_landw,fake_landw)

     call write_real_plume_array(model%geometry%lsrf*thk0, plume_lsrf_ext, 0.d0, &
          model%general%ewn, model%general%nsn, fake_landw, fake_landw)

     !must be explicit about model%temper%temp array bounds, since it has extra cells
     call write_real_plume_array( model%temper%temp(model%general%upn-1, &
                                                    1:model%general%ewn, &
                                                    1:model%general%nsn), &
          plume_t_interior, 0.d0, &
          model%general%ewn, model%general%nsn, fake_landw,fake_landw)

     call write_real_plume_array((model%numerics%sigma(model%general%upn) - &
                                  model%numerics%sigma(model%general%upn -1)) * &
                                  model%geometry%thck * thk0, &
                                  plume_ice_dz, 0.d0, &
                                  model%general%ewn, model%general%nsn, fake_landw,fake_landw)   

     call plume_logging_initialize(trim(plume_ascii_output_dir), &
          trim(plume_output_prefix), &
          plume_suppress_logging)

     if (doCrossShelfAvg) then
	  call cross_shelf_average(plume_lsrf_ext)
     end if

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
          plume_speed_steadiness_tol, &
          plume_min_spinup_time, &
	  plume_max_spinup_time, &
	  3600.d0*24, &
          .true.,&
	  plume_reached_steady, &
          plume_write_all_states, &
          plume_write_every_n, &
	  0.d0, &	
	  plume_initial_bmlt)

     if (.not. plume_reached_steady .and. .not. plume_initial_bmlt) then
	call write_log("Plume did not reach a steady state",GM_WARNING)
     end if

     call homotopy_bmlt(plume_bmelt_out,time,model%numerics%tstart*1.d0,model%geometry%lsrf,min_melt_depth)

     call write_real_ice_array(plume_bmelt_out / scale2d_f1,model%temper%bmlt, &
          model%general%ewn, model%general%nsn, fake_landw,fake_landw)
     
     call write_real_ice_array(plume_btemp_out,model%temper%temp(model%general%upn, &
                                                                 1:model%general%ewn, &
                                                                 1:model%general%nsn), &
          model%general%ewn, model%general%nsn, fake_landw,fake_landw)

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
  prev_ice_thk = 0.d0

  model%geometry%thck(1:thk_zero_margin ,:) = 0.0_dp
  model%geometry%thck(model%general%ewn-(thk_zero_margin-1):model%general%ewn,:) = 0.0_dp

  do while(time .le. model%numerics%tend  .and. &
       .not. (check_for_steady .and. is_steady))

     model%climate%acab(:,:) = climate_cfg%accumulation_rate / scale2d_f1
     model%climate%artm(:,:) = climate_cfg%artm
     model%climate%eus = climate_cfg%eus

     print *, 'beginning tstep_p1'
     call glide_tstep_p1(model,time) ! temp evolution

     print *, 'beginning tstep_p2'
     call glide_tstep_p2(model)      ! velocities, thickness advection

     print *, 'zeroing out marginal thickness values'
     ! adjust the heights in the upstream row
     model%geometry%thck(thk_zero_margin+1:model%general%ewn-thk_zero_margin, &
                         1  ) = upstream_thck(:,1) + inflow_thk_perturb(time)

     model%geometry%thck(thk_zero_margin+1:model%general%ewn-thk_zero_margin, &
                         2  ) = upstream_thck(:,1) + inflow_thk_perturb(time)

     model%geometry%thck(thk_zero_margin+1:model%general%ewn-thk_zero_margin, &
                         3  ) = upstream_thck(:,1) + inflow_thk_perturb(time)
    
     ! and therefore zero out the thck_t values
     model%geometry%thck_t(thk_zero_margin+1:model%general%ewn-thk_zero_margin, &
                           1  ) = 0.d0
     model%geometry%thck_t(thk_zero_margin+1:model%general%ewn-thk_zero_margin, &
                           2  ) = 0.d0
     model%geometry%thck_t(thk_zero_margin+1:model%general%ewn-thk_zero_margin, &
                           3  ) = 0.d0

     ! set thickness in last row equal to the next upstream row.  This is
     ! something like a calving condition, to keep the thickness from growing 
     ! for reasons I don't understand **cvg***

     model%geometry%thck(:,model%general%nsn-4) = model%geometry%thck(:,model%general%nsn-5)
     model%geometry%thck_t(:,model%general%nsn-4) = 0.d0

     ! impose dh/dx = 0 along the lateral side walls

     ! This has been commented out.  I believe it is physically incorrect since it can act 
     ! as a thickness source along the sides, which is not the boundary conditions we
     ! want when simulating a steep-walled fjord like Petermann.
!     model%geometry%thck(model%general%ewn,:) = model%geometry%thck(model%general%ewn-1,:)
!     model%geometry%thck(                1,:) = model%geometry%thck(                  2,:)

     !apply x-direction smoothing
     non_dim_smooth_diff = &
            smooth_diff*model%numerics%tinc/(model%numerics%dew**2.d0) &
	    * (1.d0/len0**2.d0)

     print *, 'non_dim_smooth_diff = ',non_dim_smooth_diff
     if (non_dim_smooth_diff > 0.5d0) then
	print *, "Diffusion violates Courant condition"
	stop 1
     end if

     smoothed_thck(thk_zero_margin+1,4:) = &
          (1.d0-2.d0*non_dim_smooth_diff)*model%geometry%thck(thk_zero_margin+1,4:) &
               +2.d0*non_dim_smooth_diff*model%geometry%thck(thk_zero_margin+2,4:)

     smoothed_thck(model%general%ewn-thk_zero_margin,4:) = &
          (1.d0-2.d0*non_dim_smooth_diff)*model%geometry%thck(model%general%ewn - &
                                                      thk_zero_margin, 4:) &
               +2.d0*non_dim_smooth_diff*model%geometry%thck(model%general%ewn - &
                                                     thk_zero_margin-1, 4:)

     do i_smooth_col=thk_zero_margin+2,model%general%ewn-thk_zero_margin-1
        smoothed_thck(i_smooth_col,4:) = &
             (1.d0 - 2.d0*non_dim_smooth_diff)*model%geometry%thck(i_smooth_col,4:) &
	           + non_dim_smooth_diff*model%geometry%thck(i_smooth_col+1,4:) &
                   +  non_dim_smooth_diff*model%geometry%thck(i_smooth_col-1,4:)
     end do
     model%velocity_hom%H_diff_t(:,4:) = &
           (1.d0/model%numerics%tinc)*(smoothed_thck(:,4:)-model%geometry%thck(:,4:))
     model%geometry%thck(:,4:) = smoothed_thck(:,4:)


     print *, 'beginning tstep_p3'
     call glide_tstep_p3(model)      ! isostasy, upper/lower surfaces


     mean_rel_thk_change = abs(sum(model%geometry%thck_t)/ &
		         (model%general%nsn*model%general%ewn))
     max_rel_thk_change =  maxval(abs(model%geometry%thck_t))

     is_steady =  check_for_steady .and. &
	  ( max_rel_thk_change <      thk_steady_tol) .and. &
          (mean_rel_thk_change < mean_thk_steady_tol)

     time = time + get_tinc(model)

     ice_ntimestep = ice_ntimestep + 1
     print '(a,e8.2,a,i12,a)', &
           'Completed ice timestep at t == ',time, &
	   '(step ', ice_ntimestep ,' )'

     if ((model%options%use_plume == USE_PLUME) .and.  (.not. plume_const_bmlt)) then	

        ! We would normally expect the plume to reach a steady state
        ! with respect to the current ice after only a few days.
        ! If necessary we run the plume over the entire ice timestep
        ! but if the basal melt rate becomes constant to the specified
        ! tolerance than we assume the plume melts the ice at the
        ! steady rate for the rest of the ice timestep and we stop
        ! timestepping the plume.

	no_plume = GLIDE_IS_LAND(model%geometry%thkmask) .or. &
		  GLIDE_IS_GROUND(model%geometry%thkmask)

        if (use_thk_zero_margin) then
	        ! The purpose of this is to force the plume to disregard
		! the laterally outermost column of ocean cells, where the
		! ice thickness has been set to zero in order to make the 
		! ice boundary conditions correct.  Thus the plume only inhabits
		! ocean cells that are beneath ice with positive thickness.
        	no_plume(1:thk_zero_margin,:) = .true.
        	no_plume(model%general%ewn-(thk_zero_margin-1):model%general%ewn,:) = .true.
        end if

	if (hide_shelf_inflow_row) then
	       no_plume(:, 1) = .true.
	end if

        call write_logical_plume_array(   no_plume, plume_land_mask, .true.,&
             model%general%ewn, model%general%nsn, fake_landw,fake_landw)

        call write_real_plume_array(model%geometry%lsrf *thk0, plume_lsrf_ext, 0.d0, &                              
             model%general%ewn, model%general%nsn, fake_landw,fake_landw)
        call write_real_plume_array( model%temper%temp(model%general%upn-1, &
                                                       1:model%general%ewn, &
                                                       1:model%general%nsn), &
             plume_t_interior, 0.d0, &
             model%general%ewn, model%general%nsn, fake_landw,fake_landw)

        call write_real_plume_array((model%numerics%sigma(model%general%upn) - &
                                     model%numerics%sigma(model%general%upn -1)) * &
                                     model%geometry%thck * thk0, &
             plume_ice_dz, 0.d0, &
             model%general%ewn, model%general%nsn, fake_landw,fake_landw)

        if (doCrossShelfAvg) then
             call cross_shelf_average(plume_lsrf_ext)
        end if

	if (time > last_bmlt_time + big_tinc) then

	last_bmlt_time = time

	if (plume_delayed_coupling) then
	    if (is_steady) then
		is_steady = .false.
	        plume_initial_bmlt = .false.
	        do_plume_coupling = .true.
            end if
        end if

        call plume_iterate(time, &
             model%numerics%tinc*1.d0, &
             plume_lsrf_ext, &
             plume_land_mask, &                      
             plume_t_interior, &
             plume_ice_dz, &
             plume_bmelt_out, &
             plume_btemp_out, &            
             plume_steadiness_tol, &
             plume_speed_steadiness_tol, &
             plume_min_subcycle_time, &
	     plume_max_subcycle_time, &
	     3600.d0*24, &
             .false., &                   !not necessarily running to steady
	     plume_reached_steady, &
	     plume_write_all_states, &
             plume_write_every_n, &
	     plume_output_frequency, &
	     plume_initial_bmlt .or. .not.(do_plume_coupling) )

	call homotopy_bmlt(plume_bmelt_out,time,model%numerics%tstart*1.d0,model%geometry%lsrf,min_melt_depth)

        call write_real_ice_array(plume_bmelt_out / scale2d_f1,model%temper%bmlt, &
             model%general%ewn, model%general%nsn, fake_landw,fake_landw)

!	where (model%geometry%lsrf > min_melt_depth/thk0)
!             model%temper%bmlt = 0.d0
!        end where

        call write_real_ice_array(plume_btemp_out,model%temper%temp(model%general%upn, &
	                                                            1:model%general%ewn, &
							            1:model%general%nsn), &
	                          model%general%ewn, model%general%nsn, fake_landw,fake_landw)
	     
	end if
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

     deallocate(plume_land_mask,no_plume,plume_lsrf_ext,plume_ice_dz, &
          plume_t_interior, plume_bmelt_out, plume_btemp_out)

  end if

  ! finalise GLIDE
  deallocate(upstream_thck, prev_ice_thk)
  deallocate(smoothed_thck)

  call glide_finalise(model)

  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_writestats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log()

contains

  function inflow_thk_perturb(time)

	! return a scaled time-varying perturbation that is to be applied
        ! to all ice entering the domain across the inflow line

	real(kind=dp),intent(in) :: time

	real(kind=dp) :: inflow_thk_perturb
	real(kind=dp),parameter :: amplitude = 0.d0 !50.d0

	inflow_thk_perturb = amplitude * sin(2.d0*3.1415926d0*time)/ thk0

  end function inflow_thk_perturb

  subroutine homotopy_bmlt(bmelt_field, time, starttime,lsrf,min_melt_depth)

	real(kind=dp),dimension(:,:),intent(inout),target :: bmelt_field
	real(kind=dp),intent(in) :: time, starttime
 	real(kind=dp),intent(in) :: min_melt_depth
	real(kind=dp),dimension(:,:),intent(in) :: lsrf

	real(kind=dp),pointer,dimension(:,:) :: ice_bmlt
	character(len=128) :: msg 
	real(kind=dp) :: current_min_depth

	if (time .ge. starttime) then
		bmelt_field = bmelt_field * min(1.d0, (plume_homotopy_frac + \
        	                                (1.d0-plume_homotopy_frac) * \
                                        (time-starttime)/plume_homotopy_ramp))
	else
		print *, "Expected time > starttime"
	        stop 1
	end if
 
	current_min_depth = max(0.d0, &
                                (1.0-(time-starttime)/plume_homotopy_ramp))*&
	                    min_melt_depth/thk0


	ice_bmlt => bmelt_field(2+1:2+size(lsrf,1),2+1:2+size(lsrf,2))

        where (lsrf > current_min_depth)
             ice_bmlt = 0.d0
	end where

	write( msg,* ) 'min melt depth ', current_min_depth*thk0
	call io_append_output(trim(msg))

  end subroutine homotopy_bmlt

  subroutine write_real_plume_array(a_in,a_out,padding_val, in_m,in_n,padw,pads)

    real(kind=dp),dimension(:,:),intent(in) :: a_in
    real(kind=dp),dimension(:,:),intent(out):: a_out
    real(kind=dp),               intent(in) :: padding_val
    integer,                     intent(in) :: in_m,in_n,padw,pads

    a_out = padding_val
    a_out( (1+padw):(in_m+padw), (1+pads):(in_n+pads) ) = a_in

  end subroutine write_real_plume_array

  subroutine write_logical_plume_array(a_in,a_out,padding_val, in_m,in_n,padw,pads)

    logical,dimension(:,:),intent(in) :: a_in
    logical,dimension(:,:),intent(out):: a_out
    logical,               intent(in) :: padding_val
    integer,               intent(in) :: in_m,in_n,padw,pads

    a_out = padding_val
    a_out( (1+padw) : (in_m+padw), (1+pads):(in_n+pads) ) = a_in

  end subroutine write_logical_plume_array

  subroutine write_real_ice_array(a_in,a_out, in_m,in_n,padw,pads)

    real(kind=dp),dimension(:,:),intent(in) :: a_in
    real(kind=dp),dimension(:,:),intent(out):: a_out
    integer,                     intent(in) :: in_m,in_n,padw,pads

    a_out = a_in( (1+padw):(in_m+padw), (1+pads):(in_n+pads) )

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
       call GetValue(section,'mean_thk_steady_tol',mean_thk_steady_tol)
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
    call GetValue(section, 'plume_write_every_n', plume_write_every_n)
    call GetValue(section, 'plume_output_frequency', plume_output_frequency)
    call GetValue(section, 'plume_min_spinup_time',plume_min_spinup_time)
    call GetValue(section, 'plume_max_spinup_time',plume_max_spinup_time)
    call GetValue(section, 'plume_min_subcycle_time',plume_min_subcycle_time)
    call GetValue(section, 'plume_max_subcycle_time',plume_max_subcycle_time)
    call GetValue(section, 'plume_steadiness_tol', plume_steadiness_tol)
    call GetValue(section, 'plume_speed_steadiness_tol', plume_speed_steadiness_tol)
    call GetValue(section, 'plume_imin',plume_imin)
    call GetValue(section, 'plume_imax',plume_imax)
    call GetValue(section, 'plume_kmin',plume_kmin)
    call GetValue(section, 'plume_kmax',plume_kmax)
    call GetValue(section, 'plume_const_bmlt', plume_const_bmlt)
    call GetValue(section, 'plume_initial_bmlt', plume_initial_bmlt)
    call GetValue(section, 'plume_delayed_coupling', plume_delayed_coupling)
    call GetValue(section, 'plume_homotopy_frac', plume_homotopy_frac)
    call GetValue(section, 'plume_homotopy_ramp', plume_homotopy_ramp)
    call GetValue(section, 'plume_smooth_diff', smooth_diff)
    call GetValue(section, 'plume_const_bmlt_rate', plume_const_bmlt_rate)
    call GetValue(section, 'plume_do_cross_shelf_avg', doCrossShelfAvg)

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
    write(message,*) 'plume_write_every_n:', plume_write_every_n
    call write_log(message)
    write(message,*) 'plume_output_frequency:', plume_output_frequency
    call write_log(message)
    write(message,*) 'plume_min_spinup_time:',plume_min_spinup_time
    call write_log(message)
    write(message,*) 'plume_max_spinup_time:',plume_max_spinup_time
    call write_log(message)
    write(message,*) 'plume_min_subcycle_time', plume_min_subcycle_time
    call write_log(message)
    write(message,*) 'plume_max_subcycle_time', plume_max_subcycle_time
    call write_log(message)
    write(message,*) 'imin',plume_imin,'imax',plume_imax,'kmin',plume_kmin,'kmax',plume_kmax
    call write_log(message)
    write(message,*) 'plume_const_bmlt', plume_const_bmlt, 'plume_const_bmlt_rate', plume_const_bmlt_rate
    call write_log(message)
    write(message,*) 'plume_initial_bmlt:', plume_initial_bmlt
    call write_log(message)
    write(message,*) 'plume_delayed_coupling', plume_delayed_coupling
    call write_log(message)
    write(message,*) 'plume_do_cross_shelf_avg:', doCrossShelfAvg
    call write_log(message)
    write(message,*) 'plume_smooth_diff:',smooth_diff
    call write_log(message)
    write(message,*) 'plume_homotopy_frac:',plume_homotopy_frac
    call write_log(message)
    write(message,*) 'plume_homotopy_ramp:',plume_homotopy_ramp
    call write_log(message)

  end subroutine plume_read_print_config

  subroutine cross_shelf_average(lsrf)

     real(kind=dp),dimension(:,:), intent(inout) :: lsrf

     integer :: row_k

     do row_k=plume_kmin,plume_kmax
	
	lsrf(plume_imin:plume_imax,row_k) = sum(lsrf(plume_imin:plume_imax,row_k))/real(plume_imax-plume_imin+1)
   
     end do

  end subroutine 

end program shelf_driver
