module plume

  ! 2d depth-averaged plume model from model of
  ! johann jungclaus & jan backhaus, 1994
  ! (journal of geophysical research, 99, 12375) 

  ! converted to model of an ice shelf water plume using formulation of
  ! adrian jenkins and andreas bombosch, 1995
  ! (journal of geophysical research, 100, 6967)
  ! and
  ! lars smedsrud and adrian jenkins, 2004
  ! (journal of geophysical research, 109, c03025)

  ! by paul holland, 2004 onwards.  basic reference is
  ! paul holland and daniel feltham, 2006
  ! (journal of physical oceanography, 36, 2312)

  ! coupling to glimmer-cism glacier model by carl gladish
  ! 2009 - 2011

  ! model axes in code are 
  !   x/i/m - east (parallel to isobaths of ice shelf base, perp. to inflow)
  !   y/k/n - north (perp. to isobaths of ice shelf base, parallel to inflow)
  !   z     - upwards from the sea bed (defined as gldep+wcdep)

  use plume_global
  use plume_functions
  use plume_context
  use plume_io

  use omp_lib

  implicit none

  ! everything is private by default
  private

  ! except for functions mentioned here
  public :: plume_initialise,plume_runstep,plume_finalise,plume_iterate
  public :: nsteps

  ! module variables

  integer :: iwetmin,iwetmax,kwetmin,kwetmax,ikcplus
  integer :: icalcan,kcalcan,icalcen,kcalcen

  logical :: sepflag

  integer :: varoutrat
  integer :: nsteps, global_time_step_count = 0
  integer :: last_output_step = -1, last_lnot_step = -1


  !NB: all times are in seconds during exection, even though
  !    tottim and fouttim are given in days in the namelist file
  !   (all others are in seconds)
  real(kind=kdp) ::  tottim, &  ! total time to run model
       outtim, &  ! time interval between writing model state out to file
       fouttim, & ! first time at which to write out model outout
       louttim, & ! last time to write out model output
       snottim, & ! short note time interval
       lnottim, & ! long note time interval
       runtim,  & ! accumulated simulation time 
       coupletim, & ! time acculumated in current call to plume_iterate
       labtim, &     ! time interval in days used to name output files
       lastwritetim  ! time of last output (used in plume_iterate)

  real(kind=kdp) :: negdep ! total of negative depths (error diagnostic)

  ! rows in which buoyancy forcing is applied to generate gravity waves
  ! (in the case where use_periodic_forcing == .true.)
  integer,parameter :: klast=20,kfirst=2

contains

  subroutine plume_initialise(nl_filename, &
       suppress_ascii_output, &
       suppress_logging, &
       lsrf_ext,landmask_ext, &
       imin,imax,kmin,kmax)

    character(len=*),intent(in) :: nl_filename
    logical,intent(in) :: suppress_logging,suppress_ascii_output
    logical,dimension(:,:),intent(in),optional :: landmask_ext
    integer,optional,intent(in) :: imin,imax,kmin,kmax
    !NB: lsrf_ext specifies the initial height of the lower ice surface
    !    relative to sea level (therefore negative)
    real(kind=kdp),dimension(:,:),intent(in),optional :: lsrf_ext

    !local variables
    integer :: i,k,l
    integer :: itmp
    logical :: ltmp	
    real(kind=kdp) :: dtmp
    ! bpos_ext provides space to convert lsrf_ext to a surface that
    ! makes sense in plume world
    real(kind=kdp),dimension(:,:),allocatable :: bpos_ext	

    !in set parameters from namelist file
    call set_parameters(nl_filename,varoutrat,labtim, &
         tottim,outtim,fouttim,louttim,snottim,lnottim)
    ! allocate memory for global arrays, based on grid size parameters     
    call allocate_arrays()

    ! preliminary operations

    ! calculate general derived global quantities
    nsteps = int(dtswtim/dt1 + (tottim - dtswtim)/dt2)
    pi = 4.d0*atan(1.d0)
    radian = pi/180.d0
    f = 4.d0*pi*sin(phi*radian)/86164.d0 ! coriolis parameter 
    ! using a sidereal day

    dt = dt1
    gdt = grav*dt1
    fdt = f*dt1

    domain_imin = 1
    domain_imax = m_grid
    domain_kmin = 1
    domain_kmax = n_grid


    ! set up frazil parameters
    if (frazil) call frazil_calc()

    ! initialise fields, including topography and grids
    if (in_glimmer) then

       if (.not. present(imin) .or. .not. present(imax) .or. &
	    .not. present(kmin) .or. .not. present(kmax)) then
          write(*,*) "Plume initialization error: Must provide imin,imax,kmin,kmax when using plume in glimmer"
          stop
       end if

       allocate(bpos_ext(m_grid,n_grid))
       bpos_ext(:,:) = lsrf_ext(:,:) + gldep + wcdep 		

       if (present(landmask_ext)) then
          where( landmask_ext )
             !this instructs plume to define its own landmask to agree with glimmer's
             bpos_ext = 0.0
          end where
       end if

       domain_imin = imin
       domain_imax = imax
       domain_kmin = kmin
       domain_kmax = kmax

       call initialise_fields(suppress_ascii_output,iwetmin,iwetmax,kwetmin,kwetmax,bpos_ext)
       !deallocate(bpos_ext)
    else
       call initialise_fields(suppress_ascii_output,iwetmin,iwetmax,kwetmin,kwetmax)
    end if

    ! initialise wet area to surround all inflow points
    ikcplus = 2  ! set additional area in which to check wetness

    if (.not.use_min_plume_thickness) then
       icalcan = max0(domain_imin+1,iwetmin-ikcplus)
       kcalcan = max0(domain_kmin+1,kwetmin-ikcplus)
       icalcen = min0(domain_imax-1,iwetmax+ikcplus)
       kcalcen = min0(domain_kmax-1,kwetmax+ikcplus)
    else
       icalcan = domain_imin + 1
       kcalcan = domain_kmin + 1
       icalcen = domain_imax - 1
       kcalcen = domain_kmax - 1
    end if

    ! initialise time, total of negative depths, and separation and negative
    ! frazil warning counters
    runtim = 0.d0
    lastwritetim = -1.d0
    coupletim = 0.d0
    negdep = 0.d0
    sepflag = .false.
    do l=1,nice
       frzcut(l) = 0.d0
    end do

    ! initialise elapsed time and set system time and write to output
    call reset_elapsed_sys_time()

    if (.not. suppress_logging) then
       call io_append_output(' ')
       call io_output_sys_time('start system time ')
       call io_append_output(' ')
    end if

    if (.not. suppress_logging) &
       call io_write_time(get_elapsed_sys_time(), ' elapsed system time ')

  end subroutine plume_initialise


  subroutine plume_iterate(time, & 
       ice_dt, &
       lsrf, &
       landmask, &
       t_interior, &
       ice_dz, &
       bmelt_out, &
       btemp_out, &
       plume_stopping_tol, &
       plume_speed_stopping_tol, &
       min_run_time, &
       max_run_time, &
       run_plume_to_steady, &
       plume_reached_steady, &
       write_all_states, & 
       write_every_n, &
       write_frequency, &
       use_plume_initial_bmlt)

    ! NB: This subroutine is intended to be called only from a 
    ! glimmer-CISM driver that uses the plume to calculate the
    ! basal melt rate.
    ! The purpose of this subroutine is to run the plume
    ! over a period of time ice_dt, stopping if a steady-state
    ! is reached.  If run_plume_to_steady is
    ! set to .true. then the plume iterates to steady even if it
    ! takes longer than ice_dt.

    ! current simulation time in years
    real(kind=kdp),               intent(in) :: time

    ! length of external ice timestep in years
    real(kind=kdp),               intent(in) :: ice_dt 

    ! lower surface of external ice in meters
    real(kind=kdp),dimension(m_grid,n_grid),intent(in) :: lsrf

    ! thickness of lowest layer of ice in meters
    real(kind=kdp),dimension(m_grid,n_grid),intent(in) :: ice_dz

    ! landmask is equal to 1 in land or grounded ice
    logical, dimension(m_grid,n_grid), intent(in) :: landmask 

    !temperature in first interior layer of ice in Celcius
    real(kind=kdp),dimension(m_grid,n_grid),intent(in) :: t_interior 

    ! If true we will iterate until plume is steady, except 
    ! we won't exceed the max_run_time, for cases where 
    ! we expect only a statistically-steady state
    logical,                      intent(in) :: run_plume_to_steady

    ! If true, write out all plume states to file
    logical,                      intent(in) :: write_all_states

    integer,                      intent(in) :: write_every_n

    real(kind=kdp),               intent(in) :: write_frequency

    ! Steadiness tolerance.  Relative change in meltrate needed
    ! to continue plume time-stepping
    real(kind=kdp),               intent(in) :: plume_stopping_tol,plume_speed_stopping_tol

    real(kind=kdp),               intent(in) :: min_run_time !in days
    real(kind=kdp),               intent(in) :: max_run_time ! in days

    ! melt rate predicted by plume, in meters per year
    real(kind=kdp),dimension(m_grid,n_grid),intent(out) :: bmelt_out

    ! temperature of basal surface determined by the plume in Celcius
    real(kind=kdp),dimension(m_grid,n_grid),intent(out) :: btemp_out 

    ! flag to indicate that the plume reached a steady state
    logical,intent(out) :: plume_reached_steady

    logical,intent(in) :: use_plume_initial_bmlt

    !local variables

    real(kind=kdp) :: subcycling_time,ice_dt_in_sec
    real(kind=kdp),dimension(m_grid,n_grid) :: bmelt_old !in meters per year
    real(kind=kdp),dimension(m_grid,n_grid) :: rel_speed_change_array
    real(kind=kdp),dimension(m_grid,n_grid) :: speed_old,speed
    real(kind=kdp) :: max_rel_bmelt_change,max_rel_speed_change
    real(kind=kdp) :: mean_rel_bmelt_change,mean_rel_speed_change
    real(kind=kdp) :: prev_rel_change,prev_rel_speed_change
    real(kind=kdp) :: min_run_time_sec,max_run_time_sec
    character(len=512) :: log_message

    integer :: local_time_step_count !number of steps executed in this plume_iterate call
    local_time_step_count = 0
    coupletim = 0.d0



    ! convert lower surface depth into height of basal surface and interface
    ! surface, storing results in plume's global variable

    where( landmask )
       jcs = 0
       !bpos = gldep+wcdep
       bpos = 0.d0
       ipos = 0.d0
       pdep = 0.d0
    elsewhere
       jcs = 1
       bpos = lsrf + gldep + wcdep
       ipos = bpos - pdep
    end where

!    bpos = lsrf + gldep + wcdep
!    ipos = bpos - pdep

    !TODO: assign t_interior and ice_dz to globals 
    ! so that thermodynamics will pick up heat conduction into ice

    bmelt_old = bmelt
    speed_old = sqrt(su * su + sv * sv)

    subcycling_time = 0.0
    ice_dt_in_sec = ice_dt * 365.25d0*24.0d0*3600.d0
    min_run_time_sec = min_run_time * 3600.d0 * 24.d0
    max_run_time_sec = max_run_time * 3600.d0 * 24.d0

    plume_reached_steady = .false.
    prev_rel_change = 1.d0
    prev_rel_speed_change = 1.d0
    rel_speed_change_array = 0.d0

    ! while not steady

    if (.not. use_plume_initial_bmlt) then

       do while ((subcycling_time .le. min(max_run_time_sec,ice_dt_in_sec) &
                     .and. .not. run_plume_to_steady) &
           .or. &
                (run_plume_to_steady .and. &
                (subcycling_time .le. max_run_time_sec)))

          
          call plume_runstep()

          ! calculate max error
          max_rel_bmelt_change = maxval(abs(bmelt_old - bmelt)/ &
               (abs(bmelt)+epsilon(1.d0)))

          ! speed is an (m-2)*(n-2) array
          ! this step is necessary because su and sv are on different grids
          speed(2:(m_grid-1),2:(n_grid-1)) = sqrt(  &
               (5.d-1*(su(1:m_grid-2,2:n_grid-1)+ &
               su(2:m_grid-1,2:n_grid-1))) ** 2.d0 + &
               (5.d-1*(sv(2:m_grid-1,1:n_grid-2)+ &
               sv(2:m_grid-1,2:n_grid-1))) ** 2.d0 )

          if (local_time_step_count > 0) then                              

             rel_speed_change_array = abs(speed - speed_old) / &
                  (abs(speed_old)+epsilon(1.d0))

          end if

          max_rel_speed_change = maxval(rel_speed_change_array)

          where (rel_speed_change_array > plume_speed_stopping_tol)
             debug = log(rel_speed_change_array / plume_speed_stopping_tol)
          elsewhere
             debug = 0.d0
          end where

          if (subcycling_time >= min_run_time_sec .and. &
               (max_rel_bmelt_change/prev_rel_change < 1.d-1)) then
             prev_rel_change = max_rel_bmelt_change
             print '(a,e8.2,a,e8.2)', &
                  'max_rel_bmelt_change  ',max_rel_bmelt_change, &
                  '  stopping tol  ',plume_stopping_tol
          end if

          if((subcycling_time >= min_run_time_sec) .and. &
               (max_rel_speed_change/prev_rel_speed_change < 1.d-1)) then
             prev_rel_speed_change = max_rel_speed_change
             print '(a,e8.2,a,e8.2)', &
                  'max_rel_speed_change  ', max_rel_speed_change, &
                  '  stopping tol  ', plume_speed_stopping_tol
          end if

          if ((max_rel_bmelt_change < plume_stopping_tol) .and. &
               (max_rel_speed_change < plume_speed_stopping_tol)) then
             if (subcycling_time .ge. (min_run_time_sec)) then
                plume_reached_steady = .true.
                exit
             end if
          end if

          if (write_all_states .and. &
               (mod(local_time_step_count,write_every_n) == 0)) then
             call plume_netcdf_write_vars(time*3600.0d0*24.0d0*365.25d0 + &
                  subcycling_time)
          end if

          local_time_step_count = local_time_step_count + 1
          global_time_step_count = global_time_step_count + 1
          subcycling_time = subcycling_time + dt
          bmelt_old = bmelt
          speed_old = speed

       end do
    end if

    btemp_out = btemp
    bmelt_out = bmelt * (365.25d0*24.0d0*3600.0d0)

    if (.not. use_plume_initial_bmlt) then
       if (.not. plume_reached_steady) then
          call io_append_output('plume did not reach steady state')
       end if

       !we do this to make the plume output timestamps indicate the ice state
       ! with respect to which the plume was steady
       runtim = time*3600.0d0*24.0d0*365.25d0 


       write(log_message, '(a,f6.1)') 'subcycling time in days', &
            subcycling_time/(3600.0*24.0)
       call io_append_output(trim(log_message))

       if (runtim/(3600.d0*24.d0*365.25d0) - lastwritetim .ge. write_frequency) then 
          call io_write_surface_output(runtim,labtim)    
          lastwritetim = runtim/(3600.d0*24.d0*365.25d0)
       end if

    end if
  end subroutine plume_iterate

  ! **************************************************************************
  ! ******* main subroutine to carry out the current timestep ****************
  ! **************************************************************************

  subroutine plume_runstep()

    !The order in which helper routines are called is:
    ! tide
    ! momentum
    ! outflow_bound_u_v
    ! continuity
    ! outflow_bound_bmlt_entr
    ! inflow_calc
    ! subglacial_discharge
    ! update
    ! outflow_bound_interface
    ! outflow_bound_u_v
    ! scalar
    ! rho_calc
    ! outflow_bound_tsdc
    ! entrainment_correction
    ! filter_waves
    ! gwaves

    ! change timestep and recalculate basic quantities if necessary

    implicit none

    real(kind=kdp),dimension(m_grid,n_grid) :: long_waves,short_waves, as, bs
    integer :: k

    if (runtim .ge. dtswtim) then
       dt = dt2
       gdt = grav*dt2
       fdt = f*dt2
    end if

    if (global_time_step_count .ge. switch_entype_n) then

	entype = entype2
   
    end if

    ! update time in seconds

    runtim = runtim + dt
    coupletim = coupletim + dt

    ! find indices of current wetted area

    if (.not.use_min_plume_thickness) then
       icalcan = max0(domain_imin+1,iwetmin-ikcplus)
       kcalcan = max0(domain_kmin+1,kwetmin-ikcplus)
       icalcen = min0(domain_imax-1,iwetmax+ikcplus)
       kcalcen = min0(domain_kmax-1,kwetmax+ikcplus)
    end if

    call tide(icalcan,kcalcan,icalcen,kcalcen)

    ! --------------------
    ! calculate velocities
    ! --------------------
    call momentum(icalcan,kcalcan,icalcen,kcalcen)

    ! set velocity components on the open boundaries if the plume 
    ! has reached that far
    call outflow_bound_u_v(icalcan,icalcen,kcalcan,kcalcen)

    ! ---------------------------------------------------
    ! evaluate continuity equation for interface position
    ! ---------------------------------------------------

    call continuity(icalcan,kcalcan,icalcen,kcalcen)

    call outflow_bound_bmelt_entr(icalcan,icalcen,kcalcan,kcalcen)

    ! interface and scalars passed forward on open boundary
    call inflow_calc(icalcan,kcalcan,icalcen,kcalcen)

    ! mix in the prescribed subglacial discharge flux near inflow edge
!   this code is redundant now, since the sgd code has been added into 
!   the continuity and scalar subroutines
!    call subglacial_discharge(icalcan,kcalcan,icalcen,kcalcen)

    ! ----------------------------------------------------------
    ! update plume thickness, wet/dry boundaries, and velocities
    ! ----------------------------------------------------------
    call update(iwetmin,iwetmax,kwetmin,kwetmax,&
         icalcan,icalcen,kcalcan,kcalcen,negdep)


    ! update interface position and flow on open boundaries

    call outflow_bound_interface(icalcan,icalcen,kcalcan,kcalcen)


    call outflow_bound_u_v(icalcan,icalcen,kcalcan,kcalcen)


    ! ---------------------------------------------------------------
    ! solve transport equations for temperature, salinity, and frazil
    ! ---------------------------------------------------------------

    call scalar(icalcan,kcalcan,icalcen,kcalcen)             

    ! -----------------------------------
    ! calculate density and update bounds
    ! -----------------------------------

    call rho_calc(icalcan,icalcen,kcalcan,kcalcen,sepflag)


    call outflow_bound_tsdc(icalcan,icalcen,kcalcan,kcalcen)

    if (entype == 5 .or. entype == 6 .or. entype == 7) then
       call entrainment_correction(icalcan, icalcen, kcalcan, kcalcen, entype)
    end if
!    call filter_waves(pdep, jcs, short_waves, long_waves, as, bs)

!    debug2 = short_waves
    do k=kcalcan,kcalcen
       debug2(:,k) = sum(bmelt(:,k))/( m_grid - 6 ) * 3600.d0*24.d0*365.25d0
    end do


    call gwaves(icalcan, icalcen, kcalcan, kcalcen)
    ! write short step output and reset separation warning flag

    if (mod(runtim,snottim).eq.0) then
       sepflag = .false.
       !call io_write_calculated_time(runtim,int(runtim/dt))
       call io_write_calculated_time(coupletim,int(coupletim/dt))
    end if

    ! write long step output 
    if ((int(floor(runtim/lnottim)) > last_lnot_step)) then
       last_lnot_step = int(floor(runtim/lnottim))

       call io_write_long_step_output(icalcan,icalcen,kcalcan,kcalcen,&
            varoutrat,negdep)
       !call io_write_calculated_time(runtim, int(runtim/dt))
       call io_write_calculated_time(coupletim,int(coupletim/dt))
    end if

    ! write state to file
    if ((int(floor(runtim/outtim)) > last_output_step) .and. &
         ((runtim.ge.fouttim).and.(runtim.le.louttim)) .and. &
         .not.(in_glimmer)) then
       last_output_step = int(floor(runtim/outtim))

       call io_write_surface_output(runtim,labtim)
    end if

  end subroutine plume_runstep

  subroutine plume_finalise()

    call deallocate_arrays()
    call io_write_time(get_elapsed_sys_time(),' :-) total system time ')

  end subroutine plume_finalise


  ! *************************************************************************
  ! ********* private (not callable outside this module) helper subroutines *
  ! *************************************************************************

  subroutine set_parameters(nl_filename,varoutrat,labtim, &
       tottim,outtim,fouttim,louttim,snottim,lnottim)

    ! set program parameters
    ! where appropriate, read from namelist file
    ! derive some parameters at end of routine, after read from file

    implicit none

    ! parameters
    character(len = *),intent(in) :: nl_filename
    integer,intent(out) :: varoutrat
    real(kind=kdp),intent(out) :: labtim,tottim,outtim,fouttim,louttim,snottim,lnottim

    ! local variables
    integer :: l
    real(kind=kdp) :: infloain,infloein,knfloain,knfloein
    real(kind=kdp) :: cseedfix,cinffix	

    namelist /plume_nml/ &
            mixlayer &
	 ,  in_glimmer &
         ,  restart &
         ,  restart_data_filename &
         ,  frazil &
         ,  nonlin &
         ,  horturb &
         ,  entrain &
         ,  entype &
         ,  entype2 &
	 ,  switch_entype_n &
         ,  C_s &
         ,  C_n &
         ,  C_i &
         ,  basmelt &
         ,  rholinear &
         ,  thermobar &
         ,  intrace &
         ,  vardrag &
         ,  topedit &
         ,  tangle &
         ,  negfrz &
         ,  use_min_plume_thickness &
         ,  use_neutral_salinity &
         ,  plume_southern_bc &
         ,  use_periodic_forcing &
         ,  forcing_period &
         ,  periodic_forcing_amp &
         ,  tottim &
         ,  fouttim &
         ,  outtim &
         ,  labtim &
         ,  snottim &
         ,  lnottim &
         ,  dt1 &
         ,  m_grid & 
         ,  n_grid &
         ,  namb &
         ,  hx &
         ,  hy &
         ,  gldep &
         ,  ifdep &
         ,  wcdep &
         ,  plume_min_thickness &
         ,  plume_max_thickness &
         ,  gaspar_cutoff &
         ,  u_star_offset &
         ,  tidal_velocity &
         ,  entr_time_const &
         ,  detrain_time_const &
         ,  context &
         ,  bathtype &
	 ,  slope_direction &
         ,  cross_slope_wavenumber &
         ,  along_slope_deepening_exp &
         ,  channel_amplitude &
         ,  random_amplitude &
         ,  bsmoothit &
         ,  depinit &
         ,  depinffix &
         ,  meltinf &
         ,  salttop &
         ,  infloain &
         ,  infloein &
         ,  knfloain &
         ,  knfloein &
         ,  saltbot &
         ,  temptop &
         ,  tempbot &
         ,  phi &
         ,  ah &
         ,  kh &
         ,  cdb &
         ,  cl & 
         ,  ef &
         ,  tiuniform &
         ,  min_melt_depth &
         ,  nus &
         ,  nbar &
         ,  nice &
         ,  seedtype &
         ,  cseedfix &
         ,  cinffix &
         ,  gasp_m1 &
         ,  gasp_m2 &
         ,  gasp_m3 &
         ,  gasp_m4 &
         ,  gasp_m5 &
         ,  a1 &
         ,  a2 &
         ,  nk_m &
         ,  nk_n &
         ,  n_amb_ctl_pt &
         ,  amb_temp_ctl_pt &
         ,  amb_salt_ctl_pt &
         ,  amb_depth_ctl_pt &
         ,  sgd_type &
         ,  sgd_flux 
     
    ! ++++++++++++++++++
    ! set default values
    ! ++++++++++++++++++

    ! set ascii output (screen and file) properties 

    varoutrat = 1       ! output resolution varies:
    ! 1 - in both directions
    ! 2 - according to extent of plume in y-direction
    ! 3 - according to extent of plume in x-direction

    ! set switches 
    ! ------------
    use_min_plume_thickness = .false.
    use_neutral_salinity = .false.
    use_periodic_forcing = .false.
    mixlayer    = .false. ! model a mixed-layer, rather than a plume
    ! (i.e. whole domain is given initial thickness)
    in_glimmer = .false.  ! by default, not running inside glimmer (ice shelf model)
    restart= .false. ! restart from previous model dump
    restart_data_filename = '' !netcdf file to read from for a restart
    frazil      = .false. ! include frazil ice
    nonlin      = .true.  ! advection terms included
    horturb     = .false. ! horizontal diffusion included
    entrain     = .true.  ! entrainment included
    entype      = 1       ! entrainment parameterisation:
    entype2     = -1      ! switch to after switch_entype_n
    switch_entype_n = 0   ! how many steps to execute before switching entype

    ! 1 - full kochergin formulation
    ! 2 - fractional kochergin formulation (using ef)
    ! 3 - halved pedersen formulation
    ! 4 - halved and modified pedersen formulation
    ! 5 - ZM 2001 diagnosed mixed-layer thickness
    ! 6 - Gaspar 1988 CMO equations

    C_n = 0.5d0
    C_i = 20.d0
    C_s = 10.d0
    alpha = -3.87d-5
    beta = 7.86d-4

    gasp_m1 = 0.45d0
    gasp_m2 = 2.6d0
    gasp_m3 = 1.9d0
    gasp_m4 = 2.3d0
    gasp_m5 = 0.6d0
    a1 = 0.6d0
    a2 = 0.3d0

    nk_m =  0.45d0
    nk_n =  0.2d0

    n_amb_ctl_pt = 0
    amb_temp_ctl_pt = 0.d0
    amb_salt_ctl_pt = 0.d0
    amb_depth_ctl_pt = 0.d0

    sgd_type = -1         ! no discharge, by default
    sgd_flux = 0.d0       ! no flux by default

    basmelt     = .true.  ! include direct basal melting and freezing
    rholinear   = .true.  ! use linear equation of state instead of unesco
    thermobar   = .false. ! thermobaricity included (code needs checking)
    intrace     = .false. ! include tracers tracking the fate of each inflow
    vardrag     = .false. ! include spatially-varying drag coefficient
    topedit     = .false. ! apply hand-editing to read-in topography
    tangle      = .false. ! use turning angle when applying drag
    ! (probably wise to check or reprogram)
    negfrz      = .false. ! artificially keep frazil concentrations positive

    ! set simulation time parameters
    ! ------------------------------

    tottim  = 0050.0d0  ! total simulation time in days
    outtim  = 3600.0d0   ! file output frequency in seconds
    labtim  = 001.0d0   ! units in which to name files in days
    snottim =  3600.0d0   ! short note output frequency in seconds
    lnottim = 86400.0d0   ! long note output frequency in seconds

    dt1     = 0360.0d0  ! long first timestep in seconds

    ! set default grid and domain properties
    ! ------------------------------
    m_grid = 0600           ! number of east-west cells
    n_grid = 0600           ! number of north-south cells
    hx = 5000.d0       ! east-west cell dimension (m) for uniform grid
    hy = hx            ! north-south cell dimension (m) for uniform grid
    gldep = 1400.d0    ! thickness of ice shelf at grounding line
    ifdep = 285.d0     ! thickness of ice shelf at ice front (or plateau)
    wcdep = 1600.d0    ! depth of water column beneath grounding line 
    plume_min_thickness = 1.d0 ! artificially determined minimum plume thickness
    plume_max_thickness = 200.d0
    gaspar_cutoff = 100.d0
    u_star_offset = 0.d0       ! source of turbulent kinetic energy for entraining
    tidal_velocity = 0.d0      ! another idea to inject TKE 
    entr_time_const = 0.d0     ! relaxation time (in seconds) for plume thickness to attain minimum
                               ! if no value is assigned in namelist, will be set to equal timestep
    detrain_time_const = 0.d0  ! relaxation time (in seconds) for plume thickness to relax to Monin-Obukhov L
                               ! if no value is assigned in namelist, will be set to equal timestep
    plume_southern_bc = 0    !use the d/dy = 0 version of southern outflow boundary by default

    ! set topography properties
    ! -------------------------
    context = "isomip" ! topography choices to use if bathtype = 0
    ! choices are currently fris,larsen, or isomip
    ! note isomip topography incomplete

    bathtype = 01     ! ice shelf bathymetry scenario:
    ! 0 => read bathymetry from file
    ! 1 => no walls, rises northwards
    ! 2 => right-angle west wall, rises northwards
    ! 3 => 45 degree west wall, rises northwards
    ! 4 => right-angle west wall, rises || to wall
    ! 5 => 45 degree west wall, rises || to wall
    ! 6 => carlson inlet
    ! 7 => rutford ice stream 
    ! 8 => carlson inlet test 1
    ! 9 => carlson inlet test 2
    ! 10=> evans ice stream
    ! 11=> whole filchner-ronne simplified
    ! 12=> analytical ice shelf
    ! 13=> longitudinal channels 'Peterman style'

    ! For bathtype 13, the number of sinusoidal undulations 
    ! going across the shelf.  
    cross_slope_wavenumber = 0.d0
    along_slope_deepening_exp = 0.d0
    channel_amplitude = 0.d0
    random_amplitude = 0.d0

    periodic_forcing_amp = 0.d0
    forcing_period = 0.d0

    kcorn = int(135.d0*1000.d0/hx) ! dist. of corner from inflow (km)
    rad = int(35.d0*1000.d0/hx)    ! radius of rounding on corner (km)

    cweight = 4.d0         ! smoothing weight of central point 
    nweight = 1.d0         ! smoothing weight of neighbour point
    bsmoothit = 00         ! its of smoothing to apply throughout domain

    ! smooth region flags and iterations (including bsmoothit) if variable
    ! smoothing
    ! 1) fris
    smflag(01)   = .false. ! rutford ice stream
    smoothit(01) = 20
    smflag(02)   = .false. ! support force glacier
    smoothit(02) = 10
    smflag(03)   = .false. ! spare slot
    smoothit(03) = 00
    ! 2) larsen
    smflag(01)   = .false. ! choyce point
    smoothit(01) = 15
    smflag(02)   = .false. ! thuronyi bluff
    smoothit(02) = 15
    smflag(03)   = .false. ! spare slot
    smoothit(03) = 00

    ! set inflow properties and initial thickness
    ! -------------------------------------------
    ninfmin = 01             ! lowest-numbered inflow turned on
    ninfmax = 01             ! highest-numbered inflow turned on 

    depinit = 01.0d0    ! initial plume thickness if mixed-layer model
    depinffix = 01.d0   ! plume thickness assigned to inflow regions
    meltinf = 9.8d-1    ! how far along gade line inflow properties are:
    ! 1 => properties at shelf, 0 => properties of ambient

    ! inflow extent if using single inflow on simple bathymetry
    infloain = 290.d0   ! western boundary (km)
    infloein = 310.d0   ! eastern boundary (km)
    knfloain = 0.d0    ! southern boundary (km)
    knfloein = 3.d0    ! northern boundary (km)

    slope_direction = 0   ! 0 = north, 1 = east, 2 = west, 3 = south

    ! set ambient fluid properties
    ! ----------------------------
    namb = 301           ! increments in ambient water column (minimum +2)
    dzincr = (gldep+wcdep)/(namb-1)       ! depth of each increment 
    ! (so total (namb - 1)*dzincr = gldep + wcdep)
    salttop = 34.500d0   ! sea surface salinity 
    saltbot = 34.950d0   ! sea bed salinity (gldep + wcdep)
    temptop = -1.900d0   ! sea surface temperature 
    tempbot = -2.500d0   ! sea bed temperature (gldep + wcdep)

    ! set physical and geographical parameters
    ! ----------------------------------------
    grav = 9.81d0       ! gravitational constant
    phi = +00.d0     ! latitude (in degrees, negative => s. hemisphere)
    ! (set to zero to remove rotation)
    ah = 0000.d0     ! horizontal eddy viscosity
    kh = ah          ! horizontal eddy diffusivity
    cdb = 2.5d-3     ! bottom drag coefficient (background value if variable)
    cdbvar = 5.0d-3  ! bottom drag coefficient in varied areas (if variable)
    cl = 1.775d-2    ! kochergin entrainment coefficient (should be 2.75d-2)
    ef = 5.0d-1      ! kochergin entrainment factor (1=> match jung & back)
    rho0 = 1030.d0   ! reference seawater density
    rhoi = 920.d0    ! density of ice

    ! drag region flags used if using variable drag
    drflag(01) = .false.  !east of henry ice rise
    drflag(02) = .false.  !east of berkner island
    drflag(03) = .false.  !spare slot

    ! set ice-related physical parameters
    ! -----------------------------------
    lat = 335000.d0  ! latent heat of ice fusion
    c0 = 3974.d0     ! specific heat capacity of weddell seawater
    ci = 2009.d0     ! specific heat capacity of ice shelf
    nu0 = 1.95d-6    ! kinematic viscosity of seawater
    pr = 13.8d0      ! molecular prandtl number of seawater
    sc = 2432.d0     ! molecular schmidt number of seawater
    fta = -5.73d-2   ! slope of liquidus for seawater (tf decrease with s)
    ftb = 8.32d-2    ! offset of liquidus for seawater (tf at s=0, z=surface)
    ftc = -7.61d-4   ! freezing temp change with depth (tf decrease with zeta)
    tiuniform = -25.d0      ! temperature of shelf (heat conduction during melting)
    min_melt_depth = 0.d0 !minimum depth at which melting occurs (crude separation modelling)
    si = 0.d0        ! salinity of shelf (salt trapped in ice during freezing)
    nus = +1.d0      ! nusselt number: 
    ! >=0 - that constant value
    ! -1  - correct full variable (not working)
    ! -2  - correct variable (no turbulent part for large)
    ! -3  - incorrect full hammar&shen variable (not working)
    ! -4  - incorrect h&s variable (no turbulent for large)
    kt = 1.4d-7      ! molecular thermal diffusivity of seawater
    ks = 8.0d-10     ! molecular haline diffusivity of seawater
    ar = 2.0d-2      ! aspect ratio of frazil discs
    eps = 7.4d-6     ! turbulent dissipation rate
    nbar = 1.0d3     ! cap on number of crystals for secondary nucleation

    ! set frazil model parameters
    ! ---------------------------
    nice = 10        ! number of frazil size classes

    r(01)  = 0.01d-3 ! radius of frazil discs
    r(02)  = 0.05d-3    
    r(03)  = 0.15d-3    
    r(04)  = 0.3d-3    
    r(05)  = 0.4d-3    
    r(06)  = 0.5d-3    
    r(07)  = 0.6d-3    
    r(08)  = 0.8d-3    
    r(09)  = 1.0d-3    
    r(10)  = 2.0d-3    

    seedtype = 2      ! seeding strategy:
    ! 1 => seed at southernmost supercooled cells 
    !      (smedsrud & jenkins)
    ! 2 => seed newly-supercooled cells if c(i) < cseed(i) 
    !      (holland & feltham)

    cseedfix = 1.0d-7 ! seed concentration in class i
    cinffix  = 0.0d0  ! inflow concentration in class i

    ! set modelling parameters
    ! ------------------------
    edepth = 1.0d-3  ! critical min depth for entrainment
    mdepth = edepth  ! critical min depth for basal melting
    fdepth = edepth  ! critical min depth for frazil dynamics
    dcr    = edepth  ! critical min plume depth
    septol = 1.0d-2  ! smallest plume-ambient density anomaly tolerated before
                     ! plume separation is flagged
    septol = 0.0d0		 
    small = epsilon(1.d0) !1.0d-15  ! smallest number

    ! +++++++++++++++++++++
    ! read values from file
    ! +++++++++++++++++++++

    open(21,file=trim(nl_filename),status='old')

    read(21,plume_nml)

    close(21)

    if (abs(phi) < small .and. entype == 6) then
       call io_append_output('Can not run Gaspar entrainment without rotation')
       stop 1
    end if

    if (use_min_plume_thickness .and. mixlayer) then
       call io_append_output('Can not run using useminthickness = .t. &
       	                     &and mixlayer = .t.')
       stop 1
    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! some parameters dependent upon namelist-read variables
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! convert timesteps to seconds

    tottim  = tottim*24.d0*3600.d0  ! total simulation time in seconds
    labtim  = labtim*24.d0*3600.d0  ! units in which to name files
    fouttim = fouttim*24.d0*3600.d0 ! first file output in seconds
    louttim = tottim                ! last file output in seconds

    dt2     = dt1                   ! short second timestep in seconds
    dtswtim = tottim                ! timestep switch time in seconds

    !    tottim = dt1
    !    labtim = dt1
    !    outtim = dt1
    !    snottim = dt1
    !    lnottim = dt1
    !    louttim = dt1
    !    fouttim = dt1

    ! derive ideal inflow extents

    infloa = int(infloain*1000.d0/hx)! western boundary (cells)
    infloe = int(infloein*1000.d0/hx)! eastern boundary (cells)
    knfloa = int(knfloain*1000.d0/hy)! southern boundary (cells)
    knfloe = int(knfloein*1000.d0/hy)! northern boundary (cells)

    ! infloain/infloein/knfloain/knfloein are now in units of cells
    infloa = int(infloain)
    infloe = int(infloein)
    knfloa = int(knfloain)
    knfloe = int(knfloein)

    
    ! set case-specific inflow flags if using realistic bathymetry

    if (context.eq."fris") then

       inflag(01) = .true.   ! evans ice stream
       inflag(02) = .false.  ! carlson inlet
       inflag(03) = .false.  ! rutford ice stream
       inflag(04) = .false.  ! institute ice stream
       inflag(05) = .false.  ! mollereisstrom
       inflag(06) = .false.  ! foundation ice stream
       inflag(07) = .false.  ! support force glacier
       inflag(08) = .false.  ! recovery glacier
       inflag(09) = .false.  ! slessor glacier

    end if

    if (context.eq."larsen") then

       inflag(01) = .false.   ! attlee glacier
       inflag(02) = .true.    ! whole larsen b
       inflag(03) = .true.    ! whole larsen c

    end if

    if (context.eq."isomip") then

       inflag(01) = .true.    ! whole southern boundary

    end if

    ! determine ambient fluid gradients
    dzincr = (gldep+wcdep)/(namb-1)       ! depth of each increment 
    sgrad = (saltbot - salttop)/(gldep + wcdep) ! rate of s change with depth
    tgrad = (tempbot - temptop)/(gldep + wcdep) ! rate of t change with depth

    ! set frazil seeds and inflows

    do l = 1,nice     ! frazil seed population 
       cseed(l) = cseedfix
    end do
    cinftot = 0.d0    ! inflow frazil
    do l = 1,nice
       cinf(l) = cinffix
       cinftot = cinftot + cinf(l)
    end do

    if (entr_time_const == 0.d0) then
       ! no value was set in the namelist
       entr_time_const = dt1
    end if

    if (detrain_time_const == 0.d0) then
       ! no value was set in the namelist
       detrain_time_const = dt1
    end if

    if (entype2 < 0) then
	entype2 = entype
    end if
    
	
  end subroutine set_parameters

  subroutine initialise_fields(suppress_ascii_output, &
       iwetmin,iwetmax,kwetmin,kwetmax, &
       bpos_ext)

    ! read data and set all initial fields

    use plume_io

    implicit none

    logical,intent(in) :: suppress_ascii_output

    ! NB: intent(inout) means that if this subroutine doesn't change *wetmin/max 
    !     then the calling subroutine will keep the original values
    integer,intent(inout) :: iwetmin,iwetmax,kwetmin,kwetmax

    real(kind=kdp),dimension(:,:),intent(in),optional :: bpos_ext

    ! local variables
    integer :: i,k	
    real(kind=kdp),dimension(m_grid,n_grid):: zd
    real(kind=kdp),dimension(m_grid,n_grid):: depth,sambindep,tambindep,c1,c2,c3,rhoa

    ! initialise all arrays (0 is assigned to all array locations)

    utrans = 0.d0
    vtrans = 0.d0
    utransa = 0.d0
    vtransa = 0.d0
    gwave_speed = 0.d0
    gwave_crit_factor = 0.d0
    su = 0.d0
    sv = 0.0d0
    u0 = 0.d0
    v0 = 0.d0
    u0a = 0.d0
    v0a = 0.d0
    tang = 0.d0
    pdep = 0.d0
    ipos = 0.d0
    bpos = 0.d0
    separated = 0
    jcs = 0
    jcw = 0
    jcd_u = 1 !initial velocity of 0.d0 is valid
    jcd_v = 1
    !jcd_u = 0
    !jcd_v = 0
    jcd_fl = 0       
    jcd_negdep = 0
    jcd_fseed = 0
    ctot = 0.d0
    tfreeze = 0.d0
    entr = 0.d0
    sgd = 0.d0
    artf_entr_frac = 0.d0
    local_tidal_speed = 0.d0
    thk_def = 0.d0
    atemp = 0.d0
    asalt = 0.d0
    drag = 0.d0
    bmelt = 0.d0
    btemp = 0.d0
    bsalt = 0.d0
    ctempd = 0.d0
    tint = tiuniform

    if (restart) then
       
       call plume_netcdf_read_int_var(restart_data_filename, 'jcs', jcs)        
       call plume_netcdf_read_int_var(restart_data_filename, 'jcw', jcw)
       call plume_netcdf_read_int_var(restart_data_filename, 'jcd_u', jcd_u)
       call plume_netcdf_read_int_var(restart_data_filename, 'jcd_v', jcd_v)
       call plume_netcdf_read_int_var(restart_data_filename, 'jcd_fl', jcd_fl)

       call plume_netcdf_read_real_var(restart_data_filename, 'bpos', bpos)
       call plume_netcdf_read_real_var(restart_data_filename, 'ipos', ipos)
       call plume_netcdf_read_real_var(restart_data_filename, 'pdep', pdep)
       call plume_netcdf_read_real_var(restart_data_filename, 'u', utrans)
       call plume_netcdf_read_real_var(restart_data_filename, 'v', vtrans)
       call plume_netcdf_read_real_var(restart_data_filename, 'su', su)
       call plume_netcdf_read_real_var(restart_data_filename, 'sv', sv)

       call plume_netcdf_read_real_var(restart_data_filename, 'bmelt', bmelt)
       call plume_netcdf_read_real_var(restart_data_filename, 'btemp', btemp)
       call plume_netcdf_read_real_var(restart_data_filename, 'bsalt', bsalt)
       call plume_netcdf_read_real_var(restart_data_filename, 'rhop' , rhop)
       call plume_netcdf_read_real_var(restart_data_filename, 'temp' , temp)
       call plume_netcdf_read_real_var(restart_data_filename, 'salt' , salt)
       call plume_netcdf_read_real_var(restart_data_filename, 'entr' , entr)

       bmelt = bmelt / (365.25d0*24.d0*3600.d0)  ! convert the m/year quantity to m/s
       utransa = utrans
       vtransa = vtrans
       tempa = temp
       salta = salt

    end if

    ! initialise ice to zero even when frazil is off so that densities are correct
    c_ice = 0.d0
    ca_ice = 0.d0
    fmelt = 0.d0
    fppn = 0.d0
    fnuc = 0.d0

    tempinf = 0.d0
    saltinf = 0.d0
    depinf = 0.d0

    debug = 0.d0
    debug2 = 0.d0
    debug3 = 0.d0

    if (intrace) then
       intrin = 0
       intracer = 0.d0
    end if

    ! get cell dimensions
    ! -------------------

    call grid_set()

    ! get topography and inflow regions, either from data or settings
    ! ---------------------------------------------------------------


    if (in_glimmer) then
       if (restart) then
          ! do nothing
       else

          if (.not. present(bpos_ext)) then
             call io_append_output('Need to provide bpos_ext')
             stop 1
          end if
          bpos = bpos_ext
          if (.not. use_min_plume_thickness) then
             call topog_depth_inflow_set(.not. use_min_plume_thickness .and. .not. mixlayer)

          end if
        end if

    elseif (restart) then	
       !not in glimmer and doing a restart
    else

       if (bathtype.gt.0) then
          call topog_depth_inflow_set(.not. use_min_plume_thickness .and. .not. mixlayer)
       else
          call topog_read_edit()
          call topog_smooth()
          if (.not. use_min_plume_thickness) then
             if (context.eq."fris")   call inflow_set_fris()
             if (context.eq."isomip") call inflow_set_isomip() 
             if (context.eq."larsen") call inflow_set_larsen()
          end if
       end if

    end if

    ! set ambient properties
    ! ----------------------

    call set_ambient(suppress_ascii_output)

    ! set drag coefficient
    ! --------------------

    drag = cdb

    if (vardrag) call drag_set()

    ! set inflow properties according to gade (jpo 1979) lines  
    ! --------------------------------------------------------
    ! (depth of plume mid-point so that t=tf initially if meltinf = 1)

    where (depinf > 0.d0) 

       depth = wcdep + gldep - bpos + depinf/2.d0

       ! interpolate ambient values 
       tambindep = get_tamb_z(depth)
       sambindep = get_samb_z(depth)

       ! calculate gade values
       c1 = fta*(1.d0 - ci/c0)
       c2 = (ftb + ftc*depth)*(1.d0 - ci/c0) &
            - tambindep - lat/c0 + (ci/c0)*(fta*sambindep + tint)
       c3 = sambindep*(lat/c0 + (ci/c0)*(ftb + ftc*depth - tint))
       saltinf = - (c2 + dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)
       tempinf = freezing_temp_func(saltinf,depth) 

       ! calculate plume properties from proportions 
       ! of ambient and melt waters

       !saltinf = sambindep + meltinf*(saltinf - sambindep)
       !tempinf = tambindep + meltinf*(tempinf - tambindep)
       
    end where

    ! initialise whole domain to wet and set plume initial thickness if mixed-layer
    ! model

    if (mixlayer) then

       iwetmin = domain_imin
       iwetmax = domain_imax
       kwetmin = domain_kmin
       kwetmax = domain_kmax

       where (bpos > 0.d0 .and. bpos < (wcdep+gldep))
          pdep = depinit
       end where

    else if (use_min_plume_thickness) then

       ! using minimum thickness, so we will never 
       ! look at iwetmin/iwetmax/kwetmin/kwetmax
       continue

    else

       ! not using minimum thickess or mixlayer
       iwetmin = domain_imax
       iwetmax = domain_imin
       kwetmin = domain_kmax
       kwetmax = domain_kmin

       do k = domain_kmin,domain_kmax
          do i = domain_imin,domain_imax
             if (depinf(i,k).gt.0.d0) then
                iwetmin = min(iwetmin,i)
                iwetmax = max(iwetmax,i)
                kwetmin = min(kwetmin,k)
                kwetmax = max(kwetmax,k)
             end if
          end do
       end do

    end if

    ! initialise flags, depths, and scalars
    ! -------------------------------------
    if (.not. restart) then

       jcs = 1

       ipos = bpos - pdep
       where (bpos <= dcr) 
          ! set index for solid ice points (no water column)
          ipos = bpos
          jcs = 0
       end where

       if (use_min_plume_thickness) then	
          ! start off the plume thickness at the minimum everywhere
          do i=domain_imin,domain_imax
             do k=domain_kmin,domain_kmax
                if (jcs(i,k) .eq. 1) then
                   pdep(i,k) = plume_min_thickness                
                   ipos(i,k) = bpos(i,k) - pdep(i,k)	                
                end if
             end do
          end do
       else
          pdep = bpos - ipos
       end if

    
       ! set field for wet/dry points (jcw)
       where (pdep >= dcr) jcw = 1
       
    end if

    ! set fields to ambient-fluid properties for depth
    ! (used when considering newly-wet cells)

    zd = max(0.d0,wcdep + gldep - bpos)


    rhoa = get_rhoamb_z(zd)
    rhoamb = rhoa

    if (.not. restart) then
       tempa = get_tamb_z(zd)
       salta = get_samb_z(zd)
       temp = tempa
       salt = salta
       rhop = rhoa
    end if


  end subroutine initialise_fields

  subroutine grid_set()

    implicit none

    ! read grid cell data

    ! local variables
    integer :: i,k,imid,kmid

    real(kind=kdp) :: rmdx,adx,adxu
    real(kind=kdp) ::  rmdxu,rmdy,rmdyv,ady,adyv
    real(kind=kdp) :: enhance_factor

    imid=(m_grid+1)/2
    kmid=(n_grid+1)/2

    ! ----------------
    ! grid for scalars
    ! ----------------

    ! set east-west cell dimensions for uniform scalar grid
    dx = hx

    ! set north-south cell dimensions for uniform scalar grid
    dy = hy

    ! -----------------------------------------------
    ! grid for velocities and turbulence coefficients
    ! -----------------------------------------------

    ! calculate cell dimensions for velocity grid
    do k=2,n_grid
       dyv(k) = 5.0d-1*(dy(k) + dy(k-1))
       rdy(k)=1.d0/dy(k)
       rdyv(k)=1.d0/dyv(k)
    end do
    dyv(1)=dyv(2)
    rdy(1)=rdy(2)
    rdyv(1)=rdyv(2)

    do i=2,m_grid
       dxu(i) = 5.0d-1*(dx(i) + dx(i-1))
       rdx(i)= 1.d0/dx(i)
       rdxu(i)= 1.d0/dxu(i)
    end do
    dxu(1)=dxu(2)
    rdx(1)=rdx(2)
    rdxu(1)=rdxu(2)

    ! calculate coefficients for horizontal eddy viscosity
    rmdx = rdx(imid)
    rmdxu = rdxu(imid)
    adx = dt*ah*rmdx**2
    adxu = dt*ah*rmdxu**2

    do k=1,n_grid
      ahdx(:,k)  = adx/rdx
      ahdxu(:,k) = adxu/rdxu
    end do

    rmdy = rdy(kmid)
    rmdyv = rdyv(kmid)
    ady = dt*ah*rmdy**2
    adyv = dt*ah*rmdyv**2

    do i=1,m_grid
       ahdy(i,:) = ady/rdy
       ahdyv(i,:) = adyv/rdyv
    end do

    khgrid = kh

    do k=floor(n_grid*visc_enhance_location),n_grid

       enhance_factor = visc_enhance_factor*min(1.d0, &
                                                real(k-floor(n_grid*visc_enhance_location)) / &
                                                visc_enhance_smoothing_ramp)
       ahdx(:,k) = (1.d0+enhance_factor) * ahdx(:,k)
       ahdxu(:,k)= (1.d0+enhance_factor) * ahdxu(:,k)
       ahdy(:,k) = (1.d0+enhance_factor) * ahdy(:,k)
       ahdyv(:,k)= (1.d0+enhance_factor) * ahdyv(:,k)
       khgrid(:,k)=(1.d0+enhance_factor) * khgrid(:,k)

    end do

    print *, 'ahdx', ahdx(25,:)
    print *, 'khgrid', khgrid(25,:)

  end subroutine grid_set

  subroutine set_ambient(suppress_ascii_output)

    ! set ambient profiles of temperture, salinity and density
    ! density is unaffected by frazil as the ambient is ice-free

    implicit none

    logical,intent(in) :: suppress_ascii_output

    ! local variables
    integer :: i,i_ctl, i_prev_ctl_depth
    real(kind=kdp),dimension(namb) :: depth,pressure,ttt,rhopot
    real(kind=kdp) :: ctl_tgrad, ctl_sgrad
    real(kind=kdp) :: t_prev, s_prev, z_prev
    logical :: found_ctl_level 

    depth = (/ (0.d0 + i*dzincr,i=0,(namb-1)) /)

    if (n_amb_ctl_pt == 0) then

       tamb = temptop + depth*tgrad

       if (use_neutral_salinity) then

          ! figure out with salinity gradient 
          ! neutralizes the given temperature gradient
          sgrad = -alpha/beta * tgrad
          samb = saltbot + (depth-maxval(depth))*sgrad

       else

          samb = salttop + depth*sgrad

       end if

    else

       if (depth(namb) > amb_depth_ctl_pt(n_amb_ctl_pt)) then
          call io_append_output('Deepest control point is not deep enough')
          stop 1
       end if

       if (depth(namb) < amb_depth_ctl_pt(n_amb_ctl_pt-1)) then
	  call io_append_output('Deepest ambient layer should be between last two control points')
	  stop 1
       end if 

       if (depth(1) < amb_depth_ctl_pt(1)) then
          call io_append_output('Shallowest control point is not shallow enough')
          stop 1
       end if

       z_prev =  amb_depth_ctl_pt(n_amb_ctl_pt)
       t_prev =  amb_temp_ctl_pt(n_amb_ctl_pt)

       if (use_neutral_salinity) then
          s_prev = saltbot
       else
          s_prev =  amb_salt_ctl_pt(n_amb_ctl_pt)
       end if

       do i=namb,1,-1

          found_ctl_level = .false.

          do i_ctl=n_amb_ctl_pt,2,-1

             if (found_ctl_level) continue

             if (depth(i) <= amb_depth_ctl_pt(i_ctl) .and. &
                 depth(i) >= amb_depth_ctl_pt(i_ctl-1)) then
                
                found_ctl_level = .true.

                ctl_tgrad =  (amb_temp_ctl_pt(i_ctl-1)-amb_temp_ctl_pt(i_ctl)) / &
                             (amb_depth_ctl_pt(i_ctl-1)-amb_depth_ctl_pt(i_ctl))
		print *, i, ctl_tgrad

                if (use_neutral_salinity) then

                   ctl_sgrad = (-alpha/beta) * ctl_tgrad

                else

                   ctl_sgrad = (amb_salt_ctl_pt(i_ctl-1)-amb_salt_ctl_pt(i_ctl)) / &
                               (amb_depth_ctl_pt(i_ctl-1)-amb_depth_ctl_pt(i_ctl))                   
                end if
                
                samb(i) = s_prev + (depth(i)-z_prev)*ctl_sgrad
                tamb(i) = t_prev + (depth(i)-z_prev)*ctl_tgrad

                z_prev = depth(i)
                t_prev = tamb(i)
                s_prev = samb(i)        

             end if
          end do

          if (.not.(found_ctl_level)) then
             call io_append_output('did not find ctl level')
             stop 1
          end if
       end do
    end if

    ttt = 0.d0

    if (rholinear) then
       rhovf = rho_func_linear(tamb,samb)              
       rhopot = rhovf
    else
       if (thermobar) then
          pressure = depth*1.0d-1
          ttt = tinsitu_func(tamb,samb,pressure)
       else
          pressure = 0.d0    
          ttt = tamb
       end if

       rhovf = rho_func_nonlinear(ttt,samb,pressure)              
       rhopot = rho_func_nonlinear(tamb,samb,0.d0)

    end if
	
    !if (.not. suppress_ascii_output)
    call io_write_amb(namb,depth,tamb,ttt,samb,rhovf,rhopot)

  end subroutine set_ambient

  subroutine frazil_calc()

    implicit none

    ! set up various parameter arrays use in the frazil calculations

    ! local variables

    integer :: l
    real(kind=kdp) :: cdc,winew,rey,mstar

    print *, 'frazil_calc() may require checking'
    stop 1

    do l=1,nice

       ! calculate effective radii

       re(l) = r(l)*(1.5d0*ar)**(1.d0/3.d0)

       ! calculate crystal thicknesses

       thi(l) = 2.d0*r(l)*ar

       ! calculate crystal volumes

       vol(l) = pi*r(l)*r(l)*thi(l)

       ! iteratively calculate frazil rising velocity

       cdc = 10.d0
       wi(l) = 0.d0
       winew = dsqrt((4.d0*grav*ar*r(l)*(rho0-rhoi))/(rho0*cdc))
       do while (abs(wi(l)-winew).gt.small) 
          wi(l) = winew
          rey = 2.d0*winew*r(l)/nu0
          cdc = 10**(1.386d0 - 0.892d0*dlog10(rey) &
               + 0.111d0*(dlog10(rey))**2)
          winew = dsqrt((4.d0*grav*ar*r(l)*(rho0-rhoi))/(rho0*cdc))
       end do
       wi(l) = winew

       ! set the minimum of each size class to be one crystal

       cmin(l) = vol(l)

       ! set the nusselt number for each size class according to option chosen

       ! if nusselt number set to positive, use that constant value
       nuss(l) = nus

       ! correct full variable formulation
       if (nus.eq.-1.d0) then
          write(*,*) 'error: nusselt number option not coded yet'
       end if

       ! correct variable formulation with no turbulent part for large crystals
       if (nus.eq.-2.d0) then
          mstar = r(l)/(nu0**3.d0/eps)**2.5d-1
          if (mstar.lt.(1.d0/dsqrt(pr))) then
             nuss(l) = 1.d0 + mstar*1.7d-1*dsqrt(pr)
          else
             nuss(l) = 1.d0 + mstar*5.5d-1*(pr/mstar)**(1.d0/3.d0)
          end if
       end if

       ! incorrect full hammar & shen formulation
       if (nus.eq.-3.d0) then
          write(*,*) 'error: nusselt number option not coded yet'
       end if

       ! incorrect hammar & shen formulation with no turbulent part 
       ! for large crystals
       if (nus.eq.-4.d0) then
          mstar = r(l)/(nu0**3.d0/eps)**2.5d-1
          if (mstar.lt.(1.d0/dsqrt(pr))) then
             nuss(l) = 1.d0/mstar + 1.7d-1*dsqrt(pr)
          else
             nuss(l) = 1.d0/mstar + 5.5d-1*(pr/mstar)**(1.d0/3.d0)
          end if
       end if

    end do

  end subroutine frazil_calc

  !*****************************************************************
  !******** helper functions called inside plume_runstep() ********
  ! ****************************************************************


  subroutine subglacial_discharge(icalcan, kcalcan, icalcen, kcalcen)

    implicit none

    integer, intent(in) :: icalcan,kcalcan,icalcen,kcalcen

    ! local variables
    integer :: i,k,l
    real(kind=kdp):: fresh_water_column_thk, ttt, pressure
    real(kind=kdp):: discharge_area 

    real(kind=kdp),dimension(m_grid) :: fresh_water_column_thks

    print *, 'this code should not be used' 
    stop 1

    discharge_area = sum(dx(infloa+1:infloe-1))*sum(dy(knfloa+1:knfloe-1))

    if (sgd_type == 0) then

       ! uniform flux across inflow edge

       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
    
             if ((i .gt. infloa) .and. (i .lt. infloe) .and. &
                  (k .gt. knfloa) .and. (k .lt. knfloe)) then
                
                ! sgd_flux should be in km^3 per year
                fresh_water_column_thk = sgd_flux * (1.0d9)*(dt/(365.25d0*3600.d0*24.d0))/discharge_area
                salt(i,k) = salt(i,k)*pdep(i,k) / (pdep(i,k)+fresh_water_column_thk)
                salta(i,k) = salt(i,k)
                temp(i,k) = temp(i,k)*pdep(i,k) / (pdep(i,k)+fresh_water_column_thk)
                tempa(i,k) = temp(i,k)

                pdep(i,k) = pdep(i,k) + fresh_water_column_thk

                ! calculate density 
                if (rholinear) then
                   rhop(i,k) = rho_func_linear(temp(i,k),salt(i,k))
                else
                   if (thermobar) then
                      pressure = 1.0d-1*(gldep + wcdep - ipos(i,k))
                      ttt = tinsitu_func(temp(i,k),salt(i,k),pressure)
                      rhop(i,k) =  &
                           & rho_func_nonlinear(ttt,salt(i,k),pressure)
                   else
                      rhop(i,k) = &
                           & rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
                   end if
                end if

             end if

          end do
       end do
    else
       print *, 'Unknown sub-glacial discharge type', sgd_type
       stop 1

    end if

  end subroutine subglacial_discharge

  subroutine inflow_calc(icalcan,kcalcan,icalcen,kcalcen)

    implicit none

    integer, intent(in) :: icalcan,kcalcan,icalcen,kcalcen

    ! local variables
    integer :: i,k,l
    real(kind=kdp):: pressure,ttt

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
          if (depinf(i,k).gt.0.d0) then

             ! reset inflow depth

             ipos(i,k) = bpos(i,k) - depinf(i,k)         
             pdep(i,k) = depinf(i,k)

             if (pdep(i,k).ge.dcr) then           

                ! reset inflow properties 

                jcw(i,k) = 1
                salt(i,k) = saltinf(i,k)
                salta(i,k) = saltinf(i,k)
                temp(i,k) = tempinf(i,k)
                tempa(i,k) = tempinf(i,k)

                ! calculate density of inflow water fraction

                if (rholinear) then
                   rhop(i,k) = rho_func_linear(tempinf(i,k),saltinf(i,k))
                else
                   if (thermobar) then
                      pressure = 1.0d-1*(gldep + wcdep - ipos(i,k))
                      ttt = tinsitu_func(tempinf(i,k),saltinf(i,k),pressure)
                      rhop(i,k) =  &
                           & rho_func_nonlinear(ttt,saltinf(i,k),pressure)
                   else
                      rhop(i,k) = &
                           & rho_func_nonlinear(tempinf(i,k),saltinf(i,k),0.d0)
                   end if
                end if

                ! set frazil and add density of ice fraction

                if (frazil) then
                   do l = 1,nice
                      c_ice(i,k,l) = cinf(l)
                      ca_ice(i,k,l) = cinf(l)
                   end do
                   ctot(i,k) = cinftot
                   ctota(i,k) = cinftot
                   rhop(i,k) = (1.d0 - cinftot)*rhop(i,k) + cinftot*rhoi
                end if

                ! set inflow tracers

                if (intrace) then
                   do l = ninfmin,ninfmax
                      if (intrin(i,k).eq.l) then
                         intracer(i,k,l) = 1.d0
                         intracera(i,k,l) = 1.d0
                      else
                         intracer(i,k,l) = 0.d0
                         intracera(i,k,l) = 0.d0
                      end if
                   end do
                end if

             end if

          end if
       end do
    end do

  end subroutine inflow_calc

  subroutine continuity(icalcan,kcalcan,icalcen,kcalcen)      

    ! evaluate entrainment and continuity equation to find interface position

    implicit none

    integer,intent(in) :: icalcan,kcalcan,icalcen,kcalcen

    ! local variables
    integer :: i,k,l,mflag,depthflag

    real(kind=kdp),dimension(m_grid,n_grid) :: pdepc,bspeed,speed
    real(kind=kdp) :: rhopac,rhoa,delrho,rhoq,redg,tt,vmid,umid
    real(kind=kdp) :: rich,sm,arg, phase, extra_entr
    real(kind=kdp) :: dragrt,prden,scden,gambt,gambs,tfreezeb,tfreezei,c1,c2,c3
    real(kind=kdp) :: deltam(m_grid,n_grid),delta(m_grid,n_grid),iold(m_grid,n_grid)
    real(kind=kdp) :: fppntot,ucrit,ucl,ustar
    real(kind=kdp) :: amb_depth

    real(kind=kdp) :: Ap, Sp, cp1, cp3, c4
    real(kind=kdp) ::  mon_obu_stab, h_over_l, h_over_lp, lambda, delta_b
    real(kind=kdp) :: delta_b_lower,delta_b_upper,u_star,Bh,heat_flux
    real(kind=kdp) :: buoyancy_flux_heat
    real(kind=kdp) :: lp_over_l

    real(kind=kdp):: discharge_area 
    real(kind=kdp),parameter :: sgd_flux_k = 4.d0

    ! this may be used as a mechanism to generate a positive 
    ! minimum entrainment rate.  Change to a number > 0.d0
    real(kind=kdp),parameter :: bspeed_min = 0.d-2

    ! 0. preliminaries
    ! ----------------
    prden = 12.5d0*pr**(2.d0/3.d0) - 9.d0
    scden = 12.5d0*sc**(2.d0/3.d0) - 9.d0
    delta = 0.d0
    deltam = 0.d0
    bspeed = 0.d0
    speed = 0.d0
    pdepc = 0.d0

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
	  if (jcs(i,k) .ne. 1) cycle ! don't calculate anything over 
                                     ! land or grounded ice
          deltam(i,k) = 0.d0
          delta(i,k) = 0.d0
          iold(i,k) = 0.d0
          pdepc(i,k) = bpos(i,k) - ipos(i,k)

          ! find flow speed on scalar grid
          ! NB: we require that all speeds surrounding cells
          !     with positive thickness should be physically valid speeds

          if (any(jcd_u(i-1:i,k) == 0) .or. any(jcd_v(i,k-1:k) == 0)) then
              print *, 'using invalid velocity at', i,k
              stop 1
          end if
          tt = 5.0d-1*dy(k)*rdyv(k)
          vmid = tt*sv(i,k-1) + (1.d0-tt)*sv(i,k) 
          tt = 5.0d-1*dx(i-1)*rdxu(i)
          umid = tt*su(i,k) + (1.d0-tt)*su(i-1,k)
          speed(i,k) = umid**2 + vmid**2 + small
          bspeed(i,k) = sqrt(speed(i,k))
       end do
    end do

    ! 2. basal melting
    ! ----------------

    if (basmelt) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle ! skip land/ grounded ice

             if (pdepc(i,k).gt.mdepth) then         

	        !                dragrt = dsqrt(drag(i,k))
                ustar = sqrt(drag(i,k)*(bspeed(i,k)**2.d0 + &
                                        0.5d0*local_tidal_speed(i,k)**2.d0) + &
                             u_star_offset**2.d0)
                ! find turbulent exchange coefficients
                gambt = ustar/(2.12d0*dlog(ustar*pdepc(i,k)/nu0) + prden)
                gambs = ustar/(2.12d0*dlog(ustar*pdepc(i,k)/nu0) + scden)

                ! calculate freezing point of plume at shelf base,
                ! decide if melt (mflag = 1) 
                ! or freeze and calculate freezing point of ice at shelf base
                tfreezeb = fta*salt(i,k) + ftb + ftc*(gldep + wcdep -bpos(i,k))
                mflag = (1 + int(sign(1.d0,temp(i,k) - tfreezeb)))/2
                depthflag = (1 + int(sign(1.d0, gldep + wcdep - bpos(i,k) - min_melt_depth)))/2
                tfreezei = (1 - mflag)*fta*si  &
                     + ftb + ftc*(gldep + wcdep - bpos(i,k))

                ! calculate coefficients in quadratic to be solved
                c1 = lat/c0 + mflag*(ci/c0)*(tfreezei-tint(i,k))
                c2 = gambs*(lat/c0 + mflag*(ci/c0)*(tfreezeb-tint(i,k)))  &
                     &   + gambt*(tfreezei-temp(i,k))
                c3 = gambs*gambt*(tfreezeb - temp(i,k))

                ! calculate melt rate
                if (separated(i,k) == 1) then
                   bmelt(i,k) = 0.d0
                else
                  bmelt(i,k) = -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1) * depthflag
                end if

                ! calculate basal temperature and salinity
                btemp(i,k) = (gambt*tempa(i,k)+mflag*(ci/c0)*bmelt(i,k)*tint(i,k) &
                     - (lat/c0)*bmelt(i,k) )/  &
                     (gambt + mflag*(ci/c0)*bmelt(i,k))
                
                bsalt(i,k) = (btemp(i,k) - ftb  &
                     - ftc*(gldep + wcdep - bpos(i,k)))/fta

                deltam(i,k) = deltam(i,k) + bmelt(i,k)

             end if
          end do
       end do
    end if
 

    ! 1. entrainment
    ! --------------

    if (entrain) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle ! don't calculate anything over land

             if(pdepc(i,k).gt.edepth) then          

                rhopac = rhop(i,k)
                rhoa = rhoamb(i,k)
                delrho = rhoa - rhopac
                rhoq = 5.0d-1*(rhopac+rhoa)
                redg = delrho/rhoq
                rich = grav*redg*pdepc(i,k)/speed(i,k)
                rich = dmax1(5.0d-2,rich)    

		amb_depth = wcdep + gldep - ipos(i,k)
		atemp(i,k) = get_tamb_z(amb_depth)
		asalt(i,k) = get_samb_z(amb_depth)

                ! full kochergin entrainment
                if (entype.eq.1) then
                   sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
                        &             3.16d-1*rich + 3.46d-2)))
                   arg = dmax1(small,speed(i,k) + grav*redg*pdepc(i,k)/sm)
                   entr(i,k) = cl**2*dsqrt(arg)/sm                
                end if

                ! reduced kochergin entrainment matched to pedersen
                if (entype.eq.2) then
                   sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
                        &             3.16d-1*rich + 3.46d-2)))
                   arg = dmax1(small,speed(i,k) + grav*redg*pdepc(i,k)/sm)
                   entr(i,k) = ef*cl**2*(drag(i,k)/3.0d-3)*dsqrt(arg)/sm    
                end if

                ! half pedersen entrainment
                if (entype.eq.3) then
                   entr(i,k) = 3.6d-2*(bspeed(i,k)+bspeed_min)*drag(i,k)/rich
                end if

                ! half modified pedersen entrainment
                if (entype.eq.4) then
                   entr(i,k) = 3.6d-2*(bspeed(i,k)+bspeed_min)*dsin(1.0d-3)
                end if

                if (entype.eq.5) then
                   ! This is the entrainment scheme described in Zilitinkevich and
                   ! Mironv's 1996 paper where the mixed layer thickness is 
                   ! diagnosed directly using local vertical variables only.  The 
                   ! thickness is adjusted a posteriori (at the end of the current
                   ! time step), so we take entrainment to be zero here and make the
                   ! appropriate adjustments to plume depth, transport, heat, and salt
                   ! later.
                   entr(i,k) = 0.d0 
                end if

                if (entype.eq.6) then
                   ! Gaspar CMO 
                   u_star = sqrt(drag(i,k)*(bspeed(i,k)**2.d0 + &
                                            0.5d0*local_tidal_speed(i,k)**2.d0) +  &
                                 u_star_offset**2.d0)
                   lambda = u_star / f
                   delta_b = grav*(- alpha*(btemp(i,k)-temp(i,k)) &
                                    - beta*(bsalt(i,k)-salt(i,k)) )
                   Bh = bmelt(i,k)*delta_b
                   mon_obu_stab = min(pdep(i,k)*Bh/(u_star ** 3.d0), gaspar_cutoff)
!                   debug3(i,k) = mon_obu_stab
                   h_over_l = a1 + a2*max(1.d0,pdep(i,k)/(0.4d0*lambda))* &
                                      exp(mon_obu_stab)
                   h_over_lp= a1 + a2*exp(mon_obu_stab)
                   lp_over_l = h_over_l/h_over_lp
                   cp3 = (gasp_m4*(gasp_m2+gasp_m3)-(lp_over_l) * &
                         (gasp_m2+gasp_m3-gasp_m5*gasp_m3))/3.d0
                   cp1 = ((2-2*gasp_m5)*(lp_over_l)+gasp_m4)/6.d0
                   Ap = cp3*u_star**3.d0-cp1*pdep(i,k)*Bh

                   if (u_star == 0.d0) then
                      print *, 'unexpectedly found u_star == 0.d0',i,k
                      stop 1
                   end if

                   if (Ap > 0.d0) then
                      !we are entraining
                      Sp = (gasp_m2+gasp_m3)*(u_star ** 3.d0)-0.5d0*pdep(i,k)*Bh
                      c4 = 2.d0*gasp_m4/(gasp_m1*gasp_m1)                   
                      delta_b_lower = &
                           grav*(- alpha*(temp(i,k)-atemp(i,k)) &
                                 - beta *(salt(i,k)-asalt(i,k)) )

                      entr(i,k) = (1.d0/(pdep(i,k)*delta_b_lower))* & 
                            (-(0.5d0*Ap+cp1*Sp)+ &
                             sqrt((0.5d0*Ap-cp1*Sp)**2.d0 + &
                               2.d0*c4*(h_over_l*h_over_l)*Ap*Sp)) / &
                            (c4*(h_over_l)**2.d0-cp1)

                   else
                      ! we are detraining, and will adjust the 
                      ! plume thickness at the end of the timestep
                      entr(i,k) = 0.d0

                   end if

		end if

		if (entype == 7) then

		   ! Niiler-Kraus model, from Gaspar 1988
                   u_star = sqrt(drag(i,k)*(bspeed(i,k)**2.d0 + &
                                            0.5d0*local_tidal_speed(i,k)**2.d0) + &
                                 u_star_offset**2.d0) + small
                   delta_b_upper = grav*(- alpha*(btemp(i,k)-temp(i,k)) &
                                         -  beta*(bsalt(i,k)-salt(i,k)) )
                   delta_b_lower = grav*(- alpha*(temp(i,k)-atemp(i,k)) &
                                         - beta*(salt(i,k)-asalt(i,k)) )

                   !when melting at the ice interface, we are creating two
	           ! buoyancy fluxes:
	           ! 1. The flux of heat from the mixed layer into the 
	           !    interface to warm the ice to the melting point and
	           !    to supply the latent heat of fusion needed to melt 
	           !    the ice.
	           ! 2. The meltwater which is immediately mixed down into 
	           !    the mixed layer is a buoyancy flux.

		   heat_flux = (rhoi*lat + rhoi*ci*(btemp(i,k)-tint(i,k)))
		   heat_flux = (rhoi*ci*(btemp(i,k)-tint(i,k)))
	           buoyancy_flux_heat = -alpha*heat_flux/(c0*rho0)
!                   debug3(i,k) = buoyancy_flux_heat/delta_b_upper
                   Bh = bmelt(i,k)*(-buoyancy_flux_heat + &
                                    delta_b_upper)
	
		   entr(i,k) = (1.d0/(pdep(i,k)*delta_b_lower)) * &
		       (2.d0*nk_m*u_star**3.d0 &
                         -0.5d0*pdep(i,k)*((1-nk_n)*abs(Bh)+(1+nk_n)*Bh))

		   if (entr(i,k) < 0.d0) then
			entr(i,k) = 0.d0
                   end if

                end if

             end if

             deltam(i,k) = deltam(i,k) + entr(i,k)

          end do
       end do
    end if

    !-----2.5  subglacial discharge------

    if (sgd_type == 0 .or. sgd_type == 1) then

       ! uniform flux across inflow edge

       discharge_area = sum(dx(infloa+1:infloe-1))*sum(dy(knfloa+1:knfloe-1))

       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
    
             if ((i .gt. infloa) .and. (i .lt. infloe) .and. &
                  (k .gt. knfloa) .and. (k .lt. knfloe)) then
                
                ! sgd_flux should be in km^3 per year
		sgd(i,k) = sgd_flux * (1.0d9)*(1.d0/(365.25d0*3600.d0*24.d0))/discharge_area

		if (sgd_type == 1) then

		   !Note: The idea of approximating a train of spikes using a truncated cosine series
	           !      doesn't seem to work well.  Going back to a much gentler x-variation
	           sgd(i,k) = sgd(i,k)*(  1.d0  \
	                                - 0.5d0*cos(2.d0*pi*1.d0*sgd_flux_k*real(i-infloa-1)/real(infloe-1-infloa-1)) )
	 

	        end if 

		deltam(i,k) = deltam(i,k) + sgd(i,k)

	     end if

          end do
       end do
    end if


    ! 3.frazil precipitation
    ! ----------------------

    if (frazil) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle !skip land

             if (pdepc(i,k).gt.fdepth) then         

                fppntot = 0.d0
                do l=1,nice
                   if (c_ice(i,k,l).gt.0.d0) then

                      ! calculate precipitation of frazil
                      ucl = &
                           dsqrt((1.0d-1*grav*re(l)*(rho0-rhoi))/(rho0*drag(i,k)))
                      ucrit = (1.d0 - bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                      fppn(i,k,l) = &
                           dmin1(-(rhoi/rho0)*c_ice(i,k,l)*wi(l)*ucrit,0.d0)

                      fppntot = fppntot + fppn(i,k,l)
                   end if
                end do

                deltam(i,k) = deltam(i,k) + fppntot

             end if
          end do
       end do
    end if

    ! final divergence plus entrainment and melting
    ! ---------------------------------------------
    !    
    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

	  if (jcs(i,k) .ne. 1) cycle

          delta(i,k) = dt*(rdxu(i)*(utrans(i-1,k)-utrans(i,k)) &
                         + rdyv(k)*(vtrans(i,k-1)-vtrans(i,k))) &
                     + dt* deltam(i,k)

          if (use_min_plume_thickness .and. .not. (entype == 5)) then

            thk_def(i,k) = max(0.d0,plume_min_thickness - (pdepc(i,k)+delta(i,k)))
            if (thk_def(i,k) > 0.d0) then

                ! make up for thickness deficiency by increasing entrainment
                ! articicially

                entr(i,k)  = entr(i,k)  +  thk_def(i,k) / entr_time_const  
	        delta(i,k) = delta(i,k) +  thk_def(i,k) / entr_time_const * dt
	              
!                if (entype.eq.5  .or.  entype.eq.6) then
!                   ! set entrainment back to 0, since we don't want to set it here
!                   entr(i,k) = 0.d0
!                end if

            end if

             !this is for output purposes so we can see what percentage of the thicknes
             !change was due to the imposed minimum thickness
            if (entr(i,k) .ne. 0.d0) then
              artf_entr_frac(i,k) = (thk_def(i,k)/entr_time_const) / entr(i,k)
            end if
          end if

          if (use_periodic_forcing .and. &
               (k .le. klast) .and. (k .ge. kfirst)) then

             print *, 'use_periodic_forcing ought to be checked'
             stop 1
            phase = 2.d0*pi*runtim/forcing_period

            !produce a rather sharp pulse of extra entrainment each cycle
            extra_entr = periodic_forcing_amp * exp(-real((1.d0 - cos(phase))**0.25d0))

            !a pulse that joins continuously with the surrounding values
            entr(i,k) = entr(i,k) + extra_entr*cos(pi*(k-(kfirst+klast)/2.d0)/(klast-kfirst))
            delta(i,k) =  delta(i,k) +  extra_entr * dt

          end if

       end do
   end do

   ! update interface position 

   iold = ipos

   delta = min(delta, iold)
   ipos = ipos - delta

   ! check for negative predicted depth of ambient fluid 
   jcd_negdep = 0
   where (bpos < (iold - 2.d0*delta))
       jcd_negdep = 2
   elsewhere (bpos < (iold - 3.d0*delta))
       jcd_negdep = 1
   end where

  end subroutine continuity

  subroutine ZM_1996_thickness(i,k,h,h_formula, h_f, h_B)

    integer,intent(in) :: i,k
    real(kind=kdp),intent(out) :: h,h_formula,h_f,h_B

    ! local variables
    real(kind=kdp) :: speed, u_star, N, B_s
    real(kind=kdp) :: a,b,c

    !Solving the equation
    !  (fh/ (C_n*u_star) )^2 + (h * B_s / (C_s * u_star ^3)) 
    !       + (Nh/(C_i * u_star)) = 1

    ! where f is Coriolis param
    !       h is the mixed layer thickness to be determined
    !       u_star = sqrt( drag stress )
    !       B_s = buoyancy source = g*m' * (alpha*T_b + beta*S_b)
    !         N = buoyancy freq = sqrt(-g*(rho_plume-rho_amb)/(rho0*plume_thickness))

    speed = sqrt( ((su(i-1,k)+su(i,k))**2.d0)/4.d0 + ((sv(i,k-1)+sv(i,k))**2.d0)/4.d0 + small)
    u_star = sqrt(drag(i,k)*(speed**2.d0+0.5d0*local_tidal_speed(i,k)**2.d0) + &
                  u_star_offset**2.d0)

!    debug3(i,k) = u_star

    if (rhop(i,k) < rhoamb(i,k)) then
       N = sqrt(-grav*(rhop(i,k)-rhoamb(i,k))/(rho0*pdep(i,k)))
    else
       N = 0.d0
    end if

    if (bmelt(i,k) > 0.d0) then
       B_s = grav * bmelt(i,k)* (alpha*btemp(i,k) + beta*bsalt(i,k) )
    else
       B_s = 0.d0
    end if

!    print *, C_s, C_n, C_i
    h_f = C_n*u_star/f
    h_B = (B_s/(C_s*u_star**3) + N / (C_i * u_star)) ** (-1.d0)

    a = (f/(C_n)) ** 2.d0
    b = (B_s/(C_s * u_star) + (N*u_star)/(C_i))
    c = -u_star**2.d0

    a = 1.d0/(h_f*h_f)
    b = 1.d0/h_B
    c = -1.d0

    if (abs(f) > small) then
       h = (-b + sqrt(b**2.d0 - 4.d0*a*c))/(2.d0*a)
    else
       if (b > 0.d0) then
          h = -c/b
       else
          h = 0.d0
       end if

    end if
    h_formula = h
!    print *, h_f, h_B, h_formula
!    stop 1
    ! ensure that the layer has at least the minimum thickness
    h = max(h, plume_min_thickness)
!    h = min(h, 500.d0)

  end subroutine ZM_1996_thickness

  subroutine entrainment_correction(icalcan, icalcen, kcalcan, kcalcen, entype)

    implicit none

    integer, intent(inout) :: icalcan,kcalcan,icalcen,kcalcen
    integer, intent(in) :: entype

    !local variables
    integer :: i,k
    real(kind=kdp) :: ZM_h, delta_thk, pressure, ttt
    real(kind=kdp) :: utrain, vtrain, pdepu, pdepv, train_avg
    real(kind=kdp) :: h_f, h_B, h_formula
    real(kind=kdp) :: detrain_thk
    real(kind=kdp) :: Ap, Sp, cp1, cp3, c4
    real(kind=kdp) :: mon_obu_stab, h_over_l, h_over_lp, lambda, u_star, Bh
    real(kind=kdp) :: tt, vmid, umid, speed, bspeed
    real(kind=kdp) :: delta_b_upper,delta_b_lower,heat_flux,buoyancy_flux_heat
    do i=icalcan,icalcen
       do k=kcalcan, kcalcen

             ! calculate the correct mixed-layer thickness according to
             ! Zilitinkevich and Mironov 1996, and adjust the entrainment
             ! rate to produce the required thickness.

             if (jcs(i,k) == 0) cycle

             if (entype == 5) then

                call ZM_1996_thickness(i,k,ZM_h, h_formula, h_f, h_B)
                !debug(i,k)  = h_formula
                debug2(i,k) = h_B
!                debug3(i,k) = h_f
                thk_def(i,k) = ZM_h - pdep(i,k)

             else if (entype == 6) then

	         ! doing Gaspar detrainment
                 tt = 5.0d-1*dy(k)*rdyv(k)
                 vmid = tt*sv(i,k-1) + (1.d0-tt)*sv(i,k) 
                 tt = 5.0d-1*dx(i-1)*rdxu(i)
                 umid = tt*su(i,k) + (1.d0-tt)*su(i-1,k)
                 speed = umid**2 + vmid**2 + small
                 bspeed = sqrt(speed)
                 u_star = sqrt(drag(i,k)*(bspeed**2.d0 + &
                                          0.5d0*local_tidal_speed(i,k)**2.d0) + &
                               u_star_offset**2.d0)
                 lambda = u_star / f
                 delta_b_upper = grav*(- alpha*(btemp(i,k)-temp(i,k)) &
                                       - beta *(bsalt(i,k)-salt(i,k)) )
                 Bh = bmelt(i,k)*delta_b_upper
                 mon_obu_stab = min(pdep(i,k)*Bh/(u_star ** 3.d0),gaspar_cutoff)
                 h_over_l = a1 + a2*max(1.d0,pdep(i,k)/(0.4*lambda))* &
                                    exp(mon_obu_stab)
                 h_over_lp= a1 + a2*exp(mon_obu_stab)
                 cp3 = (gasp_m4*(gasp_m2+gasp_m3)-(h_over_l/h_over_lp)* &
                       (gasp_m2+gasp_m3-gasp_m5*gasp_m3))/3.d0
                 cp1 = ((2-2*gasp_m5)*(h_over_l/h_over_lp)+gasp_m4)/6.d0
                 ! Gaspar thickness comes from Ap = 0
                 if (cp1*Bh == 0.d0) then
                    !print *, 'that is surprising'
                    !stop 1
                    detrain_thk = plume_max_thickness
                 else
                    if (use_min_plume_thickness) then
                       detrain_thk = max(plume_min_thickness, &
                            (cp3*u_star**3.d0)/(cp1*Bh))
                    else
                       detrain_thk = cp3*u_star**3.d0/(cp1*Bh)
                    end if
                 end if
                   
                 ! if (detrain_thk > pdep) then we don't need to detrain,
                 ! so we force thk_def to be a negative quantity, ie it only
                 ! has an effect when we are detraining
                 thk_def(i,k) = min(0.d0, detrain_thk - pdep(i,k))
               
             else if (entype == 7) then
	         ! Niiler-Kraus model, from Gaspar 1988

                 tt = 5.0d-1*dy(k)*rdyv(k)
                 vmid = tt*sv(i,k-1) + (1.d0-tt)*sv(i,k) 
                 tt = 5.0d-1*dx(i-1)*rdxu(i)
                 umid = tt*su(i,k) + (1.d0-tt)*su(i-1,k)
                 speed = umid**2 + vmid**2 + small
	         u_star = sqrt(drag(i,k) * (speed + 0.5d0*local_tidal_speed(i,k)**2.d0) + &
                               u_star_offset**2.d0)
                 delta_b_upper = grav*(-alpha*(btemp(i,k)-temp(i,k)) &
                                       - beta*(bsalt(i,k)-salt(i,k)) )
                 delta_b_lower = grav*(-alpha*(temp(i,k)-atemp(i,k)) &
                                       - beta*(salt(i,k)-asalt(i,k)) )
                 heat_flux = (rhoi*lat + rhoi*ci*(btemp(i,k)-tint(i,k)))
	         buoyancy_flux_heat = -alpha*heat_flux/(c0*rho0)
                 Bh = bmelt(i,k)*(-buoyancy_flux_heat + &
                                   delta_b_upper)

                 if (Bh .le. 0.d0) then
                    detrain_thk = plume_max_thickness
                 else
                    detrain_thk = 2.d0*nk_m*u_star**3.d0 / Bh 
                 end if

                 debug3(i,k) = detrain_thk

                 if (use_min_plume_thickness) then
		      detrain_thk = max(plume_min_thickness, &
	                                detrain_thk)
                 end if

		! we don't want to correct if the detrainment
		! thickness exceeds the current depth
 		thk_def(i,k) = min(0.d0, detrain_thk-pdep(i,k))

             else 
                print *, 'Cannot do thickness adjustment unless entype is 5,6 or 7'
                stop 1
             end if

             ! We use the known thickness change to infer the necessary 
             ! detrainment/entrainment rate

             if (global_time_step_count > 0) then
                train(i,k) = thk_def(i,k) / detrain_time_const
             end if

             ! adjust the plume thickness and interface position, 
             delta_thk = train(i,k) * dt
             pdep(i,k) = pdep(i,k) + delta_thk
             ipos(i,k) = bpos(i,k) - pdep(i,k)

             if (delta_thk > 0.d0) then

                ! entraining case

                ! update the average temperature,salt and density
                temp(i,k) = temp(i,k) + & 
                     delta_thk * (atemp(i,k) - temp(i,k))/(pdep(i,k)+dcr)
                salt(i,k) = salt(i,k) + &
                     delta_thk * (asalt(i,k) - salt(i,k))/(pdep(i,k)+dcr)

                if (rholinear) then
                   rhop(i,k) = rho_func_linear(temp(i,k),salt(i,k))
                else
                   if (thermobar) then
                      pressure = 1.0d-1*(wcdep + gldep - ipos(i,k))
                      ttt = tinsitu_func(temp(i,k),salt(i,k),pressure)
                      rhop(i,k) =rho_func_nonlinear(ttt,salt(i,k),pressure)
                   else
                      rhop(i,k) =  &
                           rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
                   end if
                end if

                tempa(i,k) = temp(i,k)
                salta(i,k) = salt(i,k)

             else ! (delta_thk < 0.0)

                ! detraining case

                ! temp and salt and rhop do not change since we are leaving behind a 
                ! layer of fluid with the current plume water properties

                ! Neither does the plume speed change.
                continue

             end if

          end do
       end do

       ! adjust the speed and transport variables after the thickness have all been 
       ! adjusted
       do i=icalcan,(icalcen-1)
          do k=kcalcan,(kcalcen-1)

             ! entrainment/detrainment rates on the u grid and v grid
             utrain = 0.5d0*(train(i,k) + train(i+1,k)) 
             vtrain = 0.5d0*(train(i,k) + train(i,k+1))

             ! new plume depths on the u-grid and v-grid
             pdepu = 0.5d0*(pdep(i,k) + pdep(i+1,k))
             pdepv = 0.5d0*(pdep(i,k) + pdep(i,k+1))

             if (all(jcs(i:(i+1),k) == 1))  then
                if (utrain > 0.d0) then
                   ! entraining at u-grid location
                   ! in entraining case, the introduced fluid has no momentum, so we do
                   ! not change the transport variables, but we need to update the
                   ! speeds since the thickness has increased
                   su(i,k) =   utrans(i,k) / pdepu
                else if (utrain < 0.d0) then
                   ! detraining at u-grid location
                   ! the transports change, since we are losing fluid
                   utrans(i,k) = utrans(i,k) * pdepu/(pdepu+abs(utrain*dt))
                   utransa(i,k) = utrans(i,k)
                else
                   ! do nothing
                end if
             end if

             if (all(jcs(i,k:(k+1)) == 1)) then
                if (vtrain > 0.d0) then
                   ! entraining at v-grid location
                   sv(i,k) =   vtrans(i,k) / pdepv
                else if (vtrain < 0.d0) then
                   ! detraining at v-grid location
                   vtrans(i,k) = vtrans(i,k) * pdepv/(pdepv+abs(vtrain*dt))
                   vtransa(i,k) = vtrans(i,k)
                else
                   !do nothing
                end if
             end if

          end do
       end do

  end subroutine entrainment_correction

  subroutine tide(icalcan,kcalcan,icalcen,kcalcen)

    implicit none

    integer, intent(in) :: icalcan,kcalcan,icalcen,kcalcen

    ! local variables
    integer :: i,k
    
    local_tidal_speed = tidal_velocity
!    do i = icalcan,icalcen
!       do k = kcalcan,kcalcen
!	   local_tidal_speed(i,k) = tidal_velocity * (wcdep+gldep)/bpos(i,k)
!       end do
!    end do

  end subroutine tide

  subroutine momentum(icalcan,kcalcan,icalcen,kcalcen)

    ! calculates velocity components within plume

    ! all turning angle theory by alex wilchinsky
    ! (wilchinsky, feltham and holland, submitted to j. phys. oceanogr.)
    ! it is probably wise to check or reprogram this before use

    implicit none

    integer, intent(in) :: icalcan,kcalcan,icalcen,kcalcen

    ! local variables

    integer :: i,k,jcvfac
    integer(kind=kdp) :: idel,kdel,iidx,kkdy,ihilf,khilf

    real(kind=kdp) :: one,termnl,termnl2,corx
    real(kind=kdp) :: redgu,slorho,islope,detrain
    real(kind=kdp) :: pdepu,zu,zum,arfac,tt,uu
    real(kind=kdp) :: vmid,umid
    real(kind=kdp) :: speed,tbotm,rhoa,tu,salu
    real(kind=kdp) :: rhoc,rhoe,rhou,rhoq,dxx,dyy,sx,sy,sxy,r1,r2
    real(kind=kdp) :: tlate,tlatw,tlats,tlatn,hordif,cory,redgv
    real(kind=kdp) :: pdepv,zv,zvm,tv,salv,rhon,rhov
    real(kind=kdp) :: delta,ctotu,ctotv,dragu,dragv

    logical :: skip_u_calc,skip_v_calc
    logical :: olddrag,norotation,variableekman,draginmelt 
    logical :: newudrag(m_grid,n_grid)
    real(kind=kdp) :: av,ekthick,kq,costang,sintang,thickratio
    real(kind=kdp) :: ugriddrag(m_grid,n_grid)

    sintang = 0.0
    kq = 0.0

    ! set switches for turning angle model

    olddrag       = .true. ! use the old drag magnitude
    norotation    = .false. ! switch off drag rotation
    variableekman = .false. ! use variable vertical eddy visc (hence ek thick)
    draginmelt    = .false.  ! write drag array so that drag used in melting
    ! (only works with spatially constant drag coeff)

    one = 1.d0

    ! calculate ekman layer thickness for constant vertical eddy diffusivity

    if (.not.variableekman) then
       av = 5.d-4
       ekthick = dsqrt(2.d0*av/abs(f)) 
    end if

    ! start main loop

    debug = 0.d0
    debug2 = 0.d0
    debug3 = 0.d0

    train = 0.d0

    !$omp parallel default(none) &
    !$omp private(costang,thickratio, &
    !$omp         skip_u_calc,skip_v_calc, &
    !$omp         pdepv,zv,zvm,tv,salv,rhon,rhov, &
    !$omp         delta,ctotu,ctotv,dragu,dragv, &
    !$omp         termnl,termnl2,corx, &
    !$omp         redgu,slorho,islope,pdepu,zu,zum,arfac,tt,uu, &
    !$omp         vmid,umid,speed,tbotm,rhoa,tu,salu,&
    !$omp         rhoc,rhoe,rhou,rhoq,dxx,dyy,sx,sy,sxy,r1,r2, &
    !$omp         tlate,tlatw,tlats,tlatn,hordif,cory,redgv, &
    !$omp         i,k,jcvfac,idel,kdel,iidx,kkdy,ihilf,khilf) &
    !$omp shared( bpos, pdep, ipos, jcw, jcd_negdep, su,sv,drag,ugriddrag,newudrag, &
    !$omp         u0,u0a,v0,v0a,tang,salt,temp,tins,ctot, detrain, &
    !$omp         jcd_u,jcd_v,utrans,vtrans,utransa,vtransa,jcs,gwave_speed, &
    !$omp 	  icalcen,icalcan,kcalcen,kcalcan , &
    !$omp         dx,dxu,rdx,rdxu,dy,dyv,rdy,rdyv,ahdx,ahdxu,ahdy,ahdyv, &
    !$omp         wcdep,gldep,f,cdb,tangle,dcr,dt,small) &
    !$omp firstprivate ( sintang, kq, olddrag, norotation, &
    !$omp          variableekman, draginmelt, one,nonlin,fdt,gdt, &
    !$omp          av,ekthick,horturb,rhoi,frazil,thermobar,rholinear)
    !$omp do

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

          skip_u_calc = .false.
          skip_v_calc = .false.

          !*******************************************************************
          !*** u-component  (eastern - cross-slope) **************************
          !*******************************************************************

          termnl = 0.d0
          termnl2 = 0.d0
          hordif = 0.d0
          corx = 0.d0
          slorho = 0.d0
          islope = 0.d0
          redgu = 0.d0  
          detrain = 0.d0

          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! skip u-component if western or eastern cell is dry land
          ! NB: this is a no normal flow condition 
          if (any(jcs(i:i+1,k) == 0)) then
             skip_u_calc = .true.
          end if

          ! skipping u calculation if northern/southern cells are dry land
          ! would be a no-slip condition
          if (any(jcs(i:i+1,k+1) == 0) .or. any(jcs(i:i+1,k-1) == 0)) then
             skip_u_calc = .true.
          end if

	  if (i == icalcen) then
             ! the last column of u grid is not part of the problem
             !NB: this is only relevant in the case of an open boundary
             ! on the eastern side of the domain
	      skip_u_calc = .true.
          end if

          ! plume thickness on the u-grid
          pdepu = 5.0d-1*(pdep(i,k) + pdep(i+1,k))
          pdepu = 5.0d-1*(bpos(i,k)-ipos(i,k) + &
                          bpos(i+1,k)-ipos(i+1,k)) ! to get the new pdep?
          
          islope =  ipos(i+1,k) - ipos(i,k)  ! NB: this is the new ipos

          ! final wet/dry logic
          if ((jcw(i,k) == 0).and.(jcw(i+1,k) == 0)) then
             ! western cell and eastern neighbour are dry
             skip_u_calc = .true.
          else if ((jcw(i,k).le.0).and.(islope.gt.0.d0)) then
             ! current cell is dry and interface slopes up to east
             skip_u_calc = .true.
          else if ((jcw(i+1,k).le.0).and.(islope.le.0.d0)) then
             ! eastern neighbour is dry and interface slopes up to the west
             skip_u_calc = .true.
          endif

          if (skip_u_calc) then
      !       print *, 'skipped u ', i,k
          else
             ! control of negative depth effects by enhanced friction
             jcvfac = 0
             arfac = 1.d0
             jcvfac = max0(jcd_negdep(i,k),jcd_negdep(i+1,k))
             if (jcvfac .ge. 1) arfac = 75.d0*jcvfac  

             ! velocity components on the u-grid
             tt = 5.0d-1 
             uu = 5.0d-1*dy(k)*rdyv(k)
             if (any(jcd_v(i:i+1,k-1:k)==0)) then
                print *, 'using invalid sv at i,k', i,k
                stop 1
             end if
             vmid = tt*uu       *(sv(i,k-1) + sv(i+1,k-1)) +  &
                    tt*(1.d0-uu)*(sv(i,k)   + sv(i+1,k))
             umid = su(i,k)
           
             speed = umid**2+vmid**2 + small

             ! drag coefficient on the u-grid
             dragu = 5.0d-1*(drag(i,k) + drag(i+1,k))

             ! 1)calculation of friction parameter
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

             tbotm = 1.d0/(1.d0+arfac*dragu*dsqrt(speed)*dt/(pdepu+dcr))

             ! turning angle theory

             if (tangle) then 
	     
		if (f == 0.d0) then
		   print *, "Can not use tangle code when f = 0"
		   stop 1
		end if

                ! set up stuff for rewriting drag if necessary

                if (draginmelt) then
                   dragu = cdb
                   ugriddrag(i,k) = 5.0d-1*(drag(i,k) + drag(i+1,k))
                   newudrag(i,k) = .false.
                endif

                ! find ekman layer thickness if necessary and thickness ratio

                if (variableekman) then
                   av = abs(0.16d0*0.4d0*dragu*speed / (f*2.71828d0))
                   ekthick = dsqrt(2.d0*av/abs(f) )
                end if
                thickratio = 2.d0*pdepu/ekthick

                ! calculate velocities at base of ekman layer (eq. 32)

                kq = (ekthick/av)*dragu   &
                     * dsqrt(speed*((u0a(i,k) + one)**2 + v0a(i,k)**2))
                u0(i,k) = - (one + kq)*kq / ((one + kq)**2 + one) 
                v0(i,k) = - kq / ((one + kq)**2 + one)

                ! calculate geostrophic velocity magnitude (eq. 36)

                kq = thickratio / dsqrt( &
                     (thickratio + u0(i,k) - v0(i,k))**2 &
                                + (u0(i,k) + v0(i,k))**2  )  

                ! calculate final drag magnitude (eq. 42)

                kq = (kq*av/ekthick) * &
                     dsqrt( 2.d0*(u0(i,k)**2 + v0(i,k)**2) )
                !         
                if (olddrag) kq = dragu*dsqrt(speed) 
                !        
                if (draginmelt) then
                   if (speed > small) then
                      ugriddrag(i,k) = kq/dsqrt(speed)
                   else
                      ugriddrag(i,k) = cdb
                   endif
                   newudrag(i,k) = .true.
                endif
                !           
                ! calculate turning angle          
                !           
                tang(i,k) = atan2( u0(i,k)-v0(i,k) ,           -u0(i,k)-v0(i,k) ) &
                          - atan2( u0(i,k)+v0(i,k) , thickratio+u0(i,k)-v0(i,k) )
                !      
                if (norotation) tang(i,k) = 0.0 
                !         
                costang = cos(tang(i,k))
                sintang = sin(tang(i,k))
                tbotm = one / (one + arfac*kq*costang*dt / (pdepu+dcr))

             end if

             ! 2)baroclinic pressure gradient
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

             zu = wcdep + gldep - 5.0d-1*(ipos(i,k) + ipos(i+1,k))
             zum = zu - 5.0d-1*pdepu      

	     rhoa  = get_rhoamb_z(zu)

      ! density of water fraction
             salu = 5.0d-1*(salt(i,k)+salt(i+1,k))
             if (rholinear) then
                tu = 5.0d-1*(temp(i,k)+temp(i+1,k))
                rhoc = rho_func_linear(temp(i,k),salt(i,k))
                rhoe = rho_func_linear(temp(i+1,k),salt(i+1,k))
                rhou = rho_func_linear(tu,salu)
             else
                if (thermobar) then
                   tu = 5.0d-1*(tins(i,k)+tins(i+1,k))
                   rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),zum*1.0d-1)
                   rhoe = rho_func_nonlinear(temp(i+1,k),salt(i+1,k),zum*1.d-1)
                   rhou = rho_func_nonlinear(tu,salu,zu*1.0d-1)
                else
                   tu = 5.0d-1*(temp(i,k)+temp(i+1,k))
                   rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
                   rhoe = rho_func_nonlinear(temp(i+1,k),salt(i+1,k),0.d0)
                   rhou = rho_func_nonlinear(tu,salu,0.d0)
                end if
             end if

             ! effects of frazil
             if (frazil) then
                ctotu = 5.0d-1*(ctot(i,k) + ctot(i+1,k))
                rhoc = (1.d0 - ctot(i,k))*rhoc + ctot(i,k)*rhoi
                rhoe = (1.d0 - ctot(i+1,k))*rhoe + ctot(i+1,k)*rhoi
                rhou = (1.d0 - ctotu)*rhou + ctotu*rhoi  
             end if

             rhoq = 5.0d-1*(rhou + rhoa)
             redgu = (rhoa - rhou)/rhoq

	     ! final sloping isopycnal term contribution
             slorho =5.0d-1*(rhoe-rhoc)*gdt*rdx(i)*pdepu*pdepu/rhoq

             ! 3)barotropic pressure gradient
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             islope = islope*rdx(i)*gdt*redgu*pdepu

             ! 4)coriolis
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             corx = pdepu*vmid*fdt

             ! if turning angle add other part of drag here 
             ! because uses v not u

             if (tangle) then
                corx = corx + arfac*kq*sintang*dt*vmid
             end if

             ! 5&6)nonlinear terms (selective vector upstream)
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             if (nonlin) then

                !  choose upstream quarter
                idel = - int(sign(one,umid))
                kdel = - int(sign(one,vmid))
                ihilf = i + idel
                khilf = k + kdel

                !  choose the relevant upstream cell dimensions 
                iidx = i + (idel+1)/2
                kkdy = k + (kdel-1)/2
                dxx = dxu(iidx)  
                dyy = dy(kkdy)
                !  weighting coefficients
                sx = abs(umid)*dt/dxx
                sy = abs(vmid)*dt/dyy

                if ((sx.gt.one).or.(sy.gt.one)) then
                   write(*,*) 'error: u courant exceeded at i=',i,' k=',k
                   write(11,*) 'error: u courant exceeded at i=',i,' k=',k
                end if

!		if ((abs(umid) + gwave_speed(i,k)) > dxx/dt) then
!                  write(*,*) 'error: u courant (w/ grv waves)at i=',i,' k=',k
!                  write(11,*) 'error: u courant (w/ grv waves) at i=',i,' k=',k
	!	end if

                sxy = sx -sy
                r1 = (sign(one,sxy) + 1.d0)*5.0d-1
                r2 = one - r1

                ! NB: (but check this)
                ! In cases where we are calculating u(i,k) in a cell that 
                ! touches an east-west running wall, we will have ternml = 0.d0 
                ! in cases where the flow is away from the wall in the first
                ! non-zero row.

                ! termnl is the advection term -(su*u_x + sv*u_y) treated 
                ! using an upwind scheme
                termnl=(r1*sy + r2*sx) &
                      *(utransa(ihilf,khilf)*dble(jcd_u(ihilf,khilf)) &
                     +  utransa(i,k)        *dble(1 - jcd_u(ihilf,khilf))) &
                     + r1*sxy &
                      *(utransa(ihilf,k)    *dble(jcd_u(ihilf,k)) &
                     +  utransa(i,k)        *dble(1 - jcd_u(ihilf,k))) &
                     - r2*sxy &
                      *(utransa(i,khilf)    *dble(jcd_u(i,khilf)) &
                     +  utransa(i,k)        *dble(1 - jcd_u(i,khilf))) &
                     - (r2*sy + r1*sx)* utransa(i,k)

                ! termnl2 is the divergence term -u(su_x + sv_y)
                termnl2 = - utransa(i,k)*dt &
                         *((su(ihilf,k) - su(i,k))*idel/dxx  &
                         + (sv(iidx,k) - sv(iidx,k-1))*rdyv(k))

             end if
             !          
             ! 7)lateral shear stress terms (east component)
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             if (horturb) then

                tlate = ahdxu(i+1,k)*(su(i+1,k) - su(i,k))
                tlatw = ahdxu(i,k)*(su(i,k) - su(i-1,k))
                
                if (any(jcs(i:i+1,k+1) == 0)) then
                   ! boundary to the north
                   tlatn = 0.d0
                else
                   tlatn = ahdy(i,k)*(su(i,k+1) - su(i,k))
                end if

                if (any(jcs(i:i+1,k-1) == 0)) then
                   ! boundary to the south
                   tlats = 0.d0
                else
                   tlats = ahdy(i,k-1)*(su(i,k) - su(i,k-1))
                end if


                hordif = pdepu*((tlate-tlatw)*rdx(i) + (tlatn-tlats)*rdyv(k))

             end if
             

             ! final
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             utrans(i,k) = (utransa(i,k)+corx+islope+slorho+termnl+termnl2+hordif+detrain)*tbotm

          end if

         
          !*******************************************************************
          !*** v-component  (northern - up-slope) ****************************
          !*******************************************************************

          termnl = 0.d0
          termnl2 = 0.d0
          cory  = 0.d0
          slorho = 0.d0
          islope = 0.d0
          redgv = 0.d0
          detrain = 0.d0
 
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! skip v-component if southern or northern cell is dry land
          if (any(jcs(i,k:k+1) == 0)) then
	      skip_v_calc = .true.
	  end if

          ! skipping v calculation if eastern/western cells are dry land
          ! would be a no-slip condition
          if (any(jcs(i-1,k:k+1) == 0) .or. any(jcs(i+1,k:k+1) == 0)) then
             skip_v_calc = .true.
          end if

	  if (k == kcalcen) then
             ! last row of v is not part of the problem domain
             !NB: this is only relevant in the case of an open boundary
             ! on the northern side of the domain
             skip_v_calc = .true.
          end if

          islope = ipos(i,k+1) - ipos(i,k)
          ! plume thickness on the v-grid
          pdepv = 5.0d-1*(pdep(i,k) + pdep(i,k+1))
          pdepv = 5.0d-1*sum(bpos(i,k:k+1)-ipos(i,k:k+1)) ! to get the new pdep?

          ! final wet/dry logic
          if ((jcw(i,k).le.0).and.(jcw(i,k+1).le.0)) then
             ! southern cell and northern neighbour are dry
             skip_v_calc = .true.
          else if ((jcw(i,k).le.0).and.(islope.gt.0.d0)) then
	     ! southern cell is dry and interface is rising to the north
             skip_v_calc = .true.
          else if ((jcw(i,k+1).le.0).and.(islope.le.0.d0)) then
             ! northern neighbour is dry and interface is rising to the south
             skip_v_calc = .true.
          endif
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          if (.not. skip_v_calc) then

          ! control of negative depth effects by enhanced friction
          jcvfac = 0
          arfac = 1.d0
          jcvfac = max0(jcd_negdep(i,k),jcd_negdep(i,k+1))
          if (jcvfac .ge. 1) arfac = 75.d0*jcvfac  

          ! velocity components on the v-grid 
          uu = 5.0d-1*dx(i)*rdxu(i)
          tt = 5.0d-1
          umid = tt*uu       *(su(i-1,k) + su(i-1,k+1)) + &
                 tt*(1.d0-uu)*(su(i,k)   + su(i,k+1))
          vmid = sv(i,k)

          speed = umid**2+vmid**2 + small

          ! drag coefficient on the v-grid
          dragv = 5.0d-1*(drag(i,k) + drag(i,k+1))

          ! 1)calculation of friction parameter
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          tbotm = 1.d0/(1.d0+arfac*dragv*dsqrt(speed)*dt/(pdepv+dcr))

          ! turning angle theory

          if (tangle) then 

             if (draginmelt) dragv = cdb

             ! find ekman layer thickness if necessary and thickness ratio

             if (variableekman) then
                av = abs(0.16d0*0.4d0*dragv*speed / (f*2.71828d0))
                ekthick = dsqrt(-2.d0*av/f)
             end if
             thickratio = 2.d0*pdepv/ekthick

             ! calculate velocities at base of ekman layer (eq. 32)

             kq = (ekthick/av)*dragv &
                  * dsqrt(speed*((u0a(i,k) + one)**2 + v0a(i,k)**2))
             u0(i,k) = - (one + kq)*kq / ((one + kq)**2 + one) 
             v0(i,k) = - kq / ((one + kq)**2 + one)

             ! calculate geostrophic velocity magnitude (eq. 36)

             kq = thickratio / dsqrt( &
                  (thickratio + u0(i,k) - v0(i,k))**2    &
                  + (u0(i,k) + v0(i,k))**2)  

             ! calculate final drag magnitude (eq. 42)

             kq = (kq*av/ekthick) * &
                  dsqrt( 2.d0*(u0(i,k)**2 + v0(i,k)**2) )

             if (olddrag) kq = dragv*dsqrt(speed) 

             ! calculate turning angle          
             !           
             tang(i,k) = atan2( u0(i,k)-v0(i,k) , -u0(i,k)-v0(i,k) ) &
                       - atan2( u0(i,k)+v0(i,k) , thickratio+u0(i,k)-v0(i,k) )

             if (norotation) tang(i,k) = 0.0 

             costang = dcos(tang(i,k))
             sintang = dsin(tang(i,k))
             tbotm = one / (one + arfac*kq*costang*dt / (pdepv+dcr))

          end if

          ! 2)baroclinic pressure gradient
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          zv = wcdep + gldep - 5.0d-1*(ipos(i,k) + ipos(i,k+1))
          zvm = zv - 5.0d-1*pdepv     

	  rhoa = get_rhoamb_z(zv)

   ! density of water fraction
          salv = 5.0d-1*(salt(i,k)+salt(i,k+1))
          if (rholinear) then
             tv = 5.0d-1*(temp(i,k)+temp(i,k+1))
             rhoc = rho_func_linear(temp(i,k),salt(i,k))
             rhon = rho_func_linear(temp(i,k+1),salt(i,k+1))
             rhov = rho_func_linear(tv,salv)
          else
             if (thermobar) then
                tv = 5.0d-1*(tins(i,k)+tins(i,k+1))
                rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),zvm*1.0d-1)
                rhon = rho_func_nonlinear(temp(i,k+1),salt(i,k+1),zvm*1.d-1)
                rhov = rho_func_nonlinear(tv,salv,zv*1.0d-1)
             else
                tv = 5.0d-1*(temp(i,k)+temp(i,k+1))
                rhoc = rho_func_nonlinear(temp(i,k),salt(i,k),0.d0)
                rhon = rho_func_nonlinear(temp(i,k+1),salt(i,k+1),0.d0)
                rhov = rho_func_nonlinear(tv,salv,0.d0)
             end if
          end if

          ! effects of frazil
          if (frazil) then
             ctotv = 5.0d-1*(ctot(i,k) + ctot(i,k+1))
             rhoc = (1.d0 - ctot(i,k))*rhoc + ctot(i,k)*rhoi
             rhon = (1.d0 - ctot(i,k+1))*rhon + ctot(i,k+1)*rhoi
             rhov = (1.d0 - ctotv)*rhov + ctotv*rhoi  
          end if

          rhoq = 5.0d-1*(rhov + rhoa)
          redgv = (rhoa - rhov)/rhoq
          slorho = 5.0d-1*(rhon-rhoc)*gdt*rdy(k)*pdepv*pdepv/rhoq

          ! 3)barotropic pressure gradient
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          islope = islope*rdy(k)*gdt*redgv*pdepv

          ! 4)coriolis
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          cory =  - pdepv*umid*fdt

          ! if turning angle add other part of drag here because uses u not v

          if (tangle) then
             cory = cory - arfac*kq*sintang*dt*umid
          end if

          ! 5&6)nonlinear terms (selective vector upstream)
          ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if (nonlin) then

             !  choose upstream quarter
             idel = - int(sign(one,umid))
             kdel = - int(sign(one,vmid))
             ihilf = i + idel
             khilf = k + kdel
             !  choose the relevant upstream cell dimensions 
             iidx = i + (idel-1)/2
             kkdy = k + (kdel+1)/2
             ihilf = idel + i
             khilf = kdel + k
             dxx = dx(iidx)
             dyy = dyv(kkdy)
             ! weighting coefficients
             sx = abs(umid)*dt/dxx
             sy = abs(vmid)*dt/dyy
             ! test courant number
             if ((sx.gt.one).or.(sy.gt.one)) then
                write(*,*) 'error: v courant exceeded at i=',i,' k=',k
                write(11,*) 'error: v courant exceeded at i=',i,' k=',k
             end if
!	     if ((abs(vmid)+gwave_speed(i,k)) > dyy/dt) then
!                  write(*,*) 'error: v courant (w/ grv waves)at i=',i,' k=',k
!                  write(11,*) 'error: v courant (w/ grv waves) at i=',i,' k=',k
!		end if
             ! calculate nonlinear terms
             sxy = sx -sy
             r1 = (dsign(one,sxy) + 1.d0)*5.0d-1
             r2 = one - r1
             termnl=(r1*sy + r2*sx)*(vtransa(ihilf,khilf)*dble(jcd_v(ihilf,khilf)) &
                  + dble(1 - jcd_v(ihilf,khilf))*vtransa(i,k)) &
                  + r1*sxy*(vtransa(ihilf,k)*dble(jcd_v(ihilf,k)) &
                  + dble(1 - jcd_v(ihilf,k))*vtransa(i,k)) &
                  - r2*sxy*(vtransa(i,khilf)*dble(jcd_v(i,khilf)) &
                  + dble(1 - jcd_v(i,khilf))*vtransa(i,k)) &
                  - (r2*sy + r1*sx)*vtransa(i,k)

             termnl2 = - vtransa(i,k)*dt*((sv(i,khilf) - sv(i,k))*kdel/dyy  &
                                        + (su(i,kkdy) - su(i-1,kkdy))*rdxu(i))
          end if
          !          
          ! 7)lateral shear stress terms (north component)
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if (horturb) then

             if (any(jcs(i-1,k:k+1) == 0)) then
                ! boundary to the west
                tlatw = 0.d0
             else
                tlatw = ahdx(i-1,k)*(sv(i,k) - sv(i-1,k))
             end if

             if (any(jcs(i+1,k:k+1) == 0)) then
                ! boundary to the east
                tlate = 0.d0
             else
                tlate = ahdx(i,k)*(sv(i+1,k) - sv(i,k))
             end if

             tlats = ahdyv(i,k)*(sv(i,k) - sv(i,k-1))
             tlatn = ahdyv(i,k+1)*(sv(i,k+1) - sv(i,k))

             hordif = pdepv*((tlate-tlatw)*rdxu(i) + (tlatn-tlats)*rdy(k))

          end if
          
          !          
          ! final
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          vtrans(i,k) = (vtransa(i,k)+cory+islope+slorho+termnl+termnl2+hordif+detrain)*tbotm

	  end if
       end do
    end do
   
    !$omp end do 
    !$omp end parallel

    !     
    ! write new drag on scalar grid if necessary 
    !     
    if (tangle.and.draginmelt) then
       !     
       do i = icalcan+1,icalcen
          do k = kcalcan,kcalcen
             if (newudrag(i-1,k).or.newudrag(i,k)) then
                drag(i,k) = 0.5d0*(ugriddrag(i-1,k) + ugriddrag(i,k))
             end if
          end do
       end do
       !     
    endif
    !     
    !*******************************************************************
    !*** final velocities **********************************************
    !*******************************************************************

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

          pdepu = 5.0d-1*(pdep(i,k) + pdep(i+1,k)) + small
          pdepv = 5.0d-1*(pdep(i,k) + pdep(i,k+1)) + small

	  if (i .ne. icalcen)  su(i,k) = utrans(i,k)/pdepu
	  if (k .ne. kcalcen)  sv(i,k) = vtrans(i,k)/pdepv
 
         ! test for negative depths and set flow to zero if so
!          delta = dt*(rdxu(i)*(utrans(i-1,k)-utrans(i,k)) &
!                    + rdyv(k)*(vtrans(i,k-1)-vtrans(i,k))) 

!          if (bpos(i,k) < (ipos(i,k) - delta)) then
!             utrans(i,k) = 0.d0
!             su(i,k) = 0.d0
!             utrans(i-1,k) = 0.d0
!             su(i-1,k) = 0.d0
!             vtrans(i,k) = 0.d0
!             sv(i,k) = 0.d0
!             vtrans(i,k-1) = 0.d0
!             sv(i,k-1) = 0.d0
!          endif
          utransa(i,k) = utrans(i,k)
          vtransa(i,k) = vtrans(i,k)
          u0a(i,k) = u0(i,k)
          v0a(i,k) = v0(i,k)
       end do
    end do

  end subroutine momentum

  subroutine gwaves(icalcan,icalcen,kcalcan,kcalcen)

     implicit none

     ! calculate local gravity wave speed and the ratio of
     ! advection speed to gravity wave speed

    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen
    integer :: i,k

    real(kind=kdp) :: umid, vmid,speed

    do i=icalcan,icalcen
       do k=kcalcan,kcalcen
          
          !gravity wave speed is sqrt(g'*D) 
          if (jcs(i,k) == 0) cycle   !no water 

!          if (i==43 .and. k== 4) then
!             print *,i,k,'rhoamb',rhoamb(i,k),'rhop',rhop(i,k),rho0,'pdep',pdep(i,k)
!          end if
          gwave_speed(i,k) = abs(sqrt( grav*(rhoamb(i,k)& 
                                            -rhop(i,k))/rho0 &
                                     * pdep(i,k) ))
          
          umid = 5.d-1 * (su(i-1,k) + su(i,k))
          vmid = 5.d-1 * (sv(i,k-1) + sv(i,k))
          speed = sqrt(umid*umid + vmid*vmid)
          
          if ((gwave_speed(i,k) > 0.d0) .and. (speed > 0.d0)) then
             gwave_crit_factor(i,k) = log(speed/gwave_speed(i,k))
          else
             gwave_crit_factor(i,k) = 0.d0
          end if
          
       end do
    end do

  end subroutine gwaves

  subroutine outflow_bound_u_v(icalcan,icalcen,kcalcan,kcalcen)

    implicit none

    ! calculates boundary values for velocity and transports
    ! from neumann conditions. 

    ! local variables

    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen
    integer :: i,k

    if (plume_southern_bc > 0) then
!       print *, 'plume_southern_bc > 0 not fully implemented'
!       stop 1
    end if

    ! southern boundary
    if (kcalcan.le.(domain_kmin+1)) then

       jcd_u(:,domain_kmin) = jcd_u(:,domain_kmin+1)
       jcd_v(:,domain_kmin+1) = jcd_v(:,domain_kmin+2)

       select case (plume_southern_bc)

       case (0)

           ! old d/dy = 0 scheme
	   do i = domain_imin,domain_imax - 1
              su(i,domain_kmin) = su(i,domain_kmin+1)
              utrans( i,domain_kmin) = su(i,domain_kmin+1)*5.0d-1*&
                               (pdep(i,domain_kmin) + pdep(i+1,domain_kmin))
           end do
           do i = domain_imin,domain_imax
              sv(i,domain_kmin) = min(0.d0, sv(i,domain_kmin+1)) !enforce negative
              vtrans(i,domain_kmin) =  min(0.d0, sv(i,domain_kmin+1)) * &
                           5.d-1*(pdep(i,domain_kmin)+pdep(i,domain_kmin+1))
           end do

       case (1)

           do i = domain_imin,domain_imax - 1
              su(i,domain_kmin) = 2.d0*su(i,domain_kmin+1)-su(i,domain_kmin+2)
              utrans(i,domain_kmin) = 2.d0 * (su(i,domain_kmin+1)*5.0d-1*&
                               (pdep(i,domain_kmin+1) + pdep(i+1,domain_kmin+1))) &
                                           - (su(i,domain_kmin+2)*5.0d-1*&
                               (pdep(i,domain_kmin+2) + pdep(i+1,domain_kmin+2)))
	   end do
	   do i = domain_imin,domain_imax
              sv(i,domain_kmin) = 2.d0*sv(i,domain_kmin+1)-sv(i,domain_kmin+2)
              vtrans(i,domain_kmin) = sv(i,domain_kmin+1)* &
                   5.d-1*(pdep(i,domain_kmin+1)+pdep(i,domain_kmin+2))
		 	
             ! need to decide what to do here
           end do

       case (2)

           su(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            su(domain_imin:domain_imax-1, domain_kmin+3), &
                            su(domain_imin:domain_imax-1, domain_kmin+2), &
                            su(domain_imin:domain_imax-1, domain_kmin+1) )


           sv(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            sv(domain_imin:domain_imax-1, domain_kmin+3), &
                            sv(domain_imin:domain_imax-1, domain_kmin+2), &
                            sv(domain_imin:domain_imax-1, domain_kmin+1) )

            utrans(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                        utrans(domain_imin:domain_imax-1, domain_kmin+3), &
                        utrans(domain_imin:domain_imax-1, domain_kmin+2), &
                        utrans(domain_imin:domain_imax-1, domain_kmin+1) )

            vtrans(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                        vtrans(domain_imin:domain_imax-1, domain_kmin+3), &
                        vtrans(domain_imin:domain_imax-1, domain_kmin+2), &
                        vtrans(domain_imin:domain_imax-1, domain_kmin+1) )

       case (3)

           call absorbing_south_boundary_u(su)
           call absorbing_south_boundary_v(sv)
           call absorbing_south_boundary_u(utrans)
           call absorbing_south_boundary_v(vtrans)

       end select

    end if

    ! northern boundary
    if (kcalcen.ge.(domain_kmax - 1)) then

       do i = domain_imin,domain_imax - 1
          jcd_u(:,domain_kmax) = jcd_u(:,domain_kmax -1)
          su(i,domain_kmax) = su(i,domain_kmax-1)
          utrans(i,domain_kmax) = su(i,domain_kmax-1)* 5.0d-1*  &
                    (pdep(i,domain_kmax) + pdep(i+1,domain_kmax))

       end do
       do i = domain_imin,domain_imax
           jcd_v(i,domain_kmax-1) = jcd_v(i,domain_kmax-2)
           sv(i,domain_kmax-1) = max(0.d0, sv(i,domain_kmax-2)) !enforce positive
           vtrans(i,domain_kmax-1) = sv(i,domain_kmax-1)* 5.d-1 * &
                    (pdep(i,domain_kmax-1) + pdep(i,domain_kmax))
                                 !note -1/-2 since north edge of v is -1
       end do
    end if

    ! western boundary
    if (icalcan.le.(domain_imin+1)) then
       do k = domain_kmin,domain_kmax-1
          jcd_u(domain_imin,k) = jcd_u(domain_imin+1,k)
          su(domain_imin,k) = su(domain_imin+1,k)
          utrans(domain_imin,k) = su(domain_imin+1,k)*&
                   5.d-1*(pdep(domain_imin+1,k)+pdep(domain_imin,k))
       end do
       do k = domain_kmin,domain_kmax-1
          jcd_v(domain_imin,k) = jcd_v(domain_imin+1,k)
          sv(domain_imin,k) = sv(domain_imin+1,k)
          vtrans(domain_imin,k) = sv(domain_imin,k) * &
                 5.0d-1 *(pdep(domain_imin,k)+pdep(domain_imin,k+1))
       end do
    end if

    ! eastern boundary
    if (icalcen.ge.(domain_imax - 1)) then
       do k = domain_kmin, domain_kmax - 1
          jcd_u(domain_imax-1,k) = jcd_u(domain_imax-2,k)
          su(domain_imax-1,k) = su(domain_imax-2,k)
          utrans(domain_imax-1,k) = su(domain_imax-1,k)* &
                         5.d-1*(pdep(domain_imax-2,k)+pdep(domain_imax-1,k))
       end do
       do k = domain_kmin,domain_kmax-1
          jcd_v(domain_imax,k) = jcd_v(domain_imax-1,k)
          sv(domain_imax,k) = sv(domain_imax-1,k)
          vtrans(domain_imax,k) = sv(domain_imax-1,k)*5.0d-1*&
                             (pdep(domain_imax,k)+pdep(domain_imax,k+1))
       end do
    end if

  end subroutine outflow_bound_u_v

  subroutine outflow_bound_interface(icalcan,icalcen,kcalcan,kcalcen)

    ! calculates boundary values for interface depth
    ! from neumann conditions. 

    implicit none

    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen

    ! local variables
    integer :: i,k
    real(kind=kdp):: pdepold

    if (plume_southern_bc > 0) then
!       print *, "plume_southern_bc > 0 is not fully implemented"
!       stop 1
    end if

    ! southern boundary
    if(kcalcan.le.(domain_kmin+1)) then

       select case(plume_southern_bc)

       case (0)
          ! old d/dy = 0 version
          jcw(domain_imin:domain_imax,domain_kmin) = &
               jcw(domain_imin:domain_imax,domain_kmin+1)
          jcd_fl(domain_imin:domain_imax,domain_kmin) = &
               jcd_fl(domain_imin:domain_imax,domain_kmin+1)
          pdep(domain_imin:domain_imax,domain_kmin) = &
               pdep(domain_imin:domain_imax,domain_kmin+1)
          ipos(domain_imin:domain_imax,domain_kmin) = &
               bpos(domain_imin:domain_imax,domain_kmin+1) - &
               pdep(domain_imin:domain_imax,domain_kmin)

       case (1)

       do i = domain_imin,domain_imax

          !First-order extrapolation version
          if (jcs(i,domain_kmin) == 1) then
          jcw(i,domain_kmin) = jcw(i,domain_kmin + 1)
          jcd_fl(i,domain_kmin) = jcd_fl(i,domain_kmin + 1)
          pdep(i,domain_kmin) = 2.d0*pdep(i,domain_kmin + 1) - pdep(i,domain_kmin+2)
          ipos(i,domain_kmin) = 2.d0*(bpos(i,domain_kmin + 1)-pdep(i,domain_kmin + 1))&
                               -  (bpos(i,domain_kmin + 2) - pdep(i,domain_kmin + 2))
          end if

       end do

       case (2)

       jcw(i,domain_kmin) = jcw(i,domain_kmin + 1)
       jcd_fl(i,domain_kmin) = jcd_fl(i,domain_kmin + 1)

       pdep(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                        pdep(domain_imin:domain_imax-1, domain_kmin+3), &
                        pdep(domain_imin:domain_imax-1, domain_kmin+2), &
                        pdep(domain_imin:domain_imax-1, domain_kmin+1) )

       ipos(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                        ipos(domain_imin:domain_imax-1, domain_kmin+3), &
                        ipos(domain_imin:domain_imax-1, domain_kmin+2), &
                        ipos(domain_imin:domain_imax-1, domain_kmin+1) )

       case (3)

            call absorbing_south_boundary_scalar(pdep)
            ipos(domain_imin:domain_imax, domain_kmin) = &
                       bpos(domain_imin:domain_imax, domain_kmin) - &
                       pdep(domain_imin:domain_imax, domain_kmin)

       end select
       
    end if

   ! northern boundary
    if(kcalcen.ge.(domain_kmax-1)) then
       do i = domain_imin,domain_imax

           jcw(i,domain_kmax) = jcw(i,domain_kmax - 1)
           jcd_fl(i,domain_kmax) = jcd_fl(i,domain_kmax - 1)
           pdep(i,domain_kmax) = pdep(i,domain_kmax - 1)
           ipos(i,domain_kmax) = bpos(i,domain_kmax) - pdep(i,domain_kmax)    

       end do
    end if

    ! western boundary
    if (icalcan.le.(domain_imin+1)) then
       do k = domain_kmin,domain_kmax
          jcw(domain_imin,k) = jcw(domain_imin+1,k)
          jcd_fl(domain_imin,k) = jcd_fl(domain_imin+1,k)
          ipos(domain_imin,k) = ipos(domain_imin+1,k)
          pdep(domain_imin,k) = pdep(domain_imin+1,k)
       end do
    end if

    ! eastern boundary
    if (icalcen.ge.(domain_imax-1)) then
       do k = domain_kmin,domain_kmax
          jcw(domain_imax,k) = jcw(domain_imax-1,k)
          jcd_fl(domain_imax,k) = jcd_fl(domain_imax-1,k)
          ipos(domain_imax,k) = ipos(domain_imax-1,k)
          pdep(domain_imax,k) = pdep(domain_imax-1,k)
       end do
    end if

  end subroutine outflow_bound_interface

  subroutine outflow_bound_bmelt_entr(icalcan, icalcen, kcalcan, kcalcen)

     implicit none

     integer, intent(in) :: icalcan, icalcen, kcalcan , kcalcen

     integer :: i,k,l

    ! northern boudnary

    if (kcalcen .ge. (domain_kmax-1)) then

        bmelt(:,domain_kmax) = bmelt(:,domain_kmax-1)
        entr(:,domain_kmax) = entr(:, domain_kmax-1)

    end if

    ! southern boundary

    if(kcalcan.le.(domain_kmin+1)) then

       select case(plume_southern_bc)

       case (0)
           bmelt(:,domain_kmin) = bmelt(:,domain_kmin+1)
           entr(:,domain_kmin) = entr(:,domain_kmin+1)
        
       case (1)

           do i = domain_imin,domain_imax
               ! using a first-order extrapolation as the 
               ! 'absorbing boundary condition'
               bmelt(i,domain_kmin) = 2.d0 * bmelt(i,domain_kmin+1) &
                                      - 1.d0*bmelt(i,domain_kmin+2)
               entr(i,domain_kmin) = 2.d0 * entr(i,domain_kmin+1) &
                                     - 1.d0*entr(i,domain_kmin+2)
           end do   

       case (2)

       bmelt(domain_imin:domain_imax,domain_kmin) = extrap3( &
                            bmelt(domain_imin:domain_imax,domain_kmin +3), &
                            bmelt(domain_imin:domain_imax,domain_kmin +2), &
                            bmelt(domain_imin:domain_imax,domain_kmin +1))


       entr(domain_imin:domain_imax,domain_kmin) = extrap3( &
                            entr(domain_imin:domain_imax,domain_kmin +3), &
                            entr(domain_imin:domain_imax,domain_kmin +2), &
                            entr(domain_imin:domain_imax,domain_kmin +1))
 
       case (3)

           call absorbing_south_boundary_scalar(bmelt)
           call absorbing_south_boundary_scalar(entr)
   
       end select

    end if

  end subroutine outflow_bound_bmelt_entr

  subroutine outflow_bound_tsdc(icalcan,icalcen,kcalcan,kcalcen)

    ! calculates boundary values for temperature, salinity, frazil,
    ! density, and
    ! inflow tracers from neumann conditions. 

    implicit none

    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen

    ! local variables
    integer :: i,k,l

    ! southern boundary

    if(kcalcan.le.(domain_kmin+1)) then

      select case(plume_southern_bc)

      case (0)

         temp(:,domain_kmin) = temp(:,domain_kmin+1)
         tempa(:,domain_kmin) = tempa(:,domain_kmin+1)
         salt(:,domain_kmin) = salt(:,domain_kmin+1)
         salta(:,domain_kmin) = salta(:,domain_kmin+1)
         rhop(:,domain_kmin) = rhop(:,domain_kmin+1)
         tfreeze(:,domain_kmin) = tfreeze(:,domain_kmin+1)

      case (1)

       do i = domain_imin,domain_imax
          ! using a first-order extrapolation as the 
          ! 'absorbing boundary condition'
          temp(i,domain_kmin) = 2.d0*  temp(i,domain_kmin+1) &
                                     - temp(i,domain_kmin+2)      
          tempa(i,domain_kmin) = 2.d0*tempa(i,domain_kmin+1) &
                                    - tempa(i,domain_kmin+2)
          salt(i,domain_kmin) = 2.d0*  salt(i,domain_kmin+1) &
                                  -    salt(i,domain_kmin+2)
          salta(i,domain_kmin) = 2.d0*salta(i,domain_kmin+1) &
                                    - salta(i,domain_kmin+2)    
          rhop(i,domain_kmin) = 2.d0* rhop(i,domain_kmin+1) &
                                    - rhop(i,domain_kmin+2)
          tfreeze(i,domain_kmin) = 2.d0*tfreeze(i,domain_kmin+1) &
                                       -tfreeze(i,domain_kmin+2) 

       end do

       case (2)

         temp(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            temp(domain_imin:domain_imax-1, domain_kmin+3), &
                            temp(domain_imin:domain_imax-1, domain_kmin+2), &
                            temp(domain_imin:domain_imax-1, domain_kmin+1) )
         tempa(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            tempa(domain_imin:domain_imax-1, domain_kmin+3), &
                            tempa(domain_imin:domain_imax-1, domain_kmin+2), &
                            tempa(domain_imin:domain_imax-1, domain_kmin+1) )

         salt(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            salt(domain_imin:domain_imax-1, domain_kmin+3), &
                            salt(domain_imin:domain_imax-1, domain_kmin+2), &
                            salt(domain_imin:domain_imax-1, domain_kmin+1) )

         salta(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            salta(domain_imin:domain_imax-1, domain_kmin+3), &
                            salta(domain_imin:domain_imax-1, domain_kmin+2), &
                            salta(domain_imin:domain_imax-1, domain_kmin+1) )

         rhop(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            rhop(domain_imin:domain_imax-1, domain_kmin+3), &
                            rhop(domain_imin:domain_imax-1, domain_kmin+2), &
                            rhop(domain_imin:domain_imax-1, domain_kmin+1) )

         tfreeze(domain_imin:domain_imax-1, domain_kmin) = extrap3( &
                            tfreeze(domain_imin:domain_imax-1, domain_kmin+3),&
                            tfreeze(domain_imin:domain_imax-1, domain_kmin+2),&
                            tfreeze(domain_imin:domain_imax-1, domain_kmin+1) )
      
       case (3)
        
           call absorbing_south_boundary_scalar(temp)
           call absorbing_south_boundary_scalar(tempa)
           call absorbing_south_boundary_scalar(salt)
           call absorbing_south_boundary_scalar(salta)
           call absorbing_south_boundary_scalar(rhop) 
           call absorbing_south_boundary_scalar(tfreeze)

       end select

       if (frazil) then
          do i = domain_imin,domain_imax
             ctot(i,domain_kmin) = ctot(i,domain_kmin+1)                          
             do l = 1,nice
		c_ice(i,domain_kmin,l) = c_ice(i,domain_kmin+1,l)
		ca_ice(i,domain_kmin,l) = ca_ice(i,domain_kmin+1,l)
             end do
          end do
       end if

       if (intrace) then
          do i = domain_imin,domain_imax
             do l = ninfmin,ninfmax
		intracer(i,domain_kmin,l) = intracer(i,domain_kmin+1,l)
		intracera(i,domain_kmin,l) = intracera(i,domain_kmin+1,l)
             end do
          end do
       end if

    end if

    ! northern boundary

    if(kcalcen.ge.domain_kmax-1) then

       do i = domain_imin,domain_imax
          temp(i,domain_kmax) = temp(i,domain_kmax-1)
          tempa(i,domain_kmax) = tempa(i,domain_kmax-1)        
          salt(i,domain_kmax) = salt(i,domain_kmax-1) 
          salta(i,domain_kmax) = salta(i,domain_kmax-1) 
          rhop(i,domain_kmax) = rhop(i,domain_kmax-1)    
          tfreeze(i,domain_kmax) = tfreeze(i,domain_kmax-1)          
       end do

       if (frazil) then
          do i = domain_kmin,domain_kmax
             ctot(i,domain_kmax) = ctot(i,domain_kmax-1)                          
             do l = 1,nice
		c_ice(i,domain_kmax,l) = c_ice(i,domain_kmax-1,l)
		ca_ice(i,domain_kmax,l) = ca_ice(i,domain_kmax-1,l)
             end do
          end do
       end if

       if (intrace) then
          do i = domain_imin,domain_imax
             do l = ninfmin,ninfmax
		intracer(i,domain_kmax,l) = intracer(i,domain_kmax-1,l)
		intracera(i,domain_kmax,l) = intracera(i,domain_kmax-1,l)
             end do
          end do
       end if

    end if

    ! western boundary

    if (icalcan.le.(domain_imin+1)) then

       do k = domain_kmin,domain_kmax
          temp(domain_imin,k) = temp(domain_imin+1,k)
          tempa(domain_imin,k) = tempa(domain_imin+1,k)
          salt(domain_imin,k) = salt(domain_imin+1,k)
          salta(domain_imin,k) = salta(domain_imin+1,k)
          rhop(domain_imin,k) = rhop(domain_imin+1,k)
          tfreeze(domain_imin,k) = tfreeze(domain_imin+1,k)
       end do

       if (frazil) then
          do k = domain_kmin,domain_kmax
             ctot(domain_imin,k) = ctot(domain_imin+1,k)                          
             do l = 1,nice
		c_ice(domain_imin,k,l) = c_ice(domain_imin+1,k,l)
		ca_ice(domain_imin,k,l) = ca_ice(domain_imin+1,k,l)
             end do
          end do
       end if

       if (intrace) then
          do k = domain_kmin,domain_kmax
             do l = ninfmin,ninfmax
		intracer(domain_imin,k,l) = intracer(domain_imin+1,k,l)
		intracera(domain_imin,k,l) = intracera(domain_imin+1,k,l)
             end do
          end do
       end if

    end if

    ! eastern boundary

    if (icalcen.ge.domain_imax-1) then

       do k = domain_kmin,domain_kmax
          temp(domain_imax,k) = temp(domain_imax-1,k)                        
          tempa(domain_imax,k) = tempa(domain_imax-1,k)                          
          salt(domain_imax,k) = salt(domain_imax-1,k)                          
          salta(domain_imax,k) = salta(domain_imax-1,k)                          
          rhop(domain_imax,k) = rhop(domain_imax-1,k)                          
          tfreeze(domain_imax,k) = tfreeze(domain_imax-1,k)                          
       end do

       if (frazil) then
          do k = domain_kmin,domain_kmax
             ctot(domain_imax,k) = ctot(domain_imax-1,k)                          
             do l = 1,nice
		c_ice(domain_imax,k,l) = c_ice(domain_imax-1,k,l)
		ca_ice(domain_imax,k,l) = ca_ice(domain_imax-1,k,l)
             end do
          end do
       end if

       if (intrace) then
          do k = domain_kmin,domain_kmax
             do l = ninfmin,ninfmax
		intracer(domain_imax,k,l) = intracer(domain_imax-1,k,l)
		intracera(domain_imax,k,l) = intracera(domain_imax-1,k,l)
             end do
          end do
       end if

    end if

  end subroutine outflow_bound_tsdc


  subroutine absorbing_south_boundary_scalar(data_array)

    implicit none
  
    real(kind=kdp),intent(inout),dimension(:,:) :: data_array
    integer :: m1,m2,n1,n2
    real(kind=kdp) :: hx,hy
    
    m1 = domain_imin
    m2 = domain_imax
    n1 = domain_kmin
    hx = dx(n1+1)
    hy = dy(n1)

    data_array(m1+1:m2-1, n1) = &
            (data_array(m1+1:m2-1,n1)*hx*hy + &
                  abs(min(sv(m1+1:m2-1,n1+1),0.d0)) * &
                  data_array(m1+1:m2-1,n1+1)*dt*hx - &
                  5.d-1*(su(m1:m2-2,n1+1)+su(m1+1:m2-1,n1+1))* &
                     (data_array(m1+2:m2,n1+1) - &
                      data_array(m1:m2-2,n1+1))*dt*hy) / &
             (hx*hy + dt*hx*abs(min(sv(m1+1:m2-1,n1+1),0.d0)))

  end subroutine absorbing_south_boundary_scalar

  subroutine absorbing_south_boundary_u(data_array)
  
    real(kind=kdp),intent(inout),dimension(:,:) :: data_array
    integer :: m1,m2,n1,n2

    real(kind=kdp) :: hx,hy
    
    hx = dx(n1+1)
    hy = dy(n1)

    m1 = domain_imin
    m2 = domain_imax
    n1 = domain_kmin

    data_array(m1+1:m2-1, n1) = &
            (data_array(m1+1:m2-1,n1)*hx*hy + &
                  abs(min( 5.d-1*(sv(m1+2:m2,n1+1)+sv(m1+1:m2-1,n1+1)),0.d0)) * &
                  data_array(m1+1:m2-1,n1+1)*dt*hx - &
                  su(m1+1:m2-1,n1+1)* &
                     (data_array(m1+2:m2,n1+1) - &
                      data_array(m1:m2-2,n1+1))*dt*hy) / &
             (hx*hy + dt*hx*abs(min( 5.d-1*(sv(m1+2:m2,n1+1)+sv(m1+1:m2-1,n1+1)),0.d0)))

  end subroutine absorbing_south_boundary_u

  subroutine absorbing_south_boundary_v(data_array)
  
    real(kind=kdp),intent(inout),dimension(:,:) :: data_array
    integer :: m1,m2,n1,n2

    real(kind=kdp) :: hx,hy
    
    hx = dx(n1+1)
    hy = dy(n1)
    
    m1 = domain_imin
    m2 = domain_imax
    n1 = domain_kmin

    data_array(m1+1:m2-1, n1) = &
   (data_array(m1+1:m2-1, n1)*hx*hy + &
             abs(min(sv(m1+1:m2-1,n1+1),0.d0)) * & !require a south flowing velocity  
             data_array(m1+1:m2-1,n1+1)*dt*hx - &
        5.d-1*(su(m1:m2-2,n1+1)+su(m1+1:m2-1,n1+1))* &
             (data_array(m1+2:m2,n1+1) - &
              data_array(m1:m2-2,n1+1)) * dt*hy) / &
             (hx*hy + dt*hx*abs(min(sv(m1+1:m2-1,n1+1),0.d0)))

  end subroutine absorbing_south_boundary_v

  subroutine update(iwetmin,iwetmax,kwetmin,kwetmax, &
       icalcan,icalcen,kcalcan,kcalcen,negdep)

    ! update depths, depth-averaged velocities (not depth-integrated),
    ! and wetted area

    implicit none

    integer,intent(inout) :: iwetmin,iwetmax,kwetmin,kwetmax
    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen   
    real(kind=kdp),intent(inout) :: negdep

    ! local variables

    integer :: i,k
    real(kind=kdp) :: iconti,error,pdepold,pdepc,pdepe
    real(kind=kdp) :: dunew,pdepn,dvnew,zd

    ! update interface position

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

          if (jcs(i,k) .ne. 1) cycle

          ! interface position corrected for negative depths at new timestep
          iconti = ipos(i,k)
          ipos(i,k) = max(0.d0, min(bpos(i,k),ipos(i,k)))

          ! report any errors
          error = abs(iconti - ipos(i,k))
          negdep = negdep+error               
          if (error.gt.5.0d-2) then
             write(*,*) 'error: negative depth ',error,' at i=', &
                  i,' k=',k, 'error =', error
             write(11,*) 'error: negative depth ',error,' at i=', &
                  i,' k=',k, 'error =', error
             print *, 'pdep', pdep(i,k),'ipos',iconti
             stop 1
          end if

          ! plume thickness and wetted area (main update of pdep)
          jcd_fl(i,k) = 0      
          pdepold = pdep(i,k)
          jcw(i,k) = 0
          pdep(i,k) = bpos(i,k) - ipos(i,k)

          ! find newly wet area
          if (pdepold .eq. 0.d0 .and. pdep(i,k) .gt. 0.d0) then
             jcd_fl(i,k) = 1 
!             print *, i,k
          else


          end if

          ! find wetted area boundaries (main update of jcw)
          if (pdep(i,k).ge.dcr) then
             jcw(i,k) = 1
             if (.not.(use_min_plume_thickness)) then
                if (i.lt.iwetmin) then
                   iwetmin = i
                else 
                   if (i.gt.iwetmax) iwetmax = i
                end if
                if (k.lt.kwetmin) then
                   kwetmin = k
                else 
                   if (k.gt.kwetmax) then
                      kwetmax = k
                   end if
                end if
             end if
          end if

       end do
    end do

    !print *,'post surface loop'
    ! update velocities
    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

         pdepc = bpos(i,k) - ipos(i,k)


	 if ((i .ne. icalcen) .and. any(jcs(i:(i+1),k) == 1)) then 
           ! u-component
           ! skip icalcen because u(icalcen) is past the edge of the valid region
           ! in the case of an open boundary on the east side
           jcd_u(i,k) = 0
           su(i,k) = 0.d0
           pdepe = bpos(i+1,k) - ipos(i+1,k)
           dunew = 5.0d-1*(pdepe + pdepc)
           if (dunew.gt.small) then
             su(i,k) = utrans(i,k)/dunew
             jcd_u(i,k) = 1
          endif
	  if (depinf(i,k) > 0.d0) then
	     su(i,k) = utrans(i,k)/depinf(i,k)
	     jcd_u(i,k) = 1
	  end if

        end if

	if ((k .ne. kcalcen) .and. any(jcs(i,k:(k+1)) == 1)) then 
          ! v-component
          ! skip kcalcen, because v(kcalcen) is past the edge of the valid region
          jcd_v(i,k) = 0
          sv(i,k) = 0.d0
          pdepn = bpos(i,k+1) - ipos(i,k+1)
          dvnew = 5.0d-1*(pdepn + pdepc)
          if (dvnew.gt.small) then
             sv(i,k) = vtrans(i,k)/dvnew
             jcd_v(i,k) = 1
          endif
	  if (depinf(i,k) > 0.d0) then
	     sv(i,k) = vtrans(i,k)/depinf(i,k)
	     jcd_v(i,k) = 1
	  end if
	end if

       end do
    end do

    !print *,'post velocity loop'
    
    
    ! interpolate ambient density field experienced

    !    rhoamb = get_rhoamb_z(gldep + wcdep - ipos)

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
          if (jcs(i,k).ne.1) cycle
          zd = gldep + wcdep - ipos(i,k)
          rhoamb(i,k) = get_rhoamb_z(zd)
       end do
    end do

    !print *,'post rho loop'
  end subroutine update


  subroutine scalar(icalcan,kcalcan,icalcen,kcalcen)          

    ! calculates transport and entrainment of all scalars, evolution of 
    ! frazil ice, basal melting and freezing, and inflow tracking

    implicit none

    integer,intent(in) :: icalcan,kcalcan,icalcen,kcalcen

    ! local variables

    integer :: i,k,l
    integer :: idel,kdel,idx,kdy,ihilf,khilf,mflag,depthflag,seedindex

    real(kind=kdp),dimension(m_grid,n_grid) :: deltat,deltas
    real(kind=kdp),dimension(m_grid,n_grid):: pdepc,pdepcp,vmid
    real(kind=kdp),dimension(m_grid,n_grid):: umid,speed,bspeed,depth
    real(kind=kdp) :: deltac(m_grid,n_grid,lice),deltatr(m_grid,n_grid,linf)
    real(kind=kdp) :: one,slon,sloe,slos,slow,sumslo
    real(kind=kdp) :: pressure,ttt
    real(kind=kdp) :: tt
    real(kind=kdp) :: dxx,dyy,sx,sy,sxy,r1,r2
    real(kind=kdp) :: dife,difw,difn,difs
    real(kind=kdp) :: dragrt,prden,scden,gambt,gambs,tfreezeb,tfreezei,c1,c2,c3
    real(kind=kdp) :: gamct,gamcs,ucrit,ucl
    real(kind=kdp) :: fmelttot,wturb,mfac1,mfac2,gi,gim1,mi,mip1
    real(kind=kdp) :: amb_depth, ustar
    !real(kind=kdp):: fresh_water_column_thk
    !real(kind=kdp):: discharge_area 
    real(kind=kdp),parameter :: sgd_temp = 0.d0
    real(kind=kdp),parameter :: sgd_salt = 0.d0
    

    one = 1.d0
    prden = 12.5d0*pr**(2.d0/3.d0) - 9.d0
    scden = 12.5d0*sc**(2.d0/3.d0) - 9.d0

    ! calculate values to be used for newly-wet cells using weights 
    ! of interface gradient (thus tempa etc is temp to be used for whole plume)

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
          if (jcd_fl(i,k) .eq. 1) then
             ! newly-wet cell
             if (use_min_plume_thickness) then
                print *, 'should not execute this jcd_fl'
                stop 1
             end if
             slon = dmax1(0.d0,ipos(i,k)-ipos(i,k+1))*jcw(i,k+1)
             sloe = dmax1(0.d0,ipos(i,k)-ipos(i+1,k))*jcw(i+1,k)
             slos = dmax1(0.d0,ipos(i,k)-ipos(i,k-1))*jcw(i,k-1)
             slow = dmax1(0.d0,ipos(i,k)-ipos(i-1,k))*jcw(i-1,k)
             sumslo = slon + sloe + slos + slow        

             if (sumslo .gt. dcr) then
                tempa(i,k) = (temp(i,k+1)*slon + temp(i+1,k)*sloe &
                     + temp(i,k-1)*slos + temp(i-1,k)*slow)/sumslo
                salta(i,k) = (salt(i,k+1)*slon + salt(i+1,k)*sloe &
                     + salt(i,k-1)*slos + salt(i-1,k)*slow)/sumslo

                if (rholinear) then
                   rhop(i,k) = rho_func_linear(tempa(i,k),salta(i,k))
                else
                   if (thermobar) then
                      pressure = 1.0d-1*(wcdep + gldep - ipos(i,k))
                      ttt = tinsitu_func(tempa(i,k),salta(i,k),pressure)
                      rhop(i,k) =rho_func_nonlinear(ttt,salta(i,k),pressure)
                   else
                      rhop(i,k) =  &
                           rho_func_nonlinear(tempa(i,k),salta(i,k),0.d0)
                   end if
                end if

                if (frazil) then
                   ctota(i,k) = 0.d0
                   do l = 1,nice
                      ca_ice(i,k,l) = (c_ice(i,k+1,l)*slon + c_ice(i+1,k,l)*sloe &
                           + c_ice(i,k-1,l)*slos + c_ice(i-1,k,l)*slow)/sumslo
                      ctota(i,k) = ctota(i,k) + ca_ice(i,k,l)
                   end do
                   rhop(i,k) = (1.d0-ctota(i,k))*rhop(i,k)+ctota(i,k)*rhoi
                end if

                if (intrace) then
                   do l = ninfmin,ninfmax
                      intracera(i,k,l) =  &
                           (intracer(i,k+1,l)*slon + intracer(i+1,k,l)*sloe  &
                           + intracer(i,k-1,l)*slos + intracer(i-1,k,l)*slow)/sumslo
                   end do
                end if

             end if
          end if

       end do
    end do

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 0)preliminaries
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

          if (jcs(i,k).ne.1) cycle

          deltat(i,k) = 0.d0
          deltas(i,k) = 0.d0

          pdepc(i,k) = bpos(i,k) - ipos(i,k)
          pdepcp(i,k) = pdepc(i,k) + dcr
          tt = 5.0d-1*dx(i)*rdxu(i)
          umid(i,k) = tt*su(i-1,k) + (1.d0-tt)*su(i,k)
          tt = 5.0d-1*dy(k)*rdyv(k)
          vmid(i,k) = tt*sv(i,k-1) + (1.d0-tt)*sv(i,k) 
          speed(i,k) = umid(i,k)**2 + vmid(i,k)**2 + small
          bspeed(i,k) = dsqrt(speed(i,k))

          ! calculate freezing temperature at plume 
          ! mid-depth (even if no frazil)
          depth(i,k) = gldep + wcdep - 5.0d-1*(ipos(i,k) + bpos(i,k))
          tfreeze(i,k) = freezing_temp_func(salta(i,k),depth(i,k))
       end do
    end do


    if (frazil) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             if (jcs(i,k).ne.1) cycle
             do l = 1,nice
                deltac(i,k,l) = 0.d0
             end do
          end do
       end do
    end if

    if (intrace) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             if (jcs(i,k) .ne. 1) cycle
             do l = ninfmin,ninfmax
                deltatr(i,k,l) = 0.d0
             end do
          end do
       end do
    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 1)entrainment
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (entrain) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k).ne.1) cycle

             if (pdepc(i,k).gt.edepth) then         

	     	! interpolate ambient values to interface and calculate
                ! quantities for entrainment
		amb_depth = wcdep + gldep - ipos(i,k)
		atemp(i,k) = get_tamb_z(amb_depth)
		asalt(i,k) = get_samb_z(amb_depth)

		if (use_periodic_forcing .and. (k .ge. kfirst) &
                                         .and. (k .le. klast)) then

                  ! artificially force the ambient seawater to have a much
                  ! lower salinity in order to make the forcing have a 
                  ! greater effect
		  asalt(i,k) = asalt(i,k) - 10.0

                end if

                if (entype .ne. 5) then

                   deltat(i,k) = deltat(i,k)  &
                        + dt*entr(i,k)*(atemp(i,k) - tempa(i,k))/pdepcp(i,k)
                   deltas(i,k) =  deltas(i,k) &
                        + dt*entr(i,k)*(asalt(i,k)-salta(i,k))/pdepcp(i,k)

                end if

             endif
          end do
       end do

    if (sgd_type == 0 .or. sgd_type == 1) then

       ! uniform flux across inflow edge
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
    
		deltat(i,k) = deltat(i,k)  &
                             + dt*sgd(i,k)*(sgd_temp - tempa(i,k))/pdepcp(i,k)
                deltas(i,k) =  deltas(i,k) &
                             + dt*sgd(i,k)*(sgd_salt - salta(i,k))/pdepcp(i,k)
 
          end do
       end do
    end if

       if (frazil) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (jcs(i,k) .ne. 1) cycle
                if (pdepc(i,k).gt.edepth) then         
                   do l=1,nice
                      print *, 'unexpected use of frazil code'
                      stop 1
                      deltac(i,k,l) = deltac(i,k,l)+ dt*entr(i,k)*(-ca_ice(i,k,l))/pdepcp(i,k)
                   end do
                end if
             end do
          end do
       end if

    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 2)advection (selective vector upstream)
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (nonlin) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle

             !  choose upstream quarter
             idel = -int(sign(one,umid(i,k)))
             kdel = -int(sign(one,vmid(i,k)))
             ihilf = i + idel
             khilf = k + kdel
             !  choose the relevant upstream cell dimensions 
             idx = i + (idel-1)/2
             kdy = k + (kdel-1)/2
             dxx = dx(idx)
             dyy = dy(kdy)
             !  weighting coefficients
             sx = abs(umid(i,k))*dt/dxx
             sy = abs(vmid(i,k))*dt/dyy
             sxy = sx - sy
             r1 = (dsign(one,sxy) + 1.)*5.0d-1
             r2 = one - r1

             if ((sx.gt.one).or.(sy.gt.one)) then
                write(*,*) 'error: tsc courant exceeded at i=',i,' k=',k
                write(11,*) 'error: tsc courant exceeded at i=',i,' k=',k
             end if

	deltat(i,k) = deltat(i,k)  &
                 + (r1*sy+r2*sx)*(tempa(ihilf,khilf)*     jcw(ihilf,khilf) &
                 +                tempa(i,k)        *(1 - jcw(ihilf,khilf)))  &
                 + r1*sxy*       (tempa(ihilf,k)*         jcw(ihilf,k) &
                 +                tempa(i,k)*          (1-jcw(ihilf,k))) &
                 - r2*sxy*       (tempa(i,khilf)*         jcw(i,khilf) &
                 +                tempa(i,k)*          (1-jcw(i,khilf))) &
                 - (r2*sy+r1*sx)* tempa(i,k)
       
	deltas(i,k) = deltas(i,k)  &
                + (r1*sy+r2*sx)*(salta(ihilf,khilf)*      jcw(ihilf,khilf) &
                +                salta(i,k)*         (1 - jcw(ihilf,khilf)))  &
                + r1*sxy*       (salta(ihilf,k)*          jcw(ihilf,k) &
                +                salta(i,k)*         (1 - jcw(ihilf,k))) &
                - r2*sxy*       (salta(i,khilf)*          jcw(i,khilf) &
                +                salta(i,k)*         ( 1- jcw(i,khilf))) &
                - (r2*sy+r1*sx)* salta(i,k)      

             if (frazil) then
                do l = 1,nice
                   deltac(i,k,l) = deltac(i,k,l) &
                        + (r1*sy+r2*sx)*(ca_ice(ihilf,khilf,l)*jcw(ihilf,khilf) &
                        + (1 - jcw(ihilf,khilf))*ca_ice(i,k,l))  &
                        + r1*sxy*(ca_ice(ihilf,k,l)*jcw(ihilf,k)  &
                        + (1-jcw(ihilf,k))*ca_ice(i,k,l)) &
                        - r2*sxy*(ca_ice(i,khilf,l)*jcw(i,khilf) &
                        + (1-jcw(i,khilf))*ca_ice(i,k,l)) &
                        - (r2*sy+r1*sx)*ca_ice(i,k,l)      
                end do
             end if

             if (intrace) then
                do l = ninfmin,ninfmax
                   deltatr(i,k,l) = deltatr(i,k,l) &
                        + (r1*sy+r2*sx)*(intracera(ihilf,khilf,l)*jcw(ihilf,khilf) &
                        + (1 - jcw(ihilf,khilf))*intracera(i,k,l)) &
                        + r1*sxy*(intracera(ihilf,k,l)*jcw(ihilf,k)  &
                        + (1-jcw(ihilf,k))*intracera(i,k,l)) &
                        - r2*sxy*(intracera(i,khilf,l)*jcw(i,khilf) &
                        + (1-jcw(i,khilf))*intracera(i,k,l)) &
                        - (r2*sy+r1*sx)*intracera(i,k,l)      
                end do
             end if

          end do
       end do

    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 3)horizontal diffusion
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (horturb) then

       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle

             dife = (tempa(i+1,k) - tempa(i,k))*jcw(i+1,k)*khgrid(i,k)*rdx(i)
             difw = (tempa(i,k) - tempa(i-1,k))*jcw(i-1,k)*khgrid(i,k)*rdx(i-1)
             difn = (tempa(i,k+1) - tempa(i,k))*jcw(i,k+1)*khgrid(i,k)*rdy(k)
             difs = (tempa(i,k) - tempa(i,k-1))*jcw(i,k-1)*khgrid(i,k)*rdy(k-1)

             deltat(i,k) = deltat(i,k) &
                  + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt 
 
             dife = (salta(i+1,k) - salta(i,k))*jcw(i+1,k)*khgrid(i,k)*rdx(i)
             difw = (salta(i,k) - salta(i-1,k))*jcw(i-1,k)*khgrid(i,k)*rdx(i-1)
             difn = (salta(i,k+1) - salta(i,k))*jcw(i,k+1)*khgrid(i,k)*rdy(k)
             difs = (salta(i,k) - salta(i,k-1))*jcw(i,k-1)*khgrid(i,k)*rdy(k-1)

             deltas(i,k) = deltas(i,k) &
                  + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt

          end do
       end do

       if (frazil) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen

                if (jcs(i,k) .ne. 1) cycle

                print *, 'unexpected use of frazil code'
                stop 1
                do l = 1,nice
                   dife = (ca_ice(i+1,k,l) - ca_ice(i,k,l))*jcw(i+1,k)*khgrid(i,k)*rdx(i)
                   difw = (ca_ice(i,k,l) - ca_ice(i-1,k,l))*jcw(i-1,k)*khgrid(i,k)*rdx(i-1)
                   difn = (ca_ice(i,k+1,l) - ca_ice(i,k,l))*jcw(i,k+1)*khgrid(i,k)*rdy(k)
                   difs = (ca_ice(i,k,l) - ca_ice(i,k-1,l))*jcw(i,k-1)*khgrid(i,k)*rdy(k-1)
                   deltac(i,k,l) = deltac(i,k,l)  &
                        + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt
                end do
             end do
          end do
       end if

       if (intrace) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen

                if (jcs(i,k) .ne. 1) cycle

                do l = ninfmin,ninfmax
                   dife =  &
                        (intracera(i+1,k,l) - intracera(i,k,l))*jcw(i+1,k)*khgrid(i,k)*rdx(i)
                   difw =  &
                        (intracera(i,k,l) - intracera(i-1,k,l))*jcw(i-1,k)*khgrid(i,k)*rdx(i-1)
                   difn =  &
                        (intracera(i,k+1,l) - intracera(i,k,l))*jcw(i,k+1)*khgrid(i,k)*rdy(k)
                   difs =  &
                        (intracera(i,k,l) - intracera(i,k-1,l))*jcw(i,k-1)*khgrid(i,k)*rdy(k-1)
                   deltatr(i,k,l) = deltatr(i,k,l)  &
                        + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt
                end do
             end do
          end do
       end if

    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 4)basal melting and freezing
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (basmelt) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle

             if (pdepc(i,k).gt.mdepth) then         
                dragrt = dsqrt(drag(i,k))
                ustar = sqrt(drag(i,k)*( bspeed(i,k)**2.d0 + &
                                         0.5d0*local_tidal_speed(i,k)**2.d0) + &
                             u_star_offset**2.d0)
                ! find turbulent exchange coefficients
                gambt = ustar/ &
                     (2.12d0*dlog(ustar*pdepc(i,k)/nu0) + prden)             
                gambs = ustar/ &
                     (2.12d0*dlog(ustar*pdepc(i,k)/nu0) + scden)

                ! calculate freezing point of plume at shelf base,
                ! decide if melt (mflag = 1) 
                ! or freeze and calculate freezing point of ice at shelf
                ! base

                tfreezeb = fta*salta(i,k) + ftb + ftc*(gldep+wcdep-bpos(i,k))
                mflag = (1 + int(sign(1.d0,tempa(i,k) - tfreezeb)))/2
                depthflag = (1 + int(sign(1.d0, gldep + wcdep - bpos(i,k) - min_melt_depth)))/2

                tfreezei = (1 - mflag)*fta*si  &
                     + ftb + ftc*(gldep + wcdep - bpos(i,k))

                ! calculate coefficients in quadratic to be solved
                c1 = lat/c0 + mflag*(ci/c0)*(tfreezei-tint(i,k))
                c2 = gambs*(lat/c0 + mflag*(ci/c0)*(tfreezeb-tint(i,k)))  &
                     + gambt*(tfreezei-tempa(i,k))
                c3 = gambs*gambt*(tfreezeb - tempa(i,k))

                ! calculate melt rate
                if (separated(i,k) == 1) then
                   bmelt(i,k) = 0.d0
                else
                  bmelt(i,k) = -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1) * depthflag
                end if
                
                ! calculate basal temperature and salinity
                btemp(i,k) = (gambt*tempa(i,k)+mflag*(ci/c0)*bmelt(i,k)*tint(i,k) &
                     - (lat/c0)*bmelt(i,k) )/  &
                     (gambt + mflag*(ci/c0)*bmelt(i,k))
                
                bsalt(i,k) = (btemp(i,k) - ftb  &
                     - ftc*(gldep + wcdep - bpos(i,k)))/fta

                ! calculate change of heat, salt, and frazil due
                ! to meltwater/freezewater flux
                ! (ice shelf has zero salinity and no frazil)

                deltat(i,k) = deltat(i,k)  &
                     + dt*bmelt(i,k)*(btemp(i,k) - tempa(i,k))/pdepcp(i,k)
                deltas(i,k) = deltas(i,k)  &
                     - dt*bmelt(i,k)*salta(i,k)/pdepcp(i,k)
              
                ! add diffusive flux of latent heat release/uptake
	        deltat(i,k) = deltat(i,k)  &
                     - dt*bmelt(i,k)*((lat/c0)+(ci/c0)*(btemp(i,k)-tint(i,k)))&
                      /pdepcp(i,k)

             end if
          end do
       end do

       if (frazil) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (pdepc(i,k).gt.mdepth) then         
                   do l=1,nice
                      deltac(i,k,l) = deltac(i,k,l)  &
                           - dt*bmelt(i,k)*ca_ice(i,k,l)/pdepcp(i,k)
                   end do
                end if
             end do
          end do
       end if

    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! 5)frazil nucleation, growth, melting, and precipitation
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (frazil) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             if (pdepc(i,k).gt.fdepth) then         
                !              
                fmelttot = 0.d0

                if (nice.eq.1) then
                   ! ----------------------------
                   ! jenkins and bombosch version
                   ! ----------------------------
                   if (ca_ice(i,k,1).ge.cmin(1)) then

                      ! a) calculate growth/melting of frazil 
                      gamct = (nuss(1)*kt)/(ar*r(1))
                      gamcs = (nuss(1)*ks)/(ar*r(1))

                      c1 = 5.0d-1*lat*r(1)/(c0*ca_ice(i,k,1)*pdepc(i,k)*gamct)
                      c2 = ftb + ftc*depth(i,k) - tempa(i,k)  &
                           + lat*gamcs/(c0*gamct)
                      c3 = (2.0d0*ctota(i,k)*pdepc(i,k)*gamcs/r(1)) &
                           *(fta*salta(i,k)+ ftb + ftc*depth(i,k) - tempa(i,k))

                      fmelt(i,k,1) = &
                           -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)
                      ctempd(i,k) = tempa(i,k) - fmelt(i,k,1)* &
                           (5.0d-1*lat*r(1)/(c0*ctota(i,k)*pdepc(i,k)*gamct))

                      ! b) calculate precipitation of frazil
                      ucl =  &
                           dsqrt((1.0d-1*grav*re(1)*(rho0-rhoi))/(rho0*drag(i,k)))
                      ucrit = (1.d0 - bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                      fppn(i,k,1) =  &
                           dmin1(-(rhoi/rho0)*ca_ice(i,k,1)*wi(1)*ucrit,0.d0)

                      ! c) calculate overall frazil source for this 
                      ! size class and 
                      !    total melt for ts source
                      deltac(i,k,1) = deltac(i,k,1)  &
                           - dt*fppn(i,k,1)*ca_ice(i,k,1)/pdepcp(i,k) &
                           + (rho0/rhoi)*dt*(fppn(i,k,1) &
                           - fmelt(i,k,1))/pdepcp(i,k)
                      fmelttot = fmelt(i,k,1)
                   end if

                else
                   ! ----------------------------
                   ! smedsrud and jenkins version
                   ! ----------------------------
                   if (ctota(i,k).ge.cmin(1)) then

                      ! a) calculate secondary nucleation of frazil

                      fnuc(i,k,1) = 0.d0
                      do l=2,nice
                         wturb = dsqrt(wi(l)*wi(l)  &
                              + ((4.d0*eps)/(15.d0*nu0))*re(l)*re(l))
                         fnuc(i,k,l) =  &
                              - (rhoi/rho0)*pdepc(i,k)*pi*dble(nbar)*wturb &
                              *(re(1)**3/re(l))*ca_ice(i,k,l)
                         fnuc(i,k,1) = fnuc(i,k,1) - fnuc(i,k,l)
                      end do

                      ! b) calculate growth/melting of frazil 
                      mfac1 = (2.d0*c0*kt/lat)*(tfreeze(i,k) - tempa(i,k))
                      mfac2 = (rhoi/rho0)*pdepc(i,k)

                      if (tempa(i,k).lt.tfreeze(i,k)) then
                         gi = mfac1*nuss(1)*ca_ice(i,k,1)/(r(1)*r(1))
                         fmelt(i,k,1) = mfac2*(vol(1)/(vol(2) - vol(1)))*gi

                         do l=2,nice-1
                            gi = mfac1*nuss(l)*ca_ice(i,k,l)/(r(l)*r(l))
                            gim1 = mfac1*nuss(l-1)*ca_ice(i,k,l-1)/(r(l-1)*r(l-1))
                            fmelt(i,k,l) = mfac2*( &
                                 (vol(l)/(vol(l+1) - vol(l)))*gi  &           
                                 - (vol(l)/(vol(l) - vol(l-1)))*gim1)         
                         end do

                         gim1 = mfac1*nuss(nice-1)*ca_ice(i,k,nice-1) &
                              /(r(nice-1)*r(nice-1))
                         fmelt(i,k,nice) = - mfac2* &
                              (vol(nice)/(vol(nice) - vol(nice-1)))*gim1
                      else
                         mi = mfac1*nuss(1)*ca_ice(i,k,1) &
                              *(1.d0/r(1) + 1.d0/thi(1))/r(1)
                         mip1 = mfac1*nuss(2)*ca_ice(i,k,2) &
                              *(1.d0/r(2) + 1.d0/thi(2))/r(2)
                         fmelt(i,k,1) = mfac2*  &
                              ((vol(1)/(vol(2) - vol(1)))*mip1 - mi)

                         do l=2,nice-1
                            mi = mfac1*nuss(l)*ca_ice(i,k,l) &
                                 *(1.d0/r(l) + 1.d0/thi(l))/r(l)
                            mip1 = mfac1*nuss(l+1)*ca_ice(i,k,l+1)* &
                                 (1.d0/r(l+1) + 1.d0/thi(l+1))/r(l+1)
                            fmelt(i,k,l) = mfac2*( &
                                 (vol(l)/(vol(l+1) - vol(l)))*mip1     &       
                                 - (vol(l)/(vol(l) - vol(l-1)))*mi)  
                         end do

                         mi = mfac1*nuss(nice)*ca_ice(i,k,nice)* &
                              (1.d0/r(nice) + 1.d0/thi(nice))/r(nice)
                         fmelt(i,k,nice) = - mfac2* &
                              (vol(nice)/(vol(nice) - vol(nice-1)))*mi
                      end if

                      ! c) calculate precipitation of frazil
                      do l=1,nice
                         ucl =  &
                              dsqrt((1.0d-1*grav*re(l)*(rho0-rhoi))/(rho0*drag(i,k)))
                         ucrit = (1.d0-bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                         fppn(i,k,l) =  &
                              dmin1(-(rhoi/rho0)*ca_ice(i,k,l)*wi(l)*ucrit,0.d0)
                      end do

                      ! d) calculate overall frazil source for this 
                      ! size class and 
                      !    total melt for ts source
                      do l=1,nice
                         deltac(i,k,l) = deltac(i,k,l) &
                              - dt*fppn(i,k,l)*ca_ice(i,k,l)/pdepcp(i,k) &
                              + (rho0/rhoi)*dt*(fppn(i,k,l) + fnuc(i,k,l)  &
                              - fmelt(i,k,l))/pdepcp(i,k)
                         fmelttot = fmelttot + fmelt(i,k,l)
                      end do
                   end if

                end if

                ! -------------
                ! both versions
                ! -------------
                ! calculate change of heat and salt due to frazil 
                ! meltwater/freezewater flux
                deltat(i,k) = deltat(i,k)  &
                     + dt*fmelttot*(tfreeze(i,k) - tempa(i,k))/pdepcp(i,k)
                deltas(i,k) = deltas(i,k)  &
                     - dt*fmelttot*salta(i,k)/pdepcp(i,k)

                ! add diffusive flux of latent heat release/uptake
                deltat(i,k) = deltat(i,k)-dt*(lat/c0)*fmelttot/pdepcp(i,k)

             end if
          end do
       end do

    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! final
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
          if (jcs(i,k) .ne. 1) cycle
	     temp(i,k) = tempa(i,k) + deltat(i,k)
             salt(i,k) = salta(i,k) + deltas(i,k)
       end do
    end do

    if (frazil) then

       ! if desired, check pre-emptively for negative frazil concentrations 
       ! and set to zero if found.  keep track of total negative concentration
       if (negfrz) then

          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                do l = 1,nice
                   if (ca_ice(i,k,l) + deltac(i,k,l).ge.0.d0) then
                      c_ice(i,k,l) = ca_ice(i,k,l) + deltac(i,k,l)
                   else
                      c_ice(i,k,l) = 0.d0
                      frzcut(l) = frzcut(l) - ca_ice(i,k,l) - deltac(i,k,l)
                   end if
                end do
             end do
          end do

       else

          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                do l = 1,nice
                   c_ice(i,k,l) = ca_ice(i,k,l) + deltac(i,k,l)
                end do
             end do
          end do

       end if

       ! total up frazil concentrations
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             ctot(i,k) = 0.d0
             do l = 1,nice
                ctot(i,k) = ctot(i,k) + c_ice(i,k,l) 
             end do
          end do
       end do

    end if

    if (intrace) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             do l = ninfmin,ninfmax
                intracer(i,k,l) = intracera(i,k,l) + deltatr(i,k,l)
             end do
          end do
       end do
    end if

    ! seed frazil ice
    ! ---------------
    if (frazil) then 
       !           
       ! smedsrud and jenkins seding strategy:
       ! (find and then seed southernmost supercooled cells          
       !           
       if (seedtype.eq.1) then
          !     
          seedindex = n_grid
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (pdepc(i,k).gt.fdepth) then
                   if (tempa(i,k).lt.tfreeze(i,k)) then
                      seedindex = min(k,seedindex)
                   end if
                end if
             end do
          end do

          if (seedindex.lt.n_grid) then
             do i = icalcan,icalcen
                if (pdepc(i,seedindex).gt.fdepth) then
                   ctot(i,seedindex) = 0.d0
                   do l=1,nice
                      c_ice(i,seedindex,l) = cseed(l)
                      ctot(i,seedindex) = ctota(i,seedindex) + cseed(l)
                   end do
                end if
             end do
          end if

       end if

       ! holland and feltham seding strategy:
       ! (if cell is newly supercooled (i.e. was superheated on last 
       ! timestep), seed
       ! with frazil if advection and diffusion have not exceeded seed
       ! population)

       if (seedtype.eq.2) then

          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (pdepc(i,k).gt.fdepth) then         
                   if (tempa(i,k).gt.tfreeze(i,k)) then
                      jcd_fseed(i,k) = 0
                   else
                      if (jcd_fseed(i,k).eq.0) then
                         ctot(i,k) = 0.d0
                         do l = 1,nice
                            if (c_ice(i,k,l).lt.cseed(l)) c_ice(i,k,l) = cseed(l)
                            ctot(i,k) = ctot(i,k) + c_ice(i,k,l)
                         end do
                         jcd_fseed(i,k) = 1     
                      end if
                   end if
                end if
             end do
          end do
          !           
       end if
       !           
    end if
    !           
    ! tidy up variables
    ! -----------------
    tempa = temp
    salta = salt
    !    do i = icalcan,icalcen
    !       do k = kcalcan,kcalcen
    !          tempa(i,k) = temp(i,k)
    !          salta(i,k) = salt(i,k)
    !       end do
    !    end do
    !     
    if (frazil) then
       ctota = ctot
       ca_ice = c_ice
    end if
    !       do i = icalcan,icalcen
    !          do k = kcalcan,kcalcen
    !             ctota(i,k) = ctot(i,k)
    !             do l = 1,nice
    !                ca_ice(i,k,l) = c_ice(i,k,l)
    !             end do
    !          end do
    !       end do
    !    end if
    !     
    if (intrace) then
       intracera = intracer
       !       do i = icalcan,icalcen
       !          do k = kcalcan,kcalcen
       !             do l = ninfmin,ninfmax
       !                intracera(i,k,l) = intracer(i,k,l)
       !             end do
       !          end do
       !       end do
    end if

  end subroutine scalar

  subroutine rho_calc(icalcan,icalcen,kcalcan,kcalcen,sepflag)

    ! calculate density and warn if plume has reached separation point

    implicit none

    ! local variables

    logical :: sepflag

    integer :: i,k,icalcan,icalcen,kcalcan,kcalcen

    real(kind=kdp) :: p,rhoatmp,depth


    ! get density of water fraction. if thermobaricity included, compute
    ! in situ temp from potential temp. for future reference
    if (rholinear) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             if (jcs(i,k) .ne. 1) cycle
             rhop(i,k) = rho_func_linear(temp(i,k),salta(i,k)) 
          end do
       end do
    else
       if (thermobar) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (jcs(i,k) .ne. 1) cycle
                p = (wcdep + gldep - ipos(i,k) - 5.0d-1*pdep(i,k))*1.0d-1
                tins(i,k) = tinsitu_func(temp(i,k),salta(i,k),p)  
                rhop(i,k) = rho_func_nonlinear(tins(i,k),salta(i,k),p) 
             end do
          end do
       else
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (jcs(i,k) .ne. 1) cycle
                rhop(i,k) = rho_func_nonlinear(temp(i,k),salta(i,k),0.d0) 
             end do
          end do
       end if
    end if

    ! add density of ice fraction if appropriate
    if (frazil) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             if (jcs(i,k) .ne. 1) cycle
             rhop(i,k) = (1.d0 - ctot(i,k))*rhop(i,k) + ctot(i,k)*rhoi
          end do
       end do
    end if

    ! calculate density of ambient fluid and test for separation
    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
          if (pdep(i,k).gt.dcr) then

             depth = wcdep + gldep - ipos(i,k) - 5.0d-1*pdep(i,k)
	     rhoatmp = get_rhoamb_z(depth)

             if ((rhoatmp + septol < rhop(i,k)).and. &
                  (.not.sepflag)) then
!                write(*,*) 'warning: plume separated at i=',i,' k=',k
!                write(11,*) 'warning: plume separated at i=',i,' k=',k
!                sepflag = .true.
                separated(i,k) = 1
             else
                separated(i,k) = 0
             end if
             !           
          end if
       end do
    end do

  end subroutine rho_calc

!  elemental real(kind=kdp) function get_rhoamb_z(z)

!    implicit none
!    real(kind=kdp),intent(in) :: z ! depth at which rhoamb is sought

!    get_rhoamb_z = rhovf(1)
    
!  end function get_rhoamb_z

   elemental real(kind=kdp) function get_rhoamb_z(z)

    implicit none

    real(kind=kdp),intent(in) :: z !depth at which rhoamb is sought

    integer :: izo,izu
    real(kind=kdp) :: difu,difo
    character(128) :: error_message

    izo = int(z/dzincr) + 1
    izu = izo + 1
    difu = dble(izo)*dzincr - z
    difo = dzincr - difu

    get_rhoamb_z = (difu*rhovf(izo) + difo*rhovf(izu))/dzincr

  end function get_rhoamb_z

  elemental real(kind=kdp) function get_tamb_z(z)
    implicit none
    real(kind=kdp),intent(in) :: z ! depth at which tamb is sought

    integer :: izo,izu
    real(kind=kdp) :: difu,difo

    izo = int(z/dzincr) + 1
    izu = izo + 1
    difu = dble(izo)*dzincr - z
    difo = dzincr - difu
    get_tamb_z = (difu*tamb(izo) + difo*tamb(izu))/dzincr

  end function get_tamb_z

  elemental real(kind=kdp) function get_samb_z(z)
    implicit none
    real(kind=kdp),intent(in) :: z

    integer :: izo,izu
    real(kind=kdp) :: difu,difo

    izo = int(z/dzincr) + 1
    izu = izo + 1
    difu = dble(izo)*dzincr - z
    difo = dzincr - difu
    get_samb_z = (difu*samb(izo) + difo*samb(izu))/dzincr	    

  end function get_samb_z

  elemental real(kind=kdp) function extrap3(y1,y2,y3)
  
    implicit none
    real(kind=kdp),intent(in) :: y1,y2,y3

    extrap3 = y1 - 3.d0*y2 + 3.d0*y3
    
  end function extrap3

  subroutine filter_waves(raw_data,ocean_mask,short_waves, long_waves, as, bs)

    implicit none

    real(kind=kdp),dimension(:,:),intent(in)  :: raw_data
    integer,dimension(:,:),intent(in)         :: ocean_mask
    real(kind=kdp),dimension(:,:),intent(out) ::short_waves, long_waves, as, bs

    ! raw_data is the field to be filtered
    ! ocean_mask is equal to 1 where the raw_data is valid, 0 otherwise

    integer :: m,n,i,k
    integer :: lower_ilim,lower_klim, upper_ilim,upper_klim
    integer :: win_imin,win_kmin,win_imax,win_kmax
    integer,parameter :: window_width = 3

    real(kind=kdp),dimension(window_width,window_width) :: data_window
    integer(kind=kdp),dimension(window_width,window_width) :: ocean_mask_window

    real(kind=kdp),dimension(window_width,window_width) :: xmat,ymat
    real(kind=kdp),dimension(3,3) :: matrix
    real(kind=kdp),dimension(3)   :: rhs
    real(kind=kdp)                :: multiplier, x0, y0

    as = 0.d0
    bs = 0.d0
    long_waves = 0.d0
    short_waves = 0.d0

    m = size(raw_data,1)
    n = size(raw_data,2)

    x0 = floor( (window_width-1) /2.d0 )
    y0 = floor( (window_width-1) /2.d0 )

    do i= 1,window_width
       xmat(i,:) = real(i) - (x0+1.d0)
       ymat(:,i) = real(i) - (y0+1.d0)
    end do

    !$omp parallel default(none) &
    !$omp private(i,k,lower_ilim,lower_klim,upper_ilim,upper_klim, &
    !$omp         win_imin,win_kmin,win_imax,win_kmax, & 
    !$omp         data_window, ocean_mask_window,matrix,rhs, &
    !$omp         multiplier) &
    !$omp shared(raw_data, ocean_mask, short_waves, long_waves, as,bs, &
    !$omp        x0,y0, xmat,ymat, &
    !$omp        domain_imin,domain_imax,domain_kmin,domain_kmax)
    !$omp do
    do i=domain_imin,domain_imax
       do k=domain_kmin,domain_kmax
         
          lower_ilim = max(domain_imin,i-int(floor(  (window_width-1)/2.d0)))
          upper_ilim = min(domain_imax,i+int(ceiling((window_width-1)/2.d0)))
          lower_klim = max(domain_kmin,k-int(floor(  (window_width-1)/2.d0)))
          upper_klim = min(domain_kmax,k+int(ceiling((window_width-1)/2.d0)))

          data_window = 0.d0
          win_imin = int((lower_ilim-i)+(x0+1.d0))
          win_imax = win_imin + (upper_ilim-lower_ilim)
          win_kmin = int((lower_klim-k)+(y0+1.d0))
          win_kmax = win_kmin + (upper_klim-lower_klim)

          data_window( win_imin:win_imax, &
                       win_kmin:win_kmax ) = &
               raw_data(lower_ilim:upper_ilim,lower_klim:upper_klim)

          ocean_mask_window = 0
          ocean_mask_window(win_imin:win_imax, &
                       win_kmin:win_kmax ) = &
               ocean_mask(lower_ilim:upper_ilim,lower_klim:upper_klim)

          ! try determining the value at (i,k) determined by a local
          ! linear approximation in the data_window

          if (ocean_mask(i,k) == 1) then

             matrix(1,1) = sum( xmat * xmat, mask=(ocean_mask_window==1))
             matrix(1,2) = sum( xmat * ymat, mask=(ocean_mask_window==1))
             matrix(1,3) = sum( xmat       , mask=(ocean_mask_window==1))
             matrix(2,2) = sum( ymat * ymat, mask=(ocean_mask_window==1))
             matrix(2,3) = sum( ymat       , mask=(ocean_mask_window==1))
             matrix(3,3) = sum(ocean_mask_window)
             matrix(2,1) = matrix(1,2)
             matrix(3,1) = matrix(1,3)
             matrix(3,2) = matrix(2,3)

             rhs(1) = sum( xmat * data_window, mask=(ocean_mask_window==1))
             rhs(2) = sum( ymat * data_window, mask=(ocean_mask_window==1))
             rhs(3) = sum( data_window       , mask=(ocean_mask_window==1))

             ! do Gaussian elim on matrix
             if (matrix(1,1) == 0.d0) then
                print *, i,k
                print *, xmat(:,3)
                print *, xmat(:,2)
                print *, xmat(:,1)
                print *, 'zero entry in 1,1'
                stop 1
             end if
             multiplier = matrix(2,1)/matrix(1,1)
             matrix(2,:) = matrix(2,:) - multiplier*matrix(1,:)
             rhs(2)      = rhs(2) -      multiplier*rhs(1)

             if (matrix(1,1) == 0.d0) then
                print *, i,k
                print *, matrix(1,:)
                print *, matrix(2,:)
                print *, matrix(3,:)
                print *, 'zero entry in 1,1'
                stop 1
             end if
             multiplier = matrix(3,1)/matrix(1,1)
             matrix(3,:) = matrix(3,:) - multiplier*matrix(1,:)
             rhs(3)      = rhs(3) -      multiplier*rhs(1)

             if (matrix(2,2) == 0.d0) then
                print *, i,k
                print *, matrix(1,:)
                print *, matrix(2,:)
                print *, matrix(3,:)
                print *, 'zero entry in 2,2'
                stop 1
             end if
             multiplier = matrix(3,2)/matrix(2,2)
             matrix(3,:) = matrix(3,:) - multiplier*matrix(2,:)
             rhs(3)      = rhs(3)      - multiplier*rhs(2)

             if (matrix(3,3) == 0.d0) then
                print *, i,k
                print *, matrix(1,:)
                print *, matrix(2,:)
                print *, matrix(3,:)
                print *, 'zero entry in 3,3'
                stop 1
             end if

             long_waves(i,k) = rhs(3) / matrix(3,3)
             bs(i,k) = (rhs(2) - long_waves(i,k)*matrix(2,3)) / matrix(2,2)
             as(i,k) = (rhs(1) - long_waves(i,k)*matrix(1,3) - &
                                    bs(i,k)*matrix(2,1)) / matrix(1,1)

          else

             long_waves(i,k) = 0.d0

          end if

       end do
    end do

    !$omp end do
    !$omp end parallel

    short_waves = 0.d0
    where(ocean_mask == 1)
       short_waves = raw_data - long_waves
    end where

  end subroutine filter_waves

end module plume

