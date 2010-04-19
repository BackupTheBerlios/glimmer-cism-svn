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

  ! model axes in code are 
  !   x/i/m - east (parallel to isobaths of ice shelf base, perp. to inflow)
  !   y/k/n - north (perp. to isobaths of ice shelf base, parallel to inflow)
  !   z     - upwards from the sea bed (defined as gldep+wcdep)


  use plume_global
  use plume_functions
  use plume_context
  use plume_io

  implicit none

  ! everything is private by default
  private

  ! except for functions mentioned here
  public :: plume_initialise,plume_runstep,plume_finalise,plume_iterate
  public :: get_nsteps

  ! private module variables

  integer :: iwetmin,iwetmax,kwetmin,kwetmax,ikcplus
  integer :: icalcan,kcalcan,icalcen,kcalcen

  logical :: sepflag

  integer :: varoutrat
  integer :: nsteps
  integer :: last_output_step = -1, last_lnot_step = -1


  real(kind=kdp) ::  tottim, &  ! total time to run model
       outtim, &  ! time interval between writing model state out to file
       fouttim, & ! first time at which to write out model outout
       louttim, & ! last time to write out model output
       snottim, & ! short note time interval
       lnottim, & ! long note time interval
       runtim,  & ! accumulated simulation time
       labtim     ! time interval in days used to name output files

  real(kind=kdp) :: negdep ! total of negative depths (error diagnostic)

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
    gdt = g*dt1
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
       deallocate(bpos_ext)
    else
       call initialise_fields(suppress_ascii_output,iwetmin,iwetmax,kwetmin,kwetmax)
    end if

    ! initialise wet area to surround all inflow points
    ikcplus = 2  ! set additional area in which to check wetness

    if (.not.use_min_plume_thickness) then
       icalcan = max0(2,iwetmin-ikcplus)
       kcalcan = max0(2,kwetmin-ikcplus)
       icalcen = min0(m_grid-1,iwetmax+ikcplus)
       kcalcen = min0(n_grid-1,kwetmax+ikcplus)
    else
       icalcan = domain_imin + 1
       kcalcan = domain_kmin + 1
       icalcen = domain_imax - 1
       kcalcen = domain_kmax - 1
    end if

    ! initialise time, total of negative depths, and separation and negative
    ! frazil warning counters
    runtim = 0.d0
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

!!! TODO: this code will not work with the new netcdf output
    ! if restarting from an old run
    ! -----------------------------
    if (restart) then

       ! read old fields
       open(31,file = trim(restart_data_filename))

       read(31,*) runtim,itmp,itmp,ltmp,ltmp,ltmp,ltmp,itmp, &
            &     itmp,itmp,dtmp,dtmp,iwetmin,iwetmax,kwetmin,kwetmax,negdep

       do i = 1,m_grid
          read(31,*) dtmp
       end do
       do k = 1,n_grid
          read(31,*) dtmp
       end do

       do i = 1,m_grid
          do k = 1,n_grid
             read(31,*) jcw(i,k),su(i,k),sv(i,k),temp(i,k),salt(i,k), &
                  pdep(i,k),bpos(i,k),ipos(i,k), &
                  rhoamb(i,k),rhop(i,k),dtmp,dtmp, &
                  entr(i,k),atemp(i,k),asalt(i,k), &
                  bmelt(i,k),btemp(i,k),bsalt(i,k), &
                  tf(i,k)
          end do
       end do

       if (vardrag) then
          do i = 1,m_grid
             do k = 1,n_grid
	        read(31,*) drag(i,k)
             end do
          end do
       end if

       if (frazil) then
          do i = 1,m_grid
             do k = 1,n_grid
	        read(31,*) ctot(i,k),dtmp
		do l = 1,nice
		   read(31,*) c(i,k,l),fmelt(i,k,l),fppn(i,k,l),fnuc(i,k,l)
                   ca(i,k,l) = c(i,k,l)
		end do
             end do
          end do
       end if

       if (intrace) then
          do i = 1,m_grid
             do k = 1,n_grid
		read(31,*) intrin(i,k)
		do l = ninfmin,ninfmax
                   read(31,*) intr(i,k,l)
                   intra(i,k,l) = intr(i,k,l)
		end do
             end do
          end do
       end if

       if (tangle) then
          do i = 1,m_grid
             do k = 1,n_grid
		read(31,*) u0(i,k),v0(i,k),u0a(i,k),v0a(i,k),tang(i,k)
             end do
          end do
       end if

       close(31)

       ! populate other arrays

       do i = 1,m_grid
          do k = 1,n_grid
             ua(i,k) = su(i,k)*5.0d-1*(pdep(i,k) + pdep(i+1,k))
             va(i,k) = sv(i,k)*5.0d-1*(pdep(i,k) + pdep(i,k+1))
             tempa(i,k) = temp(i,k)
             salta(i,k) = salt(i,k)
             ctota(i,k) = ctot(i,k)
          end do
       end do

       ! sort out timestep

       if (runtim.lt.dtswtim) then
          nsteps = int((dtswtim - runtim)/dt1 + (tottim - dtswtim)/dt2)
          dt = dt1
       else
          nsteps = int((tottim - runtim)/dt2)
          dt = dt2
       end if

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
       min_run_time, &
       run_plume_to_steady, &
       write_all_states)

    ! NB: This subroutine is intended to be called only from a 
    ! glimmer-CISM driver that uses the plume to calculate the
    ! basal melt rate.
    ! The purpose of this subroutine is to run the plume
    ! over a period of time ice_dt, stopping if a steady-state
    ! is reached.  If run_plume_to_steady is
    ! set to .true. then the plume iterates to steady even if it
    ! takes longer than ice_dt.

    real(kind=kdp),intent(in) :: time,ice_dt !in years
    real(kind=kdp),dimension(:,:),intent(in) :: lsrf,ice_dz !in meters
    logical,dimension(:,:),intent(in) :: landmask !equal to 1 in land areas or grounded ice
    real(kind=kdp),dimension(:,:),intent(in) :: t_interior !in Celcius
    real(kind=kdp),dimension(:,:),intent(out) :: bmelt_out,btemp_out !in meters per year
    logical,intent(in) :: run_plume_to_steady, write_all_states
    real(kind=kdp),intent(in) :: plume_stopping_tol ! percent change in meltrate need to continue plume time-stepping
    real(kind=kdp),intent(in) :: min_run_time !in days

    !local variables
    real(kind=kdp) :: subcycling_time,ice_dt_in_sec
    real(kind=kdp),dimension(m_grid,n_grid) :: bmelt_old !in meters per year
    real(kind=kdp) :: max_rel_bmelt_change,prev_rel_change,min_run_time_sec
    logical :: reached_steady
    character(len=512) :: log_message

    ! convert lower surface depth into height of basal surface and interface
    ! surface, storing results in plume's global variable
    bpos = lsrf + gldep + wcdep
    ipos = bpos - pdep

    where( landmask )
       jcs = 0
    elsewhere
       jcs = 1
    end where

    !TODO: assign t_interior and ice_dz to globals 
    ! so that thermodynamics will pick up heat conduction into ice

    bmelt_old = bmelt

    subcycling_time = 0.0
    ice_dt_in_sec = ice_dt * 365.25d0*24.0d0*3600.d0
    min_run_time_sec = min_run_time * 3600.d0 * 24.d0
    reached_steady = .false.
    prev_rel_change = 1.d0

    ! while not steady

    do while ((subcycling_time .le. ice_dt_in_sec) .or. run_plume_to_steady )


       call plume_runstep()

       !we need to set the bmelt along the domain boundaries 
       !because the plume doesn't calculate the rates there
       bmelt(:,domain_kmax) = bmelt(:,domain_kmax - 1) !North edge
       bmelt(:,domain_kmin) = bmelt(:,domain_kmin + 1) !South edge
       bmelt(domain_imin,:) = bmelt(domain_imin + 1,:) !West edge
       bmelt(domain_imax,:) = bmelt(domain_imax - 1,:) !East edge	

       if (write_all_states) then
          call plume_netcdf_write_vars(time*3600.0d0*24.0d0*365.25d0 + subcycling_time)
       end if

       ! calculate max error
       max_rel_bmelt_change = maxval(abs(bmelt_old - bmelt)/ &
            (abs(bmelt)+1.d-30))

       if (subcycling_time >= min_run_time_sec .and. (max_rel_bmelt_change/prev_rel_change < 1.d-1)) then
          prev_rel_change = max_rel_bmelt_change
          write(*,*) 'max_rel_bmelt_change',max_rel_bmelt_change, '   stopping tol',plume_stopping_tol
       end if

       if (max_rel_bmelt_change < plume_stopping_tol) then
          if (subcycling_time .ge. (min_run_time_sec)) then
             reached_steady = .true.
             exit
          end if
       end if

       subcycling_time = subcycling_time + dt
       bmelt_old = bmelt

    end do

    btemp_out = btemp
    bmelt_out = bmelt * (365.25d0*24.0d0*3600.0d0)

    if (.not. reached_steady) then
       call io_append_output('plume did not reach steady state')
    end if

    !we do this to make the plume output timestamps indicate the ice state
    ! with respect to which the plume was steady
    runtim = time*3600.0d0*24.0d0*365.25d0 

    write(log_message, '(a,f6.1)') 'subcycling time in days', subcycling_time/(3600.0*24.0)
    call io_append_output(trim(log_message))

    call io_write_long_step_output(icalcan,icalcen,kcalcan,kcalcen,&
         varoutrat,negdep)

    call io_write_surface_output(runtim,labtim)    


  end subroutine plume_iterate

  ! **************************************************************************
  ! ******* main subroutine to carry out the current timestep ****************
  ! **************************************************************************

  subroutine plume_runstep()

    !The order in which helper routines are called is:
    ! momentum
    ! bound_u_v
    ! continuity
    ! bmlt boundary condition
    ! inflow_calc
    ! update
    ! bound_interface
    ! bound_u_v
    ! scalar
    ! rho_calc
    ! bound_tsdc

    ! change timestep and recalculate basic quantities if necessary

    implicit none

    if (runtim.eq.dtswtim) then
       dt = dt2
       gdt = g*dt2
       fdt = f*dt2
    end if

    ! update time in seconds

    runtim = runtim + dt

    ! find indices of current wetted area

    if (.not.use_min_plume_thickness) then
       icalcan = max0(2,iwetmin-ikcplus)
       kcalcan = max0(2,kwetmin-ikcplus)
       icalcen = min0(m_grid-1,iwetmax+ikcplus)
       kcalcen = min0(n_grid-1,kwetmax+ikcplus)
    end if

    ! --------------------
    ! calculate velocities
    ! --------------------


    call momentum(icalcan,kcalcan,icalcen,kcalcen)

    ! set velocity components on the open boundaries if the plume 
    ! has reached that far

    call bound_u_v(icalcan,icalcen,kcalcan,kcalcen)

    ! ---------------------------------------------------
    ! evaluate continuity equation for interface position
    ! ---------------------------------------------------
    call continuity(icalcan,kcalcan,icalcen,kcalcen)   

    ! interface and scalars passed forward on open boundary
    call inflow_calc(icalcan,kcalcan,icalcen,kcalcen)

    ! ----------------------------------------------------------
    ! update plume thickness, wet/dry boundaries, and velocities
    ! ----------------------------------------------------------

    call update(iwetmin,iwetmax,kwetmin,kwetmax,&
         icalcan,icalcen,kcalcan,kcalcen,negdep)


    ! update interface position and flow on open boundaries

    call bound_interface(icalcan,icalcen,kcalcan,kcalcen)


    call bound_u_v(icalcan,icalcen,kcalcan,kcalcen)


    ! ---------------------------------------------------------------
    ! solve transport equations for temperature, salinity, and frazil
    ! ---------------------------------------------------------------
    call scalar(icalcan,kcalcan,icalcen,kcalcen)             


    ! -----------------------------------
    ! calculate density and update bounds
    ! -----------------------------------

    call rho_calc(icalcan,icalcen,kcalcan,kcalcen,sepflag)


    call bound_tsdc(icalcan,icalcen,kcalcan,kcalcen)


    ! write short step output and reset separation warning flag

    if (mod(runtim,snottim).eq.0) then
       sepflag = .false.
       call io_write_calculated_time(runtim,int(runtim/dt))
    end if

    ! write long step output 
    if ((int(floor(runtim/outtim)) > last_lnot_step) .and. (.not. (in_glimmer))) then
       last_lnot_step = int(floor(runtim/outtim))
       call io_write_long_step_output(icalcan,icalcen,kcalcan,kcalcen,&
            varoutrat,negdep)
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

  real function get_nsteps()
    get_nsteps = nsteps
  end function get_nsteps

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
    real(kind=kdp) :: infloain,infloein,cseedfix,cinffix	

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
         ,  basmelt &
         ,  rholinear &
         ,  thermobar &
         ,  intrace &
         ,  vardrag &
         ,  topedit &
         ,  tangle &
         ,  negfrz &
         ,  use_min_plume_thickness &
         ,  tottim &
         ,  outtim &
         ,  labtim &
         ,  snottim &
         ,  lnottim &
         ,  dt1 &
         ,  m_grid & 
         ,  n_grid &
         ,  hx &
         ,  hy &
         ,  gldep &
         ,  ifdep &
         ,  wcdep &
         ,  plume_min_thickness &
         ,  context &
         ,  bathtype &
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
         ,  nus &
         ,  nbar &
         ,  nice &
         ,  seedtype &
         ,  cseedfix &
         ,  cinffix

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
    mixlayer    = .false. ! model a mixed-layer, rather than a plume
    ! (i.e. whole domain is given initial thickness)
    in_glimmer = .false.  ! by default, not running inside glimmer (ice shelf model)
    restart     = .false. ! restart from previous model dump
    restart_data_filename = '' ! file to read from for a restart
    frazil      = .false. ! include frazil ice
    nonlin      = .true.  ! advection terms included
    horturb     = .false. ! horizontal diffusion included
    entrain     = .true.  ! entrainment included
    entype      = 1       ! entrainment parameterisation:
    ! 1 - full kochergin formulation
    ! 2 - fractional kochergin formulation (using ef)
    ! 3 - halved pedersen formulation
    ! 4 - halved and modified pedersen formulation
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
    outtim  = 050.0d0   ! file output frequency in days
    labtim  = 001.0d0   ! units in which to name files in days
    snottim = 005.0d0   ! short note output frequency in days
    lnottim = 050.0d0   ! long note output frequency in days

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
    plume_min_thickness = 00.01d0 ! artificially determined minimum plume thickness

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
    cross_slope_wavenumber = 0
    along_slope_deepening_exp = 0
    channel_amplitude = 0.d0
    random_amplitude = 0.d0

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
    knfloa = 1          ! southern boundary (cells)
    knfloe = 3          ! northern boundary (cells)

    ! set ambient fluid properties
    ! ----------------------------
    namb = 302           ! increments in ambient water column (minimum +2)
    dzincr = 10.d0       ! depth of each increment 
    ! (so total (namb - 2)*dzincr = gldep + wcdep)
    salttop = 34.500d0   ! sea surface salinity 
    saltbot = 34.950d0   ! sea bed salinity (gldep + wcdep)
    temptop = -1.900d0   ! sea surface temperature 
    tempbot = -2.500d0   ! sea bed temperature (gldep + wcdep)

    ! set physical and geographical parameters
    ! ----------------------------------------
    g = 9.81d0       ! gravitational constant
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
    small = 1.0d-15  ! smallest number

    ! +++++++++++++++++++++
    ! read values from file
    ! +++++++++++++++++++++

    open(21,file=trim(nl_filename),status='old')

    read(21,plume_nml)

    close(21)

    if (use_min_plume_thickness .and. mixlayer) then
       call io_append_output('Can not run using useminthickness = .t. &
            and mixlayer = .t.')
       stop
    end if

    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! some parameters dependent upon namelist-read variables
    ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! convert timesteps to seconds

    tottim  = tottim*24.d0*3600.d0  ! total simulation time in seconds
    outtim  = outtim*24.d0*3600.d0  ! file output frequency in seconds
    labtim  = labtim*24.d0*3600.d0  ! units in which to name files
    snottim = snottim*24.d0*3600.d0 ! short note output frequency in seconds
    lnottim = lnottim*24.d0*3600.d0 ! long note output frequency in seconds
    fouttim = outtim                ! first file output in seconds
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

  end subroutine set_parameters

  subroutine initialise_fields(suppress_ascii_output,iwetmin,iwetmax,kwetmin,kwetmax,bpos_ext)

    ! read data and set all initial fields
    implicit none

    logical,intent(in) :: suppress_ascii_output
    !intent(inout) means that if this subroutine doesn't change *wetmin/max 
    ! then the calling subroutine will keep the original values
    integer,intent(inout) :: iwetmin,iwetmax,kwetmin,kwetmax
    real(kind=kdp),dimension(:,:),intent(in),optional :: bpos_ext

    ! local variables
    integer :: i,k	
    real(kind=kdp),dimension(m_grid,n_grid):: zd
    real(kind=kdp),dimension(m_grid,n_grid):: depth,sambindep, &
         tambindep,c1,c2,c3,rhoa

    ! initialise all arrays (0 is assigned to all array locations)

    u = 0.d0
    v = 0.d0
    ua = 0.d0
    va = 0.d0
    su = 0.d0
    sv = 0.d0
    u0 = 0.d0
    v0 = 0.d0
    u0a = 0.d0
    v0a = 0.d0
    tang = 0.d0
    pdep = 0.d0
    ipos = 0.d0
    bpos = 0.d0
    jcs = 0
    jcw = 0
    jcd_u = 0
    jcd_v = 0
    jcd_fl = 0       
    jcd_negdep = 0
    jcd_fseed = 0
    ctot = 0.d0
    tf = 0.d0
    entr = 0.d0
    thk_def = 0.d0
    atemp = 0.d0
    asalt = 0.d0
    drag = 0.d0
    bmelt = 0.d0
    btemp = 0.d0
    bsalt = 0.d0
    ctempd = 0.d0
    tint = tiuniform

    ! initialise ice to zero even when frazil is off so that densities are correct
    c = 0.d0
    ca = 0.d0
    fmelt = 0.d0
    fppn = 0.d0
    fnuc = 0.d0

    tempinf = 0.d0
    saltinf = 0.d0
    depinf = 0.d0

    if (intrace) then
       intrin = 0
       intr = 0.d0
    end if

    ! get cell dimensions
    ! -------------------

    call grid_set()

    ! get topography and inflow regions, either from data or settings
    ! ---------------------------------------------------------------


    if (in_glimmer) then
       if (.not. present(bpos_ext)) then
          call io_append_output('Need to provide bpos_ext')
          stop
       end if
       bpos = bpos_ext
       !if (.not. present(tint_ext)) then
       !   call io_append_output('Need to provide tint_ext')
       !   stop
       !end if
       !int = tint_ext

    else	
       if (bathtype.gt.0) then
          call topog_depth_inflow_set(.not. use_min_plume_thickness)
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

       saltinf = sambindep + meltinf*(saltinf - sambindep)
       tempinf = tambindep + meltinf*(tempinf - tambindep)
    end where

    ! initialise whole domain to wet and set plume initial thickness if mixed-layer
    ! model

    if (mixlayer) then

       iwetmin = 1
       iwetmax = m_grid
       kwetmin = 1
       kwetmax = n_grid

       where (bpos > 0.d0 .and. bpos < (wcdep+gldep))
          pdep = depinit
       end where

    else if (use_min_plume_thickness) then

       ! using minimum thickness, so we will never 
       ! look at iwetmin/iwetmax/kwetmin/kwetmax
       continue
    else

       ! not using minimum thickess or mixlayer
       iwetmin = m_grid
       iwetmax = 1
       kwetmin = n_grid
       kwetmax = 1

       do k = 1,n_grid
          do i = 1,m_grid
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

    ! set fields to ambient-fluid properties for depth
    ! (used when considering newly-wet cells)

    zd = max(0.d0,wcdep + gldep - bpos)

    tempa = get_tamb_z(zd)
    salta = get_samb_z(zd)
    rhoa = get_rhoamb_z(zd)
    temp = tempa
    salt = salta
    rhoamb = rhoa
    rhop = rhoa

  end subroutine initialise_fields

  subroutine grid_set()

    implicit none

    ! read grid cell data

    ! local variables
    integer :: i,k,imid,kmid

    real(kind=kdp) :: rmdx,adx,adxu
    real(kind=kdp) ::  rmdxu,rmdy,rmdyv,ady,adyv

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

    ahdx = adx/rdx
    ahdxu = adxu/rdxu

    rmdy = rdy(kmid)
    rmdyv = rdyv(kmid)
    ady = dt*ah*rmdy**2
    adyv = dt*ah*rmdyv**2

    ahdy = ady/rdy
    ahdyv = adyv/rdyv

  end subroutine grid_set

  subroutine set_ambient(suppress_ascii_output)

    ! set ambient profiles of temperture, salinity and density
    ! density is unaffected by frazil as the ambient is ice-free

    implicit none

    logical,intent(in) :: suppress_ascii_output

    ! local variables
    integer :: i
    real(kind=kdp),dimension(namb) :: depth,pressure,ttt,rhopot

    depth = (/ (0.d0 + i*dzincr,i=0,(namb-1)) /)

    tamb = temptop + depth*tgrad
    samb = salttop + depth*sgrad
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

    if (.not. suppress_ascii_output)  call io_write_amb(namb,depth,tamb,ttt,samb,rhovf,rhopot)

  end subroutine set_ambient

  subroutine frazil_calc()

    implicit none

    ! set up various parameter arrays use in the frazil calculations

    ! local variables

    integer :: l
    real(kind=kdp) :: cdc,winew,rey,mstar

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
       winew = dsqrt((4.d0*g*ar*r(l)*(rho0-rhoi))/(rho0*cdc))
       do while (abs(wi(l)-winew).gt.small) 
          wi(l) = winew
          rey = 2.d0*winew*r(l)/nu0
          cdc = 10**(1.386d0 - 0.892d0*dlog10(rey) &
               + 0.111d0*(dlog10(rey))**2)
          winew = dsqrt((4.d0*g*ar*r(l)*(rho0-rhoi))/(rho0*cdc))
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

  subroutine inflow_calc(icalcan,kcalcan,icalcen,kcalcen)

    implicit none

    ! local variables

    integer :: i,k,l,icalcan,kcalcan,icalcen,kcalcen

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
                      c(i,k,l) = cinf(l)
                      ca(i,k,l) = cinf(l)
                   end do
                   ctot(i,k) = cinftot
                   ctota(i,k) = cinftot
                   rhop(i,k) = (1.d0 - cinftot)*rhop(i,k) + cinftot*rhoi
                end if

                ! set inflow tracers

                if (intrace) then
                   do l = ninfmin,ninfmax
                      if (intrin(i,k).eq.l) then
                         intr(i,k,l) = 1.d0
                         intra(i,k,l) = 1.d0
                      else
                         intr(i,k,l) = 0.d0
                         intra(i,k,l) = 0.d0
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

    integer :: icalcan,kcalcan,icalcen,kcalcen

    ! local variables
    integer :: i,k,l,mflag

    real(kind=kdp),dimension(m_grid,n_grid) :: pdepc,bspeed,speed
    real(kind=kdp) :: rhopac,rhoa,delrho,rhoq,redg,tt,vmid,umid
    real(kind=kdp) :: rich,sm,arg,delta,iold
    real(kind=kdp) :: dragrt,prden,scden,gambt,gambs,tfb,tfi,c1,c2,c3
    real(kind=kdp) :: deltam(m_grid,n_grid),fppntot,ucrit,ucl

    ! 0. preliminaries
    ! ----------------
    prden = 12.5d0*pr**(2.d0/3.d0) - 9.d0
    scden = 12.5d0*sc**(2.d0/3.d0) - 9.d0

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
	  if (jcs(i,k) .ne. 1) cycle ! don't calculate anything over grounded ice
          deltam(i,k) = 0.d0
          pdepc(i,k) = bpos(i,k) - ipos(i,k)
          ! find flow speed on scalar grid
          tt = 5.0d-1*dy(k)*rdyv(k)
          vmid = tt*sv(i,k-1) + (1.-tt)*sv(i,k) 
          tt = 5.0d-1*dx(i-1)*rdxu(i)
          umid = tt*su(i,k) + (1.-tt)*su(i-1,k)
          speed(i,k) = umid**2 + vmid**2 + small
          bspeed(i,k) = sqrt(speed(i,k))
       end do
    end do

    ! 1. entrainment
    ! --------------

    if (entrain) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle ! don't calculate anything over land/grounded ice

             if(pdepc(i,k).gt.edepth) then          

                rhopac = rhop(i,k)
                rhoa = rhoamb(i,k)
                delrho = rhoa - rhopac
                rhoq = 5.0d-1*(rhopac+rhoa)
                redg = delrho/rhoq
                rich = g*redg*pdepc(i,k)/speed(i,k)
                rich = dmax1(5.0d-2,rich)    

                ! full kochergin entrainment
                if (entype.eq.1) then
                   sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
                        &             3.16d-1*rich + 3.46d-2)))
                   arg = dmax1(small,speed(i,k) + g*redg*pdepc(i,k)/sm)
                   entr(i,k) = cl**2*dsqrt(arg)/sm                
                end if

                ! reduced kochergin entrainment matched to pedersen
                if (entype.eq.2) then
                   sm = rich/(7.25d-1*(rich + 1.86d-1 - dsqrt(rich*rich- &
                        &             3.16d-1*rich + 3.46d-2)))
                   arg = dmax1(small,speed(i,k) + g*redg*pdepc(i,k)/sm)
                   entr(i,k) = ef*cl**2*(drag(i,k)/3.0d-3)*dsqrt(arg)/sm    
                end if

                ! half pedersen entrainment
                if (entype.eq.3) then
                   entr(i,k) = 3.6d-2*bspeed(i,k)*drag(i,k)/rich
                end if

                ! half modified pedersen entrainment
                if (entype.eq.4) then
                   entr(i,k) = 3.6d-2*bspeed(i,k)*dsin(1.0d-3)
                end if

             end if

             deltam(i,k) = deltam(i,k) + entr(i,k)

          end do
       end do
    end if

    ! 2. basal melting
    ! ----------------

    if (basmelt) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen

             if (jcs(i,k) .ne. 1) cycle ! skip land/ grounded ice

             if (pdepc(i,k).gt.mdepth) then         

                dragrt = dsqrt(drag(i,k))

                ! find turbulent exchange coefficients
                gambt = (dragrt*bspeed(i,k))/ &
                     & (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + prden)
                gambs = (dragrt*bspeed(i,k))/ &
                     & (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + scden)

                ! calculate freezing point of plume at shelf base,
                ! decide if melt (mflag = 1) 
                ! or freeze and calculate freezing point of ice at shelf base
                tfb = fta*salt(i,k) + ftb + ftc*(gldep + wcdep -bpos(i,k))
                mflag = (1 + int(sign(1.d0,temp(i,k) - tfb)))/2
                tfi = (1 - mflag)*fta*si  &
                     + ftb + ftc*(gldep + wcdep - bpos(i,k))

                ! calculate coefficients in quadratic to be solved
                c1 = lat/c0 + mflag*(ci/c0)*(tfi-tint(i,k))
                c2 = gambs*(lat/c0 + mflag*(ci/c0)*(tfb-tint(i,k)))  &
                     &   + gambt*(tfi-temp(i,k))
                c3 = gambs*gambt*(tfb - temp(i,k))

                ! calculate melt rate
                bmelt(i,k) = -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)


                !! multiply by 10 if freezing
                !              if (mflag.eq.0) then
                !                bmelt(i,k) = bmelt(i,k)*10.d0
                !              end if

                deltam(i,k) = deltam(i,k) + bmelt(i,k)

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
                   if (c(i,k,l).gt.0.d0) then

                      ! calculate precipitation of frazil
                      ucl = &
                           dsqrt((1.0d-1*g*re(l)*(rho0-rhoi))/(rho0*drag(i,k)))
                      ucrit = (1.d0 - bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                      fppn(i,k,l) = &
                           dmin1(-(rhoi/rho0)*c(i,k,l)*wi(l)*ucrit,0.d0)

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

          delta = dt*(rdxu(i)*(u(i-1,k)-u(i,k)) &
               + rdyv(k)*(v(i,k-1)-v(i,k))) + &
               dt* deltam(i,k)

          if (use_min_plume_thickness) then
             thk_def(i,k) = max(0.d0,plume_min_thickness - (pdepc(i,k)+delta))
             if (thk_def(i,k) > 0.d0) then
                ! make up for thickness deficiency by increasing entrainment
                ! articicially
                entr(i,k) = entr(i,k) + thk_def(i,k)/dt
                delta = delta + thk_def(i,k)
             end if
             !this is for output purposes so we can see what percentage of the thicknes
             !change was due to the imposed minimum thickness
             thk_def(i,k) = thk_def(i,k) / (entr(i,k)*dt)

          end if

          ! update interface position (NB - this will change when ice surface moves
          ! in the future coupled model)
          iold = ipos(i,k)
          ipos(i,k) = ipos(i,k) - delta

          ! check for negative predicted depth of ambient fluid 
          jcd_negdep(i,k) = 0
          if (bpos(i,k) < (iold - 2.d0*delta)) jcd_negdep(i,k) = 2
          if (bpos(i,k) < (iold - 3.d0*delta)) jcd_negdep(i,k) = 1

       end do
    end do

  end subroutine continuity


  subroutine momentum(icalcan,kcalcan,icalcen,kcalcen)

    ! calculates velocity components within plume

    ! all turning angle theory by alex wilchinsky
    ! (wilchinsky, feltham and holland, submitted to j. phys. oceanogr.)
    ! it is probably wise to check or reprogram this before use

    implicit none

    ! local variables

    integer :: i,k,icalcan,kcalcan,icalcen,kcalcen,jcvfac
    integer(kind=kdp) :: idel,kdel,iidx,kkdy,ihilf,khilf

    real(kind=kdp) :: one,termnl,termnl2,corx
    real(kind=kdp) :: redgu,slorho,islope
    real(kind=kdp) :: pdepu,zu,zum,arfac,tt,uu,umid,vmid
    real(kind=kdp) :: speed,tbotm,rhoa,tu,salu
    real(kind=kdp) :: rhoc,rhoe,rhou,rhoq,dxx,dyy,sx,sy,sxy,r1,r2
    real(kind=kdp) :: tlate,tlatw,tlats,tlatn,hordif,cory,redgv
    real(kind=kdp) :: pdepv,zv,zvm,tv,salv,rhon,rhov
    real(kind=kdp) :: delta,ctotu,ctotv,dragu,dragv

    logical :: skip_u_calc
    logical :: olddrag,norotation,variableekman,draginmelt 
    logical :: newudrag(m_grid,n_grid)
    real(kind=kdp) :: av,ekthick,kq,costang,sintang,thickratio
    real(kind=kdp) :: ugriddrag(m_grid,n_grid)

    sintang = 0.0
    kq = 0.0

    ! set switches for turning angle model

    olddrag       = .false. ! use the old drag magnitude
    norotation    = .false. ! switch off drag rotation
    variableekman = .false. ! use variable vertical eddy visc (hence ek thick)
    draginmelt    = .true.  ! write drag array so that drag used in melting
    ! (only works with spatially constant drag coeff)

    one = 1.d0

    ! calculate ekman layer thickness for constant vertical eddy diffusivity

    if (.not.variableekman) then
       av = 5.d-4
       ekthick = dsqrt(-2.d0*av/f) 
    end if

    ! start main loop

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

          if (k == kcalcen) then
             !!write(*,*) 'North edge'
          end if

          skip_u_calc = .false.

          ! skip cell if dry land
          if (jcs(i,k).lt.1) cycle

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

          ! plume thickness on the u-grid
          pdepu = 5.0d-1*(pdep(i,k) + pdep(i+1,k))

          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! skip u-component if eastern cell is dry land
          ! NB: this is a no normal flow condition 
          if (jcs(i+1,k).ne.1) then
             skip_u_calc = .true.
          end if

          islope = ipos(i+1,k) - ipos(i,k)

          ! final wet/dry logic
          if ((jcw(i,k).le.0).and.(jcw(i+1,k).le.0)) then
             ! current cell and eastern neighbour are dry
             skip_u_calc = .true.
          else if ((jcw(i,k).le.0).and.(islope.gt.0.d0)) then
             ! current cell is dry and interface slopes up to east
             skip_u_calc = .true.
          else if ((jcw(i+1,k).le.0).and.(islope.le.0.d0)) then
             ! eastern neighbour is dry and interface slopes up to the west
             skip_u_calc = .true.
          endif

          if (.not.(skip_u_calc)) then
             ! control of negative depth effects by enhanced friction
             jcvfac = 0
             arfac = 1.d0
             jcvfac = max0(jcd_negdep(i,k),jcd_negdep(i+1,k))
             if (jcvfac .ge. 1) arfac = 75.d0*jcvfac  

             ! velocity components on the u-grid
             tt = 5.0d-1 
             uu = 5.0d-1*dy(k)*rdyv(k)
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

                ! set up stuff for rewriting drag if necessary

                if (draginmelt) then
                   dragu = cdb
                   ugriddrag(i,k) = 5.0d-1*(drag(i,k) + drag(i+1,k))
                   newudrag(i,k) = .false.
                endif

                ! find ekman layer thickness if necessary and thickness ratio

                if (variableekman) then
                   av = abs(0.16d0*0.4d0*dragu*speed / (f*2.71828d0))
                   ekthick = dsqrt(-2.d0*av/f) 
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
                tang(i,k) = atan2( u0(i,k)-v0(i,k) , -u0(i,k)-v0(i,k) ) &
                     - atan2( u0(i,k)+v0(i,k) , thickratio+u0(i,k)-v0(i,k) )
                !      
                if (norotation) tang(i,k) = 0.0 
                !         
                costang = dcos(tang(i,k))
                sintang = dsin(tang(i,k))
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

             ! if turning angle add other part of drag here because uses v not u

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

                sxy = sx -sy
                r1 = (sign(one,sxy) + 1.d0)*5.0d-1
                r2 = one - r1
                termnl=(r1*sy + r2*sx)*(ua(ihilf,khilf)*dble(jcd_u(ihilf,khilf)) &
                     + dble(1 - jcd_u(ihilf,khilf))*ua(i,k)) &
                     + r1*sxy*(ua(ihilf,k)*dble(jcd_u(ihilf,k)) &
                     + dble(1 - jcd_u(ihilf,k))*ua(i,k)) &
                     - r2*sxy*(ua(i,khilf)*dble(jcd_u(i,khilf)) &
                     + dble(1 - jcd_u(i,khilf))*ua(i,k)) &
                     - (r2*sy + r1*sx)*ua(i,k)

                termnl2 = - ua(i,k)*dt*((su(ihilf,k) - su(i,k))*idel/dxx  &
                     + (sv(iidx,k) - sv(iidx,k-1))*rdyv(k))

             end if
             !          
             ! 7)lateral shear stress terms (east component)
             ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             if (horturb) then

                tlate = ahdxu(i+1)*(su(i+1,k) - su(i,k))
                tlatw = ahdxu(i)*(su(i,k) - su(i-1,k))
                tlats = ahdy(k-1)*(su(i,k) - su(i,k-1))
                tlatn = ahdy(k)*(su(i,k+1) - su(i,k))

                hordif = pdepu*((tlate-tlatw)*rdx(i) + (tlatn-tlats)*rdyv(k))

             end if

             ! final
             !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             u(i,k) = (ua(i,k)+corx+islope+slorho+termnl+termnl2+hordif)*tbotm

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

          ! plume thickness on the v-grid
          pdepv = 5.0d-1*(pdep(i,k) + pdep(i,k+1))

          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          ! skip v-component if northern cell is dry land
          if (jcs(i,k+1).le.0) cycle

          islope = ipos(i,k+1) - ipos(i,k)

          ! final wet/dry logic
          if ((jcw(i,k).le.0).and.(jcw(i,k+1).le.0)) then
             ! current cell and northern neighbour are dry
             cycle
          else if ((jcw(i,k).le.0).and.(islope.gt.0.d0)) then
	     ! current cell is dry and interface is rising to the north
             cycle
          else if ((jcw(i,k+1).le.0).and.(islope.le.0.d0)) then
             ! northern neighbour is dry and interface is rising to the south
             cycle
          endif
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
             ! calculate nonlinear terms
             sxy = sx -sy
             r1 = (dsign(one,sxy) + 1.d0)*5.0d-1
             r2 = one - r1
             termnl=(r1*sy + r2*sx)*(va(ihilf,khilf)*dble(jcd_v(ihilf,khilf)) &
                  + dble(1 - jcd_v(ihilf,khilf))*va(i,k)) &
                  + r1*sxy*(va(ihilf,k)*dble(jcd_v(ihilf,k)) &
                  + dble(1 - jcd_v(ihilf,k))*va(i,k)) &
                  - r2*sxy*(va(i,khilf)*dble(jcd_v(i,khilf)) &
                  + dble(1 - jcd_v(i,khilf))*va(i,k)) &
                  - (r2*sy + r1*sx)*va(i,k)

             termnl2 = - va(i,k)*dt*((sv(i,khilf) - sv(i,k))*kdel/dyy  &
                  + (su(i,kkdy) - su(i-1,kkdy))*rdxu(i))
          end if
          !          
          ! 7)lateral shear stress terms (north component)
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          if (horturb) then

             tlate = ahdx(i)*(sv(i+1,k) - sv(i,k))
             tlatw = ahdx(i-1)*(sv(i,k) - sv(i-1,k))
             tlats = ahdyv(k)*(sv(i,k) - sv(i,k-1))
             tlatn = ahdyv(k+1)*(sv(i,k+1) - sv(i,k))

             hordif = pdepv*((tlate-tlatw)*rdxu(i) + (tlatn-tlats)*rdy(k))

          end if
          !          
          ! final
          !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          v(i,k) = (va(i,k)+cory+islope+slorho+termnl+termnl2+hordif)*tbotm

       end do
    end do
    !     
    ! write new drag on scalar grid if necessary 
    !     
    if (tangle.and.draginmelt) then
       !     
       do i = icalcan+1,icalcen
          do k = kcalcan,kcalcen
             if (jcs(i,k) .ne. 1) cycle
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
          if (jcs(i,k).ne.1) cycle
          pdepu = 5.0d-1*(pdep(i,k) + pdep(i+1,k)) + small
          pdepv = 5.0d-1*(pdep(i,k) + pdep(i,k+1)) + small
          su(i,k) = u(i,k)/pdepu
          sv(i,k) = v(i,k)/pdepv
          ! test for negative depths and set flow to zero if so
          delta = dt*(rdxu(i)*(u(i-1,k)-u(i,k)) &
               + rdyv(k)*(v(i,k-1)-v(i,k))) 

          if (bpos(i,k) < (ipos(i,k) - delta)) then
             u(i,k) = 0.d0
             su(i,k) = 0.d0
             u(i-1,k) = 0.d0
             su(i-1,k) = 0.d0
             v(i,k) = 0.d0
             sv(i,k) = 0.d0
             v(i,k-1) = 0.d0
             sv(i,k-1) = 0.d0
          endif
          ua(i,k) = u(i,k)
          va(i,k) = v(i,k)
          u0a(i,k) = u0(i,k)
          v0a(i,k) = v0(i,k)
       end do
    end do

  end subroutine momentum

  subroutine bound_u_v(icalcan,icalcen,kcalcan,kcalcen)

    implicit none

    ! calculates boundary values for velocity and transports
    ! from neumann conditions. 

    ! local variables

    integer :: i,k,icalcan,icalcen,kcalcan,kcalcen

    ! southern boundary
    if (kcalcan.le.(domain_kmin+1)) then
       do i = domain_imin,domain_imax - 1
          su(i,domain_kmin) = su(i,domain_kmin+1)
          u(i,domain_kmin) = su(i,domain_kmin+1)*5.0d-1*(pdep(i,domain_kmin+1) + pdep(i+1,domain_kmin+1))
          sv(i,domain_kmin) = sv(i,domain_kmin+1)
          v(i,domain_kmin) = sv(i,domain_kmin+1)*pdep(i,domain_kmin+1)
       end do
    end if

    ! northern boundary
    if (kcalcen.ge.(domain_kmax - 1)) then
       do i = domain_imin,domain_imax - 1
          su(i,domain_kmax) = su(i,domain_kmax-1)
          u(i,domain_kmax) = su(i,domain_kmax-1)*5.0d-1*(pdep(i,domain_kmax-1) + pdep(i+1,domain_kmax-1))
          sv(i,domain_kmax) = sv(i,domain_kmax-1)
          v(i,domain_kmax) = sv(i,domain_kmax-1)*pdep(i,domain_kmax-1)
       end do
    end if

    ! western boundary
    if (icalcan.le.(domain_imin+1)) then
       do k = domain_kmin,domain_kmax-1
          su(domain_imin,k) = su(domain_imin+1,k)
          u(domain_imin,k) = su(domain_imin+1,k)*pdep(domain_imin+1,k)
          sv(domain_imin,k) = sv(domain_imin+1,k)
          v(domain_imin,k) = sv(domain_imin+1,k)*5.0d-1*(pdep(domain_imin+1,k)+pdep(domain_imin+1,k+1))
       end do
    end if

    ! eastern boundary
    if (icalcen.ge.(domain_imax - 1)) then
       do k = domain_kmin, domain_kmax - 1
          su(domain_imax,k) = su(domain_imax-1,k)
          u(domain_imax,k) = su(domain_imax-1,k)*pdep(domain_imax-1,k)
          sv(domain_imax,k) = sv(domain_imax-1,k)
          v(domain_imax,k) = sv(domain_imax-1,k)*5.0d-1*(pdep(domain_imax-1,k)+ &
               pdep(domain_imax-1,k+1))
       end do
    end if

  end subroutine bound_u_v

  subroutine bound_interface(icalcan,icalcen,kcalcan,kcalcen)

    ! calculates boundary values for interface depth
    ! from neumann conditions. 

    implicit none

    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen

    ! local variables
    integer :: i,k
    real(kind=kdp):: pdepold

    ! southern boundary
    if(kcalcan.le.(domain_kmin+1)) then
       do i = domain_imin,domain_imax

          !Paul's new version
          jcw(i,domain_kmin) = jcw(i,domain_kmin + 1)
          jcd_fl(i,domain_kmin) = jcd_fl(i,domain_kmin + 1)
          pdep(i,domain_kmin) = pdep(i,domain_kmin + 1)
          ipos(i,domain_kmin) = bpos(i,domain_kmin) + pdep(i,domain_kmin)

          !          jcd_fl(i,domain_kmin) = 0       
          !          jcw(i,domain_kmin) = 0            
          !          pdepold = pdep(i,domain_kmin)

          !          pdep(i,domain_kmin) = pdep(i,domain_kmin + 1)
          !          ipos(i,domain_kmin) = bpos(i,domain_kmin) - pdep(i,domain_kmin)

          !          if (pdep(i,domain_kmin).ge.dcr) jcw(i,domain_kmin) = 1
          !          if (pdepold.lt.small .and. pdep(i,domain_kmin) .ge. small) &
          !               jcd_fl(i,domain_kmin) = 1      
       end do
    end if

    ! northern boundary
    if(kcalcen.ge.(domain_kmax-1)) then
       do i = domain_imin,domain_imax

          !Paul's version
          jcw(i,domain_kmax) = jcw(i,domain_kmax - 1)
          jcd_fl(i,domain_kmax) = jcd_fl(i,domain_kmax - 1)
          pdep(i,domain_kmax) = pdep(i,domain_kmax - 1)
          ipos(i,domain_kmax) = bpos(i,domain_kmax) - pdep(i,domain_kmax)          

          !          jcd_fl(i,domain_kmax) = 0       
          !          jcw(i,domain_kmax) = 0
          !          pdepold = pdep(i,domain_kmax)        

          !          pdep(i,domain_kmax) = pdep(i,domain_kmax - 1)
          !          ipos(i,domain_kmax) = bpos(i,domain_kmax) - pdep(i,domain_kmax)

          !          if (pdep(i,domain_kmax).ge.dcr) jcw(i,domain_kmax) = 1
          !          if (pdepold < small .and. pdep(i,domain_kmax) > small) &
          !               jcd_fl(i,domain_kmax) = 1      


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

  end subroutine bound_interface

  subroutine bound_tsdc(icalcan,icalcen,kcalcan,kcalcen)

    ! calculates boundary values for temperature, salinity, frazil,
    ! density, and
    ! inflow tracers from neumann conditions. 

    implicit none

    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen
    ! local variables

    integer :: i,k,l


    ! southern boundary

    if(kcalcan.le.(domain_kmin+1)) then

       do i = domain_imin,domain_imax
          temp(i,domain_kmin) = temp(i,domain_kmin+1)                          
          tempa(i,domain_kmin) = tempa(i,domain_kmin+1)                          
          salt(i,domain_kmin) = salt(i,domain_kmin+1)                          
          salta(i,domain_kmin) = salta(i,domain_kmin+1)                          
          rhop(i,domain_kmin) = rhop(i,domain_kmin+1)                          
          tf(i,domain_kmin) = tf(i,domain_kmin+1)                          
       end do

       if (frazil) then
          do i = domain_imin,domain_imax
             ctot(i,domain_kmin) = ctot(i,domain_kmin+1)                          
             do l = 1,nice
		c(i,domain_kmin,l) = c(i,domain_kmin+1,l)
		ca(i,domain_kmin,l) = ca(i,domain_kmin+1,l)
             end do
          end do
       end if

       if (intrace) then
          do i = domain_imin,domain_imax
             do l = ninfmin,ninfmax
		intr(i,domain_kmin,l) = intr(i,domain_kmin+1,l)
		intra(i,domain_kmin,l) = intra(i,domain_kmin+1,l)
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
          tf(i,domain_kmax) = tf(i,domain_kmax-1)                          
       end do

       if (frazil) then
          do i = domain_kmin,domain_kmax
             ctot(i,domain_kmax) = ctot(i,domain_kmax-1)                          
             do l = 1,nice
		c(i,domain_kmax,l) = c(i,domain_kmax-1,l)
		ca(i,domain_kmax,l) = ca(i,domain_kmax-1,l)
             end do
          end do
       end if

       if (intrace) then
          do i = domain_imin,domain_imax
             do l = ninfmin,ninfmax
		intr(i,domain_kmax,l) = intr(i,domain_kmax-1,l)
		intra(i,domain_kmax,l) = intra(i,domain_kmax-1,l)
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
          tf(domain_imin,k) = tf(domain_imin+1,k)                          
       end do

       if (frazil) then
          do k = domain_kmin,domain_kmax
             ctot(domain_imin,k) = ctot(domain_imin+1,k)                          
             do l = 1,nice
		c(domain_imin,k,l) = c(domain_imin+1,k,l)
		ca(domain_imin,k,l) = ca(domain_imin+1,k,l)
             end do
          end do
       end if

       if (intrace) then
          do k = domain_kmin,domain_kmax
             do l = ninfmin,ninfmax
		intr(domain_imin,k,l) = intr(domain_imin+1,k,l)
		intra(domain_imin,k,l) = intra(domain_imin+1,k,l)
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
          tf(domain_imax,k) = tf(domain_imax-1,k)                          
       end do

       if (frazil) then
          do k = domain_kmin,domain_kmax
             ctot(domain_imax,k) = ctot(domain_imax-1,k)                          
             do l = 1,nice
		c(domain_imax,k,l) = c(domain_imax-1,k,l)
		ca(domain_imax,k,l) = ca(domain_imax-1,k,l)
             end do
          end do
       end if

       if (intrace) then
          do k = domain_kmin,domain_kmax
             do l = ninfmin,ninfmax
		intr(domain_imax,k,l) = intr(domain_imax-1,k,l)
		intra(domain_imax,k,l) = intra(domain_imax-1,k,l)
             end do
          end do
       end if

    end if

  end subroutine bound_tsdc


  subroutine update(iwetmin,iwetmax,kwetmin,kwetmax, &
       icalcan,icalcen,kcalcan,kcalcen,negdep)

    ! update depths, depth-averaged velocities (not depth-integrated),
    ! and wetted area

    implicit none

    integer,intent(inout) :: iwetmin,iwetmax,kwetmin,kwetmax
    integer,intent(in) :: icalcan,icalcen,kcalcan,kcalcen   
    ! local variables

    integer :: i,k


    real(kind=kdp) :: iconti,error,pdepold,pdepc,pdepe
    real(kind=kdp) :: dunew,pdepn,dvnew,zd,negdep

    ! update interface position


    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

          if (jcs(i,k) .ne. 1) cycle

          ! interface position corrected for negative depths at new timestep
          iconti = ipos(i,k)
          ipos(i,k) = dmin1(bpos(i,k),ipos(i,k))

          ! report any errors
          error = abs(iconti - ipos(i,k))
          negdep = negdep+error               
          if (error.gt.5.0d-1) then
             write(*,*) 'error: negative depth ',error,' at i=', &
                  i,' k=',k
             write(11,*) 'error: negative depth ',error,' at i=', &
                  i,' k=',k
          end if

          sv(i,k) = 0.d0
          su(i,k) = 0.d0

          ! plume thickness and wetted area (main update of pdep)
          jcd_fl(i,k) = 0      
          pdepold = pdep(i,k)
          jcw(i,k) = 0
          pdep(i,k) = bpos(i,k) - ipos(i,k)

          ! find newly wet area
          if (pdepold .eq. 0.d0 .and. pdep(i,k) .gt. 0.d0) then
             jcd_fl(i,k) = 1 
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

    ! update velocities
    do i = icalcan,icalcen
       do k = kcalcan,kcalcen

          if (jcs(i,k).ne.1) cycle

          jcd_u(i,k) = 0
          jcd_v(i,k) = 0
          pdepc = bpos(i,k) - ipos(i,k)
          ! u-component
          pdepe = bpos(i+1,k) - ipos(i+1,k)
          dunew = 5.0d-1*(pdepe + pdepc)
          if (dunew.gt.small) then
             su(i,k) = u(i,k)/dunew
             jcd_u(i,k) = 1
          endif
          ! v-component
          pdepn = bpos(i,k+1) - ipos(i,k+1)
          dvnew = 5.0d-1*(pdepn + pdepc)
          if (dvnew.gt.small) then
             sv(i,k) = v(i,k)/dvnew
             jcd_v(i,k) = 1
          endif
       end do
    end do

    ! interpolate ambient density field experienced

    !    rhoamb = get_rhoamb_z(gldep + wcdep - ipos)

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
          if (jcs(i,k).ne.1) cycle
          zd = gldep + wcdep - ipos(i,k)
          rhoamb(i,k) = get_rhoamb_z(zd)
       end do
    end do

  end subroutine update


  subroutine scalar(icalcan,kcalcan,icalcen,kcalcen)          

    ! calculates transport and entrainment of all scalars, evolution of 
    ! frazil ice, basal melting and freezing, and inflow tracking

    implicit none

    ! local variables

    integer :: i,k,l,icalcan,kcalcan,icalcen,kcalcen
    integer :: idel,kdel,idx,kdy,ihilf,khilf,mflag,seedindex

    real(kind=kdp),dimension(m_grid,n_grid) :: deltat,deltas
    real(kind=kdp),dimension(m_grid,n_grid):: pdepc,pdepcp,vmid
    real(kind=kdp),dimension(m_grid,n_grid):: umid,speed,bspeed,depth
    real(kind=kdp) :: deltac(m_grid,n_grid,lice),deltatr(m_grid,n_grid,linf)
    real(kind=kdp) :: one,slon,sloe,slos,slow,sumslo
    real(kind=kdp) :: pressure,ttt
    real(kind=kdp) :: tt
    real(kind=kdp) :: dxx,dyy,sx,sy,sxy,r1,r2
    real(kind=kdp) :: dife,difw,difn,difs
    real(kind=kdp) :: dragrt,prden,scden,gambt,gambs,tfb,tfi,c1,c2,c3
    real(kind=kdp) :: gamct,gamcs,ucrit,ucl
    real(kind=kdp) :: fmelttot,wturb,mfac1,mfac2,gi,gim1,mi,mip1
    real(kind=kdp) :: amb_depth

    one = 1.d0
    prden = 12.5d0*pr**(2.d0/3.d0) - 9.d0
    scden = 12.5d0*sc**(2.d0/3.d0) - 9.d0

    ! calculate values to be used for newly-wet cells using weights 
    ! of interface gradient (thus tempa etc is temp to be used for whole plume)

    do i = icalcan,icalcen
       do k = kcalcan,kcalcen
          if (jcd_fl(i,k) .eq. 1) then
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
                      ca(i,k,l) = (c(i,k+1,l)*slon + c(i+1,k,l)*sloe &
                           + c(i,k-1,l)*slos + c(i-1,k,l)*slow)/sumslo
                      ctota(i,k) = ctota(i,k) + ca(i,k,l)
                   end do
                   rhop(i,k) = (1.d0-ctota(i,k))*rhop(i,k)+ctota(i,k)*rhoi
                end if

                if (intrace) then
                   do l = ninfmin,ninfmax
                      intra(i,k,l) =  &
                           (intr(i,k+1,l)*slon + intr(i+1,k,l)*sloe  &
                           + intr(i,k-1,l)*slos + intr(i-1,k,l)*slow)/sumslo
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
          tf(i,k) = freezing_temp_func(salta(i,k),depth(i,k))
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

                deltat(i,k) = deltat(i,k)  &
                     + dt*entr(i,k)*(atemp(i,k) - tempa(i,k))/pdepcp(i,k)
                deltas(i,k) = deltas(i,k)  &
                     + dt*entr(i,k)*(asalt(i,k) - salta(i,k))/pdepcp(i,k)

             endif
          end do
       end do

       if (frazil) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (jcs(i,k) .ne. 1) cycle
                if (pdepc(i,k).gt.edepth) then         
                   do l=1,nice
                      deltac(i,k,l) = deltac(i,k,l)  &
                           + dt*entr(i,k)*(-ca(i,k,l))/pdepcp(i,k)
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
                  + (r1*sy+r2*sx)*(tempa(ihilf,khilf)*jcw(ihilf,khilf) &
                  + (1 - jcw(ihilf,khilf))*tempa(i,k))  &
                  + r1*sxy*(tempa(ihilf,k)*jcw(ihilf,k) &
                  + (1-jcw(ihilf,k))*tempa(i,k)) &
                  - r2*sxy*(tempa(i,khilf)*jcw(i,khilf) &
                  + (1-jcw(i,khilf))*tempa(i,k)) &
                  - (r2*sy+r1*sx)*tempa(i,k)                         

             deltas(i,k) = deltas(i,k)  &
                  + (r1*sy+r2*sx)*(salta(ihilf,khilf)*jcw(ihilf,khilf) &
                  + (1 - jcw(ihilf,khilf))*salta(i,k))  &
                  + r1*sxy*(salta(ihilf,k)*jcw(ihilf,k) &
                  + (1-jcw(ihilf,k))*salta(i,k)) &
                  - r2*sxy*(salta(i,khilf)*jcw(i,khilf) &
                  + (1-jcw(i,khilf))*salta(i,k)) &
                  - (r2*sy+r1*sx)*salta(i,k)                         

             if (frazil) then
                do l = 1,nice
                   deltac(i,k,l) = deltac(i,k,l) &
                        + (r1*sy+r2*sx)*(ca(ihilf,khilf,l)*jcw(ihilf,khilf) &
                        + (1 - jcw(ihilf,khilf))*ca(i,k,l))  &
                        + r1*sxy*(ca(ihilf,k,l)*jcw(ihilf,k)  &
                        + (1-jcw(ihilf,k))*ca(i,k,l)) &
                        - r2*sxy*(ca(i,khilf,l)*jcw(i,khilf) &
                        + (1-jcw(i,khilf))*ca(i,k,l)) &
                        - (r2*sy+r1*sx)*ca(i,k,l)      
                end do
             end if

             if (intrace) then
                do l = ninfmin,ninfmax
                   deltatr(i,k,l) = deltatr(i,k,l) &
                        + (r1*sy+r2*sx)*(intra(ihilf,khilf,l)*jcw(ihilf,khilf) &
                        + (1 - jcw(ihilf,khilf))*intra(i,k,l)) &
                        + r1*sxy*(intra(ihilf,k,l)*jcw(ihilf,k)  &
                        + (1-jcw(ihilf,k))*intra(i,k,l)) &
                        - r2*sxy*(intra(i,khilf,l)*jcw(i,khilf) &
                        + (1-jcw(i,khilf))*intra(i,k,l)) &
                        - (r2*sy+r1*sx)*intra(i,k,l)      
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

             dife = (tempa(i+1,k) - tempa(i,k))*jcw(i+1,k)*kh*rdx(i)
             difw = (tempa(i,k) - tempa(i-1,k))*jcw(i-1,k)*kh*rdx(i-1)
             difn = (tempa(i,k+1) - tempa(i,k))*jcw(i,k+1)*kh*rdy(k)
             difs = (tempa(i,k) - tempa(i,k-1))*jcw(i,k-1)*kh*rdy(k-1)
             deltat(i,k) = deltat(i,k) &
                  + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt 

             dife = (salta(i+1,k) - salta(i,k))*jcw(i+1,k)*kh*rdx(i)
             difw = (salta(i,k) - salta(i-1,k))*jcw(i-1,k)*kh*rdx(i-1)
             difn = (salta(i,k+1) - salta(i,k))*jcw(i,k+1)*kh*rdy(k)
             difs = (salta(i,k) - salta(i,k-1))*jcw(i,k-1)*kh*rdy(k-1)
             deltas(i,k) = deltas(i,k)  &
                  + ((dife - difw)*rdxu(i) + (difn-difs)*rdyv(k))*dt

          end do
       end do

       if (frazil) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen

                if (jcs(i,k) .ne. 1) cycle

                do l = 1,nice
                   dife = (ca(i+1,k,l) - ca(i,k,l))*jcw(i+1,k)*kh*rdx(i)
                   difw = (ca(i,k,l) - ca(i-1,k,l))*jcw(i-1,k)*kh*rdx(i-1)
                   difn = (ca(i,k+1,l) - ca(i,k,l))*jcw(i,k+1)*kh*rdy(k)
                   difs = (ca(i,k,l) - ca(i,k-1,l))*jcw(i,k-1)*kh*rdy(k-1)
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
                        (intra(i+1,k,l) - intra(i,k,l))*jcw(i+1,k)*kh*rdx(i)
                   difw =  &
                        (intra(i,k,l) - intra(i-1,k,l))*jcw(i-1,k)*kh*rdx(i-1)
                   difn =  &
                        (intra(i,k+1,l) - intra(i,k,l))*jcw(i,k+1)*kh*rdy(k)
                   difs =  &
                        (intra(i,k,l) - intra(i,k-1,l))*jcw(i,k-1)*kh*rdy(k-1)
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

                ! find turbulent exchange coefficients
                gambt = (dragrt*bspeed(i,k))/ &
                     (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + prden)
                gambs = (dragrt*bspeed(i,k))/ &
                     (2.12d0*dlog((dragrt*bspeed(i,k)*pdepc(i,k))/nu0) + scden)

                ! calculate freezing point of plume at shelf base,
                ! decide if melt (mflag = 1) 
                ! or freeze and calculate freezing point of ice at shelf
                ! base.  note that heat 
                ! conduction into the ice shelf is taken into account here
                ! if melting (as removal 
                ! of warmest ice steepens gradient) but is not represented
                ! in heat source terms

                tfb = fta*salta(i,k) + ftb + ftc*(gldep+wcdep-bpos(i,k))
                mflag = (1 + int(sign(1.d0,tempa(i,k) - tfb)))/2
                tfi = (1 - mflag)*fta*si  &
                     + ftb + ftc*(gldep + wcdep - bpos(i,k))

                ! calculate coefficients in quadratic to be solved
                c1 = lat/c0 + mflag*(ci/c0)*(tfi-tint(i,k))
                c2 = gambs*(lat/c0 + mflag*(ci/c0)*(tfb-tint(i,k)))  &
                     + gambt*(tfi-tempa(i,k))
                c3 = gambs*gambt*(tfb - tempa(i,k))

                ! calculate melt rate
                bmelt(i,k) = -(c2 - dsqrt(c2*c2 - 4.d0*c1*c3))/(2.d0*c1)

                !! multiply by 10 if freezing
                !              if (mflag.eq.0) then
                !                bmelt(i,k) = bmelt(i,k)*10.d0
                !              end if

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
                     - dt*(lat/c0)*bmelt(i,k)/pdepcp(i,k)

             end if
          end do
       end do

       if (frazil) then
          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                if (pdepc(i,k).gt.mdepth) then         
                   do l=1,nice
                      deltac(i,k,l) = deltac(i,k,l)  &
                           - dt*bmelt(i,k)*ca(i,k,l)/pdepcp(i,k)
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
                   if (ca(i,k,1).ge.cmin(1)) then

                      ! a) calculate growth/melting of frazil 
                      gamct = (nuss(1)*kt)/(ar*r(1))
                      gamcs = (nuss(1)*ks)/(ar*r(1))

                      c1 = 5.0d-1*lat*r(1)/(c0*ca(i,k,1)*pdepc(i,k)*gamct)
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
                           dsqrt((1.0d-1*g*re(1)*(rho0-rhoi))/(rho0*drag(i,k)))
                      ucrit = (1.d0 - bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                      fppn(i,k,1) =  &
                           dmin1(-(rhoi/rho0)*ca(i,k,1)*wi(1)*ucrit,0.d0)

                      ! c) calculate overall frazil source for this 
                      ! size class and 
                      !    total melt for ts source
                      deltac(i,k,1) = deltac(i,k,1)  &
                           - dt*fppn(i,k,1)*ca(i,k,1)/pdepcp(i,k) &
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
                              *(re(1)**3/re(l))*ca(i,k,l)
                         fnuc(i,k,1) = fnuc(i,k,1) - fnuc(i,k,l)
                      end do

                      ! b) calculate growth/melting of frazil 
                      mfac1 = (2.d0*c0*kt/lat)*(tf(i,k) - tempa(i,k))
                      mfac2 = (rhoi/rho0)*pdepc(i,k)

                      if (tempa(i,k).lt.tf(i,k)) then
                         gi = mfac1*nuss(1)*ca(i,k,1)/(r(1)*r(1))
                         fmelt(i,k,1) = mfac2*(vol(1)/(vol(2) - vol(1)))*gi

                         do l=2,nice-1
                            gi = mfac1*nuss(l)*ca(i,k,l)/(r(l)*r(l))
                            gim1 = mfac1*nuss(l-1)*ca(i,k,l-1)/(r(l-1)*r(l-1))
                            fmelt(i,k,l) = mfac2*( &
                                 (vol(l)/(vol(l+1) - vol(l)))*gi  &           
                                 - (vol(l)/(vol(l) - vol(l-1)))*gim1)         
                         end do

                         gim1 = mfac1*nuss(nice-1)*ca(i,k,nice-1) &
                              /(r(nice-1)*r(nice-1))
                         fmelt(i,k,nice) = - mfac2* &
                              (vol(nice)/(vol(nice) - vol(nice-1)))*gim1
                      else
                         mi = mfac1*nuss(1)*ca(i,k,1) &
                              *(1.d0/r(1) + 1.d0/thi(1))/r(1)
                         mip1 = mfac1*nuss(2)*ca(i,k,2) &
                              *(1.d0/r(2) + 1.d0/thi(2))/r(2)
                         fmelt(i,k,1) = mfac2*  &
                              ((vol(1)/(vol(2) - vol(1)))*mip1 - mi)

                         do l=2,nice-1
                            mi = mfac1*nuss(l)*ca(i,k,l) &
                                 *(1.d0/r(l) + 1.d0/thi(l))/r(l)
                            mip1 = mfac1*nuss(l+1)*ca(i,k,l+1)* &
                                 (1.d0/r(l+1) + 1.d0/thi(l+1))/r(l+1)
                            fmelt(i,k,l) = mfac2*( &
                                 (vol(l)/(vol(l+1) - vol(l)))*mip1     &       
                                 - (vol(l)/(vol(l) - vol(l-1)))*mi)  
                         end do

                         mi = mfac1*nuss(nice)*ca(i,k,nice)* &
                              (1.d0/r(nice) + 1.d0/thi(nice))/r(nice)
                         fmelt(i,k,nice) = - mfac2* &
                              (vol(nice)/(vol(nice) - vol(nice-1)))*mi
                      end if

                      ! c) calculate precipitation of frazil
                      do l=1,nice
                         ucl =  &
                              dsqrt((1.0d-1*g*re(l)*(rho0-rhoi))/(rho0*drag(i,k)))
                         ucrit = (1.d0-bspeed(i,k)*bspeed(i,k)/(ucl*ucl))
                         fppn(i,k,l) =  &
                              dmin1(-(rhoi/rho0)*ca(i,k,l)*wi(l)*ucrit,0.d0)
                      end do

                      ! d) calculate overall frazil source for this 
                      ! size class and 
                      !    total melt for ts source
                      do l=1,nice
                         deltac(i,k,l) = deltac(i,k,l) &
                              - dt*fppn(i,k,l)*ca(i,k,l)/pdepcp(i,k) &
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
                     + dt*fmelttot*(tf(i,k) - tempa(i,k))/pdepcp(i,k)
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
                   if (ca(i,k,l) + deltac(i,k,l).ge.0.d0) then
                      c(i,k,l) = ca(i,k,l) + deltac(i,k,l)
                   else
                      c(i,k,l) = 0.d0
                      frzcut(l) = frzcut(l) - ca(i,k,l) - deltac(i,k,l)
                   end if
                end do
             end do
          end do

       else

          do i = icalcan,icalcen
             do k = kcalcan,kcalcen
                do l = 1,nice
                   c(i,k,l) = ca(i,k,l) + deltac(i,k,l)
                end do
             end do
          end do

       end if

       ! total up frazil concentrations
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             ctot(i,k) = 0.d0
             do l = 1,nice
                ctot(i,k) = ctot(i,k) + c(i,k,l) 
             end do
          end do
       end do

    end if

    if (intrace) then
       do i = icalcan,icalcen
          do k = kcalcan,kcalcen
             do l = ninfmin,ninfmax
                intr(i,k,l) = intra(i,k,l) + deltatr(i,k,l)
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
                   if (tempa(i,k).lt.tf(i,k)) then
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
                      c(i,seedindex,l) = cseed(l)
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
                   if (tempa(i,k).gt.tf(i,k)) then
                      jcd_fseed(i,k) = 0
                   else
                      if (jcd_fseed(i,k).eq.0) then
                         ctot(i,k) = 0.d0
                         do l = 1,nice
                            if (c(i,k,l).lt.cseed(l)) c(i,k,l) = cseed(l)
                            ctot(i,k) = ctot(i,k) + c(i,k,l)
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
       ca = c
    end if
    !       do i = icalcan,icalcen
    !          do k = kcalcan,kcalcen
    !             ctota(i,k) = ctot(i,k)
    !             do l = 1,nice
    !                ca(i,k,l) = c(i,k,l)
    !             end do
    !          end do
    !       end do
    !    end if
    !     
    if (intrace) then
       intra = intr
       !       do i = icalcan,icalcen
       !          do k = kcalcan,kcalcen
       !             do l = ninfmin,ninfmax
       !                intra(i,k,l) = intr(i,k,l)
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
                write(*,*) 'warning: plume separated at i=',i,' k=',k
                write(11,*) 'warning: plume separated at i=',i,' k=',k
                sepflag = .true.
             end if
             !           
          end if
       end do
    end do

  end subroutine rho_calc

  elemental real(kind=kdp) function get_rhoamb_z(z)

    implicit none

    real(kind=kdp),intent(in) :: z !depth at which rhoamb is sought

    integer :: izo,izu
    real(kind=kdp) :: difu,difo

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


end module plume

