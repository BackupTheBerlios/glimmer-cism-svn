! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide.f90 - part of the GLIMMER ice model                + 
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
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide
  !*FD the top-level GLIDE module

  use glide_types
  use glide_stop
  use glide_nc_custom
  use glide_io
  use lithot_io
  use lithot
  use glide_profile
  use glimmer_config
  integer, private, parameter :: dummyunit=99

contains

  subroutine glide_config(model,config)
    !*FD read glide configuration from file and print it to the log
    use glide_setup
    use isostasy
    use lithot, only : lithot_readconfig,lithot_printconfig
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file

    type(ConfigSection), pointer :: ncconfig
   
    ! read configuration file
    call glide_readconfig(model,config)
    call glide_printconfig(model)
    ! Read alternate sigma levels from config file, if necessary
    call glide_read_sigma(model,config)
    ! read isostasy configuration file
    call isos_readconfig(model%isos,config)
    call isos_printconfig(model%isos)
    ! read gthf configuration
    call lithot_readconfig(model%lithot,config)
    call lithot_printconfig(model%lithot)
    ! read mapping from config file
    ! **** Use of dew and dns here is an ugly fudge that
    ! **** allows the use of old [GLINT projection] config section
    ! **** for backwards compatibility. It will be deleted soon.
    ! **** (You have been warned!)
    ! **** N.B. Here, dew and dns are unscaled - i.e. real distances in m
    call glimmap_readconfig(model%coordinates%projection,config, &
         model%numerics%dew, &
         model%numerics%dns)

    ! netCDF I/O
    if (trim(model%funits%ncfile).eq.'') then
       ncconfig => config
    else
       call ConfigRead(process_path(model%funits%ncfile),ncconfig)
    end if
    call glimmer_nc_readparams(ncconfig,model%numerics%tstart,model%funits%out_first,model%funits%in_first)
  end subroutine glide_config

  subroutine glide_initialise(model)
    !*FD initialise GLIDE model instance
    use glide_setup
    use glimmer_ncio
    use glide_velo
    use glide_thck
    use glide_temp
    use glimmer_log
    use glimmer_mask
    use glimmer_scales
    use isostasy
    use glimmer_map_init
    use glimmer_horizcoord
    use glide_glenflow, only: glenflow_init,calcflwa
    use glide_tempFullSoln, only: init_tempFullSoln
    use glide_thckADI, only: thckADI_init
    use glide_thckCommon, only: glide_calclsrf
    implicit none
    type(glide_global_type) :: model        !*FD model instance

    character(len=100), external :: glimmer_version_char

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters
    call glide_scale_params(model)
    ! set up coordinate systems
    model%coordinates%ice_grid = horizCoord_new(0.d0, 0.d0, &
         model%numerics%dew, model%numerics%dns, &
         model%general%ewn, model%general%nsn)
    model%coordinates%velo_grid = horizCoord_new(model%numerics%dew/2.,model%numerics%dns/2., &
         model%numerics%dew,model%numerics%dns, &
         model%general%ewn-1,model%general%nsn-1)

    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! allocate arrays
    call glide_allocarr(model)

    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const

    ! set uniform basal heat flux
    model%temper%bheatflx = model%lithot%geot

    ! open all input files
    call openall_in(model%funits%in_first,model%coordinates)
    ! and read first time slice
    call glide_io_readall(model,model)
    ! Write projection info to log
    call glimmap_printproj(model%coordinates%projection)

    ! read lithot if required
    if (model%lithot%do_lithot) then
       call lithot_io_readall(model,model)
    end if

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model%isos,model%numerics%tstart)
    select case(model%options%whichrelaxed)
    case(1) ! Supplied topography is relaxed
       model%isos%relx = model%geometry%topg
    case(2) ! Supplied topography is in equilibrium
       call isos_relaxed(model%isos,model%geometry%topg,model%geometry%thck,model%climate%eus)
    end select

    ! open all output files
    call openall_out(model%funits%out_first,model%coordinates)
    ! create glide variables
    call glide_io_createall(model)

    ! initialise glide components
    call glenflow_init(model%glenflow,model%paramets%fiddle)
    call init_velo(model%velowk,model%paramets%bpar)

    call init_tempFullSoln(    &
         model%tempFullSoln,   &
         model%coordinates, &
         model%numerics%thklim,&
         model%options%periodic_ew)

    call thckADI_init(model%thckADI, &
         model%general%ewn, &
         model%general%nsn, &
         model%general%upn, &
         model%numerics%sigma, &
         model%numerics%dew, &
         model%numerics%dns, &
         model%numerics%thklim, &
         model%numerics%alpha, &
         model%options%basal_mbal, &
         model%options%periodic_ew)

    if (model%lithot%do_lithot) then
       call lithot_io_createall(model)
       call init_lithot(model%lithot,model%numerics%dt,(model%options%hotstart.eq.1))
    end if

    if (model%options%hotstart.ne.1) then
       ! initialise Glen's flow parameter A using an isothermal temperature distribution
       call calcTemp_asSurfTemp(model%temper%temp,model%climate%artm)
       ! Calculate Glenn's A --------------------------------------------------------
       call calcflwa(model%glenflow,          &
            model%temper%flwa,     &
            model%temper%temp,     &
            model%geometry%thck,   &
            model%options%whichflwa, &
            model%numerics%thklim,   &
            model%numerics%sigma) 
    end if

    ! calculate mask
    call glimmer_set_mask( &
         model%geometry%thkmask, &
         model%geometry%thck,    &
         model%geometry%topg,    &
         model%climate%eus,      &
         model%numerics%thklim)

    model%geometry%iarea = calc_iarea(model%geometry%thkmask, &
         model%geometry%thck, &
         model%numerics%dew,  &
         model%numerics%dns)

    model%geometry%ivol = calc_ivol(model%geometry%thkmask, &
         model%geometry%thck, &
         model%numerics%dew,  &
         model%numerics%dns)

    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! initialise profile
#ifdef PROFILING
    call glide_prof_init(model)
#endif
  end subroutine glide_initialise
  
  subroutine glide_tstep_p1(model,time)
    !*FD Performs first part of time-step of an ice model instance.
    !*FD calculate velocity and temperature
    use glimmer_global, only : rk
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_tempFullSoln, only: tstep_tempFullSoln
    use glimmer_mask
    use glide_glenflow, only: calcflwa
    use glimmer_log,    only: write_log, GM_FATAL
    use glimmer_utils,  only: stagvarb
    use glide_bwat,     only: calcbwat
    use glimmer_pmpt,   only: calcbpmp
    use glimmer_deriv, only : df_field_2d_staggered

    implicit none

    type(glide_global_type) :: model        !*FD model instance
    real(rk),  intent(in)   :: time         !*FD Current time in years

    ! Update internal clock
    model%numerics%time=time  
    model%temper%newtemps = .false.

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------     
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%geomderv)
#endif
    call stagvarb(model%geometry% thck, &
         model%geomderv% stagthck,&
         model%general%  ewn, &
         model%general%  nsn)

    call df_field_2d_staggered(model%geometry%usrf, &
         model%numerics%dew, model%numerics%dns, &
         model%geomderv%dusrfdew, & 
         model%geomderv%dusrfdns, &
         .false., .false.)

    call df_field_2d_staggered(model%geometry%thck, &
         model%numerics%dew, model%numerics%dns, &
         model%geomderv%dthckdew, & 
         model%geomderv%dthckdns, &
         .false., .false.)

#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%geomderv)
#endif


    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    if (model%lithot%do_lithot) then
       call calc_lithot(model%lithot,model%geometry%thkmask,model%temper%temp(model%general%upn,:,:),model%climate%artm,model%temper%bheatflx)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%temperature)
#endif
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

       ! Calculate vertical velocities
       call calcVerticalVelocity(model%velowk,model%timederivs,model%options%periodic_ew, model%options%whichwvel, &
            model%numerics%thklim, model%geometry%thck, model%geometry%usrf, model%geomderv%stagthck, &
            model%geomderv%dthckdew, model%geomderv%dthckdns, model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
            model%temper%bmlt, model%climate%acab, model%numerics%time, &
            model%velocity%uvel, model%velocity%vvel, &
            model%velocity%wvel, model%velocity%wgrd)

       ! Do ice temperature calculation
       select case(model%options%whichtemp)
       case(0)
          call calcTemp_asSurfTemp(model%temper%temp,model%climate%artm)
       case(1)
          call tstep_tempFullSoln(    &
               model%tempFullSoln,    &
               model%temper%temp,     &
               model%climate%artm,    &
               model%geometry%thck,   &
               model%geometry%usrf,   &
               model%geometry%thkmask,&
               model%geometry%topg,   &
               model%velocity%uvel,   &
               model%velocity%vvel,   &
               model%velocity%ubas,   &
               model%velocity%vbas,   &
               model%velocity%wvel,   &
               model%velocity%wgrd,   &
               model%temper%flwa,     &
               model%temper%bheatflx, &
               model%temper%bwat,     &
               model%temper%bmlt,     &
               model%numerics%dttem)

          ! Calculate basal water depth ------------------------------------------------

          call calcbwat(model%options%whichbwat,     &
               model%temper%bmlt,                    &
               model%temper%bwat,                    &
               model%geometry%thck,                  &
               model%geometry%topg,                  &
               model%temper%temp(model%general%upn,:,:), &
               is_float(model%geometry%thkmask),     &
               model%numerics%thklim,                &
               model%numerics%dttem,                 &
               model%general%ewn,                    &
               model%general%nsn,                    &
               model%numerics%dew,                   &
               model%numerics%dns,                   &
               model%options%periodic_ew,            &
               model%paramets%bwat_smooth,           &
               model%paramets%hydtim)

          ! now also calculate basal water in velocity coord system

          call stagvarb(model%temper%bwat,  &
               model%temper%stagbwat,       &
               model%general%ewn,           &
               model%general%nsn)

          ! Transform basal temperature and pressure melting point onto velocity grid -

          call stagvarb(model%temper%temp(model%general%upn,1:model%general%ewn,1:model%general%nsn), &
               model%temper%stagbtemp,      &
               model%general%ewn,           &
               model%general%nsn)
       
          ! Calculate the basal pressure melting point and stagger
          call calcbpmp(model%geometry%thck,model%temper%bpmp)

          call stagvarb(model%temper%bpmp,  &
               model%temper%stagbpmp,       &
               model%general%ewn,           &
               model%general%nsn)

       case(2)
          call calcTemp_VerticalProfile( &
               model%temper%temp,    &
               model%climate%artm,   &
               model%numerics%sigma, &
               model%geometry%thck,  &
               model%temper%bwat)
       case default
          call write_log('Unknown temperature option',GM_FATAL)
       end select

       ! Calculate Glenn's A --------------------------------------------------------
       call calcflwa(model%glenflow,          &
            model%temper%flwa,     &
            model%temper%temp,     &
            model%geometry%thck,   &
            model%options%whichflwa, &
            model%numerics%thklim,   &
            model%numerics%sigma) 
       model%temper%newtemps = .true.
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%temperature)
#endif

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 
    call calc_btrc(model%velowk,model%options%whichbtrc,model%velocity%bed_softness, &
         model%temper%stagbwat, model%temper%stagbtemp, model%temper%stagbpmp, &
         model%temper%bmlt, model%isos%relx, model%velocity%btrc)

  end subroutine glide_tstep_p1


  subroutine glide_tstep_p2(model,no_write)
    !*FD Performs second part of time-step of an ice model instance.
    !*FD write data and move ice
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glimmer_mask
    use glide_thckADI, only: thckADI_tstep, get_uflx, get_vflx
    use isostasy
    implicit none

    type(glide_global_type) :: model        !*FD model instance
    logical,optional :: no_write

    logical nw

    ! ------------------------------------------------------------------------ 
    ! write to netCDF file
    ! ------------------------------------------------------------------------ 

    if (present(no_write)) then
       nw=no_write
    else
       nw=.false.
    end if

    if (.not. nw) then
       call glide_io_writeall(model,model)
       if (model%lithot%do_lithot) then
          call lithot_io_writeall(model,model)
       end if
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_evo)
#endif
    select case(model%options%whichevol)
    case(0) ! Use precalculated uflx, vflx -----------------------------------

       call thck_nonlin_evolve(    &
            model%thckADI,         &
            model%geometry%thck,   &
            model%geometry%usrf,   &
            model%geometry%lsrf,   &
            model%geometry%topg,   &
            model%climate%acab,    &
            model%temper%bmlt,     &
            model%velocity%btrc,   &
            model%temper%flwa,     &
            model%velocity%uvel,   &
            model%velocity%vvel,   &
            model%velocity%diffu,  &
            model%velocity%ubas,   &
            model%velocity%vbas,   &
            model%climate%eus,     &
            model%temper%newtemps, &
            .true.,                &
            model%numerics%dt)

    case(1) ! Use explicit leap frog method with uflx,vflx -------------------

       call thckADI_tstep(model%thckADI, &
            model%geometry%thck, &
            model%climate%acab,  &
            model%geometry%lsrf, &
            model%geometry%usrf, &
            model%geometry%topg, &
            model%velocity%btrc, &
            model%velocity%ubas, &
            model%velocity%vbas, &
            model%temper%bmlt,   &
            model%temper%flwa,   &
            model%velocity%uvel, &
            model%velocity%vvel, &
            model%velocity%diffu,&
            model%numerics%dt,   &
            model%climate%eus,   &
            model%temper%newtemps)

    case(2) ! Use non-linear calculation that incorporates velocity calc -----

       call thck_nonlin_evolve(    &
            model%thckADI,         &
            model%geometry%thck,   &
            model%geometry%usrf,   &
            model%geometry%lsrf,   &
            model%geometry%topg,   &
            model%climate%acab,    &
            model%temper%bmlt,     &
            model%velocity%btrc,   &
            model%temper%flwa,     &
            model%velocity%uvel,   &
            model%velocity%vvel,   &
            model%velocity%diffu,  &
            model%velocity%ubas,   &
            model%velocity%vbas,   &
            model%climate%eus,     &
            model%temper%newtemps, &
            .false.,               &
            model%numerics%dt)

    end select

    ! Retrieve fluxes from thickness model
    call get_uflx(model%thckADI,model%velocity%uflx)
    call get_vflx(model%thckADI,model%velocity%vflx)

#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_evo)
#endif

    ! ------------------------------------------------------------------------
    ! get new mask
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask2)
#endif
    call glimmer_set_mask( &
         model%geometry%thkmask, &
         model%geometry%thck,    &
         model%geometry%topg,    &
         model%climate%eus,      &
         model%numerics%thklim)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask2)
#endif

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 
    call glide_marinlim(model%options%  whichmarn, &
         model%geometry% thck,      &
         model%isos% relx,      &
         model%geometry%topg,   &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         model%climate%calving)

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos_water)
#endif
    if (model%isos%do_isos) then
       if (model%numerics%time.ge.model%isos%next_calc) then
          model%isos%next_calc = model%isos%next_calc + model%isos%period
          call isos_icewaterload(model%isos,model%geometry%topg,model%geometry%thck,model%climate%eus)
       end if
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos_water)
#endif
    
    ! basal shear stress calculations
    call calc_basal_shear(model%geomderv%stagthck,model%geomderv%dusrfdew,model%geomderv%dusrfdns,model%velocity%tau_x,model%velocity%tau_y)
  end subroutine glide_tstep_p2

  subroutine glide_tstep_p3(model)
    !*FD Performs third part of time-step of an ice model instance:
    !*FD calculate isostatic adjustment and upper and lower ice surface
    use isostasy
    use glide_setup
    use glide_thckCommon, only: glide_calclsrf
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    
    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos)
#endif
    if (model%isos%do_isos) then
       call isos_isostasy(model%isos,model%geometry%topg,model%numerics%dt)
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos)
#endif

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

  end subroutine glide_tstep_p3

  !-------------------------------------------------------------------

!MH!  subroutine glide_write_mod_rst(rfile)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart_common
!MH!
!MH!#ifdef RESTARTS
!MH!    use glide_types
!MH!    use isostasy_types
!MH!#endif
!MH!
!MH!    type(restart_file) :: rfile      !*FD Open restart file 
!MH!
!MH!#ifdef RESTARTS
!MH!    call glide_types_modrsw(rfile)
!MH!    call isostasy_types_modrsw(rfile)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_write_mod_rst
!MH!
!MH!  !-------------------------------------------------------------------
!MH!
!MH!  subroutine glide_read_mod_rst(rfile)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart_common
!MH!
!MH!#ifdef RESTARTS
!MH!    use glide_types
!MH!    use isostasy_types
!MH!#endif
!MH!
!MH!    type(restart_file) :: rfile      !*FD Open restart file 
!MH!
!MH!#ifdef RESTARTS
!MH!    call glide_types_modrsr(rfile)
!MH!    call isostasy_types_modrsr(rfile)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_read_mod_rst
!MH!
!MH!  !-------------------------------------------------------------------
!MH!
!MH!  subroutine glide_write_restart(model,rfile)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart
!MH!    use glimmer_restart_common
!MH!    implicit none
!MH!
!MH!    type(glide_global_type) :: model !*FD model instance
!MH!    type(restart_file) :: rfile      !*FD Open restart file     
!MH!
!MH!#ifdef RESTARTS
!MH!    call glimmer_write_mod_rst(rfile)
!MH!    call glide_write_mod_rst(rfile)
!MH!    call rsw_glide_global_type(rfile,model)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_write_restart

!MH!  !-------------------------------------------------------------------
!MH!
!MH!  subroutine glide_read_restart(model,rfile,prefix)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart
!MH!    use glimmer_restart_common
!MH!    use glimmer_ncdf
!MH!    use glimmer_ncio
!MH!    implicit none
!MH!
!MH!    type(glide_global_type) :: model !*FD model instance
!MH!    type(restart_file) :: rfile      !*FD Open restart file 
!MH!    character(*),optional,intent(in) :: prefix !*FD prefix for new output files
!MH!
!MH!    character(40) :: pf
!MH!
!MH!    if (present(prefix)) then
!MH!       pf = prefix
!MH!    else
!MH!       pf = 'RESTART_'
!MH!    end if
!MH!
!MH!#ifdef RESTARTS
!MH!    call glimmer_read_mod_rst(rfile)
!MH!    call glide_read_mod_rst(rfile)
!MH!    call rsr_glide_global_type(rfile,model)
!MH!    call nc_repair_outpoint(model%funits%out_first)
!MH!    call nc_repair_inpoint(model%funits%in_first)
!MH!    call nc_prefix_outfiles(model%funits%out_first,trim(pf))
!MH!    call openall_out(model)
!MH!    call glide_io_createall(model)
!MH!    call glide_nc_fillall(model)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_read_restart

end module glide
