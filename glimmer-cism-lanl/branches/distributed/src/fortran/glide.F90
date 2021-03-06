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
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glide
  !*FD the top-level GLIDE module

  use glide_types
  use glide_stop
  use glide_nc_custom
  use glide_io
  use glide_lithot
  use glide_profile
  use glide_deriv
  use glimmer_config
  use glimmer_global
  integer, private, parameter :: dummyunit=99

contains

  subroutine glide_config(model,config)
    !*FD read glide configuration from file and print it to the log
    use glide_setup
    use isostasy
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

    ! read mapping from config file
    ! **** Use of dew and dns here is an ugly fudge that
    ! **** allows the use of old [GLINT projection] config section
    ! **** for backwards compatibility. It will be deleted soon.
    ! **** (You have been warned!)
    ! **** N.B. Here, dew and dns are unscaled - i.e. real distances in m

    call glimmap_readconfig(model%projection,config, &
         model%numerics%dew, &
         model%numerics%dns)

    ! netCDF I/O
    if (trim(model%funits%ncfile).eq.'') then
       ncconfig => config
    else
       call ConfigRead(process_path(model%funits%ncfile),ncconfig)
    end if
    call glimmer_nc_readparams(model,ncconfig)

  end subroutine glide_config

  subroutine glide_initialise(model)
    !*FD initialise GLIDE model instance
    use parallel
    use glide_setup
    use glimmer_ncio
    use glide_velo
    use glide_velo_higher
    use glide_thck
    use glide_temp
    use glissade_temp
    use glimmer_log
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_ground

    ! *sfp** added
    use glam_strs2, only : glam_velo_fordsiapstr_init

!whl - remap_glamutiles to be removed
    use remap_glamutils, only : horizontal_remap_init

    ! *sfp** added for summer modeling school
    use fo_upwind_advect, only : fo_upwind_advect_init

    !*mb* added 
    use glam_Basal_Proc, only : Basal_Proc_init

    implicit none
    type(glide_global_type) :: model        !*FD model instance

    integer :: i,j

    call write_log(glimmer_version)

    ! initialise scales
    call glimmer_init_scales
    call initnan !Initialize the NAN representation, hack to get smart compilers like gfortran to divide by zero

    ! scale parameters
    call glide_scale_params(model)
    ! set up coordinate systems
    ! time to change to the parallel values of ewn and nsn
    call distributed_grid(model%general%ewn,model%general%nsn)
    model%general%ice_grid = coordsystem_new(0.d0, 0.d0, &
         model%numerics%dew, model%numerics%dns, &
         model%general%ewn, model%general%nsn)
    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.,model%numerics%dns/2., &
         model%numerics%dew,model%numerics%dns, &
         model%general%ewn-1,model%general%nsn-1)

    ! allocate arrays
    call glide_allocarr(model)

    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const
  
    !Initialize boundary condition fields to be NaN everywhere
    model%geometry%marine_bc_normal = NaN

    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! set uniform basal heat flux
    model%temper%bheatflx = model%paramets%geot
    
    ! open all input files
    call openall_in(model)
    ! and read first time slice
    call glide_io_readall(model,model)
    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model)
    select case(model%options%whichrelaxed)
    case(1) ! Supplied topography is relaxed
       model%isos%relx = model%geometry%topg
    case(2) ! Supplied topography is in equilibrium
       call not_parallel(__FILE__,__LINE__)
       call isos_relaxed(model)
    end select


    ! open all output files
    call openall_out(model)
    ! create glide variables
    call glide_io_createall(model)

    ! initialise glide components
    call init_velo(model)

    if (model%options%whichtemp == TEMP_REMAP_ADV) then
       call glissade_init_temp(model)
    else
       call glide_init_temp(model)
    endif

    call init_thck(model)

    call glide_initialise_backstress(model%geometry%thck,&
                                     model%climate%backstressmap,&
                                     model%climate%backstress, &
                                     model%climate%stressin, &
                                     model%climate%stressout)
    if (model%options%gthf.gt.0) then
       call not_parallel(__FILE__,__LINE__)
       call init_lithot(model)
    end if

    if (model%options%which_ho_diagnostic == HO_DIAG_PATTYN_UNSTAGGERED .or. &
        model%options%which_ho_diagnostic == HO_DIAG_PATTYN_STAGGERED) then

        call init_velo_hom_pattyn(model)

    end if

    ! *sfp** added; initialization of Payne/Price HO dynamics subroutine ... name can change once agreed on
    if (model%options%which_ho_diagnostic == HO_DIAG_PP ) then

        call glam_velo_fordsiapstr_init(model%general%ewn,    model%general%nsn,  &
                                        model%general%upn,                        &
                                        model%numerics%dew,   model%numerics%dns, &
                                        model%numerics%sigma)
    end if

!whl - Subroutine horizontal_remap_init is not needed for the new remapping scheme.

    ! *sfp** added; initialization of LANL incremental remapping subroutine for thickness evolution
    if (model%options%whichevol== EVOL_INC_REMAP ) then

        if (model%options%whichtemp == TEMP_REMAP_ADV) then ! Use IR to advect temperature

#ifdef JEFFORIG
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       model%general%ewn,  model%general%nsn,   &
                                       model%general%upn,  model%numerics%sigma )    
#endif
           !JEFF In distributed, horizontal remapping is serialized, so it needs to be initialized with full grid size
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       global_ewn,  global_nsn,   &
                                       model%general%upn,  model%numerics%sigma )

        else  ! Use IR to transport thickness only
#ifdef JEFFORIG
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       model%general%ewn,  model%general%nsn)
#endif
           !JEFF In distributed, horizontal remapping is serialized, so it needs to be initialized with full grid size
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       global_ewn,  global_nsn)

        endif ! whichtemp

    endif 

    ! *sfp** added for summer modeling school
    if (model%options%whichevol== EVOL_FO_UPWIND ) then

        ! JEFF - Review for distributed. OK for parallel Trilinos. call not_parallel(__FILE__,__LINE__)
        call fo_upwind_advect_init( model%general%ewn, model%general%nsn )

    endif

    ! *mb* added; initialization of basal proc. module
    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
        model%options%which_bmod == BAS_PROC_FASTCALC) then
        
        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
                              model%numerics%ntem)
    end if      

    ! initialise ice age
    ! Currently the ice age is only computed for remapping transport
    ! (whichevol = 3 or 4)
    model%geometry%age(:,:,:) = 0._dp
 
    if (model%options%hotstart.ne.1) then
       ! initialise Glen's flow parameter A using an isothermal temperature distribution
       ! JEFF - Review for distributed.  Distributed is Ok for 0 method called here. 
       ! Method 1 involves derivs, therefore marked as not_parallel. 
       ! call not_parallel(__FILE__,__LINE__)
! KJE: Jeff your comments above refer to the commented line below? Bill has changed to glide_temp_driver FYI
!       call timeevoltemp(model,0)
       call glide_temp_driver(model,0)
    end if

    ! calculate mask
    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)
    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

    !calculate the normal at the marine margin
    call glide_marine_margin_normal(model%geometry%thck, model%geometry%thkmask, model%geometry%marine_bc_normal)

    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! initialise profile
#ifdef PROFILING
    call glide_prof_init(model)
#endif
    
    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)
  end subroutine glide_initialise
  
  subroutine glide_tstep_p1(model,time)
    !*FD Performs first part of time-step of an ice model instance.
    !*FD calculate velocity and temperature
    use parallel
    use glimmer_global, only : rk
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glissade_temp
    use glide_mask
    use glide_thckmask
    use glide_grids

    ! *mb* added for basal proc module  
    use glam_Basal_Proc, only : Basal_Proc_driver

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
    call geometry_derivs(model)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%geomderv)
#endif

#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask1)
#endif
    !TREY This sets local values of dom, mask, totpts, and empty
    call glide_maskthck(&
         model%geometry% thck,      &
         model%climate%  acab,      &
         .true.,                    &
         model%numerics%thklim,     &
         model%geometry% dom,       &
         model%geometry% mask,      &
         model%geometry% totpts,    &
         model%geometry% empty)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask1)
#endif

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    if (model%options%gthf.gt.0) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%temperature)
#endif
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

       if (model%options%whichtemp == TEMP_REMAP_ADV) then 

         ! Vert diffusion and strain heating only; no advection
         ! Remapping routine is used to advect temperature in glide_tstep_p2

         call glissade_temp_driver(model)

       else

         ! standard Glide driver, including temperature advection

         call glide_temp_driver(model, model%options%whichtemp)

       endif

       model%temper%newtemps = .true.

    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%temperature)
#endif

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 
    call calc_btrc(model,model%options%whichbtrc,model%velocity%btrc)

    ! ------------------------------------------------------------------------ 
    ! Calculate basal shear strength from Basal Proc module, if necessary
    ! ------------------------------------------------------------------------    
    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
        model%options%which_bmod == BAS_PROC_FASTCALC) then
        call Basal_Proc_driver (model%general%ewn,model%general%nsn,model%general%upn,       &
                                model%numerics%ntem,model%velocity_hom%uvel(model%general%upn,:,:), &
                                model%velocity_hom%vvel(model%general%upn,:,:), &
                                model%options%which_bmod,model%temper%bmlt,model%basalproc)
    end if

  end subroutine glide_tstep_p1


  subroutine glide_tstep_p2(model,no_write)
    !*FD Performs second part of time-step of an ice model instance.
    !*FD write data and move ice
    use parallel
    use glide_thck
    use glide_velo
    use glide_ground
    use glide_setup
    use glide_temp
    use glide_mask
    use isostasy

    ! *sfp** driver module/subroutines for Payne/Price HO dynamics and LANL inc. remapping for dH/dt 
    ! Modeled after similar routines in "glide_thck"
    use glam, only: inc_remap_driver

    ! *sfp** added for summer modeling school
    use fo_upwind_advect, only: fo_upwind_advect_driver

    ! *sfp* added so that stress tensor is populated w/ HO stress fields
    use stress_hom, only: glide_stress

    implicit none

    type(glide_global_type) :: model        !*FD model instance
    logical,optional :: no_write

    logical nw
    integer :: dummy

    ! ------------------------------------------------------------------------ 
    ! write to netCDF file
    ! ------------------------------------------------------------------------ 
    if (present(no_write)) then
       nw=no_write
    else
       nw=.false.
    end if 
    if (.not. nw) call glide_io_writeall(model,model)
    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_evo)
#endif
    select case(model%options%whichevol)
    case(EVOL_PSEUDO_DIFF) ! Use precalculated uflx, vflx -----------------------------------

       call not_parallel(__FILE__,__LINE__)
       call thck_lin_evolve(model,model%temper%newtemps)

    case(EVOL_ADI) ! Use explicit leap frog method with uflx,vflx -------------------

       call not_parallel(__FILE__,__LINE__)
       call stagleapthck(model,model%temper%newtemps)

    case(EVOL_DIFFUSION) ! Use non-linear calculation that incorporates velocity calc -----

       call not_parallel(__FILE__,__LINE__)
       call thck_nonlin_evolve(model,model%temper%newtemps)

    case(EVOL_INC_REMAP) ! Use incremental remapping scheme for advecting ice thickness ---
                         ! (and temperature too, if whichtemp = TEMP_REMAP_ADV)
       call inc_remap_driver( model )

       !JEFF Gathers required for glide_stress
       call distributed_gather_var(model%geomderv%dusrfdew, gathered_dusrfdew)
       call distributed_gather_var(model%geomderv%dusrfdns, gathered_dusrfdns)
       call distributed_gather_var(model%geomderv%dthckdew, gathered_dthckdew)
       call distributed_gather_var(model%geomderv%dthckdns, gathered_dthckdns)

       ! Tau is calculated in glide_stress and initialized in glide_types.
       ! These gathers effective just create variables of the correct size to hold the values.
       ! If I had a better mechanism, then I could just create the variables locally,
       ! then distribute at the end.
       call distributed_gather_var(model%velocity_hom%tau%xx, gathered_tauxx)
       call distributed_gather_var(model%velocity_hom%tau%yy, gathered_tauyy)
       call distributed_gather_var(model%velocity_hom%tau%xy, gathered_tauxy)
       call distributed_gather_var(model%velocity_hom%tau%scalar, gathered_tauscalar)
       call distributed_gather_var(model%velocity_hom%tau%xz, gathered_tauxz)
       call distributed_gather_var(model%velocity_hom%tau%yz, gathered_tauyz)

       if (main_task) then
          call glide_stress( model )       !*sfp* added for populating stress tensor w/ HO fields
       endif

       call parallel_barrier   ! Other tasks hold here until main_task completes

    ! *sfp** added for summer modeling school
    case(EVOL_FO_UPWIND) ! Use first order upwind scheme for mass transport
       call not_parallel(__FILE__,__LINE__)
       call fo_upwind_advect_driver( model )
 
    end select
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_evo)
#endif

    ! ------------------------------------------------------------------------
    ! get new mask
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask2)
#endif

#ifdef JEFFORIG
    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)
#endif
    call distributed_gather_var(model%geometry%topg, gathered_topg)
    call distributed_gather_var(model%geometry%thkmask, gathered_thkmask)

    if (main_task) then
	   call glide_set_mask(model%numerics, gathered_thck, gathered_topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        gathered_thkmask, model%geometry%iarea, model%geometry%ivol, exec_serial=.TRUE.)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask2)
#endif

    !calculate the normal at the marine margin
    call distributed_gather_var(model%geometry%marine_bc_normal, gathered_marine_bc_normal)

#ifdef JEFFORIG
    call glide_marine_margin_normal(model%geometry%thck, model%geometry%thkmask, model%geometry%marine_bc_normal)
#endif

    if (main_task) then
       call glide_marine_margin_normal(gathered_thck, gathered_thkmask, gathered_marine_bc_normal, exec_serial=.TRUE.)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes

    !calculate the grounding line flux after the mask is correct
    call distributed_gather_var(model%velocity%surfvel, gathered_surfvel)
    call distributed_gather_var(model%ground%gline_flux, gathered_gline_flux)
    call distributed_gather_var(model%velocity%ubas, gathered_ubas)
    call distributed_gather_var(model%velocity%vbas, gathered_vbas)

#ifdef JEFFORIG
    call calc_gline_flux(model%geomderv%stagthck,model%velocity%surfvel, &
                         model%geometry%thkmask,model%ground%gline_flux, model%velocity%ubas, &
                         model%velocity%vbas, model%numerics%dew)
#endif
    if (main_task) then
       call calc_gline_flux(gathered_stagthck, gathered_surfvel, &
                         gathered_thkmask, gathered_gline_flux, gathered_ubas, &
                         gathered_vbas, model%numerics%dew)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes

   ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 
#ifdef JEFFORIG
    call glide_marinlim(model%options%whichmarn, &
         model%geometry%thck,      &
         model%isos%relx,      &
         model%geometry%topg,   &
         model%temper%flwa,   &
         model%numerics%sigma,   &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         model%climate%calving,  &
         model%climate%backstress, &
         model%climate%tempanmly, &
         model%numerics%dew,    &
         model%numerics%dns, &
         model%climate%backstressmap, &
         model%climate%stressout, &
         model%climate%stressin, &
         model%ground, &
         model%general%nsn, &
         model%general%ewn, &
         model%geometry%usrf)
#endif
    call distributed_gather_var(model%isos%relx, gathered_relx)
    call distributed_gather_var(model%temper%flwa, gathered_flwa)
    call distributed_gather_var(model%climate%calving, gathered_calving)
    call distributed_gather_var(model%climate%backstress, gathered_backstress)
    call distributed_gather_var(model%geometry%usrf, gathered_usrf)
    call distributed_gather_var(model%climate%backstressmap, gathered_backstressmap)

    if (main_task) then
       call glide_marinlim(model%options%whichmarn, &
         gathered_thck,      &
         gathered_relx,      &
         gathered_topg,   &
         gathered_flwa,   &
         model%numerics%sigma,   &
         gathered_thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         gathered_calving,  &
         gathered_backstress, &
         model%climate%tempanmly, &
         model%numerics%dew,    &
         model%numerics%dns, &
         gathered_backstressmap, &
         model%climate%stressout, &
         model%climate%stressin, &
         model%ground, &
         model%general%nsn, &
         model%general%ewn, &
         gathered_usrf)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes

    !issues with ice shelf, calling it again fixes the mask
#ifdef JEFFORIG
   call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)
#endif
    if (main_task) then
	   call glide_set_mask(model%numerics, gathered_thck, gathered_topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        gathered_thkmask, model%geometry%iarea, model%geometry%ivol, exec_serial=.TRUE.)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes

#ifdef JEFFORIG
    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)
#endif
    if (main_task) then
	    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
	                            model%geometry%iarea, gathered_thkmask, &
	                            model%geometry%iareaf, model%geometry%iareag)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos_water)
#endif
    if (model%isos%do_isos) then
       !JEFF the isos_icewaterload() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       if (model%numerics%time.ge.model%isos%next_calc) then
          model%isos%next_calc = model%isos%next_calc + model%isos%period
          call isos_icewaterload(model)
          model%isos%new_load = .true.
       end if
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos_water)
#endif

    ! basal shear stress calculations
#ifdef JEFFORIG
    call calc_basal_shear(model)
#endif
    call distributed_gather_var(model%velocity%tau_x, gathered_tau_x)
    call distributed_gather_var(model%velocity%tau_y, gathered_tau_y)

    if (main_task) then
	    call calc_basal_shear(gathered_stagthck, &
	                          gathered_dusrfdew, &
	                          gathered_dusrfdns, &
	                          gathered_tau_x, &
	                          gathered_tau_y)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes

  end subroutine glide_tstep_p2

  subroutine glide_tstep_p3(model)
    !*FD Performs third part of time-step of an ice model instance:
    !*FD calculate isostatic adjustment and upper and lower ice surface
    use parallel
    use isostasy
    use glide_setup
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    
    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos)
#endif
    if (model%isos%do_isos) then
       !JEFF the isos_isostasy() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       call isos_isostasy(model)
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos)
#endif

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------
#ifdef JEFFORIG
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)
#endif
    call distributed_gather_var(model%geometry%lsrf, gathered_lsrf)

    if (main_task) then
	    call glide_calclsrf(gathered_thck, gathered_topg, model%climate%eus, gathered_lsrf)
	    gathered_usrf = max(0.d0,gathered_thck + gathered_lsrf)
    endif

    call parallel_barrier   ! Other tasks hold here until main_task completes

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

  end subroutine glide_tstep_p3

  !-------------------------------------------------------------------

  subroutine glide_write_mod_rst(rfile)

    use glimmer_log
    use glimmer_restart_common

#ifdef RESTARTS
    use glide_types
    use isostasy_types
#endif

    type(restart_file) :: rfile      !*FD Open restart file 

#ifdef RESTARTS
    call glide_types_modrsw(rfile)
    call isostasy_types_modrsw(rfile)
#else
    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
#endif

  end subroutine glide_write_mod_rst

  !-------------------------------------------------------------------

  subroutine glide_read_mod_rst(rfile)

    use glimmer_log
    use glimmer_restart_common

#ifdef RESTARTS
    use glide_types
    use isostasy_types
#endif

    type(restart_file) :: rfile      !*FD Open restart file 

#ifdef RESTARTS
    call glide_types_modrsr(rfile)
    call isostasy_types_modrsr(rfile)
#else
    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
#endif

  end subroutine glide_read_mod_rst

  !-------------------------------------------------------------------

  subroutine glide_write_restart(model,rfile)

    use glimmer_log
    use glimmer_restart
    use glimmer_restart_common
    implicit none

    type(glide_global_type) :: model !*FD model instance
    type(restart_file) :: rfile      !*FD Open restart file     

#ifdef RESTARTS
    call glimmer_write_mod_rst(rfile)
    call glide_write_mod_rst(rfile)
    call rsw_glide_global_type(rfile,model)
#else
    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
#endif

  end subroutine glide_write_restart

  !-------------------------------------------------------------------

  subroutine glide_read_restart(model,rfile,prefix)

    use glimmer_log
    use glimmer_restart
    use glimmer_restart_common
    use glimmer_ncdf
    use glimmer_ncio
    implicit none

    type(glide_global_type) :: model !*FD model instance
    type(restart_file) :: rfile      !*FD Open restart file 
    character(*),optional,intent(in) :: prefix !*FD prefix for new output files

    character(40) :: pf

    if (present(prefix)) then
       pf = prefix
    else
       pf = 'RESTART_'
    end if

#ifdef RESTARTS
    call glimmer_read_mod_rst(rfile)
    call glide_read_mod_rst(rfile)
    call rsr_glide_global_type(rfile,model)
    call nc_repair_outpoint(model%funits%out_first)
    call nc_repair_inpoint(model%funits%in_first)
    call nc_prefix_outfiles(model%funits%out_first,trim(pf))
    call openall_out(model)
    call glide_io_createall(model)
    call glide_nc_fillall(model)
#else
    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
#endif

  end subroutine glide_read_restart

end module glide
