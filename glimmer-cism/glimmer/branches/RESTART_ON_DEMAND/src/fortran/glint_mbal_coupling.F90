! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_mbal_coupling.f90 - part of the GLIMMER ice model  + 
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
#include <config.inc>
#endif

module glint_mbal_coupling

  use glint_mbal
  use glimmer_config

  !*FD Module to handle the accumulation of inputs and calculation of mass-balance

  type glint_mbc
     real(dp),dimension(:,:),pointer :: prcp_save  => null() !*FD used to accumulate precip
     real(dp),dimension(:,:),pointer :: ablt_save  => null() !*FD used to accumulate ablation
     real(dp),dimension(:,:),pointer :: acab_save  => null() !*FD used to accumulate mass-balance
     real(dp),dimension(:,:),pointer :: artm_save  => null() !*FD used to average air-temperature
     real(sp),dimension(:,:),pointer :: snowd      => null() !*FD Keeps track of snow depth
     real(sp),dimension(:,:),pointer :: siced      => null() !*FD Keeps track of superimposed ice depth 
     real(sp),dimension(:,:),pointer :: prcp       => null() !*FD Instantaneous precip
     real(sp),dimension(:,:),pointer :: ablt       => null() !*FD Instantaneous ablation
     real(sp),dimension(:,:),pointer :: acab       => null() !*FD Instantaneous mass-balance
     real(sp),dimension(:,:),pointer :: artm       => null() !*FD Instantaneous air temperature
     real(sp),dimension(:,:),pointer :: xwind      => null() !*FD Instantaneous x-wind
     real(sp),dimension(:,:),pointer :: ywind      => null() !*FD Instantaneous y-wind
     real(sp),dimension(:,:),pointer :: humidity   => null() !*FD Instantaneous humidity
     real(sp),dimension(:,:),pointer :: SWdown     => null() !*FD Instantaneous sw down
     real(sp),dimension(:,:),pointer :: LWdown     => null() !*FD Instantaneous lw down
     real(sp),dimension(:,:),pointer :: Psurf      => null() !*FD Instantaneous psurf
     real(sp),dimension(:,:),pointer :: snowd_save => null() !*FD Saves snow depth
     real(sp),dimension(:,:),pointer :: siced_save => null() !*FD Saves superimposed ice depth 
     integer :: av_count =0 !*FD Counter for averaging temperature input
     logical :: new_accum=.true.
     type(glint_mbal_params) :: mbal
  end type glint_mbc

  !MAKE_RESTART
#ifdef RESTARTS
#define RST_GLINT_MBAL_COUPLING
#include "glimmer_rst_head.inc"
#undef RST_GLINT_MBAL_COUPLING
#endif

contains

#ifdef RESTARTS
#define RST_GLINT_MBAL_COUPLING
#include "glimmer_rst_body.inc"
#undef RST_GLINT_MBAL_COUPLING
#endif

  subroutine glint_mbc_init(params,lgrid,config,whichacab,snowd,siced,nx,ny,dx)

    use glimmer_coordinates

    type(glint_mbc)  :: params
    type(coordsystem_type) :: lgrid
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    integer          :: whichacab
    real(sp),dimension(:,:),intent(in) :: snowd !*FD Initial snow-depth field
    real(sp),dimension(:,:),intent(in) :: siced !*FD Initial superimposed ice field
    integer          :: nx,ny  !*FD grid dimensions (for SMB)
    real(rk)         :: dx     !* Grid length (for SMB)

    ! Deallocate if necessary

    if (associated(params%prcp_save))  deallocate(params%prcp_save)
    if (associated(params%ablt_save))  deallocate(params%ablt_save)
    if (associated(params%acab_save))  deallocate(params%acab_save)
    if (associated(params%artm_save))  deallocate(params%artm_save)
    if (associated(params%snowd))      deallocate(params%snowd)
    if (associated(params%siced))      deallocate(params%siced)
    if (associated(params%prcp))       deallocate(params%prcp)
    if (associated(params%ablt))       deallocate(params%ablt)
    if (associated(params%acab))       deallocate(params%acab)
    if (associated(params%artm))       deallocate(params%artm)
    if (associated(params%xwind))      deallocate(params%xwind)
    if (associated(params%ywind))      deallocate(params%ywind)
    if (associated(params%humidity))   deallocate(params%humidity)
    if (associated(params%SWdown))     deallocate(params%SWdown)
    if (associated(params%LWdown))     deallocate(params%LWdown)
    if (associated(params%Psurf))      deallocate(params%Psurf)
    if (associated(params%snowd_save)) deallocate(params%snowd_save)
    if (associated(params%siced_save)) deallocate(params%siced_save)

    ! Allocate arrays and zero

    call coordsystem_allocate(lgrid,params%prcp_save);  params%prcp_save = 0.0
    call coordsystem_allocate(lgrid,params%ablt_save);  params%ablt_save = 0.0
    call coordsystem_allocate(lgrid,params%acab_save);  params%acab_save = 0.0
    call coordsystem_allocate(lgrid,params%artm_save);  params%artm_save = 0.0
    call coordsystem_allocate(lgrid,params%snowd);      params%snowd = 0.0
    call coordsystem_allocate(lgrid,params%siced);      params%siced = 0.0
    call coordsystem_allocate(lgrid,params%prcp);       params%prcp = 0.0
    call coordsystem_allocate(lgrid,params%ablt);       params%ablt = 0.0
    call coordsystem_allocate(lgrid,params%acab);       params%acab = 0.0
    call coordsystem_allocate(lgrid,params%artm);       params%artm = 0.0
    call coordsystem_allocate(lgrid,params%xwind);      params%xwind = 0.0
    call coordsystem_allocate(lgrid,params%ywind);      params%ywind = 0.0
    call coordsystem_allocate(lgrid,params%humidity);   params%humidity = 0.0
    call coordsystem_allocate(lgrid,params%SWdown);     params%SWdown = 0.0
    call coordsystem_allocate(lgrid,params%LWdown);     params%LWdown = 0.0
    call coordsystem_allocate(lgrid,params%Psurf);      params%Psurf = 0.0
    call coordsystem_allocate(lgrid,params%snowd_save); params%snowd_save = 0.0
    call coordsystem_allocate(lgrid,params%siced_save); params%siced_save = 0.0

    ! Initialise the mass-balance scheme and other components

    call glint_mbal_init(params%mbal,config,whichacab,nx,ny,dx)

    ! Copy snow and ice depths if relevant

    if (mbal_has_snow_model(params%mbal)) then
       params%snowd=snowd
       params%siced=siced
    end if

  end subroutine glint_mbc_init

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_accumulate(params,artm,arng,prcp,snowd,siced,xwind,ywind,local_orog, &
       thck,humidity,SWdown,LWdown,Psurf)

    type(glint_mbc)  :: params
    real(sp),dimension(:,:),intent(inout) :: artm      !*FD Mean air temperature (degC)
    real(sp),dimension(:,:),intent(in) :: arng         !*FD Air temperature half-range (degC)
    real(sp),dimension(:,:),intent(inout) :: prcp      !*FD Precipitation (m)
    real(sp),dimension(:,:),intent(in) :: snowd        !*FD Snow depth (m)
    real(sp),dimension(:,:),intent(in) :: siced        !*FD Superimposed ice depth (m)
    real(rk),dimension(:,:),intent(in) :: xwind        !*FD $x$-component of surface winds (m/s)
    real(rk),dimension(:,:),intent(in) :: ywind        !*FD $y$-component of surface winds (m/s)
    real(sp),dimension(:,:),intent(in) :: local_orog   !*FD Local orography (m)
    real(rk),dimension(:,:),intent(in) :: thck         !*FD Ice thickness (m)
    real(rk),dimension(:,:),intent(in) :: humidity     !*FD Relative humidity (%)
    real(rk),dimension(:,:),intent(in) :: SWdown       !*FD Downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in) :: LWdown       !*FD Downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in) :: Psurf        !*FD Surface pressure (Pa)

    real(sp),dimension(size(artm,1),size(artm,2)) :: ablt,acab

    ! Things to do the first time

    if (params%new_accum) then

       params%new_accum=.false.
       params%av_count =0

       ! Initialise 

       params%snowd=snowd
       params%siced=siced

       params%prcp_save=0.0
       params%ablt_save=0.0
       params%acab_save=0.0
       params%artm_save=0.0

    end if

    params%av_count=params%av_count+1

    ! Call mass-balance

    call glint_mbal_calc(params%mbal,artm,arng,prcp,(local_orog>0.0.or.thck>0.0),params%snowd, &
         params%siced,ablt,acab,thck,xwind,ywind,humidity,SWdown,LWdown,Psurf) 

    ! Accumulate

    params%prcp_save = params%prcp_save + prcp
    params%ablt_save = params%ablt_save + ablt
    params%acab_save = params%acab_save + acab
    params%artm_save = params%artm_save + artm

    ! Copy instantaneous fields

    params%prcp=prcp
    params%ablt=ablt
    params%acab=acab
    params%artm=artm
    params%xwind=xwind
    params%ywind=ywind
    params%humidity=humidity
    params%SWdown=SWdown
    params%LWdown=LWdown
    params%Psurf=Psurf

  end subroutine glint_accumulate

  ! +++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_get_mbal(params,artm,prcp,ablt,acab,snowd,siced)

    type(glint_mbc)  :: params
    real(sp),dimension(:,:),intent(out)   :: artm   !*FD Mean air temperature (degC)
    real(sp),dimension(:,:),intent(out)   :: prcp   !*FD Precipitation (m)
    real(sp),dimension(:,:),intent(out)   :: ablt   !*FD Ablation
    real(sp),dimension(:,:),intent(out)   :: acab   !*FD Mass-balance
    real(sp),dimension(:,:),intent(inout) :: snowd  !*FD Snow depth (m)
    real(sp),dimension(:,:),intent(inout) :: siced  !*FD Superimposed ice depth (m)

    if (.not.params%new_accum) then
       params%artm_save=params%artm_save/real(params%av_count)
    end if

    params%new_accum=.true.

    artm=params%artm_save
    prcp=params%prcp_save
    ablt=params%ablt_save
    acab=params%acab_save
    snowd=params%snowd
    siced=params%siced

    where (snowd<0.0) snowd=0.0
    where (siced<0.0) siced=0.0

  end subroutine glint_get_mbal

end module glint_mbal_coupling
