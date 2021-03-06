! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_mbal.f90 - part of the GLIMMER ice model           + 
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

module glint_mbal

  use glimmer_pdd
  use glimmer_daily_pdd
  use glimmer_global

#ifdef USE_ENMABAL
  use smb_mecons
#else
  use smb_dummy
#endif

  implicit none

  !*FD Unified wrapper for different mass-balance codes

  type glint_mbal_params
     type(glimmer_pdd_params),      pointer :: annual_pdd => null() !*FD Pointer to annual PDD params
     type(glimmer_daily_pdd_params),pointer :: daily_pdd => null()  !*FD Pointer to daily PDD params
     type(smb_params),              pointer :: smb => null()        !*FD Pointer to SMB params
     integer :: which !*FD Flag for chosen mass-balance type
     integer :: tstep !*FD Timestep of mass-balance scheme in hours
  end type glint_mbal_params

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLINT_MBAL
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLINT_MBAL
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLINT_MBAL
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLINT_MBAL
!MH!#endif

  subroutine glint_mbal_init(params,config,which,nx,ny,dxr)

    use glimmer_config
    use glimmer_log
    use glint_constants

    !*FD Initialise mass-balance schemes

    type(glint_mbal_params)      :: params !*FD parameters to be initialised
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    integer,intent(in)           :: which  !*FD selector for pdd type
    integer                      :: nx,ny  !*FD grid dimensions (for SMB)
    real(rk)                     :: dxr    !* Grid length (for SMB)

    ! Copy selector

    params%which=which

    ! Deallocate if necessary

    if (associated(params%annual_pdd)) deallocate(params%annual_pdd)
    if (associated(params%daily_pdd))  deallocate(params%daily_pdd)
    if (associated(params%smb))        deallocate(params%smb)

    ! Allocate desired type and initialise
    ! Also check we have a valid value of which

    select case(which)
    case(1)
       allocate(params%annual_pdd)
       call glimmer_pdd_init(params%annual_pdd,config)
       params%tstep=years2hours
    case(2)
       params%tstep=years2hours
    case(3)
       allocate(params%smb)
       params%tstep=6
       call SMBInitWrapper(params%smb,nx,ny,nint(dxr),params%tstep*60,'/data/ggdagw/src/smb/smb_config/online')
    case(4)
       allocate(params%daily_pdd)
       call glimmer_daily_pdd_init(params%daily_pdd,config)
       params%tstep=days2hours
    case default
       call write_log('Invalid value of whichacab',GM_FATAL,__FILE__,__LINE__)
    end select

  end subroutine glint_mbal_init

  ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_mbal_calc(params,artm,arng,prcp,landsea,snowd,siced,ablt,acab, &
       thck,U10m,V10m,humidity,SWdown,LWdown,Psurf)

    use glimmer_log

    type(glint_mbal_params)                 :: params  !*FD parameters to be initialised
    real(sp), dimension(:,:), intent(in)    :: artm    !*FD Mean air-temperature ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: arng    !*FD Temperature half-range ($^{\circ}$C)
    real(sp), dimension(:,:), intent(in)    :: prcp    !*FD Accumulated precipitation (m)
    logical,  dimension(:,:), intent(in)    :: landsea !*FD Land-sea mask (land is TRUE)
    real(sp), dimension(:,:), intent(inout) :: snowd   !*FD Snow depth (m)
    real(sp), dimension(:,:), intent(inout) :: siced   !*FD Superimposed ice depth (m)
    real(sp), dimension(:,:), intent(out)   :: ablt    !*FD Ablation (m)
    real(sp), dimension(:,:), intent(out)   :: acab    !*FD Mass-balance (m)
    real(rk), dimension(:,:), intent(in)    :: thck    !*FD Ice thickness (m)
    real(rk), dimension(:,:), intent(in)    :: U10m    !*FD Ten-metre x-wind (m/s)
    real(rk), dimension(:,:), intent(in)    :: V10m    !*FD Ten-metre y-wind (m/s)
    real(rk), dimension(:,:), intent(in)    :: humidity !*FD Relative humidity (%)
    real(rk), dimension(:,:), intent(in)    :: SWdown  !*FD Downwelling shortwave (W/m^2)
    real(rk), dimension(:,:), intent(in)    :: LWdown  !*FD Downwelling longwave (W/m^2)
    real(rk), dimension(:,:), intent(in)    :: Psurf   !*FD Surface pressure (Pa)

    real(rk),dimension(size(acab,1),size(acab,2)) :: acab_temp

    select case(params%which)
    case(1)
       call glimmer_pdd_mbal(params%annual_pdd,artm,arng,prcp,ablt,acab,landsea) 
    case(2) 
       acab = prcp
    case(3)
       ! The energy-balance model will go here...
       ! NB SLM will be thickness array...
       call SMBStepWrapper(params%smb,acab_temp,thck,real(artm,rk),real(prcp*1000.0,rk),U10m,V10m,humidity,SWdown,LWdown,Psurf)
       acab=acab_temp
       acab=acab/1000.0  ! Convert to metres
       ablt=prcp-acab    ! Construct ablation field (in m)
       ! Fix according to land-sea mask
       where (.not.landsea)
          ablt=prcp
          acab=0.0
          snowd=0.0
          siced=0.0
       end where
    case(4)
       call glimmer_daily_pdd_mbal(params%daily_pdd,artm,arng,prcp,snowd,siced,ablt,acab,landsea)
    end select

  end subroutine glint_mbal_calc

  logical function mbal_has_snow_model(params)

    type(glint_mbal_params)      :: params 

    if (params%which==4) then
       mbal_has_snow_model=.true.
    else
       mbal_has_snow_model=.false.
    end if

  end function mbal_has_snow_model

end module glint_mbal
