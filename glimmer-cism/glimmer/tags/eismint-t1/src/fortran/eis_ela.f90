! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eis_ela.f90 - part of the GLIMMER ice model              + 
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

module eis_ela
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  
  !*FD This module reproduces the ELA forcing used to drive the 
  !*FD Edinburgh Ice Sheet model. 

  use glimmer_ts
  use glimmer_global, only : fname_length

  type eis_ela_type
     !*FD Parameters for the EIS climate forcing
     real :: ela_a = 10821.                      !*FD ELA paramters
     real :: ela_b = -238.                       !*FD ELA paramters
     real :: ela_c = 1.312                       !*FD ELA paramters
     real :: zmax = 1200.                        !*FD parameters describing how the MB
     real :: bmax = 1.5                          !*FD parameters describing how the MB
                                                 !*FD varies around the ELA
     character(len=fname_length) :: fname=''     !*FD name of file containing ELA ts
     type(glimmer_tseries) :: ela_ts             !*FD ELA time series 
     real,dimension(:,:),pointer :: ela => null()!*FD ELA field
  end type eis_ela_type

  private :: ela_lat, calc_mb

contains
  subroutine eis_ela_config(config,ela)
    !*FD get ELA configuration from config file
    use glimmer_config
    implicit none
    type(eis_ela_type)           :: ela     !*FD ela data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'EIS ELA')
    if (associated(section)) then
       call GetValue(section,'ela_file',ela%fname)
       call GetValue(section,'zmax',ela%zmax)
       call GetValue(section,'bmax',ela%bmax)
       call GetValue(section,'ela_a',ela%ela_a)
       call GetValue(section,'ela_b',ela%ela_b)
       call GetValue(section,'ela_c',ela%ela_c)
    end if
  end subroutine eis_ela_config
    
  subroutine eis_ela_printconfig(ela)
    !*FD print configuration to log
    use glimmer_log
    implicit none
    type(eis_ela_type)      :: ela   !*FD ela data
    ! local variables
    character(len=100) :: message
    call write_log('EIS ELA')
    call write_log('-------')
    write(message,*) 'ELA file: ',trim(ela%fname)
    call write_log(message)
    write(message,*) 'zmax    : ',ela%zmax
    call write_log(message)
    write(message,*) 'bmax    : ',ela%bmax
    call write_log(message)
    write(message,*) 'ela A   : ',ela%ela_a
    call write_log(message)
    write(message,*) 'ela B   : ',ela%ela_B
    call write_log(message)
    write(message,*) 'ela C   : ',ela%ela_C
    call write_log(message)
    call write_log('')
  end subroutine eis_ela_printconfig

  subroutine eis_init_ela(ela,model)
    !*FD initialise ELA forcing
    use glide_types
    use physcon, only : scyr
    use paramets, only: thk0, tim0
    implicit none
    type(eis_ela_type)      :: ela   !*FD ela data
    type(glide_global_type) :: model !*FD model instance

    call glimmer_read_ts(ela%ela_ts,ela%fname)
    
    ! scale parameters
    ela%ela_ts%values = ela%ela_ts%values/thk0
    ela%zmax = ela%zmax/thk0
    ela%bmax = ela%bmax*tim0/(scyr * thk0)
    ela%ela_a = ela%ela_a/thk0
    ela%ela_b = ela%ela_b/thk0
    ela%ela_c = ela%ela_c/thk0

    ! calculate shape of mass balance field
    allocate(ela%ela(model%general%ewn,model%general%nsn))
    ela%ela = ela_lat(ela%ela_a,ela%ela_b,ela%ela_c,model%climate%lati)
  end subroutine eis_init_ela
    
  subroutine eis_massbalance(ela,model,time)
    !*FD calculate mass balance
    use glide_types
    use glimmer_global, only : rk
    implicit none
    type(eis_ela_type)        :: ela   !*FD ela data
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time

    ! local variables
    real :: ela_time = 0.

    call glimmer_ts_step(ela%ela_ts,real(time),ela_time)
    model%climate%acab = calc_mb(ela%ela+ela_time, &
         model%geometry%topg, &
         model%geometry%thck, &
         model%climate%eus, ela%zmax,ela%bmax)
  end subroutine eis_massbalance

  !*****************************************************************************
  ! private procedures
  !*****************************************************************************
  elemental function calc_mb(ela,topo,thick,eus,zmax,bmax)
    !*FD calculate mass balance
    use glimmer_global, only : dp
    implicit none
    real, intent(in) :: ela       !*FD equilibrium line altitude
    real(kind=dp), intent(in) :: topo      !*FD topography
    real(kind=dp), intent(in) :: thick     !*FD ice thickness
    real, intent(in) :: eus       !*FD eustatic sea level
    real, intent(in) :: zmax,bmax !*FD parameters describing MB variation around ELA
    real calc_mb

    ! local variables
    real z

    if (topo.ge.eus .or. thick.gt.0) then
       z = topo+thick-eus
       if (z.lt.0.) then
          calc_mb = -1.    ! ablation on ice shelf
          return
       end if
       z = z - ela
       if ((zmax.le.0.01).or.(z .ge. zmax)) then
          ! first condition allows forcing of mb=bmax.
          calc_mb = bmax
          return
       else
          calc_mb = 2.*z*bmax/zmax - z*z*bmax/(zmax*zmax)
          return
       end if
    else
       calc_mb = 0.
       return
    end if
  end function calc_mb

  elemental function ela_lat(a,b,c,lat)
    !*FD calculate ELA variation with latitude
    implicit none
    real, intent(in) :: a,b,c !*FD shape of ELA field
    real,intent(in)  :: lat   !*FD latitude
    real ela_lat

    ela_lat = a +  b * lat +  c * lat * lat
  end function ela_lat
  
end module eis_ela
